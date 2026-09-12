#!/usr/bin/env python3
"""Parse AVT Makelaars 'Er zijn nieuwe objecten gevonden' digest emails
into a normalised listing set.

Input : JSON files produced by the Gmail get_thread tool (FULL_CONTENT).
Output: househunt/data/listings.json
"""
import json, re, sys, html, glob, os, unicodedata

CARD_RE = re.compile(r'<table class="real_estate_table".*?</table>\s*</td>\s*</tr>\s*</tbody></table>', re.S)
HREF_RE = re.compile(r'href="(https://move\.nl/exchange-object/[^"]+)"')
IMG_RE = re.compile(r'<img[^>]*?src="([^"]+)"[^>]*?alt="([^"]*)"|<img[^>]*?alt="([^"]*)"[^>]*?src="([^"]+)"')
MATCH_RE = re.compile(r'Match:\s*(?:<[^>]+>\s*)*?(\d+)%', re.S)
PRICE_RE = re.compile(r'Vraagprijs:\s*&euro;|Vraagprijs:\s*€\s*([\d.]+)')
DETAIL_RE = re.compile(r'>([A-Za-zÀ-ÿ\- ]+?)\s*\|\s*(\d+)\s*m²\s*(?:/\s*(\d+)\s*m²\s*)?\|\s*(\d+)\s*kamers?(?:\s*\((\d+)\s*slaapkamers?\))?</a>')
ADDR_ANCHOR_RE = re.compile(
    r'font-weight:bold"[^>]*>\s*(.*?)\s*<span[^>]*>\s*(\d{4}\s*[A-Z]{2})\s+([^<]+?)\s*</span>', re.S)

def clean(s):
    s = html.unescape(s or '')
    s = s.replace('\r', ' ').replace('\n', ' ').replace('\t', ' ')
    return re.sub(r'\s+', ' ', s).strip()

def slug(s):
    s = unicodedata.normalize('NFKD', s.lower())
    s = ''.join(c for c in s if not unicodedata.combining(c))
    return re.sub(r'[^a-z0-9]+', '-', s).strip('-')

def parse_card(card):
    rec = {}
    m = HREF_RE.search(card)
    if m:
        rec['url'] = clean(m.group(1))

    # address: prefer the img alt ("Street 12 H, 1075 CN Amsterdam"), fall back to anchor text
    alt = None
    im = IMG_RE.search(card)
    if im:
        src, a1, a2, src2 = im.groups()
        alt = a1 or a2
        rec['image'] = clean(src or src2)
    if alt and ',' in alt:
        addr, rest = alt.split(',', 1)
        rec['address'] = clean(addr)
        pc = re.match(r'\s*(\d{4}\s*[A-Z]{2})\s+(.+)', clean(rest))
        if pc:
            rec['postcode'], rec['city'] = clean(pc.group(1)), clean(pc.group(2))
    if 'address' not in rec:
        am = ADDR_ANCHOR_RE.search(card)
        if am:
            rec['address'] = clean(am.group(1))
            rec['postcode'] = clean(am.group(2))
            rec['city'] = clean(am.group(3))
    if 'address' not in rec:
        return None

    mm = MATCH_RE.search(card)
    if mm:
        rec['match'] = int(mm.group(1))

    pm = re.search(r'(Vraagprijs|Koopsom):\s*(?:&euro;|€)\s*([\d.]+)', card)
    if pm:
        rec['price'] = int(pm.group(2).replace('.', ''))
        rec['price_label'] = pm.group(1)

    dm = DETAIL_RE.search(card)
    if dm:
        rec['type'] = clean(dm.group(1))
        rec['m2'] = int(dm.group(2))
        if dm.group(3):
            rec['plot_m2'] = int(dm.group(3))
        rec['rooms'] = int(dm.group(4))
        rec['bedrooms'] = int(dm.group(5)) if dm.group(5) else None
    return rec

def main(paths, out):
    listings = {}
    msg_count = card_count = 0
    for path in paths:
        with open(path) as f:
            blob = json.load(f)
        messages = blob.get('messages', [blob]) if isinstance(blob, dict) else blob
        for msg in messages:
            body = msg.get('htmlBody') or ''
            if 'real_estate_table' not in body:
                continue
            msg_count += 1
            date = msg.get('date', '')[:10]
            for card in CARD_RE.findall(body):
                rec = parse_card(card)
                if not rec:
                    continue
                card_count += 1
                key = slug(rec['address'] + '-' + rec.get('postcode', ''))
                rec['id'] = key
                rec['first_seen'] = date
                rec['last_seen'] = date
                if key in listings:
                    prev = listings[key]
                    prev['first_seen'] = min(prev['first_seen'], date)
                    prev['last_seen'] = max(prev['last_seen'], date)
                    # keep the richest record
                    for k, v in rec.items():
                        if k not in ('first_seen', 'last_seen') and v and not prev.get(k):
                            prev[k] = v
                else:
                    listings[key] = rec
    out_list = sorted(listings.values(), key=lambda r: (r['first_seen'], r['address']))
    with open(out, 'w') as f:
        json.dump(out_list, f, indent=1, ensure_ascii=False)
    print(f"messages with cards: {msg_count}  cards seen: {card_count}  unique listings: {len(out_list)}")
    missing = [r['address'] for r in out_list if not r.get('url') or not r.get('price') or not r.get('m2')]
    if missing:
        print(f"incomplete ({len(missing)}): {missing[:10]}")

if __name__ == '__main__':
    main(sorted(glob.glob(sys.argv[1])), sys.argv[2])
