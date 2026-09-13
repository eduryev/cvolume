#!/usr/bin/env python3
"""Fill in the fields the AVT digest emails don't carry, from the move.nl listing pages.

The listing pages expose a structured "kenmerken" table (label/value pairs grouped
under section headers). We read energy class, garden/terrace, sale status and tenure
from it. Requires the cloud environment's network policy to allow move.nl.

    python3 househunt/enrich.py [--only-missing] [--limit N]

Writes the fields back into househunt/data/listings.json, resumable: it saves after
every listing, so interrupting and re-running picks up where it left off.
"""
import json, re, html, sys, time, os, urllib.request, urllib.error

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = f'{ROOT}/househunt/data/listings.json'
UA = 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120 Safari/537.36'
DELAY = 0.7

TOKEN = re.compile(
    r'kenmerk-section-header">(.*?)</div>'
    r'|kenmerk-label">(.*?)</div><div class="kenmerk-value">(.*?)</div>')

def strip(s):
    return html.unescape(re.sub(r'<[^>]+>', '', s or '')).strip()

def kenmerken(page):
    """-> {section: {label: value}} in document order."""
    out, section = {}, '?'
    for m in TOKEN.finditer(page):
        if m.group(1):
            section = strip(m.group(1))
            out.setdefault(section, {})
        else:
            out.setdefault(section, {})[strip(m.group(2))] = strip(m.group(3))
    return out

def flat(k):
    """All label/value pairs, section-independent — move.nl files Balkon under
    'Indeling' but Tuin under 'Buitenruimte'."""
    out = {}
    for pairs in k.values():
        out.update(pairs)
    return out

def present(value):
    """A kenmerk value that means 'yes, there is one'."""
    if not value:
        return False
    v = value.lower().strip()
    return not (v.startswith('geen') or v.startswith('niet') or v in ('nee', 'onbekend', '-'))

def outdoor(k):
    """(garden, terrace). The 'Tuin' field names the TYPE of outdoor space, so
    'Zonneterras' is a terrace, not a garden, and 'Geen tuin' is neither."""
    f = flat(k)
    garden = terrace = False
    for key in ('Tuin', 'Hoofdtuin'):
        v = (f.get(key) or '').lower()
        if not present(v):
            continue
        if 'terras' in v:
            terrace = True
        elif 'tuin' in v or 'patio' in v:
            garden = True
    if present(f.get('Balkon')) or present(f.get('Dakterras')):
        terrace = True
    return garden, terrace

def fetch(url, tries=2):
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': UA})
            with urllib.request.urlopen(req, timeout=45) as r:
                return r.read().decode('utf8', 'replace')
        except urllib.error.HTTPError as e:
            return None if e.code in (404, 410) else (time.sleep(2) or None if attempt else None)
        except Exception:
            if attempt + 1 == tries:
                return None
            time.sleep(2)
    return None

def extract(page):
    k = kenmerken(page)
    energie = k.get('Energie', {})
    buiten = k.get('Buitenruimte', {})
    over = k.get('Overdracht', {})
    kad = k.get('Kadastrale gegevens', {})

    rec = {}
    cls = energie.get('Energieklasse') or energie.get('Energielabel')
    if cls:
        cls = cls.strip().upper()
        if re.fullmatch(r'A\++|[A-G]', cls):
            rec['energy'] = cls
    if not cls:                                   # fall back to the free-text description
        m = re.search(r'Energielabel\s*:?\s*(A\+*|[B-G])\b', page)
        if m:
            rec['energy'] = m.group(1).upper()

    rec['garden'], rec['terrace'] = outdoor(k)
    if rec['garden'] and buiten.get('Oppervlakte hoofdtuin'):
        rec['garden_size'] = buiten['Oppervlakte hoofdtuin'].split('(')[0].strip()

    if over.get('Status'):
        rec['status'] = over['Status']
    if over.get('Bijdrage VVE p/m'):
        rec['vve'] = over['Bijdrage VVE p/m']
    if kad.get('Eigendomssituatie'):
        rec['tenure'] = kad['Eigendomssituatie']
    rec['enriched'] = time.strftime('%Y-%m-%d')
    return rec

def main():
    only_missing = '--only-missing' in sys.argv
    limit = None
    if '--limit' in sys.argv:
        limit = int(sys.argv[sys.argv.index('--limit') + 1])

    listings = json.load(open(DATA))
    todo = [r for r in listings if not (only_missing and r.get('enriched'))]
    if limit:
        todo = todo[:limit]
    print(f'{len(todo)} listings to enrich', flush=True)

    ok = fail = 0
    for i, r in enumerate(todo, 1):
        page = fetch(r['url'])
        if not page:
            fail += 1
            r['enrich_failed'] = True
            print(f'{i:>3}/{len(todo)} FAIL {r["address"]}', flush=True)
        else:
            r.pop('enrich_failed', None)
            r.update(extract(page))
            ok += 1
            print(f'{i:>3}/{len(todo)} ok   {r["address"]:<40} '
                  f'{r.get("energy","-"):<4} G={int(r["garden"])} T={int(r["terrace"])} '
                  f'{r.get("status","")}', flush=True)
        json.dump(listings, open(DATA, 'w'), indent=1, ensure_ascii=False)
        time.sleep(DELAY)

    print(f'\ndone: {ok} enriched, {fail} failed', flush=True)

if __name__ == '__main__':
    main()
