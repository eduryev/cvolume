#!/usr/bin/env python3
"""Build docs/index.html: the AVT house-hunt board with its listing data inlined.

    python3 househunt/build.py

Re-run after parse_emails.py picks up new AVT digests.
"""
import json, os, datetime

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CUTOFF = '2026-09-10'          # anything first seen before this lands in Archived

def main():
    listings = json.load(open(f'{ROOT}/househunt/data/listings.json'))
    pc4 = json.load(open(f'{ROOT}/househunt/pc4.json'))

    slim = []
    for r in listings:
        slim.append({
            'id': r['id'],
            'address': r['address'],
            'postcode': r['postcode'],
            'url': r['url'],
            'image': r.get('image'),
            'price': r['price'],
            'm2': r['m2'],
            'rooms': r['rooms'],
            'bedrooms': r.get('bedrooms'),
            'type': r.get('type'),
            'first_seen': r['first_seen'],
            'energy': r.get('energy'),
            'garden': r.get('garden'),
            'terrace': r.get('terrace'),
            'garden_size': r.get('garden_size'),
            'status': r.get('status'),
            'tenure': r.get('tenure'),
            'enriched': r.get('enriched'),
            'defaultStage': 'backlog' if r['first_seen'] >= CUTOFF else 'archived',
        })
    slim.sort(key=lambda r: (r['first_seen'], r['address']))

    tpl = open(f'{ROOT}/househunt/template.html').read()
    html = (tpl
            .replace('__LISTINGS__', json.dumps(slim, ensure_ascii=False, separators=(',', ':')))
            .replace('__PC4__', json.dumps(pc4, ensure_ascii=False, separators=(',', ':')))
            .replace('__BUILT__', datetime.date.today().isoformat())
            .replace('__CUTOFF__', CUTOFF))

    os.makedirs(f'{ROOT}/docs/househunt', exist_ok=True)
    out = f'{ROOT}/docs/househunt/index.html'
    open(out, 'w').write(html)

    live = sum(1 for r in slim if r['defaultStage'] != 'archived')
    withE = sum(1 for r in slim if r.get('energy'))
    print(f'{out}: {len(slim)} listings ({live} backlog / {len(slim)-live} archived), '
          f'{withE} with energy class, {os.path.getsize(out)//1024} KB')

if __name__ == '__main__':
    main()
