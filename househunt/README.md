# AVT House Hunt

A kanban board + map over the apartment listings AVT Makelaars emails to
eduryev@gmail.com (cc katiazoritch@gmail.com). Published as a static page at
`docs/househunt/index.html`.

## Pipeline

```
Gmail digests  ──parse_emails.py──>  data/listings.json  ──build.py──>  docs/househunt/index.html
```

| file | what it does |
|---|---|
| `parse_emails.py` | Parses AVT "Er zijn nieuwe objecten gevonden" digests into one deduplicated record per address. Handles both `Vraagprijs:` and `Koopsom:` pricing, and the `living m² / plot m²` form used for houses. |
| `data/listings.json` | 153 unique listings from 95 digest emails (27 Aug – 11 Sep 2026). |
| `pc4.json` | Amsterdam PC4 → neighbourhood name + approximate centroid. Used to place pins instantly before real geocoding resolves. |
| `template.html` | The page. `__LISTINGS__`, `__PC4__`, `__BUILT__`, `__CUTOFF__` are substituted at build time. |
| `build.py` | Inlines the data and writes the single self-contained page. |

## Rebuilding

`parse_emails.py` reads Gmail `get_thread` JSON dumps (FULL_CONTENT format):

```sh
python3 househunt/parse_emails.py "/path/to/thread-dumps/*.txt" househunt/data/listings.json
python3 househunt/build.py
```

`CUTOFF` in `build.py` (currently `2026-09-10`) decides what starts in Archived.

## Fields

Address, move.nl link, price, m², rooms/bedrooms and photo all come from the
emails. **Energy label, garden (G) and terrace (T) are not in the emails** — they
start blank and are filled in from the board's Edit panel.

## State and sharing

Board state (column, energy, G/T, notes) lives in the viewer's `localStorage`.
"Share board" exports/merges JSON; the merge keeps whichever edit is newer per
listing. Committing that JSON as `docs/househunt/board.json` makes it the shared
starting point for both viewers.

Geocoding runs in the visitor's browser against PDOK (the Dutch national
address service) and is cached locally; until it resolves, pins sit at their
postcode centroid.
