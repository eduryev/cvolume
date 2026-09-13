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
| `enrich.py` | Reads each move.nl listing page for the fields the emails omit: energy class, garden/terrace, sale status, tenure. Resumable. |
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

Address, move.nl link, price, m², rooms/bedrooms and photo come from the emails.
Energy class, garden (G), terrace (T), sale status and tenure are **not** in the
emails and come from `enrich.py` reading the move.nl listing page.

The listing pages server-render a structured `kenmerk-label` / `kenmerk-value`
table, so no browser or JS execution is needed — plain HTTP is enough. This
requires the cloud environment's network policy to allow `move.nl`
(**Custom** network access with `move.nl` and `*.move.nl` in Allowed domains,
plus "Also include default list of common package managers"). Without it every
fetch fails and the fields stay blank.

Two traps in that markup, both of which produced silently wrong data on the
first pass:

- **`Balkon` is filed under the `Indeling` section, not `Buitenruimte`.** A
  section-scoped lookup misses every balcony. Match labels across all sections.
- **`Tuin` names the *type* of outdoor space.** `Zonneterras` is a terrace, not
  a garden; `Geen tuin` is neither. Classify the value, don't just test that the
  field is present.

Scraped values are **defaults**. Anything set in the board's Edit panel wins,
including deliberately clearing a field; the panel shows which it is displaying.

## State and sharing

Board state (column, energy, G/T, notes) is shared through **`househunt/board.json`
on the working branch**. Every open board polls it each 10 seconds and merges by
per-listing timestamp, so the most recent edit to a given listing wins — and only
for that listing. An edit also pushes within half a second rather than waiting out
the interval. `localStorage` still holds a local copy, so the board works offline
and catches up when it reconnects.

Reading needs nothing: the repo is public, and a board with no token polls
`raw.githubusercontent.com` read-only. **Saving** needs a GitHub
[fine-grained token](https://github.com/settings/personal-access-tokens) limited to
`eduryev/cvolume` with **Contents: Read and write**, pasted once into the Sync
panel. It is kept in that browser's `localStorage` only.

Concurrent writes use the file's blob SHA for optimistic concurrency: if the other
board committed first, GitHub rejects the write, and the next tick re-reads, merges
and retries. Conflicts resolve per listing, never by clobbering the whole file.

> A token in browser storage is readable by anything else served from the same
> origin (`eduryev.github.io`). Keep it scoped to this one repository, set an
> expiry, and revoke it when the hunt is over.

**Because the board commits to the working branch, anything that pushes to that
branch must `git pull --rebase` first** — see the daily routine.

Geocoding runs in the visitor's browser against PDOK (the Dutch national
address service) and is cached locally; until it resolves, pins sit at their
postcode centroid.
