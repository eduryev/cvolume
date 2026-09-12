# Daily refresh routine

Routine `AVT house hunt — daily listing refresh` (id `trig_01RMr1KEum9Kfe2th7XQ3y4g`)
runs at 06:00 UTC / 08:00 Amsterdam and re-parses any new AVT digests into the board.

**It needs the Gmail connector attached to it.** It was created from a session that
could not pass connector grants through, so as scheduled it has no Gmail tools and
will stop on step 2 without doing anything. Fix it at
claude.ai → Routines → this routine → enable the **Gmail** connector.

The prompt it runs is below, so it can be recreated from the Routines UI if that is
easier than editing it.

---

Refresh the AVT house-hunt board in the eduryev/cvolume repo with any new apartment listings.

Work on the branch `claude/apartment-kanban-dashboard-wir16n`. Read `househunt/README.md` first — it documents the whole pipeline. Do not change the page design or the build scripts unless something is actually broken.

Steps:

1. Check out the branch and read `househunt/README.md`, `househunt/parse_emails.py` and `househunt/build.py`.

2. In Gmail, search `from:info@avtmakelaars.nl Vraagprijs` and find digest emails ("Er zijn nieuwe objecten gevonden") newer than the newest `first_seen` date already present in `househunt/data/listings.json`. If there are none, stop — do not commit, do not message anyone. (If the Gmail connector is not available in this session, say so and stop.)

3. For each new thread call `get_thread` with `messageFormat: "FULL_CONTENT"`. Results too large for context are written to a file automatically; that is the good path — point `parse_emails.py` at those files rather than reading the HTML into context. If a thread comes back inline instead, hand-build an equivalent minimal JSON entry, and note that `parse_emails.py`'s card regex requires the card to end `</table></td></tr></tbody></table>`.

4. Re-run the parser over ALL saved thread dumps plus the existing data, so deduplication by address still works:
   `python3 househunt/parse_emails.py "<dump-dir>/*.txt" househunt/data/listings.json`
   Sanity-check the output: every listing needs url, price, m2, rooms and postcode, and the total must not go down. If the parser reports incomplete records, look at the raw card — AVT varies the wording (`Vraagprijs:` vs `Koopsom:`, and `living m² / plot m²` for houses) and a new variant may need a small regex fix.

5. Any new postcode area (PC4) that is not already in `househunt/pc4.json` needs an entry — neighbourhood name plus approximate lat/lng — or its pins will not place.

6. `python3 househunt/build.py`, then commit and push to the same branch. Commit message: what was added, e.g. "Add 4 listings from the 12-13 Sep AVT digests".

Leave `CUTOFF` in build.py alone — Eduard set it deliberately at 2026-09-10. New listings arrive in Backlog automatically because their first_seen is later than the cutoff.

Note: energy label, garden and terrace are deliberately blank — they are not in the emails and this cloud environment cannot reach move.nl (the egress proxy blocks it). Do not try to scrape them and do not invent them. Eduard and Katia fill those in from the board's Edit panel.

Only notify if you actually added listings, or if the parser broke and you could not fix it. Say how many listings were added and their addresses.
