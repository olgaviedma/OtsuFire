# MANIFEST — `inst/scripts/long_run/` canonical runners

These are the **canonical, versioned** long-run runner files. The same content
hashes are in `MANIFEST.csv`. Use the `sha256` (or the LF-normalized content
md5, which survives CRLF normalization) to detect divergence between a working
copy under a run's `_ORCHESTRATION/` and this canonical source.

| File | bytes | LF-content md5 | sha256 (prefix) |
|------|------:|----------------|-----------------|
| `00_common.R` | 21000 | `7982570cd615e0e290733f6d8b43e16c` | `eb70d1bc…` |
| `05_env_snapshot.R` | 3181 | `8704f9b81608cc13ae5d822b6444af6d` | `a2c2df00…` |
| `10_shared_upstream.R` | 5628 | `866cbc1aa0cc15fbb76d3a4c2172da10` | `27ef8c49…` |
| `20_full.R` | 10044 | `ab4aa604cd2bc31cc160b1bb443fba79` | `8524f447…` |
| `30_no_hotspot.R` | 9573 | `2a4e80e6ff789e5fe8ea4aeba954c949` | `f0c4f4d5…` |
| `98_dry_verify_routes.R` | 5435 | `600ca77620e1d935325a316a8ff28e4d` | `c997d53c…` |
| `MASTER_RUN.bat` | 4631 | `43ad31d9eec83d2f820e6e753ca6d5c7` | `8920e5b7…` |
| `RESUME_RUN.bat` | 5151 | `a8f3191feacd5ac205732cdd4fa79180` | `230bb1a7…` |
| `README.md` | 5912 | `42bbef3553c3229fab91d82b45ebb5ff` | `886e78fe…` |

The `reports/` subfolder holds the versioned audit reports (small `.md`/`.csv`
only). Heavy run artifacts (models, feature GPKG, RDS handoffs) are **NOT**
copied into the repo; they are referenced by absolute path + SHA256 in
`reports/RUN_MANIFEST.md`.
