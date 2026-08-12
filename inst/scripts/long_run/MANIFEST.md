# MANIFEST — `inst/scripts/long_run/` canonical runners

These are the **canonical, versioned** long-run runner files. The same content
hashes are in `MANIFEST.csv`. Use the `sha256` (or the LF-normalized content
md5, which survives CRLF normalization) to detect divergence between a working
copy under a run's `_ORCHESTRATION/` and this canonical source.

| File | bytes | LF-content md5 | sha256 (prefix) |
|------|------:|----------------|-----------------|
| `00_common.R` | 21137 | `fae2fd9ca7a3093c207594b534bf7764` | `49f5af02…` |
| `05_env_snapshot.R` | 3181 | `8704f9b81608cc13ae5d822b6444af6d` | `a2c2df00…` |
| `10_shared_upstream.R` | 5644 | `291a5ad7f5f0e5fa0d206a2037be30ec` | `62511e4c…` |
| `20_full.R` | 8167 | `015cf509cc84bbe72505e91f359c254d` | `6528afc2…` |
| `30_no_hotspot.R` | 8859 | `0e236f00b9017a7d965d34c1b7a2663d` | `e5169219…` |
| `98_dry_verify_routes.R` | 5412 | `0043ec67c8a3b52d77b6130491669885` | `fd3bbd35…` |
| `MASTER_RUN.bat` | 4631 | `43ad31d9eec83d2f820e6e753ca6d5c7` | `8920e5b7…` |
| `RESUME_RUN.bat` | 5151 | `a8f3191feacd5ac205732cdd4fa79180` | `230bb1a7…` |
| `README.md` | 5966 | `306311f878f741259bf29238873a2adb` | `7b5b92f8…` |

The `reports/` subfolder holds the versioned audit reports (small `.md`/`.csv`
only). Heavy run artifacts (models, feature GPKG, RDS handoffs) are **NOT**
copied into the repo; they are referenced by absolute path + SHA256 in
`reports/RUN_MANIFEST.md`.
