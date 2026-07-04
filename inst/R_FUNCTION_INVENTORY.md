# `R/` function inventory (index)

Maintainer index for symbols defined under **`R/`**. Split into two lists:

| Document | Contents |
|----------|----------|
| **[R_EXPORTED_AND_DOCUMENTED.md](R_EXPORTED_AND_DOCUMENTED.md)** | `NAMESPACE` exports grouped by **glmbayes** / **lmebayes** overlap (including **glmbayes** retain vs phase-out split), S3 methods split by **glmbayes** registration, internal `man/` topics, doc topics. |
| **[R_CORE_ONLY_EXPORTS.md](R_CORE_ONLY_EXPORTS.md)** | Same exports organized by function type (simulation, envelopes, two-block, block-Gibbs, …) with **lmebayes** direct/indirect notes. |
| **[R_EXPORT_REACHABILITY.md](R_EXPORT_REACHABILITY.md)** | Dead / inactive export analysis (reachability from **glmbayes** retain re-exports and **lmebayes** drivers). |
| **[R_INTERNAL_HELPERS.md](R_INTERNAL_HELPERS.md)** | `@noRd` / `@keywords internal` helpers, C++ glue, attach hooks. |

Scratch checks and one-off scripts live in `data-raw/` (not run by `test_check`).

