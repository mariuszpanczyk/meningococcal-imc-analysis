# Final release checklist

## Scientific
- [ ] Reconcile admission-to-antibiotic p/q discrepancy between HTML and Table 2.
- [ ] Approve the reconstructed permutation statistic, new seed, bootstrap stream and supplementary WY implementation.
- [ ] Compare all final manuscript values with a private run; retain an internal reconciliation log.
- [ ] Verify source measurement units, cohort eligibility and any source corrections privately.
- [ ] Confirm the stated laboratory window/event timing from clinical records; do not claim this notebook verified timestamps absent from input.
- [ ] Retain the documented population distinction: 53 in the full descriptive cohort; 51 in IMC analyses after excluding two deaths within 24 hours of admission because IMC could not be assessed. Do not code early deaths as no IMC.
- [ ] If full-cohort Table 1/Figure 1 reproduction is required, supply a separate authorized 53-person descriptive input and validate the optional branch. This does not indicate missing records in the 51-person analytical database.
- [ ] Approve the native Python Firth implementation; independent logistf comparison is recommended.

## Privacy
- [ ] Upload only the explicit public files in this archive.
- [ ] Do not upload the input CSV, original HTML, manuscript attachments, private results or logs.
- [ ] No participant rows, IDs, dates, private paths, credentials or private hashes embedded in source/Markdown.
- [ ] No saved outputs, execution counts, attachments or hidden cell metadata.
- [ ] Restore all public defaults; input paths None; upload/export disabled.
- [ ] Inspect Git history and all tracked files; .gitignore does not remove previously committed data.
- [ ] Review any separately published aggregate outputs for small-cell disclosure risk.

## Publication
- [ ] Complete verified software title and full approved author list/ORCIDs in CITATION.cff.
- [ ] Choose/add approved LICENSE and consistent Zenodo license metadata.
- [ ] Complete repository URL, manuscript identifier and version-specific DOI only when real.
- [ ] Create reviewed public GitHub repository, enable Zenodo integration, create a release and verify archived files/metadata.
- [ ] State actual data-access conditions separately; do not promise access without authorization.
