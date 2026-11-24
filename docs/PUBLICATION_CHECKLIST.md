# Decisions before public release

The code/data package is assembled for review at
<https://github.com/NickelRaney/OS_network_release>. The author requested upload
to this repository on 2026-11-24. No software release version, archival DOI or
license has been selected on behalf of the authors.
The author requested on 2026-11-23 that authorship, paper identifiers and
licensing remain unset for now. These items are deferred, not a request for
approval to continue the requested code preparation and repository upload.

- Scope now includes both original MATLAB models, following the author's
  subsequent request to package Low-LGN as well. High-cx F1 EE/EI scales are
  2.25/1.5 (author-confirmed); Low-LGN uses 1.3/1.2 from its original parameter
  notes. BG/grating LGN gains are respectively 1.15/1.20 and 1.11/1.15;
  the historical onset helper starts at 1.11/1.15 for both.
- Provide author names, order and optional ORCIDs, the paper title, journal or
  preprint information and DOI when assigned. Complete `CITATION.cff.in`, then
  rename it to `CITATION.cff`; do not upload the template as finalized metadata.
- Confirm ownership and permission to distribute the selected MATLAB sources,
  both fixed connectivity/parameter snapshots and derived historical histograms.
  Specify a code license and, if different, a data license. Add their actual
  texts and preserve any required upstream attribution.
- Identify the foundational model/MFE paper and other source references for
  numerical parameters. Source comments mention an MFE paper but do not give a
  reliable bibliographic identifier. No guessed citation has been inserted.
- The repository name/URL is now supplied. Select a release version when
  ready. After metadata review, a DOI archive such as an institutional/software
  repository can be linked if the authors choose one; no archival identifier
  has been assigned to this candidate.
- Decide whether the submission promises reproduction of additional figures.
  This package supplies fresh simulation examples and an exact frozen complex
  onset-analysis check. Original many-trial spike datasets, simple-cell onset
  analysis and figure builders are not represented as reproduced by this demo.

The unmodified experiment archive retains all source data. Additional datasets
can be selected deliberately after the paper's reproduction scope is fixed;
they should not be copied wholesale from the historical working directory.
