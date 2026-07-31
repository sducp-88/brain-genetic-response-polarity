# Data access, version, redistribution, and attribution audit

Audit date: 31 July 2026  
Scope: frozen low-target manuscript and aggregate-only reproducibility bundle

## Decision summary

The completed analyses do **not** depend on downloading controlled Synapse
participant-level files. The primary PsychAD expression inputs were the public
processed H5AD assets distributed by CZ CELLxGENE Discover. The PsychAD
Synapse DOI remains the persistent identifier for the originating dataset and
the AD Knowledge Portal acknowledgement still applies.

No third-party H5AD, XLSX, or SEA-AD object will be redistributed in the code
bundle. Only retrieval metadata, source versions, scripts, and
non-identifying aggregate model outputs are included.

## PsychAD / CELLxGENE

- Collection: `84ce6837-548d-4a1f-919f-0bc0d9a3952f`
- Collection version inspected through the official API:
  `52ba1325-53e8-4477-9ce0-f50836019a1f`
- RADC dataset ID: `5e57cd50-8e42-42d6-940d-5c1660d06864`
- RADC dataset version: `54293783-669c-410e-919d-474960f8761b`
- MSSM dataset ID: `37a17b78-4864-4a42-b67b-31c00962795a`
- MSSM dataset version: `0e853475-e298-4b09-881a-ed0b60d5a8c9`
- Source PsychAD study: `syn52160016`
- Source PsychAD data DOI: `10.7303/syn60084804`
- Descriptor: Fullard et al., *Scientific Data* 2025,
  `10.1038/s41597-025-04687-5`

Use the dataset-specific CELLxGENE citation strings, the Scientific Data
descriptor, and the PsychAD dataset DOI. The official API identifies the
assets as public but does not expose a collection-specific licence field in
the returned metadata; therefore this bundle conservatively avoids
redistributing the H5AD files.

The AD Knowledge Portal requires the general acknowledgement:

> The results published here are in whole or in part based on data obtained
> from the AD Knowledge Portal.

The descriptor's acknowledgements also identify the relevant tissue sources
and support: NIH NeuroBioBank at Mount Sinai (NIMH-75N95019C00049), Rush
Alzheimer's Disease Center (P30AG10161, P30AG72975, R01AG15819, R01AG17917,
R01AG22018, U01AG46152, and U01AG61356), and the NIMH-IRP Human Brain
Collection Core (ZIC MH002903).

## SingleBrain genetic-expression anchors

- Original article DOI: `10.1038/s41588-026-02541-x`
- Publisher correction DOI: `10.1038/s41588-026-02581-3`
- Local supplementary workbook SHA-256:
  `081bc422bd39d50f0a8dcad4f0c69ad8441f507318e20c46ef7cc19d8a834468`

The correction changes the Fig. 1a cohort label from
`ROS/MAP-Columbia` to `ROS/MAP-MIT` and fixes one citation range. The current
publisher version is therefore the reference version. The publisher workbook
is not included in the reproducibility bundle.

## SEA-AD

- Public bucket: `s3://sea-ad-single-cell-profiling`
- Registry: `https://registry.opendata.aws/allen-sea-ad-atlas/`
- Objects used: 12 multiregion donor-pseudobulk H5AD files whose names contain
  release date `2026-06-22`
- Local downloads completed between 30 and 31 July 2026
- Exact URLs, sizes, and official/local multipart ETags:
  `SEA_AD_OBJECT_MANIFEST.csv`

Required registry citation:

> Seattle Alzheimer's Disease Brain Cell Atlas (SEA-AD) was accessed on
> 31 July 2026 from https://registry.opendata.aws/allen-sea-ad-atlas/.

The registry identifies support from NIA grant U19AG060909. Allen Institute
terms permit research and other noncommercial use with proper attribution and
allow a limited set of content in scholarly publication. To avoid any
redistribution ambiguity, the public bundle contains no Allen Institute H5AD
content.

## Public-release boundary

Safe to release:

- author-owned scripts without credentials or absolute local paths;
- frozen anchor-, locus-, and group-level statistical tables;
- aggregate figures and figure QA;
- source identifiers, download instructions, checksums, and ETags.

Do not release:

- CELLxGENE, SingleBrain, or SEA-AD source files;
- donor-level pseudobulk matrices or metadata;
- any table containing donor, sample, participant, individual, subject,
  person, or specimen identifiers;
- tokens, account information, or machine-specific monitoring/download
  wrappers.

## Submission items not solved by this audit

- **Completed:** MIT License applied to author-owned code and documentation;
  third-party data remain excluded and are not relicensed.
- **Completed for the release candidate:** ordered authors, affiliations,
  contribution roles, funding, correspondence, and available ORCID
  identifiers have been added.
- Publish the prepared public repository, freeze the release commit, and mint
  an archival DOI.
- Confirm the corresponding author's local institutional determination for
  secondary analysis of de-identified public data.
- Confirm competing interests and obtain collective author approval of the
  contribution statement.
