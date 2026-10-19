# Rare-Disease Benchmark Protocol

## 1. Purpose

Evaluate the phenotype-to-disease ranking engine on independently sourced rare-disease cases with documented phenotype observations and reference diagnoses.

This benchmark is for research and software evaluation only. It does not establish clinical validity or diagnostic readiness.

## 2. Candidate Data Source

- Resource: GA4GH Phenopacket Store
- Official repository: https://github.com/monarch-initiative/phenopacket-store
- Required data: case-level HPO phenotype terms and reference diagnosis
- Dataset release, download URL, retrieval date, license, and SHA-256 checksum must be recorded before evaluation.
- Do not commit downloaded benchmark data unless its license explicitly permits it.
- Candidate release: Phenopacket Store 0.1.27, subject to verification of the release asset and license before download.

## 3. Case Eligibility

A case is eligible only when:

- It has a usable case identifier.
- It contains at least one valid, observed HPO phenotype.
- It has a documented reference diagnosis that can be mapped to the disease identifiers used by the ranking engine.
- Its source and provenance can be recorded.
- Its input features do not expose the reference diagnosis or other information that would reveal the answer to the ranker.

Record exclusions and their reasons. Do not silently discard cases.

## 4. Leakage Assessment

Assess whether the case's source publication contributed to the disease-phenotype annotations used by the ranker.

Report results separately for:

1. Cases with no identified source-publication overlap.
2. Cases with identified overlap.
3. Cases where overlap status is unknown.

Absence of identified overlap is not proof that all leakage has been eliminated. Do not claim a leak-free benchmark without evidence.

## 5. Evaluation Rules

- Freeze and record the ontology and annotation releases.
- Exclude diagnosis labels and publication-derived answer hints from ranking inputs.
- Use deterministic scoring and tie-breaking.
- For tied scores, use disease identifier as a deterministic secondary sort key. This is a software convention, not evidence that tied candidates have different biological likelihoods.
- Document every reference-diagnosis-to-disease-ID mapping.
- Record ambiguous and unmappable diagnoses with explicit exclusion reasons.
- Keep all cases from the same source publication in the same evaluation group.
- Do not tune scoring weights on the final test cohort.

## 6. Metrics

Report:

- Recall@1
- Recall@5
- Recall@10
- Recall@20
- Mean Reciprocal Rank (MRR)
- Median rank
- Number and percentage of unranked cases
- Number of eligible and excluded cases

The metric implementation and tie-handling policy must be documented before evaluation.

## 7. Robustness Analysis

After establishing the unmodified baseline, evaluate controlled changes such as:

- Phenotype dropping
- Replacement of specific terms with broader ancestors
- Addition of unrelated, nonspecific phenotype terms

Record random seeds and report these results separately from the unmodified baseline.

## 8. Reproducibility

Record:

- Dataset release and source URL
- Retrieval date and license
- SHA-256 checksums
- Ontology and annotation versions
- Ranking configuration and software commit
- Case inclusion and exclusion counts
- Leakage assessment method
- Metric definitions and tie-handling policy

## 9. Interpretation

Synthetic benchmark results validate pipeline behavior, not clinical accuracy. Results from published cases may be affected by publication bias, annotation overlap, and differences between reported cases and real-world referrals.

The system remains a research prototype and must not be used as a stand-alone diagnostic tool.
