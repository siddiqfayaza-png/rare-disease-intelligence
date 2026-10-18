# Project Scope: AI-Powered Rare Disease Diagnosis & Variant-to-Mechanism Intelligence

## 1. Executive Summary & Mission Statement

The **AI-Powered Rare Disease Diagnosis & Variant-to-Mechanism Intelligence** project is a research-grade computational biology initiative designed to accelerate the identification of causative rare genetic diseases from patient phenotypic features and genomic variants.

The mission of this project is to develop an open, modular, verifiable, and mathematically grounded framework that:
- Normalizes complex multi-system clinical presentations into standardized ontological representations.
- Employs semantic similarity and graph algorithms to rank candidate rare disorders with transparent scoring attribution.
- Establishes decoupled interfaces for genomic variant interpretation and biological mechanism exploration in future phases.
- Maintains rigorous scientific standards through leak-free evaluation methodologies and deterministic reproducibility.

---

## 2. Scope Boundaries

To prevent scope creep and maintain research rigor, explicit boundaries are enforced across all phases of development.

### In-Scope (Phase 0 and Planned Roadmap)
- **Phase 0 (Current):**
  - Project scaffolding, packaging configuration, developer environment setup.
  - Architectural blueprints and module interface contracts.
  - Data governance, licensing evaluation, and data dictionary specifications.
  - Metric definitions and leak-free evaluation protocol design.
  - Baseline project-structure test harness.
- **Phase 1 (MVP Baseline):**
  - Parsing and graph modeling of the Human Phenotype Ontology (HPO).
  - Patient phenotype term normalization and validation.
  - Information Content (IC) calculation across disease-phenotype annotation graphs.
  - Deterministic semantic similarity algorithms (e.g., Resnik, Lin, Phenomizer-style scoring).
  - Baseline evaluation on curated synthetic and literature-derived benchmark cases.
- **Subsequent Phases (Incremental Expansion):**
  - Ingestion of verified, curated biomedical evidence (Monarch, ClinVar, gnomAD).
  - Abstracted genomic variant ingestion (VCF / HGVS) and coordinate normalization.
  - Functional pathway mapping and mechanism investigation (Reactome, Gene Ontology).
  - Comprehensive reporting with provenance tracking.

### Out-of-Scope (Strictly Excluded)
- **Direct Clinical Diagnosis:** The project is not a diagnostic device or a substitute for professional clinical judgement, genetic counseling, or certified medical assays.
- **Raw Sequencing Primary Alignment:** Secondary alignment (FASTQ to BAM/CRAM) and variant calling are external upstream processes; this project begins with structured phenotype lists and normalized variant calls.
- **Opaque Black-Box Predictions:** Unexplainable end-to-end black-box models lacking term-level attribution or biological provenance are excluded.
- **Premature Web / UI Deployment:** Production frontend frameworks, cloud deployment pipelines, or mobile interfaces are out of scope during core algorithmic development.
- **Direct Electronic Health Record (EHR) Ingestion:** Direct integration with live hospital EHR systems (e.g., HL7 FHIR live endpoints) is deferred until core algorithmic validation is complete.

---

## 3. Target Personas and Use Cases

### Persona 1: Computational Geneticist / Bioinformatician
- **Role:** Designs, benchmarks, and optimizes variant prioritization and phenotype-driven disease ranking pipelines.
- **Need:** Modular Python APIs, reproducible graph algorithms, exact version pinning of ontologies and annotations, and reproducible evaluation benchmarks with Recall@k metrics.

### Persona 2: Clinical Genetics Researcher
- **Role:** Investigates undiagnosed or atypical rare disease cases through retrospective research cohorts.
- **Need:** Transparent ranking explanations showing which phenotypic terms contributed to candidate disease scores, backed by verified citations and ontology cross-references.

---

## 4. Architectural Principles

1. **Strict Modularity & Decoupling:** Each pipeline stage operates through explicit data contracts without hidden side effects. Phenotype ranking does not depend on variant interpretation, allowing independent algorithmic iteration.
2. **Deterministic Reproducibility:** Every algorithm must produce identical outputs given identical inputs and ontology versions. No hidden stochastic seeds or nondeterministic tie-breaking.
3. **Traceable Provenance:** Every candidate disease score must be decomposable into constituent phenotype matches, information content weights, and source database releases.
4. **Minimal Dependency Footprint:** Avoid bloated frameworks, premature ORMs, heavyweight agent libraries, or external cloud dependencies. Prefer clean, standard Python with targeted numerical libraries when necessary.

---

## 5. Phased Roadmap

```text
Phase 0 ──► Phase 1 ──► Phase 2 ──► Phase 3 ──► Phase 4 ──► Phase 5
Scaffold    Phenotype   Biomedical  Variant     Mechanism   Reporting &
& Specs     Ranking     Evidence    Interface   Analysis    Benchmarking
(CURRENT)   (MVP)       (Monarch)   (VCF/HGVS)  (Pathways)  (Validation)
```

- **Phase 0:** Scaffolding, repository governance, testing setup, architectural documentation (Current).
- **Phase 1:** Phenotype normalization, ontology graph ingestion, deterministic semantic similarity ranking.
- **Phase 2:** Curated knowledge integration, disease-gene association tables, license verification.
- **Phase 3:** Genomic variant data interface, coordinate validation, pathogenicity evidence linking.
- **Phase 4:** Biological pathway enrichment, variant-to-mechanism contextualization.
- **Phase 5:** Audit reporting, visualization of evidence graphs, comprehensive multi-cohort evaluation.

---

## 6. Regulatory and Safety Stance

This software is developed strictly as a research tool under the MIT License. It does not possess CE-IVD, FDA 510(k), or CLIA/CAP accreditation. Any clinical adoption or translation requires formal clinical validation, verification of diagnostic accuracy under certified regulatory frameworks, and institutional ethics approval.
