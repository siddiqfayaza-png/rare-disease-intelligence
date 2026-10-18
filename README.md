# AI-Powered Rare Disease Diagnosis & Variant-to-Mechanism Intelligence

[![Python Version](https://img.shields.io/badge/python-3.11%20%7C%203.12-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Development Status](https://img.shields.io/badge/status-Phase%200%20(Scaffolding)-orange.svg)]()

> [!WARNING]
> **RESEARCH PROTOTYPE ONLY — NOT FOR CLINICAL USE**
> This software is an exploratory research prototype developed exclusively for computational biology, bioinformatics evaluation, and scientific inquiry. It is **not** a certified medical device, has **not** been clinically validated, and **must not** be used as a diagnostic tool, for direct patient care, clinical decision support, or treatment selection.

---

## 1. The Problem: The Diagnostic Odyssey in Rare Diseases

Rare diseases collectively impact approximately 300 to 400 million individuals globally, encompassing between 7,000 and 10,000 distinct identified disorders. Despite their collective frequency, individual conditions exhibit low prevalence, often presenting with:

- **Prolonged Diagnostic Delay:** Patients frequently endure an arduous "diagnostic odyssey" lasting an average of 5 to 7+ years, undergoing numerous specialist consultations and repeated misdiagnoses.
- **Phenotypic Heterogeneity & Overlap:** Distinct genetic diseases often present with overlapping clinical signs, and single diseases can exhibit variable expressivity and incomplete penetrance.
- **Cognitive & Computational Complexity:** Clinicians must cross-reference multi-system phenotypic observations against evolving biomedical literature, ontologies, and rapidly expanding genomic databases.

Bridging structured clinical observations with computational knowledge representations is critical to reducing diagnostic friction in a transparent, verifiable manner.

---

## 2. Project Objectives

The project aims to construct a modular, auditable, and reproducible computational pipeline that:

1. **Normalizes Patient Phenotypes:** Standardizes disparate clinical signs and symptoms into ontological concepts using the [Human Phenotype Ontology (HPO)](https://hpo.jax.org/).
2. **Prioritizes Candidate Diseases:** Ranks candidate rare diseases through semantic similarity and information-theoretic graph algorithms over curated disease-phenotype knowledge bases.
3. **Decouples Variant Ingestion:** Provides a clean, abstracted interface for genomic variant data (HGVS / VCF) to prepare for downstream variant-to-mechanism prioritization without premature pipeline coupling.
4. **Maintains Complete Evidence Provenance:** Emphasizes auditable, rule-grounded associations where every scoring contribution can be traced back to its underlying ontology term or curated assertion.
5. **Enforces Rigorous Evaluation:** Implements leak-free benchmark suites evaluating retrieval performance (Recall@k, MRR) against realistic, noise-perturbed clinical vignettes.

---

## 3. Intended Research Use

This system is engineered for researchers in computational genetics, bioinformatics, and clinical informatics to:

- Benchmark phenotype-matching and semantic similarity algorithms.
- Evaluate the impact of clinical term imprecision, phenotypic noise, and missing observations on disease ranking.
- Study methods for integrating genomic variant signals with phenotypic evidence in future stages.
- Formulate reproducible, version-controlled research workflows for rare disease intelligence.

---

## 4. MVP Scope & Staged Roadmap

To ensure software quality and architectural integrity, development proceeds incrementally:

- **Phase 0: Project Scaffolding & Architectural Foundation (Current)**
  Establish workspace structure, package layout (`src/rare_disease_ai`), quality tooling (Ruff, pytest), strict `.gitignore`, dependency specifications, and core architectural documentation. Zero external dependencies, datasets, or models.
- **Phase 1: Phenotype Normalization & Ranking Baseline**
  Implement HPO graph ingestion, term normalization, Information Content (IC) calculation, and baseline semantic similarity ranking (e.g., Resnik, Lin, Phenomizer-style scoring).
- **Phase 2: Biomedical Evidence & Knowledge Integration**
  Incorporate curated gene-disease-phenotype association tables with explicit license verification, provenance tracking, and data dictionary enforcement.
- **Phase 3: Variant Data Interface**
  Implement schema validation and coordinate normalization for genomic variants (GRCh37/GRCh38), decoupled from phenotype ranking.
- **Phase 4: Biological Mechanism Investigation**
  Integrate pathway enrichment and functional annotation to contextualize prioritized variants and candidate disease mechanisms.
- **Phase 5: Reporting, Auditing, and Benchmark Validation**
  Deliver human-readable diagnostic summaries with transparent score breakdowns and comprehensive benchmark reporting.

---

## 5. Planned Modules

The codebase is organized into discrete, single-responsibility modules within the `rare_disease_ai` package:

| Module | Scope & Responsibility |
| :--- | :--- |
| `case_validation` | Validates patient input schemas (phenotype lists, demographic context, genomic coordinates). |
| `phenotype_normalization` | Maps raw terms and text descriptions to standardized HPO identifiers (`HP:XXXXXXX`). |
| `phenotype_ranking` | Scores candidate diseases using ontology graph traversal and semantic similarity metrics. |
| `variant_interface` | Ingests, normalizes, and validates genomic variants (HGVS, VCF) across reference builds. |
| `biomedical_evidence` | Interfaces with verified, curated gene-disease and variant-disease knowledge sources. |
| `mechanism_investigation`| Connects genes and variants to biological pathways (e.g., Reactome, GO) for functional insights. |
| `reporting` | Generates transparent, auditable diagnostic evidence summaries with explicit provenance. |
| `evaluation` | Computes ranking metrics (Recall@1, Recall@5, MRR) with strict data leakage prevention. |

---

## 6. Current Development Status

- **Stage:** **Phase 0 (Project Scaffolding)**
- **Active Code:** Package scaffolding (`src/rare_disease_ai/__init__.py`), packaging configuration (`pyproject.toml`), test harness (`tests/test_project_structure.py`), and documentation suite (`docs/`).
- **No external datasets have been downloaded.**
- **No external APIs or web services are called.**
- **No ML models, LLMs, or agent frameworks are integrated.**

---

## 7. Local Setup Instructions

### Prerequisites
- Python 3.11 or Python 3.12
- `pip` package manager

### 1. Clone or Open the Repository
```bash
cd rare-disease-intelligence
```

### 2. Create and Activate a Virtual Environment

**On Linux / macOS:**
```bash
python3 -m venv .venv
source .venv/bin/activate
```

**On Windows (PowerShell):**
```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
```

### 3. Install Development Dependencies
```bash
pip install -e ".[dev]"
```

### 4. Run the Test Suite
Using `pytest` (once installed in virtual environment):
```bash
pytest
```

Using Python's built-in `unittest` runner (runs with zero third-party dependencies):
```bash
# On Linux / macOS
PYTHONPATH=src python3 -m unittest discover -s tests -p "test_*.py" -v

# On Windows (PowerShell)
$env:PYTHONPATH="src"; python -m unittest discover -s tests -p "test_*.py" -v
```

### 5. Run Linting and Code Quality Checks
```bash
ruff check .
```

---

## 8. Repository Layout

```text
rare-disease-intelligence/
├── .env.example              # Template environment configuration (no secrets)
├── .gitignore                # Comprehensive ignores for caches, data, environments
├── LICENSE                   # Open-source license definition
├── README.md                 # Project overview and research scope (this file)
├── pyproject.toml            # Build metadata, minimal deps, pytest and ruff configs
├── docs/
│   ├── architecture.md       # Modular pipeline design and technical specifications
│   ├── data-sources.md       # Proposed data sources, licensing, and integration criteria
│   ├── evaluation-plan.md    # Metrics (Recall@k, MRR), leakage prevention, test cases
│   └── project-scope.md      # Mission, boundaries, personas, and phased roadmap
├── src/
│   └── rare_disease_ai/      # Core Python package root
│       └── __init__.py       # Package metadata and version definition
└── tests/
    ├── __init__.py
    └── test_project_structure.py  # Initial smoke tests for package scaffolding
```

---

## 9. Ethics, Governance, and Regulatory Stance

This software is strictly a computational research vehicle. It is not designed to replace clinical expertise, diagnostic assays, or genetic counseling. Any diagnostic hypothesis produced by future iterations of this software must be validated by qualified medical professionals using accredited clinical diagnostic procedures.
