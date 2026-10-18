# Data Sources: Candidate Knowledge Bases, Licensing, and Integration Governance

> [!IMPORTANT]
> **GOVERNANCE NOTICE: ZERO EXTERNAL DATA INGESTION IN PHASE 0**
> No external datasets, ontologies, annotations, or genomic files have been downloaded, bundled, or cached within this repository during Phase 0. The sources cataloged below represent vetted candidates for subsequent development phases. Prior to integrating any external data source, explicit verification of current license terms, access mechanisms, release versions, and schema invariants is strictly required.

---

## 1. Candidate Biomedical Data Sources

The table below summarizes the external knowledge bases proposed to power the phenotype prioritization, evidence synthesis, and variant interpretation stages of the pipeline:

| Resource | Primary Domain | Upstream Provider | Target Pipeline Stage | Proposed Access Method |
| :--- | :--- | :--- | :--- | :--- |
| **Human Phenotype Ontology (HPO)** | Phenotype terminology, hierarchy, and disease associations | The Jackson Laboratory / HPO Consortium | Phase 1 (Phenotype Ranking) | Bulk release download (`hp.obo`, `phenotype.hpoa`) |
| **Monarch Initiative** | Integrated cross-species disease-gene-phenotype knowledge graph | Monarch Initiative / Mondo Consortium | Phase 2 (Biomedical Evidence) | Bulk release graph dumps (TSV/Parquet) or Mondo OBO |
| **ClinVar** | Curated variant clinical significance & pathogenicity assertions | NCBI / National Library of Medicine | Phase 3 (Variant Interface) | Compressed bulk VCF (`clinvar.vcf.gz`) & summary TSVs |
| **gnomAD** | Population allele frequencies & gene intolerance/constraint | Broad Institute | Phase 3 (Variant Interface) | Cloud bucket downloads (AWS/GCP Open Data) |
| **Orphanet (Orphadata)** | Rare disease clinical epidemiology and classifications | INSERM / Orphanet Consortium | Phase 2 (Biomedical Evidence) | Structured XML/JSON data packs |

---

## 2. In-Depth Source Profiles & Pre-Integration Requirements

### 2.1 Human Phenotype Ontology (HPO)

- **Intended Use:**
  - Serves as the canonical clinical phenotype taxonomy.
  - Hierarchical Directed Acyclic Graph (DAG) defines ancestor-descendant relationships used for calculating semantic similarity, finding the Most Informative Common Ancestor (MICA), and computing Information Content (IC).
  - Disease-phenotype annotation table (`phenotype.hpoa`) establishes associations between HPO terms and rare disease entities (OMIM, Orphanet, DECIPHER).
- **Proposed Artifacts:**
  - `hp.obo` or `hp.json` (Ontology definition graph)
  - `phenotype.hpoa` (Curated annotations with evidence codes, frequencies, and onset modifiers)
- **Pre-Integration Verification Checklist:**
  - [ ] **License Verification:** Confirm current licensing terms (historically CC-BY 4.0; verify attribution requirements and non-commercial/commercial redistribution clauses).
  - [ ] **Access Method:** Programmatic download script with pinned release URL and automated SHA-256 checksum verification.
  - [ ] **Version Pinning:** Pin exact release date (e.g., `release-YYYY-MM-DD`) in configuration; disallow unpinned latest redirects in production/eval.
  - [ ] **Known Limitations:**
    - *Ascertainment bias:* Frequently studied diseases have more granular phenotype annotations than ultra-rare, newly described syndromes.
    - *Variable penetrance:* Not all patients with a disease display all annotated phenotypes.
    - *Annotation sparsity:* ultra-rare conditions may possess only 2-3 annotated terms.

---

### 2.2 Monarch Initiative (Mondo & Knowledge Graph)

- **Intended Use:**
  - Standardizes disease entity identifiers using the Mondo Disease Ontology, reconciling conflicting definitions between OMIM, Orphanet, and ICD-10/11.
  - Ingests cross-species gene-to-phenotype and disease-to-gene associations for evidence consolidation.
- **Proposed Artifacts:**
  - `mondo.obo` or `mondo.json`
  - Curated edge tables (e.g., disease-to-gene and gene-to-phenotype relations)
- **Pre-Integration Verification Checklist:**
  - [ ] **License Verification:** Verify CC-BY 4.0 license status across aggregated sub-components, especially where third-party proprietary ontologies might be mapped.
  - [ ] **Access Method:** Evaluate bulk TSV downloads versus Monarch GraphQL/REST API; prioritize local pinned flat files for determinism.
  - [ ] **Version Pinning:** Record exact schema version and KG release identifier.
  - [ ] **Known Limitations:**
    - Heterogeneity in edge confidence scores across distinct source databases.
    - Model organism phenotypes (e.g., mouse MGI) must be clearly distinguished from human clinical presentations.

---

### 2.3 ClinVar

- **Intended Use:**
  - Provides clinical assertions regarding the pathogenicity of genomic variants (Pathogenic, Likely Pathogenic, Variant of Uncertain Significance, Likely Benign, Benign).
  - Supplies review status (0 to 4 stars) and submitter-level evidence.
- **Proposed Artifacts:**
  - `clinvar.vcf.gz` (Monthly release aligned to GRCh38 and GRCh37)
  - `variant_summary.txt.gz`
- **Pre-Integration Verification Checklist:**
  - [ ] **License Verification:** Confirm public domain / US Government Work status; verify citation and attribution compliance.
  - [ ] **Access Method:** Automated retrieval from NCBI FTP servers with MD5 checksum verification.
  - [ ] **Version Pinning:** Lock specific monthly archive releases (e.g., `clinvar_YYYYMMDD.vcf.gz`).
  - [ ] **Known Limitations:**
    - *Assertion conflicts:* Different clinical laboratories occasionally submit conflicting classifications for the same variant (e.g., VUS vs. Likely Pathogenic).
    - *Historical classifications:* Older records may lack contemporary ACMG/AMP criteria re-evaluations.
    - *Review status variance:* Low-star (0 or 1 star) records require downstream confidence discounting.

---

### 2.4 Genome Aggregation Database (gnomAD)

- **Intended Use:**
  - Assesses population allele frequency (AF) to filter out common polymorphisms that are incompatible with rare, high-penetrance Mendelian disease models.
  - Supplies gene-level constraint metrics (e.g., pLI for loss-of-function intolerance, LOEUF scores).
- **Proposed Artifacts:**
  - Population frequency tables and site-level VCFs (or regional summary queries).
  - Gene constraint summary tables (`gnomad.v4.0.constraint_metrics.tsv`).
- **Pre-Integration Verification Checklist:**
  - [ ] **License Verification:** Confirm open-data license terms (ODC-By for gnomAD v4; verify terms regarding redistribution of derived files).
  - [ ] **Access Method:** S3/GCP Public Datasets bucket access; extract only required target fields to avoid massive multi-terabyte overhead.
  - [ ] **Version Pinning:** Choose and document target version (e.g., v4.1 for GRCh38 or v2.1.1 for GRCh37).
  - [ ] **Known Limitations:**
    - Variable sequencing coverage in GC-rich and repetitive genomic regions.
    - Unequal population representation across continental ancestries.

---

### 2.5 Supplementary Sources: Orphanet & OMIM

- **Intended Use:**
  - Clinical epidemiology (prevalence estimates, age of onset distributions, inheritance modes).
- **Pre-Integration Verification Checklist:**
  - [ ] **License Verification:**
    - *OMIM:* Requires explicit non-commercial academic licensing or commercial agreement; **direct redistribution of OMIM databases is restricted**. Integration must respect these legal terms.
    - *Orphanet / Orphadata:* Freely accessible for non-commercial research under specific data use agreements.
  - [ ] **Governance Rule:** Never bundle or mirror restricted third-party data directly into the public source code repository.

---

## 3. Data Integration Governance Protocol

Prior to the introduction of any external data into the codebase during Phase 1, Phase 2, or Phase 3, the following steps must be completed and documented:

1. **Licensing & Legal Audit:** Validate that license terms permit the intended computational usage, derivative work generation, and research dissemination.
2. **Deterministic Version Pinning:** Every external file must have a documented release version, source URI, and cryptographic hash (SHA-256).
3. **Data Dictionary Specification:** Create a formal schema definition document detailing expected fields, nullability, types, and primary keys.
4. **Local Caching Discipline:** External datasets must reside exclusively in ignored local directories (`data/external/` or `data/cache/`), never committed to git history.
5. **Ingestion Pipelines as Pure Code:** All data ingestion must be executed through reproducible, automated Python scripts capable of rebuilding local caches from verified upstream sources.
