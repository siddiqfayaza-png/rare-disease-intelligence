import json

import pytest

from rare_disease_ai.mechanism_evidence import MechanismEvidence
from rare_disease_ai.mechanism_investigation import (
    build_mechanism_index,
    investigate_gene_mechanisms,
    mechanism_report_to_dict,
)


def _annotation(
    *,
    gene_symbol: str = "BRCA1",
    annotation_id: str = "GO:0006281",
    qualifier: str | None = None,
) -> MechanismEvidence:
    return MechanismEvidence(
        gene_symbol=gene_symbol,
        annotation_id=annotation_id,
        annotation_label="DNA repair",
        annotation_type="biological_process",
        source_name="Gene Ontology Annotation (UniProtKB)",
        source_release="test-release",
        source_url="https://example.org/go-annotations",
        evidence_code="IDA",
        reference_ids=("PMID:12345678",),
        qualifier=qualifier,
        taxon_id="taxon:9606",
    )


def test_builds_case_insensitive_gene_index() -> None:
    annotations = (
        _annotation(gene_symbol="BRCA1"),
        _annotation(gene_symbol="brca1", annotation_id="GO:0003677"),
    )

    index = build_mechanism_index(annotations)

    assert len(index) == 1
    assert len(index["brca1"]) == 2


def test_retrieves_annotations_and_preserves_provenance() -> None:
    annotation = _annotation(qualifier="NOT|involved_in")
    index = build_mechanism_index((annotation,))

    report = investigate_gene_mechanisms(("BRCA1",), index)

    assert report.requested_gene_count == 1
    assert report.annotations_returned == 1
    assert report.unmatched_gene_symbols == ()
    assert report.results[0].gene_symbol == "BRCA1"
    assert report.results[0].annotations[0] == annotation
    assert report.results[0].annotations[0].qualifier == "NOT|involved_in"
    assert report.results[0].annotations[0].reference_ids == ("PMID:12345678",)


def test_lookup_is_case_insensitive() -> None:
    index = build_mechanism_index((_annotation(gene_symbol="BRCA1"),))

    report = investigate_gene_mechanisms(("brca1",), index)

    assert len(report.results) == 1
    assert report.results[0].gene_symbol == "brca1"
    assert report.annotations_returned == 1


def test_reports_genes_without_annotations() -> None:
    index = build_mechanism_index((_annotation(),))

    report = investigate_gene_mechanisms(("BRCA1", "UNKNOWN1"), index)

    assert report.requested_gene_count == 2
    assert len(report.results) == 1
    assert report.unmatched_gene_symbols == ("UNKNOWN1",)
    assert report.annotations_returned == 1


def test_collapses_duplicate_requested_genes() -> None:
    index = build_mechanism_index((_annotation(),))

    report = investigate_gene_mechanisms(
        ("BRCA1", "brca1", "BRCA1"),
        index,
    )

    assert report.requested_gene_count == 1
    assert len(report.results) == 1


def test_preserves_distinct_annotations() -> None:
    annotations = (
        _annotation(annotation_id="GO:0006281"),
        _annotation(annotation_id="GO:0003677"),
    )
    index = build_mechanism_index(annotations)

    report = investigate_gene_mechanisms(("BRCA1",), index)

    assert report.annotations_returned == 2
    assert {item.annotation_id for item in report.results[0].annotations} == {
        "GO:0006281",
        "GO:0003677",
    }


def test_serializes_report_to_json() -> None:
    index = build_mechanism_index((_annotation(),))
    report = investigate_gene_mechanisms(("BRCA1", "UNKNOWN1"), index)

    payload = mechanism_report_to_dict(report)
    encoded = json.dumps(payload)
    decoded = json.loads(encoded)

    assert decoded["report_type"] == "gene_mechanism_investigation"
    assert decoded["research_only"] is True
    assert decoded["annotations_returned"] == 1
    assert decoded["unmatched_gene_symbols"] == ["UNKNOWN1"]
    assert decoded["results"][0]["annotations"][0]["source_release"] == ("test-release")
    assert decoded["interpretation_limits"]


def test_rejects_invalid_annotation_records() -> None:
    with pytest.raises(TypeError, match="MechanismEvidence"):
        build_mechanism_index(("not an annotation",))  # type: ignore[arg-type]


@pytest.mark.parametrize("genes", [("",), ("  ",), (None,)])
def test_rejects_invalid_gene_symbols(genes: tuple[object, ...]) -> None:
    index = build_mechanism_index((_annotation(),))

    with pytest.raises(ValueError, match="non-empty string"):
        investigate_gene_mechanisms(genes, index)  # type: ignore[arg-type]
