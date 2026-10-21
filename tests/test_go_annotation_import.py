from pathlib import Path

from rare_disease_ai.go_annotation_import import read_gaf_annotations

SOURCE_URL = "https://current.geneontology.org/annotations/gaf/HUMAN-uniprot.gaf.gz"


def _gaf_row(
    *,
    symbol: str = "GENE1",
    qualifier: str = "",
    go_id: str = "GO:0006281",
    aspect: str = "P",
    taxon: str = "taxon:9606",
    references: str = "PMID:12345678",
) -> str:
    fields = [
        "UniProtKB",
        "P12345",
        symbol,
        qualifier,
        go_id,
        references,
        "IDA",
        "",
        aspect,
        "Example protein",
        "",
        "protein",
        taxon,
        "20261001",
        "UniProt",
        "",
        "",
    ]
    return "\t".join(fields)


def _write_gaf(tmp_path: Path, *rows: str) -> Path:
    path = tmp_path / "test.gaf"
    path.write_text(
        "!gaf-version: 2.2\n" + "\n".join(rows) + "\n",
        encoding="utf-8",
    )
    return path


def test_imports_annotations_and_preserves_provenance(
    tmp_path: Path,
) -> None:
    path = _write_gaf(tmp_path, _gaf_row())

    result = read_gaf_annotations(
        path,
        source_release="test-release",
        source_url=SOURCE_URL,
        term_labels={"GO:0006281": "DNA repair"},
    )

    assert result.rows_read == 1
    assert result.row_errors == ()
    assert result.rows_filtered_by_taxon == 0
    assert len(result.annotations) == 1

    annotation = result.annotations[0]
    assert annotation.gene_symbol == "GENE1"
    assert annotation.annotation_label == "DNA repair"
    assert annotation.annotation_type == "biological_process"
    assert annotation.evidence_code == "IDA"
    assert annotation.reference_ids == ("PMID:12345678",)
    assert annotation.source_release == "test-release"
    assert annotation.taxon_id == "taxon:9606"


def test_preserves_negative_qualifier(tmp_path: Path) -> None:
    path = _write_gaf(
        tmp_path,
        _gaf_row(qualifier="NOT|involved_in"),
    )

    result = read_gaf_annotations(
        path,
        source_release="test-release",
        source_url=SOURCE_URL,
        term_labels={"GO:0006281": "DNA repair"},
    )

    assert result.row_errors == ()
    assert len(result.annotations) == 1
    assert result.annotations[0].qualifier == "NOT|involved_in"


def test_filters_nonhuman_taxon(tmp_path: Path) -> None:
    path = _write_gaf(
        tmp_path,
        _gaf_row(taxon="taxon:10090"),
    )

    result = read_gaf_annotations(
        path,
        source_release="test-release",
        source_url=SOURCE_URL,
        term_labels={"GO:0006281": "DNA repair"},
    )

    assert result.annotations == ()
    assert result.row_errors == ()
    assert result.rows_filtered_by_taxon == 1


def test_filters_row_when_human_taxon_is_not_primary(
    tmp_path: Path,
) -> None:
    path = _write_gaf(
        tmp_path,
        _gaf_row(taxon="taxon:562|taxon:9606"),
    )

    result = read_gaf_annotations(
        path,
        source_release="test-release",
        source_url=SOURCE_URL,
        term_labels={"GO:0006281": "DNA repair"},
    )

    assert result.annotations == ()
    assert result.row_errors == ()
    assert result.rows_filtered_by_taxon == 1


def test_reports_missing_term_labels(tmp_path: Path) -> None:
    path = _write_gaf(tmp_path, _gaf_row())

    result = read_gaf_annotations(
        path,
        source_release="test-release",
        source_url=SOURCE_URL,
        term_labels={},
    )

    assert result.annotations == ()
    assert len(result.row_errors) == 1
    assert "no ontology label" in result.row_errors[0]


def test_reports_malformed_rows(tmp_path: Path) -> None:
    path = _write_gaf(tmp_path, "too\tfew\tcolumns")

    result = read_gaf_annotations(
        path,
        source_release="test-release",
        source_url=SOURCE_URL,
        term_labels={},
    )

    assert result.rows_read == 1
    assert result.annotations == ()
    assert len(result.row_errors) == 1
    assert "expected 17 GAF columns" in result.row_errors[0]


def test_maps_all_three_go_aspects(tmp_path: Path) -> None:
    path = _write_gaf(
        tmp_path,
        _gaf_row(go_id="GO:0006281", aspect="P"),
        _gaf_row(go_id="GO:0005515", aspect="F"),
        _gaf_row(go_id="GO:0005737", aspect="C"),
    )

    result = read_gaf_annotations(
        path,
        source_release="test-release",
        source_url=SOURCE_URL,
        term_labels={
            "GO:0006281": "DNA repair",
            "GO:0005515": "protein binding",
            "GO:0005737": "cytoplasm",
        },
    )

    assert result.row_errors == ()
    assert len(result.annotations) == 3
    assert [item.annotation_type for item in result.annotations] == [
        "biological_process",
        "molecular_function",
        "cellular_component",
    ]
