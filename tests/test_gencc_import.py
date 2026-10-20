import csv

import pytest

from rare_disease_ai.gencc_import import (
    REQUIRED_COLUMNS,
    read_gencc_csv,
)


def write_csv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=REQUIRED_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)


def valid_row(**overrides):
    row = {column: "" for column in REQUIRED_COLUMNS}
    row.update(
        {
            "sgc_id": "SGC-100001",
            "version_number": "1",
            "gene_curie": "HGNC:10896",
            "gene_symbol": "SKI",
            "disease_curie": "MONDO:0008426",
            "disease_title": "Example disease",
            "disease_original_curie": "OMIM:182212",
            "classification_curie": "GENCC:100001",
            "classification_title": "Definitive",
            "moi_curie": "HP:0000006",
            "moi_title": "Autosomal dominant",
            "submitter_curie": "GENCC:000101",
            "submitter_title": "Example submitter",
            "submitted_as_date": "2024-01-01",
            "submitted_as_pmids": "PMID:12345; PMID:67890",
        }
    )
    row.update(overrides)
    return row


def test_imports_assertion_and_preserves_both_disease_ids(tmp_path):
    path = tmp_path / "gencc.csv"
    write_csv(path, [valid_row()])

    result = read_gencc_csv(path)

    assert result.rows_read == 1
    assert result.row_errors == ()
    assert len(result.assertions) == 1
    assertion = result.assertions[0]
    assert assertion.disease_curie == "MONDO:0008426"
    assert assertion.disease_original_curie == "OMIM:182212"
    assert assertion.submission_id == "SGC-100001"
    assert assertion.pmids == ("PMID:12345", "PMID:67890")


def test_reports_bad_rows_without_discarding_valid_rows(tmp_path):
    path = tmp_path / "gencc.csv"
    write_csv(
        path,
        [
            valid_row(),
            valid_row(sgc_id="SGC-BAD", version_number="invalid"),
        ],
    )

    result = read_gencc_csv(path)

    assert result.rows_read == 2
    assert len(result.assertions) == 1
    assert len(result.row_errors) == 1
    assert "line 3" in result.row_errors[0]


def test_rejects_missing_required_columns(tmp_path):
    path = tmp_path / "gencc.csv"
    path.write_text("sgc_id,gene_symbol\nSGC-1,SKI\n", encoding="utf-8")

    with pytest.raises(ValueError, match="missing required columns"):
        read_gencc_csv(path)


def test_empty_csv_returns_empty_result(tmp_path):
    path = tmp_path / "gencc.csv"
    write_csv(path, [])

    result = read_gencc_csv(path)

    assert result.rows_read == 0
    assert result.assertions == ()
    assert result.row_errors == ()
