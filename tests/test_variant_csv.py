import pytest

from rare_disease_ai.variant_prioritization import (
    read_variant_csv,
)

CSV_HEADER = (
    "chromosome,position,reference,alternate,gene_symbol,consequence,zygosity,allele_frequency\n"
)


def test_reads_valid_variant_csv(tmp_path):
    csv_path = tmp_path / "variants.csv"
    csv_path.write_text(
        CSV_HEADER + "1,12345,A,T,GENE1,missense_variant,heterozygous,0.001\n",
        encoding="utf-8",
    )

    result = read_variant_csv(csv_path)

    assert result.rows_read == 1
    assert len(result.variants) == 1
    assert result.variants[0].gene_symbol == "GENE1"
    assert result.row_errors == ()


def test_invalid_row_is_reported_and_valid_row_is_kept(tmp_path):
    csv_path = tmp_path / "variants.csv"
    csv_path.write_text(
        CSV_HEADER
        + "1,not-a-number,A,T,GENE1,missense_variant,heterozygous,0.001\n"
        + "2,22222,C,G,GENE2,frameshift_variant,heterozygous,0.002\n",
        encoding="utf-8",
    )

    result = read_variant_csv(csv_path)

    assert result.rows_read == 2
    assert len(result.variants) == 1
    assert result.variants[0].gene_symbol == "GENE2"
    assert len(result.row_errors) == 1
    assert result.row_errors[0].startswith("Row 2:")


def test_missing_required_column_raises_error(tmp_path):
    csv_path = tmp_path / "variants.csv"
    csv_path.write_text(
        "chromosome,position,reference,alternate\n1,123,A,T\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="missing required columns"):
        read_variant_csv(csv_path)


def test_blank_frequency_is_allowed(tmp_path):
    csv_path = tmp_path / "variants.csv"
    csv_path.write_text(
        CSV_HEADER + "1,12345,A,T,GENE1,missense_variant,heterozygous,\n",
        encoding="utf-8",
    )

    result = read_variant_csv(csv_path)

    assert len(result.variants) == 1
    assert result.variants[0].allele_frequency is None


def test_invalid_frequency_is_reported(tmp_path):
    csv_path = tmp_path / "variants.csv"
    csv_path.write_text(
        CSV_HEADER + "1,12345,A,T,GENE1,missense_variant,heterozygous,not-a-number\n",
        encoding="utf-8",
    )

    result = read_variant_csv(csv_path)

    assert result.variants == ()
    assert len(result.row_errors) == 1
    assert result.row_errors[0].startswith("Row 2:")
