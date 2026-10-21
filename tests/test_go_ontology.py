from pathlib import Path

import pytest

from rare_disease_ai.go_ontology import read_go_term_labels


def test_reads_active_go_term_labels(tmp_path: Path) -> None:
    obo_file = tmp_path / "go.obo"
    obo_file.write_text(
        "[Term]\n"
        "id: GO:0006281\n"
        "name: DNA repair\n"
        "\n"
        "[Term]\n"
        "id: GO:0003677\n"
        "name: DNA binding\n",
        encoding="utf-8",
    )

    assert read_go_term_labels(obo_file) == {
        "GO:0006281": "DNA repair",
        "GO:0003677": "DNA binding",
    }


def test_skips_obsolete_terms(tmp_path: Path) -> None:
    obo_file = tmp_path / "go.obo"
    obo_file.write_text(
        "[Term]\n"
        "id: GO:0000001\n"
        "name: obsolete example\n"
        "is_obsolete: true\n"
        "\n"
        "[Term]\n"
        "id: GO:0000002\n"
        "name: active example\n",
        encoding="utf-8",
    )

    assert read_go_term_labels(obo_file) == {
        "GO:0000002": "active example",
    }


def test_ignores_typedef_and_non_go_terms(tmp_path: Path) -> None:
    obo_file = tmp_path / "go.obo"
    obo_file.write_text(
        "[Term]\n"
        "id: CL:0000000\n"
        "name: non-GO term\n"
        "\n"
        "[Typedef]\n"
        "id: part_of\n"
        "name: part of\n",
        encoding="utf-8",
    )

    assert read_go_term_labels(obo_file) == {}


def test_rejects_duplicate_active_go_term_ids(tmp_path: Path) -> None:
    obo_file = tmp_path / "go.obo"
    obo_file.write_text(
        "[Term]\n"
        "id: GO:0006281\n"
        "name: DNA repair\n"
        "\n"
        "[Term]\n"
        "id: GO:0006281\n"
        "name: duplicate DNA repair\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="Duplicate active GO term ID"):
        read_go_term_labels(obo_file)
