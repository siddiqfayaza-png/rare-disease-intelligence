"""Parse active Gene Ontology term identifiers and labels from OBO files."""

from pathlib import Path


def read_go_term_labels(path: str | Path) -> dict[str, str]:
    """Return active GO term IDs mapped to names from an OBO file."""
    labels: dict[str, str] = {}
    in_term = False
    term_id: str | None = None
    term_name: str | None = None
    is_obsolete = False

    def save_term() -> None:
        if not in_term or is_obsolete:
            return
        if term_id is None or term_name is None:
            return
        if not term_id.startswith("GO:"):
            return
        if term_id in labels:
            raise ValueError(f"Duplicate active GO term ID: {term_id}")
        labels[term_id] = term_name

    with Path(path).open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()

            if line.startswith("["):
                save_term()
                in_term = line == "[Term]"
                term_id = None
                term_name = None
                is_obsolete = False
                continue

            if not in_term or not line or line.startswith("!"):
                continue

            key, separator, value = line.partition(": ")
            if not separator:
                continue

            if key == "id":
                term_id = value.strip()
            elif key == "name":
                term_name = value.strip()
            elif key == "is_obsolete" and value.strip().lower() == "true":
                is_obsolete = True

    save_term()
    return labels
