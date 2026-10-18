"""Index HPO disease annotations for phenotype-based candidate ranking."""

from __future__ import annotations

from dataclasses import dataclass

from rare_disease_ai.hpo_annotations import HPOAnnotation


@dataclass(frozen=True)
class DiseasePhenotypeProfile:
    """Phenotype annotations associated with one disease."""

    database_id: str
    disease_name: str
    positive_hpo_ids: frozenset[str]
    negative_hpo_ids: frozenset[str]
    annotations: tuple[HPOAnnotation, ...]
    disease_name_variants: tuple[str, ...] = ()


def build_disease_phenotype_index(
    annotations: tuple[HPOAnnotation, ...],
) -> dict[str, DiseasePhenotypeProfile]:
    """Group annotations by disease and preserve disease-name variants."""
    grouped: dict[str, list[HPOAnnotation]] = {}

    for annotation in annotations:
        grouped.setdefault(annotation.database_id, []).append(annotation)

    index: dict[str, DiseasePhenotypeProfile] = {}

    for database_id, disease_annotations in grouped.items():
        # Keep original spellings, ignoring only case and repeated whitespace
        # when identifying duplicate variants.
        disease_names: dict[str, str] = {}

        for annotation in disease_annotations:
            name = " ".join(annotation.disease_name.split())

            if name:
                normalized_name = name.casefold()
                disease_names.setdefault(normalized_name, name)

        name_variants = tuple(sorted(disease_names.values(), key=str.casefold))
        disease_name = name_variants[0] if name_variants else database_id

        positive_hpo_ids = frozenset(
            annotation.hpo_id
            for annotation in disease_annotations
            if annotation.qualifier.strip().upper() != "NOT"
        )

        negative_hpo_ids = frozenset(
            annotation.hpo_id
            for annotation in disease_annotations
            if annotation.qualifier.strip().upper() == "NOT"
        )

        index[database_id] = DiseasePhenotypeProfile(
            database_id=database_id,
            disease_name=disease_name,
            positive_hpo_ids=positive_hpo_ids,
            negative_hpo_ids=negative_hpo_ids,
            annotations=tuple(disease_annotations),
            disease_name_variants=name_variants,
        )

    return index
