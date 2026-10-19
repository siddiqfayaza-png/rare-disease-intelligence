"""Build toy disease profiles for pipeline validation only."""

from rare_disease_ai.disease_phenotypes import DiseasePhenotypeProfile


def build_synthetic_disease_index() -> dict[str, DiseasePhenotypeProfile]:
    """Return invented disease profiles for controlled software tests."""
    profiles = (
        DiseasePhenotypeProfile(
            database_id="SYNTHETIC:DISEASE_001",
            disease_name="Synthetic Disease 001",
            positive_hpo_ids=frozenset({"HP:0001250", "HP:0001263"}),
            negative_hpo_ids=frozenset(),
            annotations=(),
        ),
        DiseasePhenotypeProfile(
            database_id="SYNTHETIC:DISEASE_002",
            disease_name="Synthetic Disease 002",
            positive_hpo_ids=frozenset({"HP:0001250", "HP:0001290"}),
            negative_hpo_ids=frozenset(),
            annotations=(),
        ),
        DiseasePhenotypeProfile(
            database_id="SYNTHETIC:DISEASE_003",
            disease_name="Synthetic Disease 003",
            positive_hpo_ids=frozenset({"HP:0001263", "HP:0001290"}),
            negative_hpo_ids=frozenset(),
            annotations=(),
        ),
    )

    return {profile.database_id: profile for profile in profiles}
