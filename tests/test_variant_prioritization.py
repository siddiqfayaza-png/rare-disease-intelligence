import pytest

from rare_disease_ai.variant_prioritization import (
    VariantRecord,
    prioritize_variants,
    validate_variant,
)


def make_variant(**overrides):
    values = {
        "chromosome": "1",
        "position": 12345,
        "reference": "A",
        "alternate": "T",
        "gene_symbol": "GENE1",
        "consequence": "missense_variant",
        "zygosity": "heterozygous",
        "allele_frequency": 0.001,
    }
    values.update(overrides)
    return VariantRecord(**values)


def test_valid_variant_passes_validation():
    result = validate_variant(make_variant())
    assert result.is_valid
    assert result.errors == ()


def test_invalid_position_fails_validation():
    result = validate_variant(make_variant(position=0))
    assert not result.is_valid
    assert any("Position" in error for error in result.errors)


def test_invalid_allele_frequency_fails_validation():
    result = validate_variant(make_variant(allele_frequency=1.5))
    assert not result.is_valid


def test_unsupported_consequence_fails_validation():
    result = validate_variant(make_variant(consequence="unknown_effect"))
    assert not result.is_valid


def test_higher_consequence_score_ranks_first():
    variants = (
        make_variant(
            position=200,
            gene_symbol="GENE2",
            consequence="missense_variant",
        ),
        make_variant(
            position=100,
            gene_symbol="GENE1",
            consequence="frameshift_variant",
        ),
    )

    ranked = prioritize_variants(variants)

    assert ranked[0].variant.consequence == "frameshift_variant"


def test_missing_frequency_does_not_get_rarity_bonus():
    variant = make_variant(allele_frequency=None)
    ranked = prioritize_variants((variant,))

    assert ranked[0].score == 2.0
    assert any("unavailable" in reason for reason in ranked[0].reasons)


def test_invalid_variants_are_excluded_from_ranking():
    variants = (make_variant(position=-1), make_variant(position=123))

    ranked = prioritize_variants(variants)

    assert len(ranked) == 1
    assert ranked[0].variant.position == 123


def test_invalid_frequency_threshold_raises_error():
    with pytest.raises(ValueError):
        prioritize_variants((), max_allele_frequency=2.0)
