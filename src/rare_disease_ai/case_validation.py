"""Validation models for structured patient phenotype cases."""

from dataclasses import dataclass


@dataclass(frozen=True)
class PatientCase:
    """A minimal, validated representation of a research patient case."""

    case_id: str
    phenotype_terms: tuple[str, ...]

    def __post_init__(self) -> None:
        """Validate the case fields and normalize whitespace."""
        if not isinstance(self.case_id, str) or not self.case_id.strip():
            raise ValueError("case_id must be a non-empty string.")

        if not isinstance(self.phenotype_terms, tuple):
            raise TypeError("phenotype_terms must be a tuple of strings.")

        if not self.phenotype_terms:
            raise ValueError("At least one phenotype term is required.")

        if any(not isinstance(term, str) or not term.strip() for term in self.phenotype_terms):
            raise ValueError("Each phenotype term must be a non-empty string.")

        object.__setattr__(self, "case_id", self.case_id.strip())
        object.__setattr__(
            self,
            "phenotype_terms",
            tuple(term.strip() for term in self.phenotype_terms),
        )
