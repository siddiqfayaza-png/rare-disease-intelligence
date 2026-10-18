"""Integrity checks for a parsed HPO ontology and graph."""

from dataclasses import dataclass

from rare_disease_ai.hpo_graph import HPOGraph
from rare_disease_ai.hpo_ontology import HPOOntology


@dataclass(frozen=True)
class OntologyValidationReport:
    """Summary of structural ontology integrity checks."""

    term_count: int
    missing_parent_references: tuple[tuple[str, str], ...]
    inconsistent_parent_child_links: tuple[tuple[str, str], ...]
    missing_replacement_references: tuple[tuple[str, str], ...]
    cycle_detected: bool

    @property
    def is_valid(self) -> bool:
        """Return whether all structural checks passed."""
        return (
            not self.missing_parent_references
            and not self.inconsistent_parent_child_links
            and not self.missing_replacement_references
            and not self.cycle_detected
        )


def validate_ontology(
    ontology: HPOOntology,
    graph: HPOGraph,
) -> OntologyValidationReport:
    """Check references, graph consistency, and directed cycles."""
    missing_parents: set[tuple[str, str]] = set()
    inconsistent_links: set[tuple[str, str]] = set()
    missing_replacements: set[tuple[str, str]] = set()

    for hpo_id, term in ontology.terms.items():
        for parent_id in term.parents:
            if parent_id not in ontology.terms:
                missing_parents.add((hpo_id, parent_id))
            elif hpo_id not in graph.children.get(
                parent_id, set()
            ) or parent_id not in graph.parents.get(hpo_id, set()):
                inconsistent_links.add((hpo_id, parent_id))

        for replacement_id in term.replaced_by:
            if replacement_id not in ontology.terms:
                missing_replacements.add((hpo_id, replacement_id))

    # Detect cycles using depth-first traversal with three states:
    # 0 = unseen, 1 = currently visiting, 2 = fully processed.
    state: dict[str, int] = {hpo_id: 0 for hpo_id in ontology.terms}
    cycle_detected = False

    def visit(hpo_id: str) -> None:
        nonlocal cycle_detected

        state[hpo_id] = 1

        for parent_id in graph.parents.get(hpo_id, set()):
            if state[parent_id] == 1:
                cycle_detected = True
            elif state[parent_id] == 0:
                visit(parent_id)

        state[hpo_id] = 2

    for hpo_id in ontology.terms:
        if state[hpo_id] == 0:
            visit(hpo_id)

    return OntologyValidationReport(
        term_count=len(ontology.terms),
        missing_parent_references=tuple(sorted(missing_parents)),
        inconsistent_parent_child_links=tuple(sorted(inconsistent_links)),
        missing_replacement_references=tuple(sorted(missing_replacements)),
        cycle_detected=cycle_detected,
    )
