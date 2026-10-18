"""Graph utilities for Human Phenotype Ontology relationships."""

from __future__ import annotations

from collections import defaultdict

from rare_disease_ai.hpo_ontology import HPOOntology


class HPOGraph:
    """Represent HPO parent relationships and support ancestor queries."""

    def __init__(self, ontology: HPOOntology) -> None:
        self.ontology = ontology
        self.parents: dict[str, set[str]] = {}
        self.children: dict[str, set[str]] = defaultdict(set)

        for hpo_id, term in ontology.terms.items():
            known_parents = {parent_id for parent_id in term.parents if parent_id in ontology.terms}

            self.parents[hpo_id] = known_parents

            for parent_id in known_parents:
                self.children[parent_id].add(hpo_id)

    def get_parents(self, hpo_id: str) -> tuple[str, ...]:
        """Return the direct parents of an HPO term."""
        if hpo_id not in self.ontology.terms:
            raise KeyError(f"Unknown HPO identifier: {hpo_id}")

        return tuple(sorted(self.parents[hpo_id]))

    def get_children(self, hpo_id: str) -> tuple[str, ...]:
        """Return the direct children of an HPO term."""
        if hpo_id not in self.ontology.terms:
            raise KeyError(f"Unknown HPO identifier: {hpo_id}")

        return tuple(sorted(self.children.get(hpo_id, set())))

    def get_ancestors(self, hpo_id: str) -> tuple[str, ...]:
        """Return all ancestors, excluding the term itself."""
        if hpo_id not in self.ontology.terms:
            raise KeyError(f"Unknown HPO identifier: {hpo_id}")

        visited: set[str] = set()
        pending = list(self.parents[hpo_id])

        while pending:
            current = pending.pop()

            if current in visited:
                continue

            visited.add(current)
            pending.extend(self.parents.get(current, set()) - visited)

        return tuple(sorted(visited))
