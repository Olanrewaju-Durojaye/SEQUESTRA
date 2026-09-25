from __future__ import annotations

import unittest

from sequestra.decision_rules import (
    CATALYTIC_REFERENCE,
    NON_CATALYTIC_REFERENCE,
    DecisionSettings,
    classify_catalytic_liability,
    has_catalytic_liability,
    has_stronger_affinity,
    mode_workflow,
    retain_top_fraction,
)


class DecisionRuleTests(unittest.TestCase):
    def test_top_fraction_is_ranked_and_rounded_up(self) -> None:
        candidates = [
            {"candidate_id": "b", "boltzgen_affinity_probability": 0.8},
            {"candidate_id": "a", "boltzgen_affinity_probability": 0.8},
            {"candidate_id": "c", "boltzgen_affinity_probability": 0.2},
        ]
        kept = retain_top_fraction(candidates, 0.5)
        self.assertEqual([row["candidate_id"] for row in kept], ["a", "b"])

    def test_catalytic_exclusion_requires_both_conditions(self) -> None:
        self.assertTrue(has_catalytic_liability(
            candidate_kcat=0.11, candidate_km=0.0009,
            reference_kcat=0.1, reference_km=0.001,
        ))
        self.assertFalse(has_catalytic_liability(
            candidate_kcat=0.1, candidate_km=0.0009,
            reference_kcat=0.1, reference_km=0.001,
        ))
        self.assertFalse(has_catalytic_liability(
            candidate_kcat=0.11, candidate_km=0.001,
            reference_kcat=0.1, reference_km=0.001,
        ))

    def test_affinity_comparison_is_strict(self) -> None:
        self.assertFalse(has_stronger_affinity(candidate_affinity=1.0, reference_affinity=1.0))
        self.assertTrue(has_stronger_affinity(candidate_affinity=1.01, reference_affinity=1.0))

    def test_modes_allow_deferred_shortlist_selection(self) -> None:
        DecisionSettings().validate(CATALYTIC_REFERENCE)
        DecisionSettings(0.1, catalytic_liability_enabled=True).validate(NON_CATALYTIC_REFERENCE)
        workflow = mode_workflow(NON_CATALYTIC_REFERENCE, True)
        liability = next(row for row in workflow if row["stage"] == "catalytic_liability")
        self.assertTrue(liability["optional"])
        self.assertEqual(liability["decision"], "informational_only_no_candidate_exclusion")
        self.assertIn("prediction-only", liability["interpretation"])

    def test_reference_relative_classifications(self) -> None:
        both = classify_catalytic_liability(
            candidate_kcat=2.0, candidate_km=0.5,
            reference_kcat=1.0, reference_km=1.0,
        )
        self.assertEqual(both["decision"], "exclude")
        partial = classify_catalytic_liability(
            candidate_kcat=2.0, candidate_km=1.0,
            reference_kcat=1.0, reference_km=1.0,
        )
        self.assertEqual(partial["decision"], "advance_with_partial_liability_flag")
        unavailable = classify_catalytic_liability(
            candidate_kcat=None, candidate_km=1.0,
            reference_kcat=1.0, reference_km=1.0,
        )
        self.assertIn("never_treat_as_non_catalytic", unavailable["decision"])


if __name__ == "__main__":
    unittest.main()
