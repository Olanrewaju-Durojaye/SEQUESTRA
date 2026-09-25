"""Frozen, auditable scientific decision rules for SEQUESTRA."""

from __future__ import annotations

from dataclasses import dataclass
from math import ceil, isfinite
from typing import Any, Iterable


CATALYTIC_REFERENCE = "catalytic-reference"
NON_CATALYTIC_REFERENCE = "non-catalytic-reference"
PROJECT_MODES = (CATALYTIC_REFERENCE, NON_CATALYTIC_REFERENCE)


def shortlist_recommendation(number_of_designs: int) -> dict[str, Any]:
    """Return a transparent, scale-aware suggestion; the user remains in control."""
    if number_of_designs < 1:
        raise ValueError("number_of_designs must be positive")
    if number_of_designs <= 100:
        low, high = 10, 20
    elif number_of_designs <= 1000:
        low, high = 5, 10
    else:
        low, high = 1, 5
    return {
        "minimum_percent": low,
        "maximum_percent": high,
        "rationale": (
            "Later kinetic and Boltz-2 stages are more computationally expensive; "
            "choose the final percentage only after inspecting the BoltzGen ranking."
        ),
        "binding": False,
    }


@dataclass(frozen=True)
class DecisionSettings:
    """Numerical settings whose values must be fixed before a run starts."""

    top_fraction: float | None = None
    affinity_comparison_margin: float = 0.0
    catalytic_liability_enabled: bool = False

    def validate(self, mode: str) -> None:
        if mode not in PROJECT_MODES:
            raise ValueError(f"mode must be one of: {', '.join(PROJECT_MODES)}")
        if self.top_fraction is not None and not 0 < self.top_fraction <= 1:
            raise ValueError("top_fraction must be greater than 0 and at most 1")
        if not isfinite(self.affinity_comparison_margin) or self.affinity_comparison_margin < 0:
            raise ValueError("affinity_comparison_margin must be finite and nonnegative")


def retain_top_fraction(
    candidates: Iterable[dict[str, Any]], fraction: float
) -> list[dict[str, Any]]:
    """Rank by BoltzGen affinity probability and retain ceil(N*fraction)."""
    rows = list(candidates)
    if not 0 < fraction <= 1:
        raise ValueError("fraction must be greater than 0 and at most 1")
    for row in rows:
        probability = row.get("boltzgen_affinity_probability")
        if not isinstance(probability, (int, float)) or not isfinite(probability):
            raise ValueError("each candidate requires a finite BoltzGen affinity probability")
    ranked = sorted(
        rows,
        key=lambda row: (-float(row["boltzgen_affinity_probability"]), str(row.get("candidate_id", ""))),
    )
    return ranked[: ceil(len(ranked) * fraction)]


def has_catalytic_liability(
    *, candidate_kcat: float, candidate_km: float, reference_kcat: float, reference_km: float
) -> bool:
    """Exclude only for strictly higher kcat AND strictly lower Km than reference."""
    return candidate_kcat > reference_kcat and candidate_km < reference_km


def classify_catalytic_liability(
    *,
    candidate_kcat: float | None,
    candidate_km: float | None,
    reference_kcat: float | None,
    reference_km: float | None,
    within_supported_interpretation: bool = True,
) -> dict[str, str]:
    """Classify reference-relative predictions without equating predictions to activity."""
    if not within_supported_interpretation:
        return {
            "classification": "prediction_outside_supported_interpretation",
            "detail": "kinetic prediction is outside supported interpretation",
            "decision": "manual_review",
        }
    values = (candidate_kcat, candidate_km, reference_kcat, reference_km)
    if any(value is None for value in values):
        return {
            "classification": "prediction_unavailable",
            "detail": "one or more required predictions are unavailable",
            "decision": "indeterminate_never_treat_as_non_catalytic",
        }
    if any(not isfinite(float(value)) or float(value) < 0 for value in values):
        return {
            "classification": "prediction_outside_supported_interpretation",
            "detail": "one or more predictions are non-finite or negative",
            "decision": "manual_review",
        }
    high_kcat = float(candidate_kcat) > float(reference_kcat)
    low_km = float(candidate_km) < float(reference_km)
    if high_kcat and low_km:
        return {
            "classification": "potential_catalytic_liability_signal",
            "detail": "higher_kcat_and_lower_km_than_reference",
            "decision": "exclude",
        }
    if high_kcat or low_km:
        detail = "higher_kcat_only" if high_kcat else "lower_km_only"
        return {
            "classification": "potential_catalytic_liability_signal",
            "detail": detail,
            "decision": "advance_with_partial_liability_flag",
        }
    return {
        "classification": "no_catalytic_liability_signal_detected",
        "detail": "neither_reference_relative_adverse_condition_detected",
        "decision": "advance",
    }


def has_stronger_affinity(
    *, candidate_affinity: float, reference_affinity: float, margin: float = 0.0
) -> bool:
    """Compare Boltz-2 affinity scores where larger values mean stronger affinity."""
    return candidate_affinity > reference_affinity + margin


def mode_workflow(mode: str, liability_enabled: bool = False) -> list[dict[str, Any]]:
    """Return the ordered workflow and its frozen decisions for a project mode."""
    if mode not in PROJECT_MODES:
        raise ValueError(f"unknown project mode: {mode}")
    common = [
        {"stage": "boltzgen", "action": "generate_binders"},
        {
            "stage": "boltzgen",
            "action": "rank_candidates",
            "metric": "affinity_probability_binary",
            "direction": "higher_is_better",
        },
        {
            "stage": "shortlist",
            "action": "retain_top_fraction",
            "rounding": "ceiling",
            "minimum_if_candidates_exist": 1,
        },
    ]
    liability = {
        "stage": "catalytic_liability",
        "action": "predict_kcat_and_km",
        "decision": "exclude_only_if_both_strict_reference_relative_conditions_are_true",
        "kcat_condition": "candidate_kcat > reference_kcat",
        "km_condition": "candidate_km < reference_km",
    }
    affinity = {
        "stage": "boltz2_affinity",
        "action": "benchmark_against_reference",
        "decision": "candidate_affinity > reference_affinity + comparison_margin",
        "strict_comparison": True,
    }
    confidence = {
        "stage": "structural_confidence",
        "action": "assess_and_rank_affinity_qualified_candidates",
    }
    if mode == CATALYTIC_REFERENCE:
        return common + [liability, affinity, confidence]
    optional = dict(liability)
    optional.update({
        "action": "predict_kcat_and_km_as_optional_liability_indicators",
        "decision": "informational_only_no_candidate_exclusion",
        "optional": True,
        "enabled": liability_enabled,
        "interpretation": (
            "prediction-only catalytic-liability flag; not evidence that a candidate "
            "is catalytic or non-catalytic"
        ),
    })
    return common + [optional, affinity, confidence]
