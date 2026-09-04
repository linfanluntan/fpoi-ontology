#!/usr/bin/env python3
"""
FPOI V7.5 -- PROTEAS treatment-target pairing and Rx-resolution audit.

Implements the gates that V7.4 lacked (audit items B2, C10, C13):

  G0  Rx resolution        Can a single, unambiguous D_Rx be established for the
                           archived course, reconciling RTPLAN DICOM against the
                           clinical spreadsheet? Detects margin-vs-maximum
                           prescription ambiguity (the P37 20 Gy / 35 Gy case).
  G1  Dose-scale QA        Existing V7.4 gate, reimplemented like-for-like:
                           reports P99.9/P99.9 and Dmax/Dmax alongside the
                           original P99.9/Dmax ratio.
  G2  Target identity      Does the baseline lesion coincide with a high-dose
                           island of the archived plan? Connected-component
                           analysis of the >= D_Rx isodose in registered space.
  G3  Rx plausibility      Is the baseline lesion Dmean/D95 consistent with the
                           lesion having been the prescribed target?
  G4  Course pairing       For multi-course archives (P17a/b, P23a/b): which
                           course best explains this lesion, and by what margin?

Everything is evidence-first: the script emits per-lesion measurements and a
verdict, and never rescales, reassigns, or repairs anything silently.

Design note on coordinates
--------------------------
Plan isocentre lives in DICOM patient coordinates; PROTEAS lesion masks live in
processed-MRI/BraTS space. The transform is not published, so isocentre distance
is NOT computable. G2 therefore works entirely inside the registered RTP NIfTI
grid, using the dose field's own high-dose islands as the target proxy. This is
the reason G2 is a proxy gate and must be reported as such.

Usage
-----
    python proteas_target_pairing_audit.py --manifest manifest.csv \
        --config gates.json --out results/

manifest.csv columns (one row per archive):
    archive_id, patient_id, course_id,
    rtp_nifti, rtplan_dcm, rtdose_dcm,
    baseline_mask, followup_masks, clinical_rx_gy, clinical_fractions

  followup_masks : ';'-separated paths, in acquisition order
  clinical_rx_gy : prescription from the clinical spreadsheet, blank if absent

Outputs
-------
    gate_evidence.csv      one row per (archive, mask) with every measurement
    rx_resolution.csv      one row per archive, all Rx candidates and the verdict
    course_pairing.csv     one row per (patient, lesion, candidate course)
    audit_summary.json     counts, config hash, provenance
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
import pydicom
from scipy import ndimage

# --------------------------------------------------------------------------
# Configuration-frozen gate constants.
#
# These are the analysis constants. Freeze this dict (or the --config file that
# overrides it) in a commit BEFORE inspecting results, and cite that commit in
# the manuscript instead of the word "prespecified" (audit item C9).
# --------------------------------------------------------------------------
DEFAULT_CONFIG = {
    # G0 -- Rx resolution
    "rx_agreement_tol": 0.05,          # rel. tolerance, DICOM vs clinical Rx
    # G1 -- dose-scale QA (V7.4 interval retained + sensitivity sweep)
    "scale_interval": [0.75, 1.25],
    "scale_sensitivity_intervals": [[0.90, 1.10], [0.80, 1.20], [0.70, 1.30]],
    # G2 -- target identity
    "isodose_island_min_cc": 0.02,     # drop islands < 0.02 cc as noise
    "target_overlap_min": 0.50,        # >= 50% of baseline lesion inside island
    "target_centroid_max_mm": 10.0,    # lesion centroid to island centroid
    # G3 -- prescription plausibility
    "baseline_dmean_min_frac": 0.90,   # Dmean >= 0.90 * D_Rx(margin)
    # G4 -- course pairing
    "course_margin_min": 1.5,          # winning course must beat runner-up by 1.5x
    # G5 -- longitudinal lesion correspondence (manifest-supplied prior QA)
    "correspondence_pass_values": ["PASS", "GREEN", "green_plausible"],
    "correspondence_amber_values": ["AMBER", "amber_low_overlap", "amber_centroid_shift"],
    "correspondence_fail_values": ["FAIL", "RED", "red_large_centroid_shift", "red_geometry_mismatch"],
    # numerics
    "nonfinite_policy": "flag",        # 'flag' | 'zero_outside_support'
    "negative_clip": True,
}

AMBER, PASS, FAIL, NOT_EVALUABLE = "AMBER", "PASS", "FAIL", "NOT_EVALUABLE"


# --------------------------------------------------------------------------
# Provenance
# --------------------------------------------------------------------------
def _git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], stderr=subprocess.DEVNULL
        ).decode().strip()
    except Exception:
        return "unavailable"


def _config_hash(cfg: dict) -> str:
    return hashlib.sha256(
        json.dumps(cfg, sort_keys=True).encode()
    ).hexdigest()[:16]


def _file_hash(path) -> str:
    if path is None or not Path(path).exists():
        return "missing"
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()[:16]


# --------------------------------------------------------------------------
# G0 -- Rx resolution
# --------------------------------------------------------------------------
@dataclass
class RxCandidates:
    archive_id: str
    target_prescription_dose: float | None = None
    target_maximum_dose: float | None = None
    dose_reference_type: str | None = None
    dose_reference_description: str | None = None
    n_fractions: int | None = None
    beam_dose_total: float | None = None
    clinical_rx_gy: float | None = None
    clinical_fractions: int | None = None
    implied_isodose_line: float | None = None
    dicom_rx_value: float | None = None
    dicom_rx_source: str | None = None
    fraction_match: bool | None = None
    rx_margin_gy: float | None = None
    rx_source: str = "unresolved"
    verdict: str = FAIL
    notes: list[str] = field(default_factory=list)


def resolve_rx(rtplan_path, clinical_rx, clinical_fx, cfg) -> RxCandidates:
    """Reconcile DICOM and clinical prescription statements.

    Important correction after diagnostic audit:
    - TargetPrescriptionDose is treated as prescription when present.
    - Beam-dose total is the fallback DICOM prescription source.
    - A clinical/DICOM mismatch is NOT automatically reinterpreted as
      prescription-to-isodose-line encoding.
    - DICOM-only evidence remains unresolved for normalized downstream gates.
    """
    rc = RxCandidates(
        archive_id=str(rtplan_path),
        clinical_rx_gy=float(clinical_rx) if pd.notna(clinical_rx) else None,
        clinical_fractions=int(clinical_fx) if pd.notna(clinical_fx) else None,
    )

    if rtplan_path is None or not Path(rtplan_path).exists():
        rc.notes.append("RTPLAN absent")
        return rc

    ds = pydicom.dcmread(str(rtplan_path), stop_before_pixels=True)

    for dr in getattr(ds, "DoseReferenceSequence", []) or []:
        if getattr(dr, "TargetPrescriptionDose", None) is not None:
            rc.target_prescription_dose = float(dr.TargetPrescriptionDose)
        if getattr(dr, "TargetMaximumDose", None) is not None:
            rc.target_maximum_dose = float(dr.TargetMaximumDose)
        rc.dose_reference_type = getattr(dr, "DoseReferenceStructureType", None)
        rc.dose_reference_description = getattr(dr, "DoseReferenceDescription", None)

    beam_total = 0.0
    for fg in getattr(ds, "FractionGroupSequence", []) or []:
        if getattr(fg, "NumberOfFractionsPlanned", None) is not None:
            rc.n_fractions = int(fg.NumberOfFractionsPlanned)
        for rb in getattr(fg, "ReferencedBeamSequence", []) or []:
            bd = getattr(rb, "BeamDose", None)
            if bd is not None:
                beam_total += float(bd) * (rc.n_fractions or 1)
    rc.beam_dose_total = beam_total or None

    if rc.target_prescription_dose is not None:
        rc.dicom_rx_value = rc.target_prescription_dose
        rc.dicom_rx_source = "TargetPrescriptionDose"
    elif rc.beam_dose_total is not None:
        rc.dicom_rx_value = rc.beam_dose_total
        rc.dicom_rx_source = "BeamDoseTotal"

    if rc.clinical_fractions is not None and rc.n_fractions is not None:
        rc.fraction_match = (rc.clinical_fractions == rc.n_fractions)

    tol = cfg["rx_agreement_tol"]

    if rc.clinical_rx_gy is not None and rc.dicom_rx_value is not None:
        ratio = rc.clinical_rx_gy / rc.dicom_rx_value
        dose_agree = abs(ratio - 1.0) <= tol
        fraction_agree = (rc.fraction_match is not False)

        if dose_agree and fraction_agree:
            rc.rx_margin_gy = rc.clinical_rx_gy
            rc.rx_source = f"clinical+{rc.dicom_rx_source}_agree"
            rc.verdict = PASS
            rc.notes.append(
                f"clinical and DICOM prescription agree (ratio {ratio:.3f}); "
                f"DICOM source={rc.dicom_rx_source}")
        else:
            rc.rx_margin_gy = None
            rc.rx_source = "clinical_dicom_linkage_conflict"
            rc.verdict = FAIL
            problems = []
            if not dose_agree:
                problems.append(
                    f"dose mismatch: clinical {rc.clinical_rx_gy:g} Gy vs "
                    f"DICOM {rc.dicom_rx_value:g} Gy from {rc.dicom_rx_source}")
            if not fraction_agree:
                problems.append(
                    f"fraction mismatch: clinical {rc.clinical_fractions} vs "
                    f"DICOM {rc.n_fractions}")
            rc.notes.append(
                "clinical-record-to-course linkage/prescription conflict; " +
                " | ".join(problems))
    elif rc.dicom_rx_value is not None:
        rc.rx_margin_gy = None
        rc.rx_source = f"{rc.dicom_rx_source}_only_unresolved"
        rc.verdict = AMBER
        rc.notes.append(
            f"DICOM prescription available from {rc.dicom_rx_source}, but no "
            "clinical course value is available for reconciliation; normalized "
            "downstream gates are not evaluable")
    elif rc.clinical_rx_gy is not None:
        rc.rx_margin_gy = None
        rc.rx_source = "clinical_only_unresolved"
        rc.verdict = AMBER
        rc.notes.append(
            "clinical prescription available but no DICOM prescription statement "
            "is recoverable; normalized downstream gates are not evaluable")
    else:
        rc.notes.append("no prescription statement from any source")

    if rc.target_maximum_dose and rc.target_prescription_dose:
        rc.notes.append(
            f"TargetMaximumDose/TargetPrescriptionDose = "
            f"{rc.target_maximum_dose / rc.target_prescription_dose:.3f}")

    return rc


# --------------------------------------------------------------------------
# dose field loading + G1
# --------------------------------------------------------------------------
def load_dose(rtp_nifti, cfg):
    img = nib.load(str(rtp_nifti))
    dose = np.asanyarray(img.dataobj, dtype=np.float64)
    nonfinite = int((~np.isfinite(dose)).sum())
    if cfg["nonfinite_policy"] == "zero_outside_support":
        dose = np.nan_to_num(dose, nan=0.0, posinf=0.0, neginf=0.0)
    else:
        dose = np.where(np.isfinite(dose), dose, np.nan)
    negatives = int(np.nansum(dose < 0))
    if cfg["negative_clip"]:
        dose = np.where(np.isnan(dose), np.nan, np.clip(dose, 0, None))
    zooms = np.asarray(img.header.get_zooms()[:3], dtype=float)
    return dose, img.affine, zooms, nonfinite, negatives


def scale_qa(dose, rtdose_path, cfg):
    """Dose-scale fidelity.

    Primary gate: registered Dmax / native DICOM RTDOSE Dmax.
    This avoids the support-dependence of percentile comparisons when the native
    and resampled grids differ in extent and voxel size.

    P99.9/P99.9 is retained as a diagnostic only and must not drive G1.
    """
    out = {"rtp_p999": None, "rtp_max": None,
           "native_p999": None, "native_max": None,
           "ratio_p999_over_max_v74": None,
           "ratio_p999_over_p999": None,
           "ratio_max_over_max": None,
           "scale_metric": "Dmax/Dmax",
           "scale_verdict": FAIL, "scale_notes": []}

    finite = dose[np.isfinite(dose)]
    if finite.size == 0:
        out["scale_notes"].append("registered dose has no finite voxels")
        return out
    out["rtp_p999"] = float(np.percentile(finite, 99.9))
    out["rtp_max"] = float(finite.max())

    if rtdose_path is None or not Path(rtdose_path).exists():
        out["scale_notes"].append("native RTDOSE absent")
        return out

    ds = pydicom.dcmread(str(rtdose_path))
    native = ds.pixel_array.astype(np.float64) * float(ds.DoseGridScaling)
    out["native_p999"] = float(np.percentile(native, 99.9))
    out["native_max"] = float(native.max())

    if out["native_max"] > 0:
        out["ratio_p999_over_max_v74"] = out["rtp_p999"] / out["native_max"]
        out["ratio_max_over_max"] = out["rtp_max"] / out["native_max"]
    if out["native_p999"] > 0:
        out["ratio_p999_over_p999"] = out["rtp_p999"] / out["native_p999"]

    lo, hi = cfg["scale_interval"]
    r = out["ratio_max_over_max"]
    if r is not None and lo <= r <= hi:
        out["scale_verdict"] = PASS
    elif r is not None:
        out["scale_verdict"] = FAIL
        out["scale_notes"].append(
            f"Dmax/Dmax scale ratio {r:.3f} outside [{lo}, {hi}]")

    # Sensitivity sweep applies to the same primary Dmax/Dmax metric.
    for a, b in cfg["scale_sensitivity_intervals"]:
        if r is not None:
            out[f"scale_max_pass_{a}_{b}"] = bool(a <= r <= b)

    if out["ratio_p999_over_p999"] is not None:
        out["scale_notes"].append(
            "P99.9/P99.9 retained as support-dependent diagnostic only: "
            f"{out['ratio_p999_over_p999']:.3f}")
    return out


# --------------------------------------------------------------------------
# G2 -- target identity via high-dose islands
# --------------------------------------------------------------------------
def isodose_islands(dose, zooms, rx_gy, cfg):
    """Connected components of the >= D_Rx isodose, i.e. the treated targets."""
    if rx_gy is None or not np.isfinite(rx_gy) or rx_gy <= 0:
        return []
    binary = np.isfinite(dose) & (dose >= rx_gy)
    if not binary.any():
        return []
    lab, n = ndimage.label(binary)
    vox_cc = float(np.prod(zooms)) / 1000.0  # mm^3 -> cc
    islands = []
    for i in range(1, n + 1):
        m = lab == i
        vol = m.sum() * vox_cc
        if vol < cfg["isodose_island_min_cc"]:
            continue
        cen = np.array(ndimage.center_of_mass(m)) * zooms
        islands.append({"island_id": i, "volume_cc": vol,
                        "centroid_mm": cen, "mask": m})
    islands.sort(key=lambda d: -d["volume_cc"])
    return islands


def mask_stats(dose, mask, zooms):
    vals = dose[mask & np.isfinite(dose)]
    vox_cc = float(np.prod(zooms)) / 1000.0
    if vals.size == 0:
        return {"n_vox": 0, "volume_cc": float(mask.sum()) * vox_cc,
                "dmean": np.nan, "dmedian": np.nan, "dmin": np.nan,
                "dmax": np.nan, "d95": np.nan}
    return {
        "n_vox": int(vals.size),
        "volume_cc": float(mask.sum()) * vox_cc,
        "dmean": float(vals.mean()),
        "dmedian": float(np.median(vals)),
        "dmin": float(vals.min()),
        "dmax": float(vals.max()),
        # D95 = dose received by 95% of the volume = 5th percentile
        "d95": float(np.percentile(vals, 5)),
    }


def target_identity(dose, lesion, zooms, rx_gy, cfg):
    out = {"n_islands": 0, "best_island_volume_cc": None,
           "lesion_in_island_frac": None, "island_centroid_dist_mm": None,
           "identity_verdict": FAIL, "identity_notes": []}

    # Always materialize sensitivity fields so non-evaluable archives remain
    # explicit in downstream tables rather than silently disappearing.
    for t in cfg.get("target_overlap_sensitivity", []):
        out[f"overlap_pass_{t}"] = None
    for t in cfg.get("target_centroid_sensitivity_mm", []):
        out[f"centroid_pass_{t}mm"] = None

    # G2 is not evaluable without a resolved prescription threshold. This is
    # semantically distinct from a genuine pairing failure in an evaluable
    # dose field.
    if rx_gy is None or not np.isfinite(rx_gy) or rx_gy <= 0:
        out["identity_verdict"] = NOT_EVALUABLE
        out["identity_notes"].append(
            "G2 not evaluable because D_Rx is unresolved; no treatment-target "
            "correspondence failure is inferred")
        return out

    islands = isodose_islands(dose, zooms, rx_gy, cfg)
    out["n_islands"] = len(islands)
    if not islands:
        out["identity_notes"].append(
            "resolved D_Rx produced no >= D_Rx isodose island in the registered dose field")
        return out
    if lesion.sum() == 0:
        out["identity_notes"].append("empty lesion mask")
        return out

    lcen = np.array(ndimage.center_of_mass(lesion)) * zooms
    best, best_frac = None, -1.0
    for isl in islands:
        frac = float((lesion & isl["mask"]).sum()) / float(lesion.sum())
        if frac > best_frac:
            best, best_frac = isl, frac

    out["best_island_volume_cc"] = best["volume_cc"]
    out["lesion_in_island_frac"] = best_frac
    out["island_centroid_dist_mm"] = float(
        np.linalg.norm(lcen - best["centroid_mm"]))

    ok_overlap = best_frac >= cfg["target_overlap_min"]
    ok_dist = out["island_centroid_dist_mm"] <= cfg["target_centroid_max_mm"]

    # Post-diagnostic sensitivity outputs; these are descriptive only and do not
    # alter the configured gate verdict.
    for t in cfg.get("target_overlap_sensitivity", []):
        out[f"overlap_pass_{t}"] = bool(best_frac >= float(t))
    for t in cfg.get("target_centroid_sensitivity_mm", []):
        out[f"centroid_pass_{t}mm"] = bool(
            out["island_centroid_dist_mm"] <= float(t)
        )

    if ok_overlap and ok_dist:
        out["identity_verdict"] = PASS
    elif ok_overlap or ok_dist:
        out["identity_verdict"] = AMBER
        out["identity_notes"].append(
            "partial evidence: overlap and centroid criteria disagree")
    else:
        out["identity_verdict"] = FAIL
        out["identity_notes"].append(
            f"baseline lesion does not coincide with any treated target "
            f"(overlap {best_frac:.2f}, centroid distance "
            f"{out['island_centroid_dist_mm']:.1f} mm); lesion-to-plan "
            f"pairing is not established")
    return out


# --------------------------------------------------------------------------
# G3 -- prescription plausibility
# --------------------------------------------------------------------------
def rx_plausibility(stats, rx_gy, cfg):
    """Baseline prescription plausibility for dataset-provided tumor masks.

    Operational gate: Dmean / D_Rx only.
    D95 is retained as a diagnostic because these MRI tumor segmentations are not
    planning PTVs and a PTV-style D95 coverage threshold is not calibrated here.
    """
    out = {"dmean_over_rx": None, "d95_over_rx": None,
           "plausibility_metric": "Dmean/D_Rx",
           "plausibility_verdict": FAIL, "plausibility_notes": []}

    if rx_gy is None or not np.isfinite(rx_gy) or rx_gy <= 0:
        out["plausibility_verdict"] = NOT_EVALUABLE
        out["plausibility_notes"].append(
            "G3 not evaluable because D_Rx is unresolved")
        return out
    if not np.isfinite(stats["dmean"]):
        out["plausibility_notes"].append("no finite dose in lesion")
        return out

    out["dmean_over_rx"] = stats["dmean"] / rx_gy
    out["d95_over_rx"] = stats["d95"] / rx_gy

    threshold = cfg["baseline_dmean_min_frac"]
    if out["dmean_over_rx"] >= threshold:
        out["plausibility_verdict"] = PASS
    else:
        out["plausibility_verdict"] = FAIL
        out["plausibility_notes"].append(
            f"baseline lesion Dmean is {out['dmean_over_rx']:.3f} x D_Rx, "
            f"below operational Dmean threshold {threshold:.2f}")

    # D95 is diagnostic only.
    if np.isfinite(out["d95_over_rx"]):
        out["plausibility_notes"].append(
            f"D95/D_Rx={out['d95_over_rx']:.3f} (diagnostic only; not gated)")

    for t in cfg.get("dmean_sensitivity_thresholds", []):
        out[f"dmean_pass_{t}"] = bool(out["dmean_over_rx"] >= t)

    return out


# --------------------------------------------------------------------------
# G4 -- true cross-course pairing across archives of the same patient
# --------------------------------------------------------------------------
def _same_grid(img_a, img_b, atol=1e-4):
    return (img_a.shape[:3] == img_b.shape[:3]
            and np.allclose(img_a.affine, img_b.affine, atol=atol, rtol=0))


def _g5_from_manifest(value, cfg):
    """Map an externally computed longitudinal correspondence verdict to G5.

    Baseline observations do not need a longitudinal same-lesion gate. Follow-up
    verdicts may be supplied from the prior correspondence QA as a semicolon-
    separated manifest field aligned with followup_masks. Missing values are
    AMBER (unresolved), never silently PASS.
    """
    if value is None or pd.isna(value) or str(value).strip() == "":
        return AMBER, "longitudinal correspondence verdict not supplied"
    v = str(value).strip()
    if v in cfg["correspondence_pass_values"]:
        return PASS, f"source correspondence verdict: {v}"
    if v in cfg["correspondence_amber_values"]:
        return AMBER, f"source correspondence verdict: {v}"
    if v in cfg["correspondence_fail_values"]:
        return FAIL, f"source correspondence verdict: {v}"
    return AMBER, f"unrecognized correspondence verdict retained unresolved: {v}"


def temporal_order_gate(days_from_baseline, cfg):
    """Require follow-up timing to satisfy the configured minimum interval."""
    if days_from_baseline is None or pd.isna(days_from_baseline):
        return AMBER, "temporal ordering unresolved: days_from_baseline missing"
    d = float(days_from_baseline)
    min_days = float(cfg.get("temporal_order_min_days", 0))
    if d < min_days:
        return FAIL, (
            f"temporal ordering failure: days_from_baseline={d:g} "
            f"is below configured minimum {min_days:g}"
        )
    return PASS, (
        f"temporal ordering valid: days_from_baseline={d:g} "
        f">= configured minimum {min_days:g}"
    )


def combine_g5(spatial_v, temporal_v):
    """Composite longitudinal correspondence verdict."""
    if FAIL in (spatial_v, temporal_v):
        return FAIL
    if AMBER in (spatial_v, temporal_v):
        return AMBER
    return PASS


def course_pairing_cross(man, cache, cfg):
    """Score each fixed baseline lesion against every candidate course.

    This corrects the earlier implementation, which compared each course's own
    lesion-to-own-dose score and therefore did not perform a true cross-course
    test. Cross-scoring is attempted only when the lesion mask and candidate
    registered-dose NIfTI are demonstrably on the same grid (shape + affine).
    If candidate course spaces are not comparable, G4 is AMBER/unresolved.
    """
    recs = []
    archive_verdict = {}
    if "patient_id" not in man.columns:
        return pd.DataFrame(recs), archive_verdict

    for pid, pg in man.groupby("patient_id"):
        aids = [str(x) for x in pg["archive_id"].tolist()]
        if len(aids) < 2:
            for aid in aids:
                archive_verdict[aid] = (PASS, "single archived course; G4 not applicable")
            continue

        for _, lesion_row in pg.iterrows():
            lesion_aid = str(lesion_row["archive_id"])
            lesion_path = lesion_row.get("baseline_mask")
            try:
                lesion_img = nib.load(str(lesion_path))
                lesion = np.asanyarray(lesion_img.dataobj) > 0
            except Exception as e:
                archive_verdict[lesion_aid] = (FAIL, f"baseline lesion unreadable: {e}")
                continue

            scored = []
            for _, course_row in pg.iterrows():
                cand_aid = str(course_row["archive_id"])
                c = cache.get(cand_aid)
                if not c:
                    scored.append({"candidate_course": cand_aid, "comparable_grid": False,
                                   "score": np.nan, "identity_verdict": AMBER,
                                   "notes": "candidate course missing from cache"})
                    continue
                try:
                    cand_img = nib.load(str(course_row["rtp_nifti"]))
                except Exception as e:
                    scored.append({"candidate_course": cand_aid, "comparable_grid": False,
                                   "score": np.nan, "identity_verdict": AMBER,
                                   "notes": f"candidate RTP unreadable: {e}"})
                    continue
                if not _same_grid(lesion_img, cand_img):
                    scored.append({"candidate_course": cand_aid, "comparable_grid": False,
                                   "score": np.nan, "identity_verdict": AMBER,
                                   "notes": "baseline lesion and candidate course are not on the same registered grid"})
                    continue

                rx = c["rx"].rx_margin_gy
                if rx is None or not np.isfinite(rx) or rx <= 0:
                    scored.append({"candidate_course": cand_aid, "comparable_grid": True,
                                   "score": np.nan, "identity_verdict": AMBER,
                                   "notes": "candidate course D_Rx unresolved"})
                    continue

                ti = target_identity(c["dose"], lesion, c["zooms"], rx, cfg)
                frac = ti.get("lesion_in_island_frac")
                dist = ti.get("island_centroid_dist_mm")
                # Score is primarily overlap; distance softly penalizes ties.
                score = np.nan
                if frac is not None and np.isfinite(frac):
                    penalty = 1.0
                    if dist is not None and np.isfinite(dist):
                        penalty = 1.0 / (1.0 + dist / max(cfg["target_centroid_max_mm"], 1e-6))
                    score = float(frac * penalty)
                scored.append({"candidate_course": cand_aid, "comparable_grid": True,
                               "score": score, "identity_verdict": ti["identity_verdict"],
                               "lesion_in_island_frac": frac,
                               "island_centroid_dist_mm": dist,
                               "notes": " | ".join(ti.get("identity_notes", []))})

            finite = [x for x in scored if x["comparable_grid"] and np.isfinite(x["score"])]
            finite.sort(key=lambda x: x["score"], reverse=True)
            if len(finite) < 2:
                verdict, reason = AMBER, "candidate courses cannot be cross-scored on a common resolved grid"
                best = runner = None
                margin = np.nan
            else:
                best, runner = finite[0], finite[1]
                margin = (best["score"] / runner["score"]
                          if runner["score"] > 0 else np.inf)
                if (best["candidate_course"] == lesion_aid
                        and best["identity_verdict"] == PASS
                        and margin >= cfg["course_margin_min"]):
                    verdict, reason = PASS, "own course is the uniquely best-supported pairing"
                elif (best["candidate_course"] != lesion_aid
                      and best["identity_verdict"] == PASS
                      and margin >= cfg["course_margin_min"]):
                    verdict, reason = FAIL, "another archived course better explains this baseline lesion"
                else:
                    verdict, reason = AMBER, "course pairing is not uniquely resolved"

            archive_verdict[lesion_aid] = (verdict, reason)
            for sc in scored:
                recs.append({
                    "patient_id": pid,
                    "lesion_archive_id": lesion_aid,
                    **sc,
                    "best_course": None if best is None else best["candidate_course"],
                    "runner_up_course": None if runner is None else runner["candidate_course"],
                    "separation_margin": margin,
                    "pairing_verdict": verdict,
                    "pairing_notes": reason,
                })

    return pd.DataFrame(recs), archive_verdict


# --------------------------------------------------------------------------
# gate composition
# --------------------------------------------------------------------------
def compose(rx_v, scale_v, identity_v, plaus_v, course_v, correspondence_v):
    """A normalized FailurePattern may be instantiated only if every gate passes.

    NOT_EVALUABLE is retained as a separate per-gate state so that missing
    prerequisites do not inflate empirical FAIL counts. It still withholds the
    interpretation.
    """
    verdicts = [rx_v, scale_v, identity_v, plaus_v, course_v, correspondence_v]
    if FAIL in verdicts:
        return FAIL, "observation retained; FailurePattern not instantiable"
    if AMBER in verdicts or NOT_EVALUABLE in verdicts:
        return AMBER, "exploratory/unresolved evidence only; FailurePattern withheld"
    return PASS, "all gates satisfied; FailurePattern instantiable"


def run(manifest_path, config_path, outdir):
    cfg = dict(DEFAULT_CONFIG)
    if config_path:
        cfg.update(json.loads(Path(config_path).read_text()))
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    man = pd.read_csv(manifest_path)
    # Optional field: semicolon-separated verdicts aligned to followup_masks.
    if "followup_correspondence_verdicts" not in man.columns:
        man["followup_correspondence_verdicts"] = ""

    evidence, rx_rows, cache = [], [], {}

    # First pass: course-level G0/G1 and baseline G2/G3.
    for _, row in man.iterrows():
        aid = str(row["archive_id"])
        rc = resolve_rx(row.get("rtplan_dcm"), row.get("clinical_rx_gy"),
                        row.get("clinical_fractions"), cfg)
        rc.archive_id = aid
        rx_rows.append({**asdict(rc), "notes": " | ".join(rc.notes)})

        dose, affine, zooms, nonfin, negs = load_dose(row["rtp_nifti"], cfg)
        sq = scale_qa(dose, row.get("rtdose_dcm"), cfg)
        cache[aid] = {"rx": rc, "dose": dose, "affine": affine, "zooms": zooms,
                      "scale": sq, "nonfin": nonfin, "negs": negs}

        bpath = row.get("baseline_mask")
        bimg = nib.load(str(bpath))
        bmask = np.asanyarray(bimg.dataobj) > 0
        geom_ok = (bmask.shape == dose.shape
                   and np.allclose(bimg.affine, affine, atol=1e-4, rtol=0))
        if geom_ok:
            bst = mask_stats(dose, bmask, zooms)
            bti = target_identity(dose, bmask, zooms, rc.rx_margin_gy, cfg)
            bpl = rx_plausibility(bst, rc.rx_margin_gy, cfg)
        else:
            bst = {k: np.nan for k in ["n_vox", "volume_cc", "dmean", "dmedian", "dmin", "dmax", "d95"]}
            bti = {"identity_verdict": FAIL, "identity_notes": ["baseline geometry mismatch"]}
            bpl = {"plausibility_verdict": FAIL, "plausibility_notes": ["baseline geometry mismatch"]}
        cache[aid].update({"baseline_path": bpath, "baseline_geom": geom_ok,
                           "baseline_stats": bst, "baseline_identity": bti,
                           "baseline_plausibility": bpl})

    # G4 requires the fixed baseline lesion to be scored against every candidate course.
    pairing_df, pairing_by_archive = course_pairing_cross(man, cache, cfg)
    for aid in man["archive_id"].astype(str):
        pairing_by_archive.setdefault(aid, (PASS, "single archived course; G4 not applicable"))

    # Second pass: emit baseline + follow-ups, inheriting baseline G2/G3 and course G4.
    for _, row in man.iterrows():
        aid = str(row["archive_id"])
        c = cache[aid]
        rc, sq = c["rx"], c["scale"]
        g4_v, g4_note = pairing_by_archive[aid]
        base_g2 = c["baseline_identity"]["identity_verdict"]
        base_g3 = c["baseline_plausibility"]["plausibility_verdict"]

        raw_f = row.get("followup_masks")
        fups = "" if (raw_f is None or pd.isna(raw_f)) else str(raw_f)
        fup_paths = [x.strip() for x in fups.split(";") if x.strip()]

        raw_i = row.get("followup_indices")
        idxs = "" if (raw_i is None or pd.isna(raw_i)) else str(raw_i)
        fup_indices = [int(float(x.strip())) for x in idxs.split(";") if x.strip()]
        if not fup_indices:
            fup_indices = []
            for i, p in enumerate(fup_paths):
                m = re.search(r"(?:fu|followup[_-]?)(\\d+)", Path(p).name, re.I)
                fup_indices.append(int(m.group(1)) if m else i + 1)

        raw_c = row.get("followup_correspondence_verdicts")
        corr = "" if (raw_c is None or pd.isna(raw_c)) else str(raw_c)
        corr_vals = [x.strip() for x in corr.split(";")] if corr else []

        raw_d = row.get("followup_days_from_baseline")
        days = "" if (raw_d is None or pd.isna(raw_d)) else str(raw_d)
        day_vals = [float(x.strip()) if x.strip() else np.nan for x in days.split(";")] if days else []

        masks = [("baseline", row["baseline_mask"], PASS, PASS, PASS,
                  "baseline: G5 not applicable", "baseline: temporal gate not applicable")]
        for i, p in enumerate(fup_paths):
            val = corr_vals[i] if i < len(corr_vals) else ""
            spatial_v, spatial_note = _g5_from_manifest(val, cfg)
            d = day_vals[i] if i < len(day_vals) else np.nan
            temporal_v, temporal_note = temporal_order_gate(d, cfg)
            g5_v = combine_g5(spatial_v, temporal_v)
            idx = fup_indices[i] if i < len(fup_indices) else i + 1
            masks.append((f"followup_{idx}", p, spatial_v, temporal_v, g5_v,
                          spatial_note, temporal_note))

        for role, mpath, g5_spatial_v, g5_temporal_v, g5_v, g5_note, temporal_note in masks:
            mimg = nib.load(str(mpath))
            m = np.asanyarray(mimg.dataobj) > 0
            geom_ok = (m.shape == c["dose"].shape
                       and np.allclose(mimg.affine, c["affine"], atol=1e-4, rtol=0))
            st = mask_stats(c["dose"], m, c["zooms"]) if geom_ok else {
                k: np.nan for k in ["n_vox", "volume_cc", "dmean", "dmedian", "dmin", "dmax", "d95"]}

            # G2/G3 are baseline-course gates and are inherited by follow-ups.
            # NOT_EVALUABLE remains distinct from empirical FAIL.
            ti = c["baseline_identity"]
            pl = c["baseline_plausibility"]
            overall, reason = compose(rc.verdict, sq["scale_verdict"],
                                      base_g2, base_g3, g4_v, g5_v)
            notes = (sq["scale_notes"] + ti.get("identity_notes", [])
                     + pl.get("plausibility_notes", [])
                     + [g4_note, g5_note, temporal_note])
            evidence.append({
                "archive_id": aid,
                "patient_id": row.get("patient_id"),
                "course_id": row.get("course_id"),
                "mask_role": role,
                "mask_path": mpath,
                "mask_sha256_16": _file_hash(mpath),
                "geometry_match": geom_ok,
                "nonfinite_dose_voxels": c["nonfin"],
                "negative_dose_voxels": c["negs"],
                "rx_margin_gy": rc.rx_margin_gy,
                "rx_source": rc.rx_source,
                "gate_G0_rx": rc.verdict,
                "gate_G1_scale": sq["scale_verdict"],
                "gate_G2_target_identity": base_g2,
                "gate_G3_rx_plausibility": base_g3,
                "gate_G4_course_pairing": g4_v,
                "gate_G5_spatial_correspondence": g5_spatial_v,
                "gate_G5_temporal_order": g5_temporal_v,
                "gate_G5_lesion_correspondence": g5_v,
                "overall_verdict": overall,
                "overall_reason": reason,
                **{k: v for k, v in st.items()},
                **{
                    (("archive_" + k) if k.startswith("scale_max_pass_") else k): v
                    for k, v in sq.items()
                    if k != "scale_notes"
                },
                **{
                    ("baseline_" + k): v
                    for k, v in ti.items()
                    if k != "identity_notes"
                },
                **{
                    ("baseline_" + k): v
                    for k, v in pl.items()
                    if k != "plausibility_notes"
                },
                "notes": " | ".join([str(x) for x in notes if x]),
            })

    ev = pd.DataFrame(evidence)
    ev.to_csv(outdir / "gate_evidence.csv", index=False)
    pd.DataFrame(rx_rows).to_csv(outdir / "rx_resolution.csv", index=False)
    pairing_df.to_csv(outdir / "course_pairing.csv", index=False)

    gate_cols = ["gate_G0_rx", "gate_G1_scale", "gate_G2_target_identity",
                 "gate_G3_rx_plausibility", "gate_G4_course_pairing",
                 "gate_G5_spatial_correspondence", "gate_G5_temporal_order",
                 "gate_G5_lesion_correspondence", "overall_verdict"]
    follow = ev[ev.mask_role != "baseline"]

    # Archive-level final status is the honest independent treatment-course unit:
    # five gates are course/baseline level and only G5 varies by observation.
    archive_overall = (
        ev.groupby("archive_id")["overall_verdict"]
          .apply(lambda s: PASS if (s == PASS).all()
                 else (FAIL if (s == FAIL).any() else AMBER))
          .to_dict()
    )

    summary = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "git_commit": _git_commit(),
        "config": cfg,
        "config_file_role": "authoritative_postdiagnostic_clean_gate_set",
        "config_sha256_16": _config_hash(cfg),
        "n_archives": int(man.shape[0]),
        "n_patients": int(man["patient_id"].nunique()) if "patient_id" in man else None,
        "n_masks": int(ev.shape[0]),
        "n_baseline_masks": int((ev.mask_role == "baseline").sum()),
        "n_followup_masks": int((ev.mask_role != "baseline").sum()),
        "gate_counts": {g: ev[g].value_counts().to_dict() for g in gate_cols},
        "n_followup_failure_patterns_instantiable": int((follow.overall_verdict == PASS).sum()),
        "archive_overall_verdicts": archive_overall,
        "archive_overall_counts": pd.Series(archive_overall).value_counts().to_dict(),
        "note": (
            "Five gates are evaluated at archive/course level and inherited by "
            "follow-ups; G5 is observation-level. Therefore report both archive-level "
            "and observation-level results."
        )
    }
    (outdir / "audit_summary.json").write_text(json.dumps(summary, indent=2))

    print(json.dumps(summary["gate_counts"], indent=2))
    print(f"\nFollow-up FailurePatterns instantiable: "
          f"{summary['n_followup_failure_patterns_instantiable']} / {summary['n_followup_masks']}")
    return summary

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--config", default=None)
    ap.add_argument("--out", default="results")
    a = ap.parse_args(argv)
    run(a.manifest, a.config, a.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
