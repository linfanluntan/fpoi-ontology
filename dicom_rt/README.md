# DICOM-RT interoperability and object-ablation experiment (C8)

The five study packages are not redistributed because their provenance and redistribution status were not independently established.

## Intact-package screen

The packaged DICOM series show one natural incomplete control and four packages containing all required object classes:

- Center 1: CT + RTDOSE only — incomplete.
- Centers 2–5: CT + RTSTRUCT + RTPLAN + RTDOSE — required object set complete; the manuscript's prior reference-chain audit classified these four as complete.

## Systematic object-ablation controls — completed

For each of the four intact packages (Centers 2–5), one required DICOM object class was removed at a time: CT, RTSTRUCT, RTPLAN, or RTDOSE. This yields **16 deterministic negative controls**.

**Result: 16/16 ablations are rejected by the required-object completeness criterion.** Each failure is surfaced before reference-chain evaluation because at least one mandatory object class is absent. Center 1 remains the natural incomplete-package control.

This experiment tests failure-state surfacing, not diagnostic performance. The repository contains the ablation procedure/result table but no study DICOM objects.
