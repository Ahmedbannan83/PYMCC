"""CyberKnife absolute-dose beam-quality helper.

Implements the workflow in Accuray Physics Essentials Guide, Chapter 2,
Absolute LINAC Dose Calibration, pp. 2-156--2-158. It is a calculation aid
for a qualified medical physicist and is not an independent clinical protocol.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
import re
from typing import Optional

import numpy as np
import pandas as pd


# BJR Supplement 25, Tables 5.2.1--5.5.1: PDD (%) at SSD=100 cm.
# Only the field sizes/depths needed for the Accuray 60-mm conversion are kept.
_BJR25 = pd.DataFrame([
    (4.0, 6.0, 59.7, 30.9), (4.0, 7.0, 60.7, 31.8), (4.0, 10.0, 63.0, 34.1),
    (5.0, 6.0, 62.0, 33.1), (5.0, 7.0, 62.9, 34.0), (5.0, 10.0, 65.0, 36.2),
    (6.0, 6.0, 64.9, 36.4), (6.0, 7.0, 65.7, 37.2), (6.0, 10.0, 67.5, 39.3),
    (8.0, 6.0, 69.1, 40.9), (8.0, 7.0, 69.7, 41.6), (8.0, 10.0, 71.0, 43.4),
], columns=["energy_mv", "field_cm", "pdd100", "pdd200"])


@dataclass(frozen=True)
class BeamQualityResult:
    measured_pdd100: float
    measured_pdd200: float
    measurement_ssd_mm: float
    cone_diameter_mm_at_800_sad: float
    equivalent_square_mm_at_1000_ssd: float
    bjr25_inferred_energy_mv: float
    equivalent_pdd100_10x10: float
    equivalent_pdd200_10x10: float
    tpr20_10: float
    ionization_chamber: Optional[str]
    kq_trs398_table16: Optional[float]
    warnings: tuple[str, ...]


def _field_interpolated(energy_mv: float, field_cm: float, column: str) -> float:
    table = _BJR25.loc[_BJR25.energy_mv == energy_mv].sort_values("field_cm")
    return float(np.interp(field_cm, table.field_cm, table[column]))


def _infer_energy(measured_pdd100: float, equivalent_square_cm: float) -> float:
    energies = np.array(sorted(_BJR25.energy_mv.unique()), dtype=float)
    reference = np.array([
        _field_interpolated(energy, equivalent_square_cm, "pdd100")
        for energy in energies
    ])
    if not reference[0] <= measured_pdd100 <= reference[-1]:
        raise ValueError(
            f"PDD100={measured_pdd100:.2f}% is outside the BJR-25 4-8 MV "
            f"reference range ({reference[0]:.2f}-{reference[-1]:.2f}%) for "
            f"a {equivalent_square_cm:.2f} cm square field."
        )
    return float(np.interp(measured_pdd100, reference, energies))


def _interpolate_energy(energy_mv: float, field_cm: float, column: str) -> float:
    energies = np.array(sorted(_BJR25.energy_mv.unique()), dtype=float)
    values = np.array([_field_interpolated(e, field_cm, column) for e in energies])
    return float(np.interp(energy_mv, energies, values))


def _read_ssd_mm(mcc_path: str | Path) -> Optional[float]:
    pattern = re.compile(r"^\s*SSD\s*=\s*([0-9.]+)")
    with open(mcc_path, encoding="utf-8") as source:
        for line in source:
            match = pattern.match(line)
            if match:
                return float(match.group(1))
    return None


def _kq_from_table16(chamber: str, tpr20_10: float, table16_csv: str | Path) -> tuple[str, float]:
    data = pd.read_csv(table16_csv)
    models = data.ionization_chamber.drop_duplicates().tolist()
    normalized = chamber.casefold().replace(" ", "")
    candidates = [model for model in models if model.casefold().replace(" ", "") in normalized
                  or normalized in model.casefold().replace(" ", "")]
    if len(candidates) != 1:
        raise ValueError(
            f"Ionization chamber '{chamber}' was not uniquely matched in TRS-398 Table 16. "
            f"Use one of: {', '.join(models)}"
        )
    selected = data.loc[data.ionization_chamber == candidates[0]].sort_values("TPR20_10")
    if not selected.TPR20_10.min() <= tpr20_10 <= selected.TPR20_10.max():
        raise ValueError("TPR20,10 is outside the Table 16 interpolation range 0.56-0.82.")
    return candidates[0], float(np.interp(tpr20_10, selected.TPR20_10, selected.kQ))


def calculate_beam_quality(
    pdd100: float,
    pdd200: float,
    ionization_chamber: Optional[str] = None,
    measurement_ssd_mm: float = 1000.0,
    fixed_cone_mm_at_800_sad: float = 60.0,
    table16_csv: str | Path | None = None,
) -> BeamQualityResult:
    """Convert a 60-mm-cone PDD measurement to TPR20,10 and optional TRS-398 kQ.

    PDD100 and PDD200 are values (%) at 100 and 200 mm, respectively, normalized
    to dmax. The prescribed Accuray measurement setup is SSD=1000 mm.
    """
    if not (0 < pdd200 < pdd100 < 100):
        raise ValueError("Expected 0 < PDD200 < PDD100 < 100 (percent).")
    warnings: list[str] = []
    if measurement_ssd_mm != 1000.0:
        warnings.append(
            "The Accuray procedure specifies SSD=1000 mm; this result is not a "
            "protocol-compliant absolute-calibration beam-quality determination."
        )
    equivalent_square_mm = 1.125 * fixed_cone_mm_at_800_sad
    equivalent_square_cm = equivalent_square_mm / 10.0
    energy = _infer_energy(pdd100, equivalent_square_cm)

    # BJR-25 field-size correction, preserving the measured beam quality.
    delta_100 = _interpolate_energy(energy, 10.0, "pdd100") - _interpolate_energy(energy, equivalent_square_cm, "pdd100")
    delta_200 = _interpolate_energy(energy, 10.0, "pdd200") - _interpolate_energy(energy, equivalent_square_cm, "pdd200")
    pdd100_10 = pdd100 + delta_100
    pdd200_10 = pdd200 + delta_200

    # TG-51 / existing pymcc Q-index relation for the 10x10 equivalent beam.
    tpr20_10 = 1.2661 * (pdd200_10 / pdd100_10) - 0.0595

    matched_chamber: Optional[str] = None
    kq: Optional[float] = None
    if ionization_chamber:
        if table16_csv is None:
            table16_csv = Path(__file__).with_name("data") / "trs398_table16_kq.csv"
        matched_chamber, kq = _kq_from_table16(ionization_chamber, tpr20_10, table16_csv)

    return BeamQualityResult(
        measured_pdd100=float(pdd100), measured_pdd200=float(pdd200),
        measurement_ssd_mm=float(measurement_ssd_mm),
        cone_diameter_mm_at_800_sad=float(fixed_cone_mm_at_800_sad),
        equivalent_square_mm_at_1000_ssd=equivalent_square_mm,
        bjr25_inferred_energy_mv=energy,
        equivalent_pdd100_10x10=pdd100_10,
        equivalent_pdd200_10x10=pdd200_10,
        tpr20_10=tpr20_10, ionization_chamber=matched_chamber,
        kq_trs398_table16=kq, warnings=tuple(warnings),
    )


def calculate_from_mcc(
    mcc_path: str | Path,
    ionization_chamber: Optional[str] = None,
    table16_csv: str | Path | None = None,
) -> BeamQualityResult:
    """Read the first photon PDD in an MCC file and calculate beam quality."""
    from pymcc.readmcc import read_file

    pdds = [curve for curve in read_file(str(mcc_path)) if curve.curve_type == "PDD" and curve.modality == "X"]
    if len(pdds) != 1:
        raise ValueError(f"Expected exactly one photon PDD curve; found {len(pdds)}.")
    pdd = pdds[0]
    return calculate_beam_quality(
        pdd100=float(pdd.dose_100()), pdd200=float(pdd.dose_200()),
        ionization_chamber=ionization_chamber,
        measurement_ssd_mm=_read_ssd_mm(mcc_path) or 1000.0,
        table16_csv=table16_csv,
    )


def result_as_dict(result: BeamQualityResult) -> dict:
    """Return a serializable dictionary for scripts, QATrack+ or JSON output."""
    return asdict(result)
