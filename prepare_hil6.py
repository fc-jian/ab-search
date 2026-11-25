#!/usr/bin/env python3
"""Prepare exact assay lookups and out-of-fold probabilities for AVIDa-hIL6."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sqlite3
import tempfile
from contextlib import closing
from pathlib import Path

import numpy as np
import pandas as pd

from prepare_patents import digest, download


RELEASE = "cd783d573e015f83364e42287810ad597850aeed"
DATASET = "https://huggingface.co/datasets/COGNANO/AVIDa-hIL6"
FILES = {
    "AVIDa-hIL6.csv": ("avida_hil6.csv", "917d90f111d759edc3440f577458d42ef89ff0315b186d6c266c120143b20cde"),
    "antigen_sequences.csv": ("avida_antigens.csv", "8f023581932e626ba1e2da4f7f861b998b8fcb7516aedb785d831f8c36984e68"),
}
ASSAY = "AVIDa-hIL6 published phage-display enrichment label"


def panel_predictions(profile: np.ndarray, calibration: np.ndarray):
    """Rows are known VHHs. Both matrices exclude every test antigen column."""
    measured = (profile >= 0).sum(axis=1)
    positive = (profile == 1).sum(axis=1)
    rates = positive / np.maximum(measured, 1)
    # Bin each prior response profile by hit rate and number of observed assays.
    bins = np.searchsorted([0, 0.25, 0.5, 0.75, 1], rates)
    bins[rates == 1] = 5
    bins += 6 * (measured >= 10)
    n = (calibration >= 0).sum(axis=1) * (measured >= 3)
    k = (calibration == 1).sum(axis=1) * (measured >= 3)
    # Exclude the query VHH from its calibration bin, including all its replicates.
    count = np.bincount(bins, weights=n)[bins] - n
    success = np.bincount(bins, weights=k)[bins] - k
    usable = (measured >= 3) & (count >= 30)
    probs = np.full(len(measured), np.nan)
    probs[usable] = (success[usable] + 1) / (count[usable] + 2)
    # The baseline also excludes the evaluated VHH and uses calibration constructs only.
    baseline = (k.sum() - k + 1) / (n.sum() - n + 2)
    return probs, np.where(usable, count, 0).astype(int), baseline


def build_panel(reference_dir: Path):
    inputs = []
    for filename, expected in FILES.values():
        path = reference_dir / filename
        if digest(path) != expected:
            raise ValueError(f"Unexpected AVIDa-hIL6 release file: {path}")
        inputs.append(path)
    measurements = pd.read_csv(inputs[0], usecols=["Ag_label", "VHH_sequence", "label"])
    antigens = pd.read_csv(inputs[1]).set_index("Ag_label")["Ag_sequence"]
    # The key is the complete released antigen construct, including its tag.
    matrix = measurements.pivot(index="VHH_sequence", columns="Ag_label", values="label")
    names = sorted(antigens.index, key=lambda k: hashlib.sha256(antigens[k].encode()).hexdigest())
    observations = matrix[names].fillna(-1).to_numpy(dtype=np.int8)
    folds = np.arange(len(names)) % 5
    output = reference_dir / "avida_hil6.sqlite"
    fd, filename = tempfile.mkstemp(dir=output.parent, suffix=".sqlite.tmp")
    os.close(fd)
    staging = Path(filename)
    reports = []
    try:
        with closing(sqlite3.connect(staging)) as db, db:
            db.executescript("""
                CREATE TABLE antigens(id INTEGER PRIMARY KEY, name TEXT, sequence TEXT UNIQUE);
                CREATE TABLE antibodies(id INTEGER PRIMARY KEY, sequence TEXT);
                CREATE TABLE assays(antigen INTEGER, antibody INTEGER, outcome INTEGER,
                                    probability REAL, calibration_n INTEGER,
                                    PRIMARY KEY(antigen,antibody)) WITHOUT ROWID;
                CREATE TABLE validation(name TEXT PRIMARY KEY, fold INTEGER, measured INTEGER, evaluated INTEGER,
                                        positives INTEGER, brier REAL, baseline_brier REAL);
                CREATE TABLE release(revision TEXT PRIMARY KEY);
            """)
            db.execute("INSERT INTO release VALUES(?)", (RELEASE,))
            db.executemany("INSERT INTO antibodies VALUES(?,?)", enumerate(matrix.index))
            db.executemany("INSERT INTO antigens VALUES(?,?,?)", ((i, s, antigens[s]) for i, s in enumerate(names)))
            for fold in range(5):
                fit = (folds == (fold + 1) % 5)
                test = (folds == fold)
                probs, counts, baseline = panel_predictions(observations[:, ~(fit | test)], observations[:, fit])
                for i in np.flatnonzero(test):
                    values = observations[:, i]
                    observed = values >= 0
                    evaluated = observed & np.isfinite(probs)
                    residual = values[evaluated] - probs[evaluated]
                    baseline_residual = values[evaluated] - baseline[evaluated]
                    report = (names[i], fold, int(observed.sum()), int(evaluated.sum()),
                              int((values[evaluated] == 1).sum()), float(np.mean(residual**2)),
                              float(np.mean(baseline_residual**2)))
                    db.execute("INSERT INTO validation VALUES(?,?,?,?,?,?,?)", report)
                    reports.append(report)
                    db.executemany("INSERT INTO assays VALUES(?,?,?,?,?)", (
                        (int(i), int(j), int(values[j]), float(probs[j]) if evaluated[j] else None, int(counts[j]))
                        for j in np.flatnonzero(observed)
                    ))
        staging.replace(output)
    finally:
        staging.unlink(missing_ok=True)
    evaluated = sum(r[3] for r in reports)
    summary = {"observed": len(measurements), "evaluated": evaluated, "constructs": len(names),
               "vhh": len(matrix), "measured_positive": int(measurements["label"].sum()),
               "evaluated_positive": sum(r[4] for r in reports),
               "brier": sum(r[3] * r[5] for r in reports) / evaluated,
               "baseline_brier": sum(r[3] * r[6] for r in reports) / evaluated}
    pd.DataFrame(reports, columns=["construct", "test_fold", "measured", "evaluated", "evaluated_positive",
                                  "brier", "baseline_brier"]).sort_values("construct").to_csv(
        reference_dir / "avida_validation.csv", index=False)
    print(json.dumps(summary, indent=2))


def assay_evidence(db, query, other_rows):
    antigen = db.execute("SELECT id,name FROM antigens WHERE sequence=?", (query,)).fetchone()
    if not antigen:
        return []
    others = {r["heavy"] for r in other_rows if not r["light"]}
    rows = []
    records = db.execute("""SELECT b.id,b.sequence,s.outcome,s.probability,s.calibration_n
                            FROM assays s JOIN antibodies b ON b.id=s.antibody WHERE s.antigen=?""", (antigen[0],))
    for antibody_id, heavy, bound, prob, count in records:
        if not bound and heavy not in others:
            continue
        rows.append(dict(source="AVIDa-hIL6", record_id=f"{antigen[1]}:{antibody_id}", document_id="AVIDa-hIL6",
                         antibody="", heavy=heavy, light="", target=antigen[1], organism="Homo sapiens",
                         match="exact", alignment_identity=1, method=ASSAY, origin="assay_construct",
                         chain_basis="assayed_VHH_construct", antibody_format="VHH",
                         tier="exact_assay" if bound else "review", reason="" if bound else "observed_negative",
                         binding_probability=prob, calibration_n=count, calibration_assay=ASSAY,
                         probability_status="panel_calibrated" if prob is not None else "insufficient_profile",
                         citation=f"{DATASET}/tree/{RELEASE}"))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-dir", type=Path, default=Path("reference"))
    parser.add_argument("--download", action="store_true")
    args = parser.parse_args()
    if args.download:
        for remote, (local, _) in FILES.items():
            download(f"{DATASET}/resolve/{RELEASE}/{remote}", args.reference_dir / local)
    build_panel(args.reference_dir)


if __name__ == "__main__":
    main()
