"""
mesa.py — MESA multimodal analysis.
Delegates all model classes (MESA_modality, MESA, MESA_CV) directly to
the installed mesa package. Only orchestration logic lives here.
"""

import os
import sys
import pickle
import time

import numpy as np
import pandas as pd
from sklearn.base import clone as sk_clone
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import LeaveOneOut
from sklearn.svm import SVC

try:
    from xgboost import XGBClassifier
    _HAS_XGB = True
except ImportError:
    _HAS_XGB = False

# ── import from installed mesa package ───────────────────────────────────────
from mesa import MESA_modality, MESA, MESA_CV

# modality_performance stays as our own implementation (uses feature_preprocessing)
from analysis.modality_performance import modality_performance


def _disp(msg):
    print(f"@{time.asctime()}\t{msg}", file=sys.stderr)


def _load_matrices(paths):
    """Load feature matrices with positional, unique in-memory feature labels.

    Targeted WPS matrices can legitimately retain overlapping windows with the
    same ``chr:start-end`` identifier.  Those labels are useful metadata in
    TSVs but become DataFrame column names after transpose, where recent
    scikit-learn/Narwhals versions require uniqueness.  Positional labels keep
    every row/value and avoid changing the scientific feature set.
    """
    matrices = []
    for filepath in paths:
        matrix = pd.read_csv(filepath, sep="\t", index_col=0).T
        matrix.columns = pd.RangeIndex(matrix.shape[1])
        matrices.append(matrix)
    return matrices


# ── classifier registry ───────────────────────────────────────────────────────

def _build_clf_dist():
    d = {
        1: RandomForestClassifier(random_state=0, n_jobs=-1),
        2: LogisticRegression(random_state=0, max_iter=1000, n_jobs=-1),
        3: SVC(random_state=0, probability=True),
    }
    if _HAS_XGB:
        d[4] = XGBClassifier(random_state=0, n_jobs=-1, eval_metric="logloss")
    return d


CLF_DIST = _build_clf_dist()


# ── clf resolution from performance table ─────────────────────────────────────

def _resolve_clf_per_modality(args, n_modalities, performance=None, perf_tsv=None):
    """Determine best classifier per modality from performance table or args.clf."""
    import re
    names = list(args.modality)

    def _from_df(df):
        clfs = []
        for name in names:
            if name in df.index:
                raw    = df.loc[name, "best_classifier, idx"]
                digits = re.findall(r'\d+', str(raw))
                pos    = int(digits[-1]) if digits else 0
                clf_id = args.clf[pos] if pos < len(args.clf) else args.clf[0]
                clfs.append(sk_clone(CLF_DIST.get(clf_id, CLF_DIST[args.clf[0]])))
            else:
                clfs.append(sk_clone(CLF_DIST[args.clf[0]]))
        return clfs

    if performance is not None:
        return _from_df(performance[0])
    if perf_tsv and os.path.exists(perf_tsv):
        return _from_df(pd.read_csv(perf_tsv, sep="\t", index_col=0))
    if len(args.clf) == n_modalities:
        return [sk_clone(CLF_DIST[c]) for c in args.clf]
    return [sk_clone(CLF_DIST[args.clf[0]])] * n_modalities


# ── single-modality LOOCV ─────────────────────────────────────────────────────

def _loocv_single(name, X, clf, top_n):
    """
    Run LOOCV for one modality using MESA_CV with LeaveOneOut.
    Exactly matches original mesa_loocv.py behavior:
      - MESA_CV clones MESA_modality each fold
      - each fold: transform_predict_proba (pipeline.transform → predictor_.predict_proba)
      - Parallel execution across folds
    """
    _disp(f"LOOCV [{name}] n={X.shape[0]} feat={X.shape[1]}")
    mod = MESA_modality(classifier=clf, top_n=top_n)
    cv  = MESA_CV(modality=mod, cv=LeaveOneOut())
    cv.fit(X, y=None)   # y passed via fit signature below
    return cv


def _run_single_loocv(name, X, clf, y, top_n):
    """
    Run LOOCV for one modality. Returns (proba, y_true) arrays in sample order.
    LeaveOneOut is sequential so cv_result[i] = fold for sample i.
    """
    _disp(f"LOOCV [{name}] n={X.shape[0]} feat={X.shape[1]}")
    mod = MESA_modality(classifier=clf, top_n=top_n)
    cv  = MESA_CV(modality=mod, cv=LeaveOneOut())
    cv.fit(X, y)
    proba  = np.array([cv._fold_prediction(f[0])[0] for f in cv.cv_result])
    y_true = np.array([f[1][0] for f in cv.cv_result])
    return proba, y_true


# ── multimodal LOOCV ──────────────────────────────────────────────────────────

def _run_multimodal_loocv(X_list, clfs, y, top_n):
    """
    Manual LOO loop for multimodal MESA stacking.
    Each fold: fresh MESA object fit on train, predict on test.
    MESA.fit uses internal OOF CV for meta-learner (no leakage).
    Matches original mesa_loocv.py _run_multimodal_loocv exactly.
    """
    n = X_list[0].shape[0]
    _disp(f"LOOCV [Multimodal] n={n} modalities={len(X_list)}")

    proba_arr  = np.empty(n)
    y_true_arr = np.empty(n)

    for i, (tr, te) in enumerate(LeaveOneOut().split(X_list[0])):
        base_modalities = [
            MESA_modality(classifier=sk_clone(clf), top_n=top_n)
            for clf in clfs
        ]
        mesa_model = MESA(
            modalities=base_modalities,
            meta_estimator=LogisticRegression(random_state=0, max_iter=1000),
        )
        mesa_model.fit([X.iloc[tr] for X in X_list], y[tr])
        proba_arr[te[0]]  = mesa_model.predict_proba(
            [X.iloc[te] for X in X_list]
        )[0, 1]
        y_true_arr[te[0]] = y[te[0]]

        if (i + 1) % 10 == 0 or (i + 1) == n:
            _disp(f"  [Multimodal] fold {i+1}/{n}")

    return proba_arr, y_true_arr


# ── public API ────────────────────────────────────────────────────────────────

def run_modality_performance(args):
    """Evaluate per-modality classification performance and save TSV."""
    os.makedirs(args.output_dir, exist_ok=True)
    matrices   = _load_matrices(args.infile)
    sample_ids = matrices[0].index.tolist()
    label_df   = pd.read_table(args.label, header=None, index_col=0)
    label_df.columns = ["label"]
    y    = label_df.reindex(sample_ids)["label"].values.reshape(-1)
    clfs = [sk_clone(CLF_DIST[c]) for c in args.clf]

    perf, _, _ = modality_performance(
        args.modality, matrices, y, clfs,
        feature_size=args.size,
        subset=args.subset,
        repeat=args.repeat,
        n_jobs=args.cores,
    )
    out = os.path.join(args.output_dir, "modality_performance.tsv")
    perf.to_csv(out, sep="\t")
    _disp(f"Performance saved → {out}")
    return perf, args.modality, clfs


def run_mesa_model(args, performance=None):
    """Train full MESA stacking model and save pickle."""
    os.makedirs(args.output_dir, exist_ok=True)
    matrices   = _load_matrices(args.infile)
    sample_ids = matrices[0].index.tolist()
    label_df   = pd.read_table(args.label, header=None, index_col=0)
    label_df.columns = ["label"]
    y    = label_df.reindex(sample_ids)["label"].values
    clfs = _resolve_clf_per_modality(
        args, len(matrices), performance=performance,
        perf_tsv=getattr(args, "perf_tsv", None),
    )
    mods  = [MESA_modality(classifier=clf, top_n=args.size) for clf in clfs]
    model = MESA(modalities=mods, meta_estimator=LogisticRegression(
        random_state=0, max_iter=1000
    ))
    model.fit(matrices, y)
    out = os.path.join(args.output_dir, "MESA_model.pkl")
    pickle.dump(model, open(out, "wb"))
    _disp(f"MESA model saved → {out}")
    return model


def run_mesa_loocv(args, performance=None):
    """
    Run LOOCV for all single modalities and multimodal MESA.
    Single-modality: MESA_CV(LeaveOneOut) — matches original mesa_loocv.py.
    Multimodal: manual LOO loop with MESA.fit — matches original mesa_loocv.py.
    """
    os.makedirs(args.output_dir, exist_ok=True)
    matrices   = _load_matrices(args.infile)
    sample_ids = matrices[0].index.tolist()
    # align matrices to same sample order
    matrices   = [m.reindex(sample_ids) for m in matrices]

    label_df = pd.read_table(args.label, header=None, index_col=0)
    label_df.columns = ["label"]
    y = label_df.reindex(sample_ids)["label"].values.astype(int)
    if np.any(pd.isnull(y)):
        missing = [s for s in sample_ids if s not in label_df.index]
        sys.exit(
            f"[mesa] ERROR: {len(missing)} sample(s) not in label.tsv: {missing[:5]}"
        )

    for name, mat, fp in zip(args.modality, matrices, args.infile):
        _disp(f"[mesa] {name}: {mat.shape[0]} samples × {mat.shape[1]} features — {fp}")

    clfs = _resolve_clf_per_modality(
        args, len(matrices), performance=performance,
        perf_tsv=getattr(args, "perf_tsv", None),
    )

    # ── single-modality LOOCV ─────────────────────────────────────────────────
    pred_dict  = {}
    y_true_ref = None

    for name, X, clf in zip(args.modality, matrices, clfs):
        proba, y_true = _run_single_loocv(name, X, clf, y, top_n=args.size)
        if y_true_ref is None:
            y_true_ref = y_true
        pred_dict[name] = proba
        auc = roc_auc_score(y_true, proba)
        _disp(f"  [{name}] AUC={auc:.4f}")

    # ── multimodal LOOCV ──────────────────────────────────────────────────────
    proba_multi, y_true_multi = _run_multimodal_loocv(
        matrices, clfs, y, top_n=args.size
    )
    auc_multi = roc_auc_score(y_true_multi, proba_multi)
    _disp(f"  [Multimodal] AUC={auc_multi:.4f}")
    pred_dict["Multimodal"] = proba_multi

    # ── assemble output ───────────────────────────────────────────────────────
    pred_df = pd.DataFrame(pred_dict, index=sample_ids)
    pred_df.insert(0, "y_true", y_true_ref.astype(int))
    pred_df.index.name = "sample_id"

    out_path = os.path.join(args.output_dir, "loocv_predictions.tsv")
    pred_df.to_csv(out_path, sep="\t")
    _disp(f"Predictions saved → {out_path}")
    return pred_df
