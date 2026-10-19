
# -*- coding: utf-8 -*-
"""
Run PaCMAP and UMATO on all *_TimeDR_forDR.mat files produced by the MATLAB export script.

Each input .mat is expected to contain (per session):
    X                   : [nPoints x nNeurons]  (per-neuron z-scored, DR-ready)
    X_raw               : [nPoints x nNeurons]  (raw window-summed; unused here)
    labels              : [nPoints]  Visual_Stim for each (trial, window) row
    trial_ids           : [nPoints]  trial index (1..nTrials) for each row
    time_window_idx     : [nPoints]  window index (1..nWin) for each row
    time_window_center  : [nPoints]  time (s) of window center for each row
    Vstim_all           : [nTrials]  Visual_Stim per trial

    DBname, experimentID, behavCond, prefOriField, mu_cols_out, sig_cols_out, etc.

For each file, this script will:
    - load X and metadata
    - compute a PaCMAP embedding (n_components dims)
    - compute a UMATO embedding (n_components dims)
    - save results in <basename>_PaCMAP_UMATO.mat, containing:
        Y_pacmap, Y_umato
        labels, trial_ids, time_window_idx, time_window_center
        Vstim_all, stim_values
        DBname, experimentID, behavCond
        n_components, hub_num_used
    Plotting will be done back in MATLAB.
"""

import os
import glob
import sys

import numpy as np
import scipy.io as sio
import h5py

import pacmap      # YingfanWang/PaCMAP
import umato       # hyungkwonko/umato


# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------

# Root directory holding D*/*_TimeDR_forDR.mat (the output of the 01 pull step).
# Pass it as the first command-line argument, or set it here.
ROOT_DIR = None

N_COMPONENTS = 3      # dimensionality of embedding; MATLAB can choose which dims to plot
DEFAULT_HUB_NUM = 50  # UMATO: number of hub points (will be clipped to n_samples)


# ---------------------------------------------------------------------
# Helpers: loading MATLAB .mat (v7 and v7.3)
# ---------------------------------------------------------------------

def load_mat_any_version(path):
    """
    Load a MATLAB .mat file, handling v7 and v7.3 (HDF5) formats.

    Returns a dict of numpy arrays. For v7.3, this uses h5py.
    """
    try:
        # SciPy loadmat works for v7 (and older), but not for v7.3
        data = sio.loadmat(path)
        # Strip MATLAB's __header__/__version__/__globals__ keys
        return {k: v for k, v in data.items() if not k.startswith("__")}
    except NotImplementedError:
        # Probably a v7.3 (HDF5) file
        out = {}
        with h5py.File(path, "r") as f:
            for k in f.keys():
                out[k] = f[k][()]
        return out


def ensure_2d_float(X):
    """
    Ensure X is a 2D float64 numpy array [n_samples x n_features].

    The MATLAB export should already save X as [nPoints x nNeurons].
    We keep a small safety check in case orientation is swapped.
    """
    X = np.array(X)
    X=X.T
    if X.ndim != 2:
        raise ValueError(f"Expected 2D array for X, got shape {X.shape}")
    return X.astype(np.float64, copy=False)


def ensure_1d(arr):
    """Ensure array is 1D."""
    return np.array(arr).ravel()


def matlab_str_to_py(x):
    """Convert a MATLAB char array or cell to a Python string, if present."""
    if x is None:
        return None
    arr = np.array(x)
    # MATLAB char arrays often come as 2D arrays of single chars
    if arr.dtype.kind in ("U", "S"):
        return "".join(arr.ravel().tolist())
    # MATLAB cellstr via loadmat can come as object arrays
    if arr.dtype == object:
        # try to join first cell as string
        try:
            return str(arr.squeeze().tolist())
        except Exception:
            return str(arr)
    if arr.dtype.kind in ("i", "f"):
        return str(arr.item())
    try:
        return str(arr.squeeze().tolist())
    except Exception:
        return str(arr)


# ---------------------------------------------------------------------
# Core DR for one file
# ---------------------------------------------------------------------

def run_dr_for_file(mat_path, n_components=N_COMPONENTS, hub_num=DEFAULT_HUB_NUM, seed=0):
    print(f"\nProcessing {mat_path}  (seed={seed})")
    data = load_mat_any_version(mat_path)

    # Basic checks
    required_vars = ["X", "labels", "trial_ids", "time_window_idx", "time_window_center"]
    missing = [v for v in required_vars if v not in data]
    if missing:
        raise KeyError(f"{mat_path} is missing variables: {missing}")

    # Core data
    X = ensure_2d_float(data["X"])                    # [n_points x n_neurons]
    labels = ensure_1d(data["labels"])
    trial_ids = ensure_1d(data["trial_ids"])
    time_window_idx = ensure_1d(data["time_window_idx"])
    time_window_center = ensure_1d(data["time_window_center"])

    n_samples, n_features = X.shape
    print(f"  X shape: {X.shape}, labels: {labels.shape}")

    # Sanity checks on lengths
    for name, arr in [
        ("labels", labels),
        ("trial_ids", trial_ids),
        ("time_window_idx", time_window_idx),
        ("time_window_center", time_window_center)
    ]:
        if arr.shape[0] != n_samples:
            raise ValueError(
                f"{mat_path}: {name} length {arr.shape[0]} "
                f"does not match X rows {n_samples}"
            )

    # Vstim_all is per trial; optional but useful for MATLAB
    Vstim_all = data.get("Vstim_all", None)
    if Vstim_all is not None:
        Vstim_all = ensure_1d(Vstim_all)

    stim_values = np.unique(labels)

    # Metadata
    DBname       = matlab_str_to_py(data.get("DBname", None))
    experimentID = matlab_str_to_py(data.get("experimentID", None))
    behavCond    = matlab_str_to_py(data.get("behavCond", None))

    dirname = os.path.dirname(mat_path)
    base = os.path.splitext(os.path.basename(mat_path))[0]
    # The MATLAB files end with "_TimeDR_forDR.mat"; drop that suffix in output names
    if base.endswith("_TimeDR_forDR"):
        base_out = base[:-len("_TimeDR_forDR")]
    else:
        base_out = base

    # ---------------- PaCMAP ----------------
    pac = pacmap.PaCMAP(
        n_components=n_components,
        n_neighbors=10,
        MN_ratio=0.5,
        FP_ratio=2.0,
        random_state=seed,
    )
    Y_pac = pac.fit_transform(X, init="pca")   # [n_samples x n_components]

    # ---------------- UMATO ----------------
    hub_num_eff = int(min(hub_num, n_samples))
    uma = umato.UMATO(hub_num=hub_num_eff, random_state=seed)
    Y_uma = uma.fit_transform(X)              # [n_samples x n_components]

    # ---------------- Save embeddings as .mat ----------------
    out_mat = os.path.join(dirname, base_out + "_PaCMAP_UMATO.mat")

    mdict = {
        "Y_pacmap": Y_pac,
        "Y_umato":  Y_uma,
        "labels":   labels,
        "trial_ids": trial_ids,
        "time_window_idx": time_window_idx,
        "time_window_center": time_window_center,
        "stim_values": stim_values,
        "n_components": np.array([[n_components]]),
        "hub_num_used": np.array([[hub_num_eff]]),
    }

    if Vstim_all is not None:
        mdict["Vstim_all"] = Vstim_all

    if DBname is not None:
        mdict["DBname"] = np.array([[DBname]], dtype=object)
    if experimentID is not None:
        mdict["experimentID"] = np.array([[experimentID]], dtype=object)
    if behavCond is not None:
        mdict["behavCond"] = np.array([[behavCond]], dtype=object)

    sio.savemat(out_mat, mdict)
    print(f"  Saved embeddings to {out_mat}")


# ---------------------------------------------------------------------
# Main: run over all *_TimeDR_forDR.mat files
# ---------------------------------------------------------------------

def main(root_dir=None, seed=0):
    if root_dir is None:
        root_dir = ROOT_DIR
    if root_dir is None:
        print("Usage: python 02_make_pacmap_umato.py <root_dir> [seed]")
        return

    pattern = os.path.join(root_dir, "D*", "*_TimeDR_forDR.mat")
    files = sorted(glob.glob(pattern))

    if not files:
        print(f"No *_TimeDR_forDR.mat files found under {root_dir}")
        return

    print(f"Found {len(files)} files under {root_dir}  (seed={seed})")
    for f in files:
        try:
            run_dr_for_file(f, seed=seed)
        except Exception as e:
            print(f"ERROR processing {f}: {e}")


if __name__ == "__main__":
    # Usage: python 02_make_pacmap_umato.py <root_dir> [seed]
    #   root_dir : folder containing D*/*_TimeDR_forDR.mat
    #   seed     : PaCMAP/UMATO random_state (default 0; main figure = 0)
    root = sys.argv[1] if len(sys.argv) > 1 else None
    sd   = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    main(root, sd)
