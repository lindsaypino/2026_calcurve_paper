"""Load figures-of-merit CSVs and config; small helpers shared by figure scripts."""
import importlib.util
import os
import re
import numpy as np
import pandas as pd
import yaml

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def load_config(path=None):
    """Load config.yaml and expand ${...} references within the raw: block."""
    path = path or os.path.join(_REPO, "config.yaml")
    with open(path) as f:
        cfg = yaml.safe_load(f)
    raw = cfg.get("raw", {})
    for _ in range(3):  # resolve nested ${...} references
        for k, v in list(raw.items()):
            if isinstance(v, str):
                raw[k] = re.sub(r"\$\{(\w+)\}", lambda m: str(raw.get(m.group(1), m.group(0))), v)
    return cfg


def repo_path(*parts):
    return os.path.join(_REPO, *parts)


def raw(key, cfg=None):
    """Resolve a raw-data path from config (e.g. raw('exploris_elib'))."""
    cfg = cfg or load_config()
    return cfg["raw"][key]


def load_tool(cfg=None):
    """Import the pinned calculate-loq.py (submodule) as a module."""
    cfg = cfg or load_config()
    path = repo_path(cfg["tool"]["calculate_loq"])
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"tool not found at {path}; run `git submodule update --init --recursive`")
    spec = importlib.util.spec_from_file_location("calcloq", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def load_fom(group, name, cfg=None):
    """Load a committed figuresofmerit CSV, e.g. load_fom('main', 'bruker_60spd')."""
    cfg = cfg or load_config()
    return pd.read_csv(repo_path(cfg["fom"][group], f"{name}.csv"))


def finite(series):
    """Finite, positive values of a FOM column (LOD/LOQ/ULOQ)."""
    v = pd.to_numeric(series, errors="coerce").to_numpy(float)
    return v[np.isfinite(v)]


def finite_count(series):
    return int(np.isfinite(pd.to_numeric(series, errors="coerce")).sum())
