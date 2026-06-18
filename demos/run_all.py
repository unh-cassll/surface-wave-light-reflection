"""
Run the seapol demo gallery end to end.

Each demo executes in its own subprocess with output captured to
output/logs/<name>.log; demos whose data requirements are not met
(ASIT trees, cached FM98 table, beta library) are skipped with the
reason printed.  Demos run in dependency order: table-building demos
precede table-loading ones.

    python run_all.py                 # everything available
    python run_all.py --list          # show the registry and status
    python run_all.py --only mc       # substring filter
"""

import argparse
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).parent
OUT = HERE / "output"
LOGS = OUT / "logs"

TABLE = OUT / "fm98_table_deep.npz"
EPSS = Path("/home/nathanlaxague/Dropbox/Professional/Github/E-PSS_paper/_data")
STATS = EPSS / "ASIT2019_wave_spectra_stats_timeseries_empirical_gain.nc"
ENV = EPSS / "ASIT2019_supporting_environmental_observations.nc"
BETA_LIB = OUT / "asit_beta_library"


def _check_table():
    return TABLE.exists(), f"needs cached FM98 table {TABLE.name} " \
        "(run demo_fm98_capillaries first)"


def _check_asit():
    ok = STATS.exists() and ENV.exists()
    return ok, f"needs ASIT stats/env data under {EPSS}"


def _check_beta_lib():
    ok = BETA_LIB.is_dir() and any(BETA_LIB.glob("*.npz"))
    return ok, f"needs beta library {BETA_LIB}"


# (script, requirement checks, runtime hint); dependency order.  Checks
# are evaluated just before each demo runs, so artifacts produced by
# earlier demos (table cache, k-f reduction) count.
REGISTRY = [
    ("demo_render_panel.py", [], "fast"),
    ("demo_time_evolution.py", [], "medium"),
    ("demo_mc_reflectance.py", [], "medium"),
    ("demo_mc_water_body.py", [], "medium"),
    ("demo_glint_foam_current.py", [], "medium"),
    ("demo_polarimetric_reconstruction.py", [], "fast"),
    ("demo_polarimetric_polarized.py", [], "fast"),
    ("demo_sky_water_gallery.py", [], "medium"),
    ("demo_near_surface_scattering.py", [],
     "medium; builds upwelling table if missing"),
    ("demo_color_scenes.py", [],
     "slow; builds per-band tables if missing"),
    ("demo_stokes_panels.py", [_check_table],
     "medium; GPU when available, reuses the deep FM98 table"),
    ("demo_fm98_crest_stokes.py", [], "medium"),
    ("demo_fm98_capillaries.py", [], "slow; builds FM98 table if missing"),
    ("demo_fm98_3d_placement.py", [], "medium; builds table if missing"),
    ("demo_slope_statistics.py", [_check_table], "slow"),
    ("demo_kw_spectrum.py", [_check_table], "slow"),
    ("demo_full_pipeline.py",
     [_check_table, _check_asit, _check_beta_lib], "slow"),
]


def run_demo(script: str, verbose: bool) -> tuple[str, float]:
    """Run one demo; returns (status, seconds)."""
    log_path = LOGS / (Path(script).stem + ".log")
    t0 = time.time()
    if verbose:
        proc = subprocess.run([sys.executable, script], cwd=HERE)
    else:
        with open(log_path, "w") as log:
            proc = subprocess.run([sys.executable, script], cwd=HERE,
                                  stdout=log, stderr=subprocess.STDOUT)
    dt = time.time() - t0
    if proc.returncode != 0 and not verbose:
        tail = log_path.read_text().splitlines()[-15:]
        print(f"      --- tail of {log_path.name} ---")
        for line in tail:
            print(f"      {line}")
    return ("PASS" if proc.returncode == 0 else "FAIL"), dt


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--list", action="store_true",
                    help="list demos with requirement status and exit")
    ap.add_argument("--only", action="append", default=[],
                    metavar="SUBSTR",
                    help="run only demos whose name contains SUBSTR "
                         "(repeatable)")
    ap.add_argument("--fail-fast", action="store_true",
                    help="stop at the first failure")
    ap.add_argument("--verbose", action="store_true",
                    help="stream demo output instead of logging")
    args = ap.parse_args()

    registry = list(REGISTRY)
    if args.only:
        registry = [r for r in registry
                    if any(s in r[0] for s in args.only)]
    if not registry:
        print("no demos match the filter")
        return 1

    if args.list:
        for script, checks, hint in registry:
            missing = [msg for chk in checks
                       for ok, msg in [chk()] if not ok]
            status = "ready" if not missing else f"SKIP: {missing[0]}"
            print(f"{script:32s} [{hint}] {status}")
        return 0

    LOGS.mkdir(parents=True, exist_ok=True)
    results = []
    for script, checks, hint in registry:
        missing = [msg for chk in checks for ok, msg in [chk()] if not ok]
        if missing:
            print(f"SKIP  {script:32s} {missing[0]}")
            results.append((script, "SKIP", 0.0))
            continue
        print(f"RUN   {script:32s} [{hint}]", flush=True)
        status, dt = run_demo(script, args.verbose)
        print(f"{status:5s} {script:32s} {dt:7.1f} s")
        results.append((script, status, dt))
        if status == "FAIL" and args.fail_fast:
            break

    n = {s: sum(1 for r in results if r[1] == s)
         for s in ("PASS", "FAIL", "SKIP")}
    total = sum(r[2] for r in results)
    print(f"\n{n['PASS']} passed, {n['FAIL']} failed, {n['SKIP']} skipped "
          f"in {total:.0f} s; logs in {LOGS}")
    for script, status, _ in results:
        if status == "FAIL":
            log = LOGS / (Path(script).stem + ".log")
            print(f"  FAILED: {script} (see {log})")
    return 1 if n["FAIL"] else 0


if __name__ == "__main__":
    sys.exit(main())
