
import pyabc
import numpy as np
import subprocess
import json
from pathlib import Path
import tempfile
from itertools import repeat
from concurrent.futures import ProcessPoolExecutor

def run_hev_model(seed, alpha, ld_map, ldr_map, r_script):

    tmp = tempfile.NamedTemporaryFile(suffix=".json", delete=False)
    target = Path(tmp.name)
    tmp.close()

    cmd = ["Rscript", str(r_script), f"seed={seed}", f"ldr={ldr_map}", f"ld={ld_map}", f"alpha={alpha}", f"target={target}"]

    subprocess.run(cmd, check=True)

    with target.open() as f:
        res = json.load(f)

    target.unlink()
    return seed, res


if __name__ == '__main__':

    db_path = f"../1_estimation/output/HEV_dual.db"
    history = pyabc.History("sqlite:///" + db_path)
    df, w = history.get_distribution(0)
    x, y, PDF = pyabc.visualization.kde.kde_2d(df, w, x="ld", y="ldr", xmin=0, xmax=2, numx=2001, ymin=0, ymax=5, numy=5001)

    imax = np.argmax(PDF)
    iy, ix = np.unravel_index(imax, PDF.shape)
    ld_map = x[iy, ix]
    ldr_map = y[iy, ix]

    R_SCRIPT = Path("HEV_NW_Model_Sim.R")
    NUM_RUNS = 100
    NUM_CORES = 48
    NUM_SCENARIOS = 11

    alpha_values = np.linspace(0.0, 1.0, NUM_SCENARIOS)
    all_results = {}

    for alpha in alpha_values:
        print(f"\n=== Scenario alpha = {alpha:.3f} ===")
        scenario_results = {}

        with ProcessPoolExecutor(max_workers=NUM_CORES) as pool:
            futures = pool.map(run_hev_model, range(1, NUM_RUNS + 1),
                               repeat(alpha), repeat(ld_map), repeat(ldr_map), repeat(R_SCRIPT))

            for seed, out in futures:
                print(f"  Completed seed {seed}")
                scenario_results[seed] = out

        all_results[f"{alpha:.3f}"] = scenario_results

    OUTPUT_FILE = "hev_sim_results.json"
    with open(OUTPUT_FILE, "w") as f:
        json.dump(all_results, f, indent=2)