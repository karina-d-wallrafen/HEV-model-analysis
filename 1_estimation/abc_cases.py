
import pyabc
import pyabc.external
import os
from numpy import loadtxt
import json
import scipy
import numpy as np

if __name__ == '__main__':

    def model(pars):
        dir = "hev-model-r"
        myModel = pyabc.external.ExternalModel(
            executable="Rscript", file=f"{dir}/HEV_NW_Model_Case_Calib.R"
        )
        mf = myModel(pars)
        with open(mf.get("loc"), "r") as f:
            res = json.load(f)
        cases_nl = res["case_result_nl"]
        cases_l = res["case_result_l"]
        seed = res["seed"]
        cases_all = res["case_result"]
        os.remove(mf.get("loc"))
        return {"cases_nl": cases_nl, "cases_l": cases_l, "seed": seed, "cases_all": cases_all}

    prior = pyabc.Distribution(ld=pyabc.RV("uniform", 0, 2))

    def distance1(x, _x0_ignored):
        return np.mean([scipy.stats.ks_2samp(x["cases_nl"], yk).statistic for yk in obs_sets_nl])

    def distance2(x, _x0_ignored):
        return np.mean([scipy.stats.ks_2samp(x["cases_l"], yk).statistic for yk in obs_sets_l])

    distance = pyabc.distance.AggregatedDistance(
        [distance1, distance2],
    )

    sampler = pyabc.sampler.MulticoreEvalParallelSampler(n_procs=48)

    abc = pyabc.ABCSMC(model, prior, distance, population_size=200, sampler=sampler)

    db_path = os.path.join('output', "HEV_cases.db")

    obs_cases_nl_raw = loadtxt("obs/observed_cases_nl", delimiter=" ", unpack=False).astype(int)
    obs_cases_l_raw = loadtxt("obs/observed_cases_l", delimiter=" ", unpack=False).astype(int)
    p = 1.0 / 1.83
    K = 200
    rng = np.random.default_rng(123)
    obs_sets_nl = [rng.binomial(n=obs_cases_nl_raw, p=p) for _ in range(K)]
    obs_sets_l = [rng.binomial(n=obs_cases_l_raw, p=p) for _ in range(K)]
    np.savez("output/thinned_cases.npz", obs_sets_nl=obs_sets_nl, obs_sets_l=obs_sets_l)

    abc.new("sqlite:///" + db_path, {"cases_nl": obs_sets_nl[0], "cases_l": obs_sets_l[0]})

    history = abc.run(minimum_epsilon=0, max_nr_populations=20)

