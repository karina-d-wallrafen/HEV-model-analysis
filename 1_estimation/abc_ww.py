
import pyabc
import pyabc.external
import os
from numpy import loadtxt
import json
import scipy

if __name__ == '__main__':

    def model(pars):
        dir = "hev-model-r"
        myModel = pyabc.external.ExternalModel(
            executable="Rscript", file=f"{dir}/HEV_NW_Model_WW_Calib.R"
        )
        mf = myModel(pars)
        with open(mf.get("loc"), "r") as f:
            res = json.load(f)
        GG_nl = res["result_GG_nl"]
        GG_l = res["result_GG_l"]
        S_nl = res["result_S_nl"]
        S_l = res["result_S_l"]
        seed = res["seed"]
        GG_all = res["result_GG"]
        S_all = res["result_S"]
        os.remove(mf.get("loc"))
        return {"GG_nl": GG_nl, "GG_l": GG_l, "S_nl": S_nl, "S_l": S_l, "seed": seed, "GG_all": GG_all, "S_all": S_all}

    prior = pyabc.Distribution(ld=pyabc.RV("uniform", 0, 2),
                               scaleG=pyabc.RV("uniform", 0, 90000000),
                               scaleS=pyabc.RV("uniform", 0, 90000000),
                               noise1G=pyabc.RV("uniform", 0, 5), noise1S=pyabc.RV("uniform", 0, 5),
                               noise2G=pyabc.RV("uniform", 0, 200), noise2S=pyabc.RV("uniform", 0, 200))

    def distance1(x, x0):
        return scipy.stats.kstest(x["GG_nl"], x0["GG_nl"])[0]

    def distance2(x, x0):
        return scipy.stats.kstest(x["GG_l"], x0["GG_l"])[0]

    def distance3(x, x0):
        return scipy.stats.kstest(x["S_nl"], x0["S_nl"])[0]

    def distance4(x, x0):
        return scipy.stats.kstest(x["S_l"], x0["S_l"])[0]

    distance = pyabc.distance.AggregatedDistance(
        [distance1, distance2, distance3, distance4],
    )

    # transition = pyabc.transition.LocalTransition(k_fraction=0.25)

    sampler = pyabc.sampler.MulticoreEvalParallelSampler(n_procs=48)

    abc = pyabc.ABCSMC(model, prior, distance, population_size=200, sampler=sampler)

    db_path = os.path.join('output', "HEV_ww.db")
    obs_GG_nl = loadtxt("obs/observed_GG_nl", delimiter=" ", unpack=False)
    obs_GG_l = loadtxt("obs/observed_GG_l", delimiter=" ", unpack=False)
    obs_S_nl = loadtxt("obs/observed_S_nl", delimiter=" ", unpack=False)
    obs_S_l = loadtxt("obs/observed_S_l", delimiter=" ", unpack=False)
    abc.new("sqlite:///" + db_path, {"GG_nl": obs_GG_nl, "GG_l": obs_GG_l, "S_nl": obs_S_nl, "S_l": obs_S_l})

    history = abc.run(minimum_epsilon=0, max_nr_populations=20)

