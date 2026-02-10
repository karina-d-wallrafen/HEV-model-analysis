
import os
import numpy as np
import pandas as pd

if __name__ == '__main__':

    os.makedirs("obs", exist_ok=True)

    # Cases
    cases = pd.read_csv(f"../0_data/cases.csv", sep=",", index_col=False)

    cases = cases.fillna(0)
    cases = cases.melt(id_vars="Week", value_vars=[c for c in cases.columns if c not in ["Week"]],
                   var_name="Year", value_name="HEV")
    cases["Year"] = cases["Year"].astype(int)
    cases = pd.concat([cases, pd.DataFrame([{"Week": 53, "Year": 2020, "HEV": 0}])], ignore_index=True)

    cases["Category"] = "Non-Emergency"
    cases.loc[(cases["Year"] == 2020) & (cases["Week"].between(12, 25)), "Category"] = "Emergency"
    cases.loc[(cases["Year"] == 2020) & (cases["Week"] >= 50), "Category"] = "Emergency"
    cases.loc[(cases["Year"] == 2021) & (cases["Week"] <= 22), "Category"] = "Emergency"
    cases.loc[(cases["Year"] == 2021) & (cases["Week"] >= 45), "Category"] = "Emergency"
    cases.loc[(cases["Year"] == 2022) & (cases["Week"] <= 19), "Category"] = "Emergency"

    HEV_l = cases.loc[cases["Category"] == "Emergency", "HEV"]
    HEV_nl = cases.loc[cases["Category"] == "Non-Emergency", "HEV"]
    HEV_l.to_csv("obs/observed_cases_l", index=False, header=False)
    HEV_nl.to_csv("obs/observed_cases_nl", index=False, header=False)

    # Wastewater
    ww = pd.read_csv(f"../0_data/ww.csv", sep=",", index_col=False)

    ww = ww.loc[:, ["Sample ID", "Name", "Date", "Comment WW", "HEV per 5 µL", "PMMoV (Cp)"]]
    ww = ww.map(lambda v: pd.NA if isinstance(v, str) and '-' in v else v)
    ww = ww.replace({"": pd.NA})
    ww = ww.replace({"n.d.": 0})
    ww = ww.dropna(subset=["PMMoV (Cp)"])
    exclude_ids = {"2095", "2136", "2099", "2138", "2105", "2131", "2119", "2137", "2122", "2139", "2107", "2132",
                   "2134", "2295"}
    ww = ww[~ww["Sample ID"].astype(str).isin(exclude_ids)]
    ww.loc[ww["Name"].str.lower() == "schmidbartlanger", "Name"] = "Schmidbartlanger"
    ww = ww[ww["Name"].isin(["Schmidbartlanger", "24h Probe"])]

    s = ww["Date"].astype(str).str.strip()
    ww["Date"] = pd.to_datetime(s, format="%m/%d/%y", errors="coerce")
    mask = ww["Date"].isna()
    ww.loc[mask, "Date"] = pd.to_datetime(s[mask], format="%m/%d/%Y", errors="coerce")

    for col in ["HEV per 5 µL", "PMMoV (Cp)"]:
        ww[col] = pd.to_numeric(ww[col].astype(str).str.replace(",", "", regex=False), errors="coerce")
    ww["ConcF"] = np.where(ww["Comment WW"] == "frozen", 38 * 2 / 0.06, 38 / 0.06)
    ww["HEV_L"] = ww["HEV per 5 µL"] * 1000 / (0.005 * ww["ConcF"])
    ww["PMMoV_1µL"] = 1.972 ** (39.65 - ww["PMMoV (Cp)"])
    ww["PMMoV_L"] = ww["PMMoV_1µL"] * 1000 / (0.001 * ww["ConcF"])
    ww["PMMoV_mean"] = np.where(ww["Name"] == "24h Probe", 6_500_000,
                                 np.where(ww["Comment WW"] == "fresh", 15_000_000, 8_300_000))
    ww["Norm"] = ww["PMMoV_mean"] / ww["PMMoV_L"]
    ww["HEV_c"] = ww["HEV_L"] * ww["Norm"]

    ww["Category"] = "Non-Emergency"
    mask1 = (ww["Date"] >= pd.to_datetime("12/09/2020", format="%m/%d/%Y")) & (
                ww["Date"] <= pd.to_datetime("06/06/2021", format="%m/%d/%Y"))
    mask2 = (ww["Date"] >= pd.to_datetime("11/11/2021", format="%m/%d/%Y")) & (
                ww["Date"] <= pd.to_datetime("05/11/2022", format="%m/%d/%Y"))
    ww.loc[mask1 | mask2, "Category"] = "Emergency"

    ww_G_l = round(ww.loc[(ww["Category"] == "Emergency") & (ww["Name"] == "24h Probe"), "HEV_c"], 2)
    ww_G_nl = round(ww.loc[(ww["Category"] == "Non-Emergency") & (ww["Name"] == "24h Probe"), "HEV_c"], 2)
    ww_S_l = round(ww.loc[(ww["Category"] == "Emergency") & (ww["Name"] == "Schmidbartlanger"), "HEV_c"], 2)
    ww_S_nl = round(ww.loc[(ww["Category"] == "Non-Emergency") & (ww["Name"] == "Schmidbartlanger"), "HEV_c"], 2)
    ww_G_l.to_csv("obs/observed_GG_l", index=False, header=False)
    ww_G_nl.to_csv("obs/observed_GG_nl", index=False, header=False)
    ww_S_l.to_csv("obs/observed_S_l", index=False, header=False)
    ww_S_nl.to_csv("obs/observed_S_nl", index=False, header=False)



