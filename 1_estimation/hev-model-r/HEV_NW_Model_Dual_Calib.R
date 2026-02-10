
suppressWarnings(suppressMessages(library(EpiModel)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(jsonlite)))

args = commandArgs(trailingOnly=TRUE)
pars = list()
names_ = list()
for (arg in args) {
  s = strsplit(arg, "=")
  names_ = c(names_, s[[1]][1])
  pars = c(pars, s[[1]][2])
}
names(pars) = names_

myModel <- function(pars) {
  
  seed <- round(as.numeric(Sys.time()))
  set.seed(seed)
  
  n <- 1000000
  subpop.prop <- 67000 / 1000000
  n.subpop <- n * subpop.prop

  source("../0_model_files/demo_params.R")

  load("../0_model_files/est-1m.RDS")
  source("../0_model_files/household_assignment.R", local = TRUE)

  num.by.hh <- persons.by.hh %>% group_by(hh) %>% summarize(num = n())
  num.by.hh$dist <- abs(cumsum(num.by.hh$num) - n.subpop)
  max.hh <- num.by.hh[which.min(num.by.hh$dist),]$hh
  persons.by.hh$subpop <- ifelse(persons.by.hh$hh <= max.hh, 1, 0)
  est[[1]]$newnetwork <- set.vertex.attribute(est[[1]]$newnetwork, "subpop", persons.by.hh$subpop)

  HEV.cases.clin <- c(53, 60) # survstat, 2022-2023, Munich SK + LK, Ref def = Yes 
  HEV.cases.noclin <- c(65, 71) # survstat, 2022-2023, Munich SK + LK, Ref def = No 
  HEV.inc <- 5.18 / 1000 * 1.83 * 10^6 # Faber incidence times Munich SK + LK pop from 2022: https://www.citypopulation.de/de/germany/admin/bayern/09184__münchen/
  report.prob.a <- round(mean(HEV.cases.noclin / (HEV.inc * 0.99)), 3) # for non-lockdown periods
  report.prob.c <- round(mean(HEV.cases.clin / (HEV.inc * 0.01)), 3) # for non-lockdown periods

  source("../0_model_files/helper_funs.R")
  source("../0_model_files/HEV_modules.R")
  param <- param.net(hh.pairs = hhPairs,
                     act.rate = c(1, 1),
                     hepE3.pp.inf.prob = c(0, 0),
                     hepE3.prop.clin = 0.01,
                     hepE3.shed.start.rate = 1 / 21,
                     hepE3.shed.stop.a.rate = 1 / 35,
                     hepE3.shed.stop.c.rate = 1 / 16,
                     hepE3.sym.start.rate = 1 / 19,
                     hepE3.sym.stop.rate = 1 / 21,
                     hepE3.rs.rate = 0.019 / 365,
                     veg.prop = c(0.90, 0.02, 0.08),
                     subpop.prop = subpop.prop,
                     mr = mr_pp_pd,
                     br = sum(mr_pp_pd * age.weights / sum(age.weights)),
                     hepE3.a.report.prob = report.prob.a,
                     hepE3.c.report.prob = report.prob.c,
                     ld.rep.scale = as.numeric(pars$ldr),
                     ld.scale = as.numeric(pars$ld),
                     ld.starts = c(76, 344, 681),
                     ld.stops = c(168, 523, 862),
                     hepE3.pork.inf.prob = 2e-05)

  init <- init.net(hepE3.r.num = 0.168 * n)

  control <- control.net(
    nsteps = 1461, # 01 Jan 2020 - 31 Dec 2023
    nsims = 1,
    ncores = 1,
    initialize.FUN = initialise,
    infection.FUN = infect,
    recovery.FUN = progress,
    dx.FUN = dx,
    prevalence.FUN = prev,
    module.order = c("infection.FUN", "recovery.FUN", "dx.FUN", "prevalence.FUN"),
    resimulate.network = TRUE,
    tergmLite = TRUE,
    verbose = FALSE,
    save.nwstats = FALSE,
    cumulative.edgelist = FALSE,
    truncate.el.cuml = Inf
  )

  sim <- netsim(est, param, init, control)
  
  # Gut Grosslappen 
  result_GG <- sim$epi$hepE3.total.shed.num$sim1
  zero_inflation_GG <- rbinom(length(result_GG), 1, 1 / (1 / as.numeric(pars$noise2G) * result_GG + 1))
  result_GG <- ifelse(zero_inflation_GG == 1, 0, result_GG)
  result_GG <- result_GG / n * as.numeric(pars$scaleG)
  result_GG <- exp(log(result_GG) + as.numeric(pars$noise1G) * rnorm(length(result_GG)))
  days_G_nl <- c(868,926,987,994,1001,1008,1015,1022,1029,1037,1046,1050,1057,1067,1071,1078,1085,1092,
                 1099,1120,1127,1134,1141,1148,1155,1162,1169,1176,1183,1190,1197,1204,1211,1218,1225,1232,
                 1239,1246,1260,1267,1274,1281,1288,1296,1302,1312,1324,1330,1337,1344,1351,1371,1379,1393,
                 1407,1414,1421)
  days_G_l <- c(688,693,700,707,714,721,728,735,749,763,770,780,784,791,798,812,819,850)
  result_GG_nl <- result_GG[days_G_nl]
  result_GG_l <- result_GG[days_G_l]
  
  # Schmidbartlanger
  result_S <- sim$epi$hepE3.subpop.shed.num$sim1
  zero_inflation_S <- rbinom(length(result_S), 1, 1 / (1 / as.numeric(pars$noise2S) * result_S + 1))
  result_S <- ifelse(zero_inflation_S == 1, 0, result_S)
  result_S <- result_S / n.subpop * as.numeric(pars$scaleS)
  result_S <- exp(log(result_S) + as.numeric(pars$noise1S) * rnorm(length(result_S)))
  days_S_nl <- c(189,203,210,219,231,275,301,329,525,532,560,595,602,609,616,623,637,645,
                 651,658,665,679,868,926,934,945,955,962,969,980,987,994,1001,1004,1008,1015,
                 1022,1037,1046,1050,1057,1067,1071,1078,1085,1092,1099,1106,1113,1120,1127,1134,1141,1148,
                 1155,1162,1169,1183,1190,1197,1204,1211,1218,1225,1232,1239,1246,1253,1260,1267,1274,1281,
                 1288,1296,1302,1428)
  days_S_l <- c(350,364,378,385,392,399,406,413,427,434,455,511,689,693,707,714,721,728,735,742,749,756,763,
                770,780,784,791,798,805,812,819,826,833,840,854)
  result_S_nl <- result_S[days_S_nl]
  result_S_l <- result_S[days_S_l]
  
  # Cases
  case_result <- sim$epi$nRepHE3$sim1
  case_result[is.na(case_result)] <- 0
  
  df = data.frame(day = 1:1461, date = seq(from = as.Date("2020-01-01"), to = as.Date("2023-12-31"), by = "day"))
  df$week <- substr(ISOweek::date2ISOweek(df$date), 1, 8)
  l_idx <- unlist(Map(`:`, sim$param$ld.starts, sim$param$ld.stops))
  df$cat <- "non-lockdown"
  df[l_idx, ]$cat <- "lockdown"
  df$dx <- case_result
  
  df <- df %>% group_by(week) %>% summarise(dx = sum(dx, na.rm = TRUE), cat = if_else(any(cat == "lockdown"), "lockdown", "non-lockdown"))
  
  case_result_l <- df[df$cat == "lockdown", ]$dx
  case_result_nl <- df[df$cat == "non-lockdown", ]$dx
  
  results <- list(result_GG_nl = result_GG_nl, result_GG_l = result_GG_l, 
                  result_S_nl = result_S_nl, result_S_l = result_S_l,
                  case_result_nl = case_result_nl, case_result_l = case_result_l, 
                  seed = seed, result_GG_all = result_GG, result_S_all = result_S,
                  result_cases_all = case_result)
  
  write_json(results, path = pars$target, auto_unbox = TRUE)
}

myModel(pars)