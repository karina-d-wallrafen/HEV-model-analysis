
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

  HEV.cases.clin <- c(39, 52, 53, 60) # survstat, 2020-2023, Munich SK + LK, Ref def = Yes 
  HEV.cases.noclin <- c(60, 74, 65, 71) # survstat, 2020-2023, Munich SK + LK, Ref def = No 
  HEV.inc <- 5.18 / 1000 * 1.83 * 10^6 # Faber incidence times Munich SK + LK pop from 2022: https://www.citypopulation.de/de/germany/admin/bayern/09184__münchen/
  report.prob.a <- round(mean(HEV.cases.noclin / (HEV.inc * 0.99)), 3)
  report.prob.c <- round(mean(HEV.cases.clin / (HEV.inc * 0.01)), 3)

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
                     #report.disp.k = as.numeric(pars$k),
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
  
  results <- list(case_result_nl = case_result_nl, case_result_l = case_result_l, 
                  seed = seed, case_result = case_result)
  
  write_json(results, path = pars$target, auto_unbox = TRUE)
}

myModel(pars)