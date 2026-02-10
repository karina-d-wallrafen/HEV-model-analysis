
infect <- function(dat, at) {

  hepE3.status <- get_attr(dat, "hepE3.status")
  veg <- get_attr(dat, "veg")
  
  ## Hepatitis E (Gen3) Pork-Borne Transmission
  hepE3.pork.inf.prob <- get_param(dat, "hepE3.pork.inf.prob") 
  ld.starts <- get_param(dat, "ld.starts")
  ld.stops <- get_param(dat, "ld.stops")
  ld.scale <- get_param(dat, "ld.scale")
  
  for (i in seq_along(ld.starts)) {
    if (at >= ld.starts[i] & at <= ld.stops[i]) {
      hepE3.pork.inf.prob <- hepE3.pork.inf.prob * ld.scale
    }
  }

  elig <- which(veg == 0 & hepE3.status == "s")

  transmit.pork <- rbinom(length(elig), 1, hepE3.pork.inf.prob)

  idsNewInfE3 <- elig[which(transmit.pork == 1)]
  nInfE3.f <- length(idsNewInfE3)
  
  if (nInfE3.f > 0) {
    dat <- EpiModel::set_attr(dat, "hepE3.status", "e", idsNewInfE3)
    dat <- EpiModel::set_attr(dat, "hepE3.infTime", at, idsNewInfE3)
    
    hepE3.prop.clin <- get_param(dat, "hepE3.prop.clin")
    dat <- EpiModel::set_attr(dat, "hepE3.clin", rbinom(nInfE3.f, 1, hepE3.prop.clin), idsNewInfE3)
  }
  
  ## Hepatitis E (Gen3) Person-to-Person Transmission
  Layers <- c("community", "household")
  nInfE3.pp <- c(0, 0)

  if (sum(get_param(dat, "hepE3.pp.inf.prob")) > 0) {
    idsInfE3 <- which(hepE3.status %in% c("a", "ip", "ic", "ih"))
    
    if (length(idsInfE3) > 0) {
      
      for (layer in seq_along(Layers)) {
        
        if (get_param(dat, "hepE3.pp.inf.prob")[layer] > 0){
        
          del <- get_discordant_edgelist(dat = dat, infs = idsInfE3, status = "hepE3.status", sus.status = "s", network = Layers[layer])
        
          if (nrow(del) > 0) {
            
            del$transProb <- get_param(dat, "hepE3.pp.inf.prob")[layer]
            del$actRate <- get_param(dat, "act.rate")[layer]
            
            del$finalProb <- 1 - (1 - del$transProb)^del$actRate
            transmit <- rbinom(nrow(del), 1, del$finalProb)
            del <- del[which(transmit == 1), , drop = FALSE]
            
            idsNewInfE3 <- unique(del$tail)
            nInfE3.pp[layer] <- length(idsNewInfE3)
            
            if (nInfE3.pp[layer] > 0) {
              dat <- EpiModel::set_attr(dat, "hepE3.status", "e", idsNewInfE3)
              dat <- EpiModel::set_attr(dat, "hepE3.infTime", at, idsNewInfE3)
              
              hepE3.prop.clin <- get_param(dat, "hepE3.prop.clin")
              dat <- EpiModel::set_attr(dat, "hepE3.clin", rbinom(nInfE3.pp[layer], 1, hepE3.prop.clin), idsNewInfE3)
            }
          }
        }
      }
    }
  }

  dat <- set_epi(dat, "nInfE3.f", at, nInfE3.f)
  dat <- set_epi(dat, "nInfE3.pp", at, sum(nInfE3.pp))
  
  return(dat)
}

resim <- function(dat, at) {

  dat <- edges_correct(dat, at)
  
  # dat.updates <- NVL(get_control(dat, "dat.updates"), function(dat, ...) dat)
  # dat <- dat.updates(dat = dat, at = at, network = 0L)
  
  for (network in seq_len(dat$num.nw)) {
    dat <- simulate_dat(dat = dat, at = at, network = network)
    # dat <- dat.updates(dat = dat, at = at, network = network)
  }
  
  dat <- summary_nets(dat, at)
  
  # for (network in seq_len(dat$num.nw)) {
  #   dat <- update_cumulative_edgelist(dat = dat, network = network, truncate = get_control(dat, "truncate.el.cuml"))
  # }
  
  return(dat)
}

progress <- function(dat, at) {

  hepE3.status <- get_attr(dat, "hepE3.status")
  hepE3.clin <- get_attr(dat, "hepE3.clin")
  hepE3.infTime <- get_attr(dat, "hepE3.infTime")
  hepE3.symTime <- get_attr(dat, "hepE3.symTime")
  hepE3.recTime <- get_attr(dat, "hepE3.recTime")
  hepE3.shedTime <- get_attr(dat, "hepE3.shedTime")
  hepE3.stopTime <- get_attr(dat, "hepE3.stopTime")
  
  hepE3.shed.start.rate <- get_param(dat, "hepE3.shed.start.rate")
  hepE3.shed.stop.a.rate <- get_param(dat, "hepE3.shed.stop.a.rate")
  hepE3.shed.stop.c.rate <- get_param(dat, "hepE3.shed.stop.c.rate")
  hepE3.sym.start.rate <- get_param(dat, "hepE3.sym.start.rate")
  hepE3.sym.stop.rate <- get_param(dat, "hepE3.sym.stop.rate")
  hepE3.rs.rate <- get_param(dat, "hepE3.rs.rate")
  
  ## Asymptomatic disease progression
  ids.asympt <- which(hepE3.clin == 0)
  asympt.status <- hepE3.status[ids.asympt]
  
  # E to A
  ids.hepE3.Ea <- which(asympt.status == "e" & hepE3.infTime[ids.asympt] < at)
  if (length(ids.hepE3.Ea) > 0) {
    vec.new.A <- rbinom(length(ids.hepE3.Ea), 1, hepE3.shed.start.rate)
    if (sum(vec.new.A) > 0) {
      ids.new.A <- ids.asympt[ids.hepE3.Ea[which(vec.new.A == 1)]]
      hepE3.status[ids.new.A] <- "a"
      hepE3.shedTime[ids.new.A] <- at
    }
  }
  
  # A to R
  ids.hepE3.A <- which(asympt.status == "a")
  if (length(ids.hepE3.A) > 0) {
    vec.new.R <- rbinom(length(ids.hepE3.A), 1, hepE3.shed.stop.a.rate)
    if (sum(vec.new.R) > 0) {
      ids.new.R <- ids.asympt[ids.hepE3.A[which(vec.new.R == 1)]]
      hepE3.status[ids.new.R] <- "r"
      hepE3.stopTime[ids.new.R] <- at
    }
  }
  
  ## Symptomatic disease progression
  ids.clin <- which(hepE3.clin == 1)
  clin.status <- hepE3.status[ids.clin]
  
  # E to Ip
  ids.hepE3.Ec <- which(clin.status == "e" & hepE3.infTime[ids.clin] < at)
  if (length(ids.hepE3.Ec) > 0) {
    vec.new.Ip <- rbinom(length(ids.hepE3.Ec), 1, hepE3.shed.start.rate)
    if (sum(vec.new.Ip) > 0) {
      ids.new.Ip <- ids.clin[ids.hepE3.Ec[which(vec.new.Ip == 1)]]
      hepE3.status[ids.new.Ip] <- "ip"
      hepE3.shedTime[ids.new.Ip] <- at
    }
  }
  
  # Ip to Ic
  ids.hepE3.Ip <- which(clin.status == "ip")
  if (length(ids.hepE3.Ip) > 0) {
    vec.new.Ic <- rbinom(length(ids.hepE3.Ip), 1, hepE3.sym.start.rate)
    if (sum(vec.new.Ic) > 0) {
      ids.new.Ic <- ids.clin[ids.hepE3.Ip[which(vec.new.Ic == 1)]]
      hepE3.status[ids.new.Ic] <- "ic"
      hepE3.symTime[ids.new.Ic] <- at
      
      hepE3.stopTime[ids.new.Ic] <- at + round(rexp(length(ids.new.Ic), rate = hepE3.shed.stop.c.rate))
      hepE3.recTime[ids.new.Ic] <- at + round(rexp(length(ids.new.Ic), rate = hepE3.sym.stop.rate))
    }
  }
  
  # Ic to Is
  ids.new.Is <- which(hepE3.stopTime[ids.clin] == at & hepE3.recTime[ids.clin] > at)
  if (length(ids.new.Is) > 0) {
    hepE3.status[ids.clin[ids.new.Is]] <- "is"
  }
  
  # Ic to Ih
  ids.new.Ih <- which(hepE3.recTime[ids.clin] == at & hepE3.stopTime[ids.clin] > at)
  if (length(ids.new.Ih) > 0) {
    hepE3.status[ids.clin[ids.new.Ih]] <- "ih"
  }
  
  # I to R
  ids.new.R <- which(hepE3.status[ids.clin] %in% c("is", "ih", "ic") & hepE3.recTime[ids.clin] <= at & hepE3.stopTime[ids.clin] <= at)
  if (length(ids.new.R) > 0) {
    hepE3.status[ids.clin[ids.new.R]] <- "r"
  }
  
  # R to S
  ids.hepE3.R <- which(hepE3.status == "r" & (is.na(hepE3.recTime) | hepE3.recTime < at) & (is.na(hepE3.stopTime) | hepE3.stopTime < at))
  if (length(ids.hepE3.R) > 0) {
    vec.new.S <- rbinom(length(ids.hepE3.R), 1, hepE3.rs.rate)
    if (sum(vec.new.S) > 0) {
      ids.new.S <- ids.hepE3.R[which(vec.new.S == 1)]
      hepE3.status[ids.new.S] <- "s"
      hepE3.symTime[ids.new.S] <- NA
      hepE3.recTime[ids.new.S] <- NA
      hepE3.shedTime[ids.new.S] <- NA
      hepE3.stopTime[ids.new.S] <- NA
      
      dat <- EpiModel::set_attr(dat, "hepE3.infTime", NA, ids.new.S)
      dat <- EpiModel::set_attr(dat, "hepE3.clin", NA, ids.new.S)
      
      # dat <- EpiModel::set_attr(dat, "hepE3.dx", NA, ids.new.S)
      # dat <- EpiModel::set_attr(dat, "hepE3.dxTime", NA, ids.new.S)
    }
  }
  
  dat <- EpiModel::set_attr(dat, "hepE3.status", hepE3.status)
  dat <- EpiModel::set_attr(dat, "hepE3.symTime", hepE3.symTime)
  dat <- EpiModel::set_attr(dat, "hepE3.recTime", hepE3.recTime)
  dat <- EpiModel::set_attr(dat, "hepE3.shedTime", hepE3.shedTime)
  dat <- EpiModel::set_attr(dat, "hepE3.stopTime", hepE3.stopTime)
  
  # record.ids <- which(hepE3.stopTime == at)
  # dat <- record_attr_history(dat, at, "clin", record.ids, hepE3.clin[record.ids])
  # dat <- record_attr_history(dat, at, "infTime", record.ids, dat$attr$hepE3.infTime[record.ids])
  # dat <- record_attr_history(dat, at, "shedTime", record.ids, hepE3.shedTime[record.ids])
  # dat <- record_attr_history(dat, at, "symTime", record.ids, hepE3.symTime[record.ids])
  # dat <- record_attr_history(dat, at, "stopTime", record.ids, hepE3.stopTime[record.ids])
  # dat <- record_attr_history(dat, at, "recTime", record.ids, hepE3.recTime[record.ids])
  
  return(dat)
  
}

initialise <- function(x, param, init, control, s) {

  dat <- create_dat_object(param, init, control)
  dat <- init_nets(dat, x)
  dat <- sim_nets_t1(dat)
  dat <- summary_nets(dat, at = 1L)

  dat <- init_status(dat)
  dat <- prev(dat, at = 1)
  
  return(dat)
}

init_status <- function(dat) {

  hepE3.r.num <- get_init(dat, "hepE3.r.num")
  
  active <- get_attr(dat, "active")
  num <- sum(active)
  
  hepE3.status <- rep("s", num)
  hepE3.status[sample(1:num, size = hepE3.r.num)] <- "r"
  dat <- EpiModel::set_attr(dat, "hepE3.status", hepE3.status)
  
  # Infection Time and related attributes
  dat <- EpiModel::set_attr(dat, "hepE3.infTime", rep(NA, num))
  dat <- EpiModel::set_attr(dat, "hepE3.clin", rep(NA, num))
  dat <- EpiModel::set_attr(dat, "hepE3.symTime", rep(NA, num))
  dat <- EpiModel::set_attr(dat, "hepE3.recTime", rep(NA, num))
  dat <- EpiModel::set_attr(dat, "hepE3.shedTime", rep(NA, num))
  dat <- EpiModel::set_attr(dat, "hepE3.stopTime", rep(NA, num))
  # dat <- EpiModel::set_attr(dat, "hepE3.dx", rep(NA, num))
  # dat <- EpiModel::set_attr(dat, "hepE3.dxTime", rep(NA, num))
  
  return(dat)
}

prev <- function(dat, at) {
  
  hepE3.status <- EpiModel::get_attr(dat, "hepE3.status")
  subpop <- EpiModel::get_attr(dat, "subpop")
  
  # dat <- set_epi(dat, "num", at, length(hepE3.status))
  
  # dat <- set_epi(dat, "hepE3.s.num", at, sum(hepE3.status == "s"))
  # dat <- set_epi(dat, "hepE3.e.num", at, sum(hepE3.status == "e"))
  # dat <- set_epi(dat, "hepE3.a.num", at, sum(hepE3.status == "a"))
  # dat <- set_epi(dat, "hepE3.ip.num", at, sum(hepE3.status == "ip"))
  # dat <- set_epi(dat, "hepE3.ic.num", at, sum(hepE3.status == "ic"))
  # dat <- set_epi(dat, "hepE3.is.num", at, sum(hepE3.status == "is"))
  # dat <- set_epi(dat, "hepE3.ih.num", at, sum(hepE3.status == "ih"))
  # dat <- set_epi(dat, "hepE3.r.num", at, sum(hepE3.status == "r"))
  
  dat <- set_epi(dat, "hepE3.total.shed.num", at, sum(hepE3.status %in% c("a", "ip", "ic", "ih")))
  dat <- set_epi(dat, "hepE3.subpop.shed.num", at, sum(subpop == 1 & hepE3.status %in% c("a", "ip", "ic", "ih")))

  return(dat)
}

aging <- function(dat, at) {

  age <- get_attr(dat, "age")

  new.age <- age + 1 / 365
  age.breaks <- c(0, 18, 65, 200)
  new.age.grp <- cut(new.age, age.breaks, labels = FALSE, right = FALSE)

  dat <- EpiModel::set_attr(dat, "age", new.age)
  dat <- EpiModel::set_attr(dat, "age.grp", new.age.grp)

  return(dat)
}

death <- function(dat, at) {

  age <- get_attr(dat, "age")
  mr <- get_param(dat, "mr")

  # General deaths
  ageCatElig <- pmin(floor(age), 100) + 1
  rates <- mr[ageCatElig]
  isDep <- runif(length(age)) < rates

  if (sum(isDep) > 0) {
    idsDep <- which(isDep)
    dat <- depart_nodes(dat, idsDep)
    dat$param$hh.pairs <- delete_vertices(dat$param$hh.pairs, idsDep)
  }

  dat <- set_epi(dat, "total.deaths", at, sum(isDep))

  return(dat)
}

birth <- function(dat, at) {

  br <- get_param(dat, "br")
  num <- dat$epi$num[1]
  nNew <- rpois(1, br * num)

  if (nNew > 0) {
    
    veg.prop <- get_param(dat, "veg.prop")
    subpop.prop <- get_param(dat, "subpop.prop")

    dat <- append_core_attr(dat, at, nNew)

    dat <- append_attr(dat, "age", 0, nNew)
    dat <- append_attr(dat, "age.grp", 1, nNew)
    dat <- append_attr(dat, "hepE3.status", "s", nNew)
    dat <- append_attr(dat, "hepE3.clin", NA, nNew)
    dat <- append_attr(dat, "hepE3.infTime", NA, nNew)
    dat <- append_attr(dat, "hepE3.symTime", NA, nNew)
    dat <- append_attr(dat, "hepE3.recTime", NA, nNew)
    dat <- append_attr(dat, "hepE3.shedTime", NA, nNew)
    dat <- append_attr(dat, "hepE3.stopTime", NA, nNew)
    # dat <- append_attr(dat, "hepE3.dx", NA, nNew)
    # dat <- append_attr(dat, "hepE3.dxTime", NA, nNew)
    dat <- append_attr(dat, "veg", sample(0:2, size = nNew, prob = veg.prop, replace = TRUE), nNew)

    # Assign each new node to a household with an adult
    age.grp <- get_attr(dat, "age.grp")
    household <- get_attr(dat, "household")
    newHH <- sample(unique(household[which(age.grp == 2)]), nNew, replace = TRUE)
    dat <- append_attr(dat, "household", newHH, nNew)
    
    # Update household edgelist
    heads <- cbind((length(household) + 1):(length(household) + nNew), newHH)
    tails <- cbind(which(household %in% newHH), household[which(household %in% newHH)])
    new.edges <- merge(heads, tails, by.x = 2, by.y = 2)[, 2:3]
    dat$param$hh.pairs <- rbind(dat$param$hh.pairs, setNames(new.edges, c(".head", ".tail")))
    attr(dat$param$hh.pairs, "n") <- attr(dat$param$hh.pairs, "n") + nNew
    
    # Set subpop based on household
    subpop.new <- get_attr(dat, "subpop", match(newHH, household))
    dat <- append_attr(dat, "subpop", subpop.new, nNew)

    dat <- arrive_nodes(dat, nNew)
  }

  return(dat)
}

dx <- function(dat, at) {

  hepE3.status <- get_attr(dat, "hepE3.status")
  hepE3.clin <- get_attr(dat, "hepE3.clin")
  
  hepE3.shedTime <- get_attr(dat, "hepE3.shedTime")
  hepE3.stopTime <- get_attr(dat, "hepE3.stopTime")
  hepE3.recTime <- get_attr(dat, "hepE3.recTime")
  
  hepE3.c.report.prob <- get_param(dat, "hepE3.c.report.prob")
  hepE3.a.report.prob <- get_param(dat, "hepE3.a.report.prob")
  
  ld.starts <- get_param(dat, "ld.starts")
  ld.stops <- get_param(dat, "ld.stops")
  ld.rep.scale <- get_param(dat, "ld.rep.scale")
  
  for (i in seq_along(ld.starts)) {
    if (at >= ld.starts[i] & at <= ld.stops[i]) {
      hepE3.c.report.prob <- pmin(hepE3.c.report.prob * ld.rep.scale, 1)
      hepE3.a.report.prob <- pmin(hepE3.a.report.prob * ld.rep.scale, 1)
    }
  }
  
  # Assign dx status and time -- upon recovery
  ids.newHE3 <- which(pmax(hepE3.stopTime, hepE3.recTime, na.rm = TRUE) == at)
  if (length(ids.newHE3) > 0) {
    dx.prob <- ifelse(hepE3.clin[ids.newHE3] == 1, hepE3.c.report.prob, hepE3.a.report.prob)
    vec.new.HE3 <- runif(length(ids.newHE3)) < dx.prob
    
    ids.newHE3.dx <- ids.newHE3[vec.new.HE3]
    if (length(ids.newHE3.dx) > 0){
      min.time <- hepE3.shedTime[ids.newHE3.dx]
      max.time <- pmax(hepE3.recTime[ids.newHE3.dx], hepE3.stopTime[ids.newHE3.dx], na.rm = TRUE)
      dx.Times <- round(runif(length(ids.newHE3.dx), min = min.time, max = max.time))
      
      for (i in seq_along(dx.Times)) {
        if (dx.Times[i] <= get_control(dat, "nsteps")) {
          nRepHE3 <- get_epi(dat, "nRepHE3", at = dx.Times[i], override.null.error = TRUE)
          nRepHE3.new <- ifelse(is.null(nRepHE3), 1, ifelse(is.na(nRepHE3), 1, nRepHE3 + 1))
          dat <- set_epi(dat, "nRepHE3", at = dx.Times[i],  nRepHE3.new)
        }
      }
    }
    
    # dat <- EpiModel::set_attr(dat, "hepE3.dx", ifelse(vec.new.HE3, 1, 0), ids.newHE3)
    # dat <- EpiModel::set_attr(dat, "hepE3.dxTime", dx.Times, ids.newHE3.dx)
  }
 
  if (at == get_control(dat, "nsteps") && is.null(get_epi(dat, "nRepHE3", override.null.error = TRUE))) {
    dat <- set_epi(dat, "nRepHE3", at = at,  0)
  }
  
  return(dat)
}
