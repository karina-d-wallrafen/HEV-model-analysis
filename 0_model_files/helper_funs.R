as_tibble_edgelist <- function(el) {
  
  if (nrow(el) > 0) {
    t_el <- tibble::tibble(head = el[, 1], tail = el[, 2])
  } else {
    t_el <- tibble::tibble(head = integer(0), tail = integer(0))
  }
  
  return(t_el)
}

get_discordant_edgelist <- function(dat, infs, status, sus.status, network) {

  if (network == "community") {
    d_el <- get_edgelist(dat, 1)
    d_el <- as_tibble_edgelist(d_el)
  } else if (network == "household") {
    d_el <- tibble::tibble(head = dat$param$hh.pairs$.head, tail = dat$param$hh.pairs$.tail)
  } else {
    stop("There is no ", network, " network to get an edgelist from")
  }

  d_el_ordered <- d_el[d_el$head %in% infs, , drop = FALSE]
  
  d_el_rev <- d_el[d_el$tail %in% infs, , drop = FALSE] |> dplyr::select(head = "tail", tail = "head")
  
  d_el <- dplyr::bind_rows(d_el_ordered, d_el_rev)
  
  tail.status <- get_attr(dat, status, d_el$tail)
  
  d_el <- d_el[which(tail.status %in% sus.status), , drop = FALSE]
  
  return(d_el)
  
}

update_cumulative_edgelist <- function(dat, network, truncate = 0) {
  
  if (!get_control(dat, "cumulative.edgelist")) {
    return(dat)
  }
  
  el <- get_edgelist(dat, network)
  el_cuml <- get_cumulative_edgelist(dat, network)
  
  el <- tibble::tibble(
    head = get_unique_ids(dat, el[, 1]),
    tail = get_unique_ids(dat, el[, 2]),
    current = TRUE
  )
  
  el_cuml <- dplyr::full_join(el_cuml, el, by = c("head", "tail"))
  
  at <- get_current_timestep(dat)
  
  new_edges <- is.na(el_cuml[["start"]])
  if (any(new_edges)) {
    el_cuml[new_edges, ][["start"]] <- at
  }
  
  terminated_edges <- is.na(el_cuml[["current"]]) & is.na(el_cuml[["stop"]])
  if (any(terminated_edges)) {
    el_cuml[terminated_edges, ][["stop"]] <- at - 1
  }
  
  if (truncate != Inf) {
    rel.age <- at - el_cuml[["stop"]]
    rel.age <- ifelse(is.na(rel.age), 0, rel.age)
    el_cuml <- el_cuml[rel.age <= truncate, ]
  }
  
  dat[["el.cuml"]][[network]] <- el_cuml[, c("head", "tail", "start", "stop")]
  
  return(dat)
}