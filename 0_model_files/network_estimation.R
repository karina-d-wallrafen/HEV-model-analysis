
# Setup -------------------------------------------------------------------

suppressMessages(library(EpiModel)) 
suppressMessages(library(socialmixr))
suppressMessages(library(matrixStats))
suppressMessages(library(dplyr))

n <- 1000000
nw <- network_initialize(n)

source("demo_params.R")

# Nodal Attributes -----------------------------------------------------------

# Age
age <- sample(0:100, size = n, prob = (age.weights / sum(age.weights)), replace = TRUE)
age.noise <- runif(n)
age <- age + age.noise
nw <- set_vertex_attribute(nw, attrname = "age", value = age)

age.breaks <- c(0, 18, 65, 101)
age.grp <- cut(age, age.breaks, labels = FALSE, right = FALSE)
nw <- set_vertex_attribute(nw, attrname = "age.grp", value = age.grp)

# Diet
veg <- sample(0:2, size = n, prob = c(0.90, 0.02, 0.08), replace = TRUE) # 0 = omnivore, 1 = pescatarian, 2 = vegetarian or vegan
nw <- set_vertex_attribute(nw, attrname = "veg", value = veg)

# NW Estimation -----------------------------------------------------------

md.cc <- suppressMessages(as.numeric(contact_matrix(polymod, countries = "Germany", age.limits = 0, filter = list(cnt_home = 0, phys_contact = 1))$matrix))

formation.cc <- ~edges
target.stats.cc <- md.cc * n / 2
coef.diss.cc <- dissolution_coefs(dissolution = ~offset(edges), duration = 1)

# set initial conditions such that exp(x) / (1 + exp(x)) = target # of edges / (n * (n-1) / 2)
est.cc <- netest(nw, formation.cc, target.stats.cc, coef.diss.cc, set.control.ergm = control.ergm(MCMLE.maxit = 600, init = -12.6804, MCMLE.density.guard = exp(10)))
dx.cc <- netdx(est.cc, nsims = 1000, dynamic = FALSE, nwstats.formula = ~edges)
print(dx.cc)

est <- list(est.cc)
est <- lapply(est, trim_netest)

save(est, file = "est-1m.RDS")
