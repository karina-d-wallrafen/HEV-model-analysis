
# Setup -------------------------------------------------------------------

# Load packages
# library(EpiModel)
# suppressMessages(library(dplyr))
# suppressMessages(library(tidyr))

# Retrieve attributes
age.grp <- est[[1]]$newnetwork$attr$age.grp
# n <- est[[1]]$newnetwork$gal$n

# Set number of households
hh.size <- 1.8
n.hh <- round(n / hh.size)
ids <- 1:n
hh.ids <- 1:n.hh

# Create empty data frames to track hh assignment
persons.by.hh <- data.frame(ids = ids, age.grp = age.grp, hh = NA)

# Set rules
prop.hh.with.child <- 0.172
prop.single.parent <- 0.152 
prop.one.person <- 0.549 

# Household Assignment ----------------------------------------------------

# Determine which households will have a member under 18 
hh.u18 <- sample(x = hh.ids, size = round(prop.hh.with.child * n.hh))

# Assign children under 18 to the selected households
persons.by.hh[persons.by.hh$ids %in% which(age.grp == 1)[1:length(hh.u18)], ]$hh <- hh.u18 
children.hh.assign <- sample(x = hh.u18, size = length(which(age.grp == 1)) - length(hh.u18), replace = TRUE)
persons.by.hh[persons.by.hh$ids %in% which(age.grp == 1)[(length(hh.u18) + 1):length(which(age.grp == 1))], ]$hh <- children.hh.assign 

# Determine which hh with children will be single-parent households
hh.single.parent <- sample(hh.u18, size = prop.single.parent * length(hh.u18), replace = FALSE)
hh.2p.parent <- setdiff(hh.u18, hh.single.parent)

# Determine which households will have only a single (adult) person
hh.18p <- setdiff(hh.ids, hh.u18)
hh.1person <- sample(hh.18p, size = round(prop.one.person * n.hh))

# Assign one 18-to-64 aged person to each household with a child
persons.by.hh[persons.by.hh$ids %in% which(age.grp == 2)[1:length(hh.u18)], ]$hh <- hh.u18 

# Assign one 18+ person to each household without a child
persons.by.hh[which(persons.by.hh$ids %in% which(age.grp %in% c(2, 3)) & is.na(persons.by.hh$hh))[1:length(hh.18p)], ]$hh <- hh.18p

# Add one or more adult to each 2-parent household and each household without a child that is not a one-person household
hh.add.adult <- setdiff(hh.ids, c(hh.single.parent, hh.1person))
persons.by.hh[which(persons.by.hh$ids %in% which(age.grp %in% c(2, 3)) & is.na(persons.by.hh$hh))[1:length(hh.add.adult)], ]$hh <- hh.add.adult

# Distribute remaining adults among households that are not single-parent or single-person households
adults.hh.assign <- sample(x = hh.add.adult, size = sum(is.na(persons.by.hh$hh)), replace = TRUE)
persons.by.hh[which(persons.by.hh$ids %in% which(age.grp %in% c(2, 3)) & is.na(persons.by.hh$hh)), ]$hh <- adults.hh.assign

# Check Household Assignment ----------------------------------------------

# # Rule 1: Expected proportion of households have at least one child
# hh.check1 <- persons.by.hh %>% group_by(age.grp) %>% summarize(num.hh = n_distinct(hh))
# hh.check1$pct.hh <- round(hh.check1$num.hh / length(unique(persons.by.hh$hh)), 3) * 100
# 
# # Rule 2: Average household size is as expected
# hh.check2 <- persons.by.hh %>% group_by(hh) %>% summarize(num.person = n())
# mean.hh.check2 <- round(mean(hh.check2$num.person), 1)
# 
# # Rule 3. Every household with a child must also have at least one adult in the 18-to-64 age range
# hh.check3 <- persons.by.hh %>% group_by(hh, age.grp) %>% summarize(num.person = n())
# hh.check3 <- spread(hh.check3, key = age.grp, value = num.person)
# hh.check3 <- hh.check3[!is.na(hh.check3$`1`) & !is.na(hh.check3$`2`), ]
# prop.child.w.adult <- sum(hh.check3$`1`) / length(which(age.grp == 1))
# 
# # Rule 4. Expected proportion of households are single-person households
# hh.check4 <- persons.by.hh %>% group_by(hh) %>% summarize(num.person = n())
# one.person.hhs <- round(nrow(hh.check4 %>% filter(num.person == 1)) / n.hh, 3)
# 
# # Rule 5. Expected proportion of households with children are single-adult households
# hh.check5 <- persons.by.hh %>% filter(hh %in% hh.u18) %>% filter(age.grp %in% c(2, 3)) %>% group_by(hh) %>% summarize(num.adults = n())
# single.parent.hh <- round(nrow(hh.check5[hh.check5$num.adults == 1, ]) / nrow(hh.check5), 3)

# Prepare household network ----------------------------------------------

# Set household attribute
est[[1]]$newnetwork <- set.vertex.attribute(est[[1]]$newnetwork, "household", persons.by.hh$hh)

# Save household edge list
hhPairs <- merge(persons.by.hh, persons.by.hh, by = "hh")
hhPairs <- subset(hhPairs, (ids.x < ids.y))
hhPairs <- hhPairs[, c(2, 4)]
hhPairs <- hhPairs[order(hhPairs[,1], hhPairs[,2]), ]
names(hhPairs) <- c(".head", ".tail")
attr(hhPairs, "n") <- n
