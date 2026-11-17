targstep <- 05
targbeg <- -60 - targstep/2
targend <- 100

targage <- seq(targbeg, targend, by=targstep)

test_age <- seq(targbeg, targend, by=0.5)
test_binnum <- as.integer(ceiling((test_age-targbeg-(targstep/2))/targstep))+1
cbind(test_age, test_binnum, targage[test_binnum], 1950 - test_age, 1950 - targage[test_binnum])
