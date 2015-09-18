h <- 6.626075e-34
c <- 299792458
em <- 9.109e-31
pm <- 1.673e-27
ec <- -1.602e-19
e0 <- 8.854e-12 # permittivity of free space
R <- 10973731.5
eV = 1.602e-19
Na <- 6.022136e23
k <- (em*ec^4)/(8*h^2*e0^2)
bohrRadius <- 5.29e-11

nMol <- function(mass_g, atomic) mass_g / atomic 
atomicWeight <- function(weight, abundance) {
  sum(weight*abundance)
}
