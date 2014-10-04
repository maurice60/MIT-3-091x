h <- 6.626075e-34
c <- 299792458
em <- 9.109e-31
pm <- 1.673e-27
ec <- -1.602e-19
e0 <- 8.854e-12
R <- 10973731.5
eV = 1.602e-19
Na <- 6.022136e23
k <- (em*ec^4)/(8*h^2*e0^2)

# week2
# Energy of a photon given wavelength
Ephot = function(x){h*c/x}
# Wavelength of a photon given energy
wavelength <- function(x){h*c/x}
# Energy of n-th energy level of element z
En <- function(z, n){-(z^2/n^2)*k}
# Wavenumber of transition from i-th to f-th energy level
wavenumber <- function(z, i, f){k*z^2*(1/f^2-1/i^2)/(h*c)}
# deBroglie wavelength of particle mass m with given energy
deBroglie <- function(m, energy){h/(m*sqrt(2*energy/m))}
# Energy of particle in box, mass m, length l, energy level n
Ebox <- function(m, l, n){(h^2*n^2)/(8*m*l^2)}
# Length of box, transition from n1 to n2 with associated energy
lenBox <- function(m, n1, n2, energy){sqrt((h^2*(n1^2-n2^2))/(8*m*energy))}

# week 3
# Velocity of an ejected photon wl - wavelength of incident photon, be - binding energy
vElec = function(wl, be){sqrt((2*(Ephot(wl) - be))/em)}
# average valence electron energy nv - vector of number of electrons per shell, ie - vector of ionisation energy
AVEE <- function(nv, ie){sum(nv*ie)/sum(nv)}
# Coulombic interaction force
Cf <- function(q1, q2, r){q1*q2/(4*pi*e0*r*r)}
# Coulombic interaction energy
Ce <- function(q1, q2, r){q1*q2/(4*pi*e0*r)}
# Bond energy of diatomic bond, c charge, r1, r2 atomic radii, n from Born experiment
be <- function(c, r1, r2, n){(-(c*ec)^2/(4*pi*e0*(r1+r2)))*(1-1/n)}
# Cohesive energy of a crystal. qp, qm charge on positive and negative ion. M Madelung constant, r0 distance
# between the centre of two atoms. n from Born experiment.
Ecry <- function(qp, qm, M, r0, n){-(Na*qp*qm*ec^2*M/(4*pi*e0*r0))*(1-1/n)}
# Ionic distance given the lattice energy
ionD <- function(qp, qm, M, ecr, n){(Na*qp*qm*ec^2*M/(4*pi*e0*ecr))*(1-1/n)}

# Week 4
# Polar covalent bond energy Ea Eb homgeneous covalent bond energy, Xa, Xb electronegativity
Ep <- function(Ea, Eb, Xa, Xb){sqrt(Ea*Eb)+96.3*(Xa-Xb)^2}
# Percent Ionic Character Xa, Xb electronegativity
Ipc <- function(Xa, Xb){100*(1-exp(-(Xa-Xb)^2/4))}
