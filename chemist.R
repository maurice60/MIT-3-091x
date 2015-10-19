source('constant.R')
# week2
# Coulomb's Law
fCoulomb <- function(q1, q2, r) {
  (q1*q2)/(4*pi*e0*r^2)
}
# Coulomb Potential energy
ePotentialCharge <- function(q1, q2, r) {
  (q1*q2)/(4*pi*e0*r)
}
# Bohr radius energy level n, atomic number z
bRadius <- function(n, Z = 1) n^2*h^2*e0/(pi*Z*em*ec^2)
# Energy of n-th energy level of element z
En <- function(n, z=1) -(z^2/n^2)*k
# Energy diff between initial i and final f energe levels
deltaE <- function(i, f, Z = 1) k*(Z^2)*(1/f^2-1/i^2)
# Velocity of electron energy level n
vLevel <- function(n) n*h/(2*pi*em*bRadius(n))
# Velocity of electron energy level n hit by photon wavelength lambda
vEject <- function(lambda, n) {
  sqrt(2*(Ephot.lambda(lambda) + En(n))/em)
}
# Energy of a photon given wavelength
Ephot.lambda = function(lambda){h*c/lambda}
# Energy of a photon given frequency
Ephot.nu = function(nu) h*nu
# Wavelength of a photon given energy
wavelength <- function(x) h*c/x
# Wavenumber of a photon given energy
wavenumber <- function(x) 1/wavelength(x)
# deBroglie wavelength of particle mass m with given energy
deBroglie <- function(m, energy) {
  p <- sqrt(2*m*energy)
  h/p
}
# Energy required for wavelength lambda
deBroglieEnergy <- function(m, lambda) {
  p <- h/lambda
  p^2/(2*m)
}

# Energy of particle in box, mass m, length l, energy level n
Ebox <- function(m, l, n) (h^2*n^2)/(8*m*l^2)
# Length of box, transition from n1 to n2 with associated energy
lenBox <- function(m, n1, n2, energy) sqrt((h^2*(n1^2-n2^2))/(8*m*energy))

# week 3
# Velocity of an ejected electron wl - wavelength of incident photon, be - binding energy
vElec = function(wl, be) {
  sqrt((2*(Ephot.lambda(wl) - be))/em)}
# average valence electron energy nv - vector of number of electrons per shell, ie - vector of ionisation energy
AVEE <- function(nv, ie){sum(nv*ie)/sum(nv)}
# # Coulombic interaction force
# Cf <- function(q1, q2, r){q1*q2/(4*pi*e0*r*r)}
# # Coulombic interaction energy
# Ce <- function(q1, q2, r){q1*q2/(4*pi*e0*r)}
# Bond energy of diatomic bond, c charge, r1, r2 atomic radii, n from Born experiment
bond.energy <- function(c, r1, r2, n){(-(c*ec)^2/(4*pi*e0*(r1+r2)))*(1-1/n)}
# Cohesive energy of a crystal. qp, qm charge on positive and negative ion. M Madelung constant, r0 distance
# between the centre of two atoms. n from Born experiment.
Ecry <- function(qp, qm, M, r0, n){-(Na*qp*qm*ec^2*M/(4*pi*e0*r0))*(1-1/n)}
nCalc <- function(qp, qm, M, r0, E) 1 / ((4*pi*e0*r0*E)/(Na*qp*qm*ec^2*M) + 1)
# Ionic distance given the lattice energy
ionD <- function(qp, qm, M, ecr, n){-(Na*qp*qm*ec^2*M/(4*pi*e0*ecr))*(1-1/n)}

atRadius <- function(n, m) n^2*h^2*e0/(pi*m*ec^2)

# Week 4
# Polar covalent bond energy Ea Eb homgeneous covalent bond energy, Xa, Xb electronegativity
covalent.energy <- function(Ea, Eb, Xa, Xb) sqrt(Ea*Eb)+96.3*(Xa-Xb)^2
Ea <- function(Et, Eb, Xa, Xb) (Et - 96.3*(Xa-Xb)^2)^2 / Eb
# Percent Ionic Character Xa, Xb electronegativity
ionic.character.pauling <- function(Xa, Xb){100*(1-exp(-(Xa-Xb)^2/4))}
# percent ionic character expermiental bond energy, charge q separation r
ionic.character <- function(pExp, r, q) 100*pExp/(r*q)
# Dipole moment
dipole.moment <- function(r1, r2, e1, e2){abs((r1+r2)*(e1-e2))}

# Week 5
activation.energy <- function(k.1, T.1, k.2, T.2) k.boltz * log(k.1/k.2) / ((1/T.2) - (1/T.1))
k.0 <- function(k.r, E.A, T.r) k.r / exp(-E.A/(k.boltz*T.r))
rate.constant <- function(k_0, E_A, T_K) k_0 * exp(-E_A/(k.boltz*T_K)) 

# Week 6

first.order <- function(C0, t, k) C0*exp(-k*t)
t.elapsed <- function(C, C0, k) (-1/k)*log(C/C0)
# Half life
k.half <- function(t.half) log(2)/t.half

# Second order reactions
rate.2 <- function(k, C.a, C.b) k*C.a*C.b
second.order <- function(C0, t, k) {
  x <- i/C0 + k*t
  1/x
}
t.2.elapsed <- function(C, C0, k) ((1/C)-(1/C0))/k 
# Second order half life
half.life <- function(k, C0) 1/(k*C0) 


