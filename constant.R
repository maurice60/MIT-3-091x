source('maths.R')
elements <- read.table('AtomicWeights.txt')
names(elements) <- c('Z', 'Symbol', 'Name', 'W')
h <- 6.626075e-34
c <- 299792458
em <- 9.109e-31
pm <- 1.673e-27
ec <- -1.602e-19
e0 <- 8.854e-12 # permittivity of free space
R <- 10973731.5
eV = -ec
Na <- 6.022136e23
k <- (em*ec^4)/(8*h^2*e0^2)
bohrRadius <- 5.29e-11

nMol <- function(mass_g, atomic) mass_g / atomic

atomicWeightAve <- function(weight, abundance) {
  sum(weight*abundance)
}

atomicWeight <- function(symbol = '', Z = 0) {
  ifelse(Z == 0, elements[elements$Symbol==symbol,'W'], elements[elements$Z == Z,'W'])
}

gramMolecularWeight <- function(form) {
  gmw <- function(element, number) {
    Z <- sapply(element, atomicWeight)
    sum(Z*number)
  }
  upper <- c('A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z')
  lower <- c('a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z')
  digit <- c('0','1','2','3','4','5','6','7','8','9')
  ele <- c()
  dig <- c()
  sts <- unlist(strsplit(form, ""))
  nch <- 0
  for (ch in sts) {
    if (ch %in% upper) {
      nch <- nch + 1
      ele[nch] <- ch
      dig[nch] = 0
    }
    if (ch %in% lower) {
      cap <- ele[nch]
      ele[nch] <- paste(cap, ch, sep = "")}
    if (ch %in% digit) {
      dig[nch] <- dig[nch]*10 + as.integer(ch)
    }
  }
  dig <- replace(dig, dig==0, 1)
  gmw(ele, dig)
}
