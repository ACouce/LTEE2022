# .. This script fits several common distributions to beneficial tails of DFE 

# .. load a useful library (to use fitdist(), which implements MLE and GoF analyses for all the desired distributions) 
library(fitdistrplus)

# .. shift values to have a zero class for the exponential fit
r <- R_beneficial$fitted1-thben
m <- T_beneficial$fitted1-thben
k <- F_beneficial$fitted1-thben

# .. maximum-likelihood fitting of Anc data to all the desired distributions
Rexp <- fitdist(r, "exp")
Rgam <- fitdist(r, "gamma")
Rlno <- fitdist(r, "lnorm")
Rwei <- fitdist(r, "weibull")
Rnor <- fitdist(r, "norm")
Rlog <- fitdist(r, "logis")

# .. print p-values
print(c("pexp",ks.test(r, "pexp", Rexp$estimate)$p.value))
print(c("pweibull",ks.test(r, "pweibull", shape=Rwei$estimate[1], scale=Rwei$estimate[2])$p.value))
print(c("pgamma",ks.test(r, "pgamma", shape=Rgam$estimate[1], rate=Rgam$estimate[2])$p.value))
print(c("plnorm",ks.test(r, "plnorm", meanlog=Rlno$estimate[1], sdlog=Rlno$estimate[2])$p.value))
print(c("pnorm",ks.test(r, "pnorm", mean=Rnor$estimate[1], sd=Rnor$estimate[2])$p.value))
print(c("plogis",ks.test(r, "plogis", location=Rlog$estimate[1], scale=Rlog$estimate[2] )$p.value))

# .. print goodness-of-fit statistics
print(gofstat(list(Rexp,Rgam, Rlno, Rwei, Rnor, Rlog), fitnames = c("Exponential", "Gamma", "LogNorm", "Weibull", "Normal", "Logistic")))


# .. maximum-likelihood fitting of 2K data to all the desired distributions
Mexp <- fitdist(m, "exp")
Mgam <- fitdist(m, "gamma")
Mlno <- fitdist(m, "lnorm")
Mwei <- fitdist(m, "weibull")
Mnor <- fitdist(m, "norm")
Mlog <- fitdist(m, "logis")

# .. print p-values
print(c("pexp",ks.test(m, "pexp", Mexp$estimate)$p.value))
print(c("pweibull",ks.test(m, "pweibull", shape=Mwei$estimate[1], scale=Mwei$estimate[2])$p.value))
print(c("pgamma",ks.test(m, "pgamma", shape=Mgam$estimate[1], rate=Mgam$estimate[2])$p.value))
print(c("plnorm",ks.test(m, "plnorm", meanlog=Mlno$estimate[1], sdlog=Mlno$estimate[2])$p.value))
print(c("pnorm",ks.test(m, "pnorm", mean=Mnor$estimate[1], sd=Mnor$estimate[2])$p.value))
print(c("plogis",ks.test(m, "plogis", location=Mlog$estimate[1], scale=Mlog$estimate[2] )$p.value))

# .. print goodness-of-fit statistics
print(gofstat(list(Mexp,Mgam, Mlno, Mwei, Mnor, Mlog), fitnames = c("Exponential", "Gamma", "LogNorm", "Weibull", "Normal", "Logistic")))

# .. maximum-likelihood fitting of 15K data to all the desired distributions
Kexp <- fitdist(k, "exp")
Kgam <- fitdist(k, "gamma")
Klno <- fitdist(k, "lnorm")
Kwei <- fitdist(k, "weibull")
Knor <- fitdist(k, "norm")
Klog <- fitdist(k, "logis")

# .. print p-values
print(c("pexp",ks.test(k, "pexp", Kexp$estimate)$p.value))
print(c("pweibull",ks.test(k, "pweibull", shape=Kwei$estimate[1], scale=Kwei$estimate[2])$p.value))
print(c("pgamma",ks.test(k, "pgamma", shape=Kgam$estimate[1], rate=Kgam$estimate[2])$p.value))
print(c("plnorm",ks.test(k, "plnorm", meanlog=Klno$estimate[1], sdlog=Klno$estimate[2])$p.value))
print(c("pnorm",ks.test(k, "pnorm", mean=Knor$estimate[1], sd=Knor$estimate[2])$p.value))
print(c("plogis",ks.test(k, "plogis", location=Klog$estimate[1], scale=Klog$estimate[2] )$p.value))

# .. print goodness-of-fit statistics
print(gofstat(list(Kexp,Kgam, Klno, Kwei, Knor, Klog), fitnames = c("Exponential", "Gamma", "LogNorm", "Weibull", "Normal", "Logistic")))
