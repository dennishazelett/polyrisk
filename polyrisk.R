setwd("~/Desktop/projects/polyrisk")
library("rstan")
library("ggplot2")
library("ggridges")
library("tidyr")
library("dplyr")
library("doParallel")
registerDoParallel(cores = 16)
source("https://bioconductor.org/biocLite.R")

genotype_prob <- function(maf) {
  #hardy-weinberg model
  p <- maf
  q <- 1-p
  x <- c(q^2, 2*p*q, p^2)
  x
} 



pt_genotypes <- function(num_pts, snp_ids, my_p, effect_size) {
  # This function simulates a population of patients from pop. variant stats
  # num_pts = number of patients to simulate
  # snp_ids
  # my_p = minor allele prob of snps
  # effect_size = odds ratio (OR)
  genotype_mat <- matrix(NA, length(snp_ids), num_pts)
  for (i in 1:length(snp_ids)) {
    genotype_mat[i,] <- sample(c(0,1,2), size = num_pts, prob = genotype_prob(my_p[i]), replace = T)
  }
  genotype_mat <- data.frame(genotype_mat, maf=my_p, row.names = snp_ids, RR = effect_size)
  genotype_mat
}


# parameter block for simulations
##################
num_pts <- 500000
num_snps <- 20
maf_sd <- 0.3 # standard dev for normal dist of minor allele freq
maf <- signif( abs(rnorm(num_snps, 0, maf_sd)), digits=2)
effect_size <- 1 + abs(rnorm(num_snps, 0, 0.4))
baserate <- 1/100
##################

pts <- pt_genotypes(num_pts, paste("rs", 1:num_snps, sep = ''), 
                    maf,
                    effect_size)


ptodds <- foreach (i=1:num_pts, .combine = c) %dopar% {sum(log(pts$RR ^ pts[,i]))}
ptprobs <- ptodds * baserate
plot(density(ptprobs))
# generate cases 'n ctrls from polygenic risk data
#########################
sick <- numeric(num_pts)
set.seed(4211)
for (i in 1:num_pts) {
  sick[i] <- sample(c(0,1), size=1, prob = c(1-ptprobs[i], ptprobs[i]))
}
cases <- which(sick==TRUE)

length(cases)

ctrls <- sample(which(sick==FALSE), size = length(cases), replace = FALSE)
#ctrls <- which(sick==FALSE)
ptrisk <- data.frame(trtgroup = c(rep("ctrls", length(ctrls)), rep("cases", length(cases))),
                     risk = c(ptprobs[ctrls], ptprobs[cases]))
head(ptrisk)
#########################

OR <- rowSums(pts[,cases])*length(ctrls)/(length(cases)*rowSums(pts[,ctrls]))
with(ptrisk, plot(risk ~ as.factor(trtgroup)))

res <- data.frame(OR=OR, RR=pts$RR)
ggplot(data=res, aes(x=OR, y=RR)) + 
  geom_point() + xlim(1,2) + ylim(1,2)

# violin plot compare distributin relative risk of cases controls
ggplot(data=ptrisk, aes(x=trtgroup, y=risk/baserate, colour=trtgroup, fill=trtgroup)) + geom_violin() 
ggplot(data=ptrisk, aes(y=trtgroup, x=risk/baserate, colour=trtgroup, fill=trtgroup)) + geom_density_ridges() 


sim_data1 <- list(nsnps = num_snps,
                  npts  = 2*length(cases),
                  genotypes = pts[,c(cases, ctrls)],
                  sick  = sick[c(cases, ctrls)])

additive_model = "
functions {
  real foo_lpdf(real x) {
    return -2*log(x);
  }
}
data {
  int<lower=1> nsnps;
  int<lower=1> npts;
  matrix<lower=0, upper=2>[nsnps, npts] genotypes;
  int<lower=0, upper=1> sick[npts];
}
parameters{
  real<lower=1, upper=3> risk_ratio[nsnps];
  real<lower=0, upper=1> baserate; 
}
//transformed parameters {
//  real<lower=0> ptrisk[npts];
//  row_vector[nsnps] snprisk[npts];
//  for (i in 1:nsnps) {
//    for (j in 1:npts) {
//      snprisk[j, i] = risk_ratio[i] ^ genotypes[i, j];
//    }
//  }
//  for (j in 1:npts) {
//    ptrisk[j] = log_sum_exp(snprisk[j]);
//  }
//}
model {
  // prior
  for (i in 1:nsnps) {
    risk_ratio[i] ~ foo();
  }
  baserate ~ beta(1, 100);
  // likelihood
  for (i in 1:npts) {
    //sick[i] ~ bernoulli(ptrisk[i] * baserate);
    for (j in 1:nsnps) {
      sick[i] ~ bernoulli(risk_ratio[j] ^ genotypes[j, i] * baserate);
    }
  }
}
"



#plot(density(rbeta(10000, 1, 100)), xlim=c(0,1))
inits1 <- list(list(odds_ratio=rep(1.1, 10), baserate=0.001))

inits4 <- list(list(odds_ratio=rep(1.1, 10), baserate=0.001),
               list(odds_ratio=rep(1.1, 10), baserate=0.001),
               list(odds_ratio=rep(1.1, 10), baserate=0.001),
               list(odds_ratio=rep(1.1, 10), baserate=0.001))

fit1 <- stan_model(model_code = additive_model)

gc()


a <- cbind(pts[,c("maf", "OR")], pred_RR=colMeans(foo$risk_ratio))
a <- mutate(a, err = OR-pred_RR)
a$joint_model <- colMeans(foo1$risk_ratio)
a <- mutate(a, err_joint = OR-joint_model)
plot(density(a$err_joint))

f1 <- stan(fit = fit1, 
           model_code = additive_model, 
           data = sim_data1,
           init = inits4,
           iter = 10000,
           warmup = 5000,
           thin = 50,
           chains = 4,
           control = list(adapt_delta=0.94,
                          max_treedepth=12))


foo1 <- rstan::extract(f1)
#bar <- get_inits(f1)
summary(foo)
save(foo, foo1, pts, ptodds, file="joint_posterior.rda")

plot(x=colMeans(foo$risk_ratio), y=pts$RR, col="red")#, xlim=c(1,3), ylim=c(1,3))
points(x=colMeans(foo$risk_ratio), y=pts$RR, col="blue")
plot(density(foo$baserate))
(1.1^-2)^10


rrs <- gather(as.data.frame(foo$risk_ratio))
ggplot(data=rrs, aes(x=value, y=key)) + 
  geom_density_ridges() 


################### OVARIAN CANCER
ovc_snps <- c("rs7651446", "rs3814113", "rs10088218", "rs9303542", "rs8170", "rs2072590", "rs34289250", 
              "rs11782652", "rs1243180", "rs11907546", "rs2046210", "rs10069690")



BRCA1
BRCA2
RAD51C
RAD51D
BARD1
BRIP1
PALB2
FANCM
CHEK2
ATM

ovc_or   <- c(1.59, 1.2658228, 1.2987013, 1.14, 1.19, 1.14, 7.95, 1.24, 1.1, 1.1111112, 1.28, 1.14)
ovc_maf  <- c(0.05, 0.32, 0.13, 0.27, 0.19, 0.68, 0.0089, 0.07, 0.31, .35903, 0.08, 0.26)
baserate <- 0.013

ovcpts <- pt_genotypes(num_pts, ovc_snps, ovc_maf, effect_size = ovc_or)
ovcodds <- exp(colSums(log(ovcpts$OR ^ ovcpts[,1:num_pts])))

popscale.factor <- baserate / mean(ovcodds)
ovcprobs <- ovcodds * popscale.factor
plot(density(ovcprobs))
valid.patients <- which(ovcprobs >=1 )


plot(density(ovcprobs))
which(ovcprobs.age>1)


#FUN=function(x) {baserate*prod(pts$OR^x)}
# generate cases 'n ctrls from polygenic risk data
#########################
set.seed(2000)
ovcsick <- numeric(length(valid.patients))
for (i in 1:length(valid.patients)) {
  ovcsick[i] <- sample(c(0,1), size=1, prob = c(1-ovcprobs.age[i], ovcprobs.age[i]))
}
cases <- which(ovcsick==TRUE)
length(cases)
ctrls <- sample(which(ovcsick==FALSE), size = length(cases), replace = FALSE)
ptrisk <- melt(data.frame(ctrls=ovcprobs.age[ctrls], cases=ovcprobs.age[cases]))
names(ptrisk) <- c("trtgroup", "risk")
#########################

# violin plot compare distributin relative risk of cases controls
ggplot(data=ptrisk, aes(x=trtgroup, y=risk/baserate, colour=trtgroup, fill=trtgroup)) + geom_violin() 
ggplot(data=data.frame(or=ovcprobs.age[cases][order(ovcprobs.age[cases])]/ovcprobs.age[ctrls][order(ovcprobs.age[ctrls])]), aes(x=or)) + geom_density(fill="grey", colour="grey")
ggplot(data=data.frame(cases=ovcprobs.age[cases][order(ovcprobs.age[cases])]/baserate, ctrls=ovcprobs.age[ctrls][order(ovcprobs.age[ctrls])]/baserate), aes(x=ctrls, y=cases)) + 
  geom_point() + #xlim(0, 20) + ylim(0, 200) + 
  geom_abline(intercept = 0, slope = 1, colour="red") +
  geom_hline(yintercept=quantile(ovcprobs.age[cases]/baserate, probs=c(.5,.95,.99)), colour="grey", lty=2) +
  geom_vline(xintercept=quantile(ovcprobs.age[ctrls]/baserate, probs=c(.5,.95,.99)), colour="grey", lty=2)

caseOR <- ovcprobs[cases][order(ovcprobs[cases])]/baserate
ctrlOR <- ovcprobs[ctrls][order(ovcprobs[ctrls])]/baserate
# % of cases OR > 5
length(which(caseOR>5))
length(which(caseOR>5))/length(caseOR)
# slope
mean(caseOR/ctrlOR)
save(pts, yrs, file = "pts.RDa")


################### BRCA Muts
prevalence <- 0.013

##                  1              2             3             4              5            6             7              8
ovc_snps <- c("rs56318008",  "rs58722170", "rs17041869", "rs2072590",   "rs2165109", "rs711830",   "rs752590",    "rs112071820", # A
              "rs7651446",   "rs9870207",  "rs13113999", "rs17329882",  "rs4691139", "rs10069690", "rs555025179", "rs6456822",   # B 
              "rs7705526",   "rs10088218", "rs11782652", "rs150293538", "rs9886651", "rs1413299",  "rs200182588", "rs320203",    # C
              "rs3814113",   "rs635634",   "rs1192691",  "rs1243180",   "rs7902587", "rs7937840",  "rs7953249",   "rs8037137",   # D 
              "rs199661266", "rs11651755", "rs183211",   "rs757210",    "rs9303542", "rs8098244",  "rs1469713",   "rs4808075",   # E
              "rs688187",    "rs6005807",  "brca1",      "brca2",       "brip1",     "palb2",      "rad51c",      "rad51d",      # F 
              "bard1")#,       "chek2",      "atm",        "nbn",         "tp53",      "rad50",      "fam174a",     "mre11a")      # G

##              1     2     3     4     5     6     7     8
ovc_or   <- c(1.11, 1.08, 1.06, 1.34, 1.09, 1.14,  1.3, 1.29, # A
              1.59, 1.19,  1.2, 1.23, 1.09, 1.14, 1.51, 1.18, # B
              1.07, 1.24, 2.19, 1.23, 1.08, 1.21, 1.53, 1.04, # C
              1.29, 1.11, 1.11,  1.1, 1.29, 1.05, 1.08, 1.07, # D 
              1.09,  1.2, 1.11, 1.11, 1.14, 1.19, 1.04, 1.18, # E
              1.33, 1.17,   29, 12.7,  6.4,  4.4,  3.4, 10.9, # F 
              4.2)#,   0.4,  2.4,  2.3,  2.9,  0.7,  1.9,  0.8) # G

##                1         2         3         4         5         6         7         8
ovc_maf  <- c(0.796725, 0.151757, 0.874002, 0.817492, 0.740016, 0.8185,   0.763778, 0.73722,  # A
              0.939896, 0.487021, 0.226438, 0.848243, 0.363219, 0.652356, 0.54173,  0.4623,   # B 
              0.321486, 0.086661, 0.945487, 0.992612, 0.725839, 0.580072, 0.4868,   0.176717, # C
              0.44389,  0.139776, 0.758387, 0.159545, 0.91893,  0.857628, 0.56889,  0.282348, # D 
              0.8676,   0.467252, 0.287784, 0.3622,   0.684904, 0.173922, 0.478035, 0.154752, # F 
              0.358427, 0.886981, 0.0031,   0.0041,   0.0017,   0.0011,   0.0011,   0.0004,   # F
              0.0005)#,   0.0082,   0.0022,   0.0014,   0.0011,   0.0024,   0.0008,   0.0007)   # G


ovcpts <- pt_genotypes(num_pts, ovc_snps, ovc_maf, effect_size = ovc_or)
ovcodds <- exp(colSums(log(ovcpts$OR ^ ovcpts[,1:num_pts])))
gc()

brca.carrier <- rep("NON-BRCA", num_pts)
brca.carrier[which(ovcpts[14,1:1e5]>0)] <- "BRCA2"
brca.carrier[which(ovcpts[13,1:1e5]>0)] <- "BRCA1" # assign brca1 last in case of overwriting double mutation status

set.seed(24)
age.scale <- function(age, onset.age, risk.age) {
  # age is the patient age
  # onset.age is the age at which the process begins to accelerate
  # risk.age is the age at which lifetime risk accumulates to its maximum
  1/risk.age * age^((age-onset.age)/(risk.age-onset.age))
}

ovcprobs.brca <- ovcodds 
plot(density(ovcodds), xlim=c(0,1e5))
plot(density(log10(ovcodds)))
abline(v = quantile(log10(ovcodds), probs = 0.99), lty=2, col="grey")
10^(mean(log10(ovcodds)))
?abline
popscale.factor <- prevalence / mean(ovcprobs.brca)
popscale.factor <- prevalence / exp(mean(log(ovcprobs.brca)))
popscale.factor <- prevalence / 4000
valid.patients <- which(ovcprobs.brca*popscale.factor < 1)

nvalid <- length(valid.patients)

ages <- c(rep(20, nvalid), 
          rep(30, nvalid),
          rep(40, nvalid),
          rep(50, nvalid),
          rep(60, nvalid),
          rep(70, nvalid),
          rep(80, nvalid))

agerisk <- data.frame(age      = as.factor(ages), 
                      risk     = rep(ovcprobs.brca[valid.patients], 7),
                      brca     = as.factor(rep(brca.carrier[valid.patients], 7)),
#                     risk.age = rep(ovcprobs.brca[valid.patients], 7) * age.scale(ages, 20, 80) * popscale.factor) 
                      risk.age = rep(ovcprobs.brca[valid.patients], 7) * inv_logit(ages, shift = 8, slope = 1/7) * popscale.factor)

gc()


ggplot(agerisk, aes(x = age, y = risk.age)) + geom_boxplot(aes(fill = brca), notch = TRUE) + labs(y="relative risk (RR)") 

# generate cases 'n ctrls from polygenic risk data
#########################
set.seed(29)
ages <- sample(40:80, size=num_pts, replace = TRUE)
probs <- ovcprobs.brca * inv_logit(ages, shift = 8, slope = 1/7)
probs <- probs[valid.patients]

ovcsick.brca <- matrix(nrow = nvalid, ncol = 100)

for (i in 1:nvalid) {
  ovcsick.brca[i,] <- sample(c(0,1), size=10, prob = c(1-popscale.factor*probs[i], popscale.factor*probs[i]), replace = TRUE)
}

pick.one <- sample(10, 1)

cases <- which(ovcsick.brca[,pick.one]==TRUE)
length(cases)
length(cases)/nvalid
probs.scale <- probs * popscale.factor

ctrls <- sample(which(ovcsick.brca[,pick.one]==FALSE), size = length(cases), replace = FALSE)
ptrisk <- melt(data.frame(ctrls=probs.scale[ctrls], cases=probs.scale[cases]))
names(ptrisk) <- c("trtgroup", "risk")
#########################


length(probs.scale)
mean(probs.scale[(brca.carrier=="BRCA1")[valid.patients]==TRUE])
mean(probs.scale[(brca.carrier=="BRCA2")[valid.patients]==TRUE])
mean(probs.scale[(brca.carrier=="NON-BRCA")[valid.patients]==TRUE])
# of cases that are BRCA carriers (see double negative?)
sum(cases %in% which((brca.carrier=="NON-BRCA")[valid.patients]==FALSE))
# fraction of cases that are BRCA1+
sum(cases %in% which((brca.carrier=="BRCA1")[valid.patients]==TRUE))/length(cases)
sum(cases %in% which((brca.carrier=="BRCA2")[valid.patients]==TRUE))/length(cases)
sum(cases %in% which((brca.carrier=="NON-BRCA")[valid.patients]==TRUE))/length(cases)

# violin plot compare distributin relative risk of cases controls
ggplot(data=ptrisk, aes(x=trtgroup, y=risk/baserate, colour=trtgroup, fill=trtgroup)) + geom_violin() + labs(x = "patient population", y = "Relative Risk (RR)")
#ggplot(data=data.frame(or=ovcprobs.brca.age[cases][order(ovcprobs.brca.age[cases])]/ovcprobs.brca.age[ctrls][order(ovcprobs.brca.age[ctrls])]), aes(x=or)) + geom_density(fill="grey", colour="grey")
ggplot(data=data.frame(cases=probs.scale[cases][order(probs.scale[cases])]/prevalence, 
                       ctrls=probs.scale[ctrls][order(probs.scale[ctrls])]/prevalence), aes(x=ctrls, y=cases)) + 
  geom_point() + #xlim(0, 20) + ylim(0, 200) + 
  geom_abline(intercept = 0, slope = 1, colour="red") +
  geom_hline(yintercept=quantile(probs.scale[cases]/prevalence, probs=c(.5,.95,.99)), colour="grey", lty=2) +
  geom_vline(xintercept=quantile(probs.scale[ctrls]/prevalence, probs=c(.5,.95,.99)), colour="grey", lty=2)


caseOR <- probs.scale[cases][order(probs.scale[cases])]/prevalence
ctrlOR <- probs.scale[ctrls][order(probs.scale[ctrls])]/prevalence
# % of cases OR > 5
length(which(caseOR>5))
length(which(caseOR>5))/length(caseOR)
# slope
mean(caseOR/ctrlOR)

# prevalence
length(cases)/nvalid
# fraction of brca1 patients that are cases
sum(cases %in% which(brca.carrier[valid.patients]=="BRCA1"))/sum(brca.carrier[valid.patients]=="BRCA1")
# fraction of BRCA2 patients that are cases
sum(cases %in% which(brca.carrier[valid.patients]=="BRCA2"))/sum(brca.carrier[valid.patients]=="BRCA2")
# fraction of cases that are brca1+
sum(cases %in% which(brca.carrier[valid.patients]=="BRCA1"))/length(cases)
# fraction of cases that are non-BRCA
sum(cases %in% which(brca.carrier[valid.patients]=="NON-BRCA"))/length(cases)


########## Tissue specific compartmentalization of risk
##########
#
# Suppose risk is subdivided in different tissues/organs
# Suppose tissue of origin (TO) and second site (SS1, SS2, ...)
# Suppose the SS risk is dependent on TO, how does that look mathematically?
# TO effective (TO_eff) = 
#
# SS_eff <- inv_logit(TO) * SS
#
# Questions, would data still fit independent risk model, but with larger variances for TO? 
# Or something different?
#
# I need to generate weights for every SNP in every tissue/organ system
# To keep it simple, let's suppose there are 4 categories, and two of them have signal associated with some SNPs, but are not truly involved

# tissue weights : (from funciVAR overlaps, just making these up right now)

#          1     2     3     4     5     6     7     8     9    10    11    12
tw1 <- c(0.00, 0.23, 0.42, 0.56, 0.00, 0.00, 0.02, 0.00, 0.00, 0.55, 0.00, 0.00)
tw2 <- c(0.00, 0.51, 0.72, 0.90, 0.00, 0.10, 0.00, 0.00, 0.00, 0.34, 0.00, 0.00)
tw3 <- c(0.41, 0.00, 0.00, 0.12, 0.53, 0.62, 0.00, 0.23, 0.05, 0.69, 0.81, 0.22)
tw4 <- c(0.01, 0.00, 0.00, 0.05, 0.00, 0.00, 0.00, 0.05, 0.00, 0.41, 0.00, 0.03)

# causal tissue (hidden; used for simulation)
ct1 <- NULL
ct2 <- c(2, 3, 4)
ct3 <- c(1, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)
ct4 <- NULL

library(boot)
plot(x=1:100, y=inv_logit(1:100, shift = 8, slope = 1/7), type='l')

inv_logit <- function(x, scale = 1, slope = 1, shift = 0) {
  x <- scale * exp(x*slope - shift)/(1+exp(x * slope - shift))
  for (i in 1:length(x)) {
    if (is.na(x[i])) {x[i] <- 1}
  }
  x
}
TO[872]
plot(x=-1:100, y=inv_logit(-1:100, slope= 1/1.2, shift = 5), type="l")
inv_logit(100000, shift = 8)
ovcpts[13,1e5+2]

#####
SS <- exp(colSums(log(ovcpts$OR[ct2] ^ ovcpts[,1:num_pts])))
#plot(density(SS))
TO <- exp(colSums(log(ovcpts$OR[ct3] ^ ovcpts[,1:num_pts])))
#plot(density(TO), xlim=c(0,100000))
SS_eff <- SS^inv_logit(TO, slope = 1, shift = 8)
#OR.mix <- TO^inv_logit(SS_eff, slope=1/1.2, shift=5)
OR.mix <- SS_eff * TO

#OR.flat <- SS * TO
plot(density(SS_eff), col="red", lty=2)
lines(density(SS))


baserate <- 0.013
popscale.factor <- baserate/exp(mean(log(OR.mix)))

#popscale.factor <- baserate/median(OR.flat)
valid.patients <- which(popscale.factor * OR.mix < 1)
#valid.patients <- which(popscale.factor * OR.flat < 1)
ovcprobs.w8d.age <- popscale.factor * OR.mix * age.scale(ages, 20, 80)
#ovcprobs.w8d.age <- popscale.factor * OR.flat * age.scale(ages, 20, 80)

ovcprobs.w8d.age <- ovcprobs.w8d.age[valid.patients]
plot(density(ovcprobs.w8d.age))
which(ovcprobs.w8d.age>1)

#########################
ovcsick.w8d <- numeric(length(valid.patients))

for (i in 1:length(valid.patients)) {
  ovcsick.w8d[i] <- sample(c(0,1), size=1, prob = c(1-ovcprobs.w8d.age[i], ovcprobs.w8d.age[i]), replace = T)
}


cases <- which(ovcsick.w8d==TRUE)
length(cases)
ctrls <- sample(which(ovcsick.w8d==FALSE), size = length(cases), replace = FALSE)
ptrisk <- melt(data.frame(ctrls=ovcprobs.w8d.age[ctrls], cases=ovcprobs.w8d.age[cases]))
names(ptrisk) <- c("trtgroup", "risk")
#########################

#ggplot(data=ptrisk, aes(x=trtgroup, y=risk/baserate, colour=trtgroup, fill=trtgroup)) + geom_violin() 
#ggplot(data=data.frame(or=ovcprobs.brca.age[cases][order(ovcprobs.brca.age[cases])]/ovcprobs.brca.age[ctrls][order(ovcprobs.brca.age[ctrls])]), aes(x=or)) + geom_density(fill="grey", colour="grey")
ggplot(data=data.frame(cases=ovcprobs.w8d.age[cases][order(ovcprobs.w8d.age[cases])]/baserate, ctrls=ovcprobs.w8d.age[ctrls][order(ovcprobs.w8d.age[ctrls])]/baserate), aes(x=ctrls, y=cases)) + 
  geom_point() + #xlim(0, 20) + ylim(0, 200) + 
  geom_abline(intercept = 0, slope = 1, colour="red") +
  geom_hline(yintercept=quantile(ovcprobs.w8d.age[cases]/baserate, probs=c(.5,.95,.99)), colour="grey", lty=2) +
  geom_vline(xintercept=quantile(ovcprobs.w8d.age[ctrls]/baserate, probs=c(.5,.95,.99)), colour="grey", lty=2)


