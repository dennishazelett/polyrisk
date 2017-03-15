setwd("~/Desktop/projects/polyrisk")
library("rstan")
library("ggplot2")
library("reshape2")



genotype_prob <- function(maf) {
  #hardy-weinberg model
  p <- maf
  q <- 1-p
  x <- c(q^2, 2*p*q, p^2)
  x
} 



pt_genotypes <- function(num_pts, snp_ids, my_p, effect_size) {
  genotype_mat <- matrix(NA, length(snp_ids), num_pts)
  for (i in 1:length(snp_ids)) {
    genotype_mat[i,] <- sample(c(0,1,2), size = num_pts, prob = genotype_prob(my_p[i]), replace = T)
  }
  genotype_mat <- data.frame(genotype_mat, maf=my_p, row.names = snp_ids, OR = effect_size)
  genotype_mat
}

#plot(density(rgamma(1000,shape = .3,rate = 1)))
#plot(density(1+abs(rnorm(num_pts, 0, 0.2))))


# parameter block
##################
num_pts <- 100000
num_snps <- 25
maf_sd <- 0.3 # standard dev for normal dist of minor allele freq
effect_size <- 1 + abs(rnorm(num_snps, 0, 0.2))
baserate <- 0.001
ages <- rnorm(num_pts, 65, 3)
#ages <- rep(65, num_pts)
hist(ages)
##################

pts <- pt_genotypes(num_pts, paste(rep("rs", num_snps), sample(5e8, size=num_snps, replace = F), sep=""), 
                    signif( abs(rnorm(num_snps, 0, maf_sd)), digits=2),
                    effect_size = effect_size)

mean(pts$OR)
mean(pts$maf)


ptodds <- exp(colSums(log(pts$OR ^ pts[,1:num_pts])))
ptprobs <- ptodds * baserate


which(ptprobs.age > 1)
mean(ptprobs)
ptprobs.age <- ptprobs/65 * ages
plot(density(ptprobs.age))

gc()

# generate cases 'n ctrls from polygenic risk data
#########################
sick <- numeric(num_pts)
set.seed(1234)
for (i in 1:num_pts) {
  sick[i] <- sample(c(0,1), size=1, prob = c(1-ptprobs.age[i], ptprobs.age[i]))
}
cases <- which(sick==TRUE)
length(cases)
ctrls <- sample(which(sick==FALSE), size = length(cases), replace = FALSE)
ptrisk <- melt(data.frame(ctrls=ptprobs[ctrls], cases=ptprobs[cases]))
names(ptrisk) <- c("trtgroup", "risk")
#########################

# violin plot compare distributin relative risk of cases controls
ggplot(data=ptrisk, aes(x=trtgroup, y=risk/baserate, colour=trtgroup, fill=trtgroup)) + geom_violin() 
#ggplot(data=data.frame(or=ptprobs[cases][order(ptprobs[cases])]/ptprobs[ctrls][order(ptprobs[ctrls])]), aes(x=or)) + geom_density(fill="grey", colour="grey")
pdf("mod_af_mod_eff_age.pdf")
ggplot(data=data.frame(cases=ptprobs[cases][order(ptprobs[cases])]/baserate, ctrls=ptprobs[ctrls][order(ptprobs[ctrls])]/baserate), aes(x=ctrls, y=cases)) + 
  geom_point() + #xlim(0, 100) + ylim(0, 100) + 
  geom_abline(intercept = 0, slope = 1, colour="red")
dev.off()

?sample

caseOR <- ptprobs[cases][order(ptprobs[cases])]/baserate
ctrlOR <- ptprobs[ctrls][order(ptprobs[ctrls])]/baserate
# % of cases OR > 5
length(which(caseOR>5))
length(which(caseOR>5))/length(caseOR)
# slope
mean(caseOR/ctrlOR)


################### OVARIAN CANCER
ovc_snps <- c("rs7651446", "rs3814113", "rs10088218", "rs9303542", "rs8170", "rs2072590", "rs34289250", 
              "rs11782652", "rs1243180", "rs11907546", "rs2046210", "rs10069690")
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
brca1.prob <- 0.025
brca2.prob <- 0.02
brca1.OR <- 20
brca2.OR <- 4
baserate <- 0.013

ovc_snps <- c("rs7651446", "rs3814113", "rs10088218", "rs9303542", "rs8170", "rs2072590", "rs34289250", 
              "rs11782652", "rs1243180", "rs11907546", "rs2046210", "rs10069690", "brca1", "brca2")
ovc_or   <- c(1.59, 1.2658228, 1.2987013, 1.14, 1.19, 1.14, 7.95, 1.24, 1.1, 1.1111112, 1.28, 1.14, brca1.OR, brca2.OR)
ovc_maf  <- c(0.05, 0.32, 0.13, 0.27, 0.19, 0.68, 0.0089, 0.07, 0.31, .35903, 0.08, 0.26, brca1.prob, brca2.prob)
prevalence <- 0.013

ovcpts <- pt_genotypes(num_pts, ovc_snps, ovc_maf, effect_size = ovc_or)
ovcodds <- exp(colSums(log(ovcpts$OR ^ ovcpts[,1:num_pts])))
gc()

ovcpts[,1:5]
#brca1.carrier <- sample(c(0,1), size=num_pts, prob=c(1-brca1.prob, brca1.prob), replace=TRUE)
#brca2.carrier <- sample(c(0,1), size=num_pts, prob=c(1-brca2.prob, brca2.prob), replace=TRUE)
brca1.carrier <- rep(FALSE, num_pts)
brca1.carrier[which(ovcpts[13,1:1e5]>0)] <- TRUE
length(brca1.carrier)
sum(brca1.carrier)
brca2.carrier <- rep(FALSE, num_pts)
brca2.carrier[which(ovcpts[14,1:1e5]>0)] <- TRUE
length(brca2.carrier)
sum(brca2.carrier)
brca.carrier <- rep("NON-BRCA", num_pts)
brca.carrier[which(ovcpts[13,1:1e5]>0)] <- "BRCA1"
brca.carrier[which(ovcpts[14,1:1e5]>0)] <- "BRCA2"

#set.seed(24)


ovcprobs.brca <- ovcodds 
popscale.factor <- prevalence / mean(ovcprobs.brca)
valid.patients <- which(ovcprobs.brca*popscale.factor < 1)
#f <- function (age) {1/age * age^((age-30)/80)}

nvalid <- length(valid.patients)

ages <- c(rep(20, nvalid), 
          rep(30, nvalid),
          rep(40, nvalid),
          rep(50, nvalid),
          rep(60, nvalid),
          rep(70, nvalid),
          rep(80, nvalid))




agerisk <- data.frame(age = as.factor(ages), 
                      risk = c(ovcprobs.brca.age.20, 
                               ovcprobs.brca.age.30,
                               ovcprobs.brca.age.40, 
                               ovcprobs.brca.age.50, 
                               ovcprobs.brca.age.60, 
                               ovcprobs.brca.age.70, 
                               ovcprobs.brca.age.80),
                      brca = as.factor(rep(brca.carrier[valid.patients], 7)))

gc()


dim(agerisk)

ggplot(agerisk, aes(x = age, y = risk)) + geom_boxplot(aes(fill = brca), notch = TRUE) + labs(y="log10 relative risk (RR)") #+ scale_y_log10()



# generate cases 'n ctrls from polygenic risk data
#########################
probs <- ovcprobs.brca.age

ovcsick.brca <- matrix(nrow = nvalid, ncol = 100)
dim(ovcsick.brca)

for (i in 1:nvalid) {
  ovcsick.brca[i,] <- sample(c(0,1), size=10, prob = c(1-popscale.factor*probs[i], popscale.factor*probs[i]), replace = TRUE)
}

pick.one <- sample(10, 1)

cases <- which(ovcsick.brca[,pick.one]==TRUE)
length(cases)
length(cases)/nvalid

ctrls <- sample(which(ovcsick.brca[,pick.one]==FALSE), size = length(cases), replace = FALSE)
ptrisk <- melt(data.frame(ctrls=probs[ctrls], cases=probs[cases]))
names(ptrisk) <- c("trtgroup", "risk")
#########################

probs.scale <- probs * popscale.factor
length(probs.scale)
mean(probs.scale[brca1.carrier[valid.patients]==T])
mean(probs.scale[brca2.carrier[valid.patients]==T])
mean(probs.scale[brca1.carrier[valid.patients]==F])
sum(cases %in% which(brca1.carrier[valid.patients]==TRUE))

# violin plot compare distributin relative risk of cases controls
ggplot(data=ptrisk, aes(x=trtgroup, y=risk/baserate, colour=trtgroup, fill=trtgroup)) + geom_violin() + labs(x = "patient population", y = "Relative Risk (RR)")
#ggplot(data=data.frame(or=ovcprobs.brca.age[cases][order(ovcprobs.brca.age[cases])]/ovcprobs.brca.age[ctrls][order(ovcprobs.brca.age[ctrls])]), aes(x=or)) + geom_density(fill="grey", colour="grey")
ggplot(data=data.frame(cases=probs[cases][order(probs[cases])]/baserate, ctrls=probs[ctrls][order(probs[ctrls])]/baserate), aes(x=ctrls, y=cases)) + 
  geom_point() + #xlim(0, 20) + ylim(0, 200) + 
  geom_abline(intercept = 0, slope = 1, colour="red") +
  geom_hline(yintercept=quantile(probs[cases]/baserate, probs=c(.5,.95,.99)), colour="grey", lty=2) +
  geom_vline(xintercept=quantile(probs[ctrls]/baserate, probs=c(.5,.95,.99)), colour="grey", lty=2)


caseOR <- probs[cases][order(probs[cases])]/baserate
ctrlOR <- probs[ctrls][order(probs[ctrls])]/baserate
# % of cases OR > 5
length(which(caseOR>5))
length(which(caseOR>5))/length(caseOR)
# slope
mean(caseOR/ctrlOR)

# prevalence
length(cases)/nvalid
# percent of brca1+ patients that are cases
sum(cases %in% which(brca1.carrier[valid.patients]==TRUE))/sum(brca1.carrier[valid.patients])
# percent of cases that are brca1+
sum(cases %in% which(brca1.carrier[valid.patients]==TRUE))/length(cases)
# percent of non-BRCA1 patients that are cases
sum(cases %in% which(brca1.carrier[valid.patients]==FALSE))/sum(brca1.carrier[valid.patients]==FALSE)
# percent of BRCA2 patients that are cases
sum(cases %in% which(brca2.carrier[valid.patients]==TRUE))/sum(brca2.carrier[valid.patients])


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
plot(x=1:100, y=inv.logit(1:100))

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


