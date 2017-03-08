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
maf_sd <- 0.1 # standard dev for normal dist of minor allele freq
effect_size <- 1 + abs(rnorm(num_snps, 0, 0.4))
baserate <- 0.01
yrs <- rnorm(num_pts, 65, 5)
##################

pts <- pt_genotypes(num_pts, paste(rep("rs", num_snps), sample(5e8, size=num_snps, replace = F), sep=""), 
                    signif( abs(rnorm(num_snps, 0, maf_sd)), digits=2),
                    effect_size = effect_size)

mean(pts$OR)
mean(pts$maf)

ptprobs <- apply(pts[,1:num_pts], MARGIN = 2, FUN=function(x) {baserate*prod(pts$OR^x)})

gc()

# generate cases 'n ctrls from polygenic risk data
#########################
sick <- numeric(num_pts)
for (i in 1:num_pts) {
  sick[i] <- sample(c(0,1), size=1, prob = c(1-ptprobs[i], ptprobs[i]))
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
ggplot(data=data.frame(cases=ptprobs[cases][order(ptprobs[cases])]/baserate, ctrls=ptprobs[ctrls][order(ptprobs[ctrls])]/baserate), aes(x=ctrls, y=cases)) + 
  geom_point() + #xlim(0, 100) + ylim(0, 100) + 
  geom_abline(intercept = 0, slope = 1, colour="red")



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

pts <- pt_genotypes(num_pts, ovc_snps, ovc_maf, effect_size = ovc_or)
ptprobs <- apply(pts[,1:num_pts], MARGIN = 2, FUN=function(x) {baserate*prod(pts$OR^x)})
gc()


#FUN=function(x) {baserate*prod(pts$OR^x)}
# generate cases 'n ctrls from polygenic risk data
#########################
sick <- numeric(num_pts)
for (i in 1:num_pts) {
  sick[i] <- sample(c(0,1), size=1, prob = c(1-ptprobs[i], ptprobs[i]))
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
png("ovarian_12_snps.png", w=4,h=4,u="in")
ggplot(data=data.frame(cases=ptprobs[cases][order(ptprobs[cases])]/baserate, ctrls=ptprobs[ctrls][order(ptprobs[ctrls])]/baserate), aes(x=ctrls, y=cases)) + 
  geom_point() + #xlim(0, 100) + ylim(0, 100) + 
  geom_abline(intercept = 0, slope = 1, colour="red")
dev.off()

(yrs-30) * pts$OR/30

caseOR <- ptprobs[cases][order(ptprobs[cases])]/baserate
ctrlOR <- ptprobs[ctrls][order(ptprobs[ctrls])]/baserate
# % of cases OR > 5
length(which(caseOR>5))
length(which(caseOR>5))/length(caseOR)
# slope
mean(caseOR/ctrlOR)
save(pts, yrs, file = "pts.RDa")
