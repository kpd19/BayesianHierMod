# original:
# gamma distribution
# flat priorts for alpha and beta, three MCMC chains, an adaptation period of 10,000 iterations
# 50,000 sampling, thinning by 50, for a total of 1000 samples from the posterior

# from text
# assumed vague priorts for v and C
# an experimentally rerived prior for the ratio parameter, rho
# explicitely allowed for error in the cadaver densities P0 
# used WAIC to compare the ability of the different models to explain the dat

library('rjags')
library('coda')
library('tidyverse')
library('ggmcmc')
library('polspline')


Data <- read.csv("DFTM_VData_2018.csv")
# Data <- Data[Data$strain!="LOV",]
MNPV <- c("TMB", "KLP", "BAN", "NMX")
# SNPV <- c("LOV", "COL", "CUB", "LST")
SNPV <- c("COL", "CUB", "LST")
strains <- c(MNPV, SNPV)
# strains = unique(Data$strain)
trees <- unique(Data$tree.sp)
Data$P0 <- Data$exact_density/0.15 # converts it to cadavers per m^2 based on Joe's average

Data <- Data[Data$strain !="LOV",]

Data$isolate <- NA
for (i in 1:length(MNPV)){
  Data[Data$strain==MNPV[i],]$isolate <- i
}

for (i in 1:length(SNPV)){
  Data[Data$strain==SNPV[i],]$isolate <- i
}

n.recov <- Data$alive # number of cats that survived the infection
n.inf <- Data$virus # number infected 
Virus <- Data$strain # the virus isolate
Morphotype <- Data$capsid # morphotype of virus isolate
Tree <- Data$tree.sp # tree species
Isolate <- Data$isolate # numeric value of the isolate for each morphotype, very hacky way to code

Data$strain <- factor(Data$strain)
Strain <- Data$strain

P <- Data$exact_density # pathogen density
N.obs <- dim(Data)[1] # number of observations
morph.index <- c(1,2,2,1,1,2,2) # index for how the isolates are ordered by morphotype, hacky


# my attempt to use initial conditions to help with the likelihood issue

# Iisolate_eff <- array(0, dim =c(7,2))
# M_eff_mean <- array(0.001,dim=c(2,2))
# n.sigma <- array(1,dim=c(2,2))
# for(i in 1:length(unique(Strain))){
#   for(j in 1:2){
#     Iisolate_eff[i,j] <- 0.001
#   }
# }
# initsList <- list(M_eff_mean = M_eff_mean,
#                n.sigma = n.sigma,
#                Iisolate_eff = Iisolate_eff)
# 

jags <- jags.model('newMod_01_Tree.bug',
                   data = list('n.recov' = n.recov,
                               'n.inf' = n.inf,
                               'Isolate' = Strain,
                               'Tree' = Tree, 
                               'P' = P,
                               'N.obs' = N.obs,
                               'len_iso' = length(unique(Strain)),
                               'morph.index' = morph.index),
                   # inits = initsList,
                   n.chains = 3,
                   n.adapt = 1000)

# I'm pretty sure this isn't working because it's going negative 
# and it can't go negative in dbinom

# need coda to do convergence diagnostics
samples <- coda.samples(jags,
                             c('M_eff_mean',
                               'n.sigma',
                               'isolate_est',
                               'isolate_eff',
                               'prob',
                               'loglik'),
                             100000)


# post_samples(filter(samples, Parameter =  ))

summary(samples)
coda_stats <- summary(samples)$statistics

stats_means <- as.numeric(coda_stats[,1])
stats_SD <- as.numeric(coda_stats[,2])

morpho_estimates <- stats_means[1:4]
isolate_estimates <- stats_means[(5+14):(5+14+13)]

virus_strains <- sort(unique(Strain))
v <- c("MNPV", "SNPV", "MNPV", "SNPV", as.character(virus_strains), as.character(virus_strains))
t <- c("DO", "DO", "GR", "GR", rep(c("DO", "GR"),each = 7))

model_output <- data.frame(v,t, c(morpho_estimates,isolate_estimates)) #, c(M_stats_SD,S_stats_SD))
colnames(model_output) <- c("strain", "tree_sp", "nu_bar")
model_output_Heir <- model_output

