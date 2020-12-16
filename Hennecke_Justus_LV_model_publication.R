#---------------------------------------------------------------------------------#
#                                                                                 #
# Lokta Volterra model to explore the effects of beta diversity of BEF            #
# relationships. Original code written by Patrick L. Thompson                     #
# and modified by F. van der Plas and Justus Hennecke                             #
# (e-mail: justus.hennecke@idiv.de)                                               #
#                                                                                 #
#                                                                                 #
#---------------------------------------------------------------------------------#


### libraries required
library(vegan)
library(Rmisc)
library(labdsv)
library(viridis)
library(reshape2)
library(ggplot2)
library(lm.beta)
library(patchwork)
library(svglite)

### part 1: the Lotka Volterra model ----------------------------------------------

BEF_beta <- function(alpha = 5, beta = 0.5, dispersal = 0.01, patches = 10,
                     env_niche = 0.2, heterogeneity = 0.5, show_plot = F,
                     max_alpha = 0.125, min_alpha=0.05, growth.rate = 2, dispersal.amount = 0.01,
                     tmax = 1000, stochasticity = 0.01, threshold = 0.05){
  
  gamma <- round(alpha/(1-beta))    # gamma diversity: total number of species
  if(is.infinite(gamma)){
    gamma <- alpha * patches
  }
  
  r <- rep(growth.rate, patches*gamma)
  a <- matrix(runif(n=(gamma*gamma), min=min_alpha, max=max_alpha), nrow=gamma, ncol=gamma) 
  diag(a) <- rep(0.1, gamma)

  env <- seq(0.5 - 0.5 * heterogeneity, 0.5 + 0.5 * heterogeneity,
                  length.out = patches) 

  env_opt <- runif(gamma, 0, 1) 
  
  if(gamma == (alpha * patches)) {
    n0 <- matrix(0, patches,gamma)
    for(s in 1:gamma){
      n0[rep(1:patches,each = alpha)[s],s] <- 1
    }
  }else{
    repeat{
      n0 <- matrix(0, patches, gamma)
      for(p in 1:patches){
        n0[p, sample(x=1:gamma, size=alpha, replace=F)] <- 1
      }
      not_pres <- which(colSums(n0)==0)
      for(np in not_pres){
        duplicates <- which(colSums(n0)>1)
        if(length(duplicates)>1){
          change_species <- sample(duplicates,1)
        }else{
          change_species <- duplicates
        }
        select_patch<-sample(which(n0[,change_species]==1),1)
        n0[select_patch,change_species] <- 0
        n0[select_patch,np] <- 1
      }
      if(sum(colSums(n0)>0)==gamma & sum(rowSums(n0)==alpha)==patches) break
    }
  }
  
  # function to count number of elements with positive value
  positive <- function(x){length(which(x > 0))}             
  
  alpha.realized.0 <- mean(apply(n0, 1, positive))          # initial alpha diversity
  beta.J.realized.0 <- mean(vegdist(n0, method="jaccard"))  # initial beta diversity (jaccard)
  beta.BC.realized.0 <- mean(vegdist(n0, method="bray"))    # initial beta diversity (Bray Curtis)
  
  n <- n0
  nsave <- array(NA, dim = c(patches, gamma, tmax))
  for(t in 1:tmax){
    nsave[,,t] <- n
    
    # Lotka Volterra equation
    nt <- r * exp(-((env - rep(env_opt, each = patches)) / 
                     env_niche)^2) * n - n * (matrix(n, ncol=nrow(a)) %*% a)

    # replace negative abundances by zero
    nt[which(nt < 0)] <- 0 
    
    # Stochasticity: actual biomass slightly lower/higher
    nt <- nt + rnorm(length(nt), 0, stochasticity*nt)
    
    # species with biomass below threshold get extinct
    nt[which(nt<threshold)] <- 0
    
    # randomly decide (exponential) which species emigrate from which patch
    emi.vector <- rexp(prod(dim(nt)),
                       rate = 1/as.vector(dispersal.amount * nt))
    
    # ensure that no more biomass can emigrate than there is in plot
    emi.vector[which(emi.vector > nt & nt > 0)] <-
      nt[which(emi.vector > nt & nt > 0)]
    
    # put in matrix
    emi.matrix <- matrix(emi.vector, nrow=dim(nt)[1])
    
    # immigration matrix: species migrate to random other plots
    immi.matrix <- matrix(rep(NA, prod(dim(nt))), nrow=nrow(nt))
    for(i in 1:ncol(emi.matrix)){
      immi.matrix[,i] <- emi.matrix[sample(c(1:nrow(emi.matrix)),
                                           nrow(emi.matrix), replace=F),i]
    }
    
    nt[which(nt<0)] <- 0
    emi.matrix[is.na(emi.matrix)] <- 0
    immi.matrix[is.na(immi.matrix)] <- 0
    
    # recalculate nt, corrected for emi- and immigration
    nt <- nt - emi.matrix + immi.matrix
    nt[which(nt<0)] <- 0
    n <- nt
  }
  
  if(show_plot==T & gamma>1){
    matplot(1:tmax, y=t(nsave[3,,]), type = 'l', col=viridis(100)[1+env_opt*99],
            lty=1, cex=2)
  }
  total.biomass <- sum(n)                                        #metacommunity biomass
  lv.df <- data.frame(alpha, beta, heterogeneity, total.biomass)
  lv.df$disp <- dispersal.amount                                 # dispersal level
  lv.df$alpha0 <- alpha                                          # initial alpha diversity                          
  lv.df$beta0 <- beta.J.realized.0                               # initial beta diversity
  lv.df$growth <- growth.rate                                    # intrinsic growth rate of all species
  lv.df$stochasticity <- stochasticity                           # level of stochasticity
  lv.df$min_alpha <- min_alpha                                   # minimum competition coefficient
  lv.df$max_alpha <- max_alpha                                   # maximum competition coefficient
  lv.df$niche <- env_niche                                       # niche width
  lv.df$patches <- patches                                       # number of patches in a metacommunity
  lv.df$tmax <- tmax                                             # number of timesteps
  lv.df$threshold <- threshold                                   # extinction threshold
  lv.df$alpha.realized.0 <- alpha.realized.0                     # initial alpha diversity
  lv.df$beta.J.realized.0 <- beta.J.realized.0                   # initial beta diversity (jaccard)
  lv.df$beta.BC.realized.0 <- beta.BC.realized.0                 # initial beta diversity (bray curtis)
  lv.df$gamma.0 <- gamma                                         # initial gamma diversity
  lv.df$alpha.realized.1 <- mean(rowSums(n>0))                   # final realized alpha diversity
  lv.df$gamma.1 <- sum(colSums(n)>0)                             # final realized gamma diversity
  distances.BC <- vegdist(n, method="bray")
  distances.BC[which(is.na(distances.BC))] <- 0
  lv.df$beta.BC.realized.1 <- mean(distances.BC)                 # final realized beta diversity (bray curtis)
  distances.J <- vegdist(n, method="jaccard")
  distances.J[which(is.na(distances.J))] <- 0
  lv.df$beta.J.realized.1 <- mean(distances.J)                   # final realized beta diversity (jaccard) 
  
  #Calculation of ommunity weighted mean of environmental optimum
  rel.abun <- n / rowSums(n) 
  env.opt.matrix <- matrix(rep(env_opt, each=patches), nrow=patches)
  env.opt.CWM <- rowSums(rel.abun * env.opt.matrix, na.rm=T)
  env.opt.CWM[which(env.opt.CWM==0)] <- NA
  
  #Calculation of environmental matching index
  env.fit <- 1-(mean(abs(env.opt.CWM-env),na.rm=T))
  lv.df$env.fit <- env.fit                                        # 'environmental matching'
  
  #Calculation of species sorting indices
  if (length(env.opt.CWM[!is.na(env.opt.CWM)])>2){
    sorting.correlation <- cor.test(env.opt.CWM, env)
    lv.df$sorting.correlation.p.new <- sorting.correlation$p.value
    lv.df$sorting.correlation.value.new <- sorting.correlation$estimate
  }else{
    sorting.correlation <- NA
    lv.df$sorting.correlation.value.new <- NA
    lv.df$sorting.correlation.p.new <- NA
  }
  
  n.complete <- n[,which(colSums(n)>0)]
  clust <- cut(env, 2, labels=FALSE)
  clust.complete <- clust[which(rowSums(n)>0)]
  if(is.matrix(n.complete)){
    n.complete.new <- n.complete[which(rowSums(n)>0),]
    if(length(unique(clust))==1){
      lv.df$species.sorting1 <- NA
      lv.df$species.sorting2 <- NA
      lv.df$species.sorting3 <- NA
      lv.df$species.sorting4 <- NA
      lv.df$species.sorting5 <- NA
    }else{
      lv.df$species.sorting1 <- mean(indval(n.complete.new,
                                            clustering=clust.complete)$indcls, na.rm=T)
      lv.df$species.sorting2 <- mean(apply(indval(n.complete.new,
                                                  clustering=clust.complete)$relabu, 1, max), na.rm=T)
      lv.df$species.sorting3 <- 1-mean(indval(n.complete.new,
                                              clustering=clust.complete)$pval, na.rm=T)
      lv.df$species.sorting4 <- 1-min(indval(n.complete.new,
                                             clustering=clust.complete)$pval, na.rm=T)
      lv.df$species.sorting5 <- length(which(indval(n.complete.new,
                                                    clustering=clust.complete)$pval<0.05))
    }
  }else{
    lv.df$species.sorting1 <- NA
    lv.df$species.sorting2 <- NA
    lv.df$species.sorting3 <- NA
    lv.df$species.sorting4 <- NA
    lv.df$species.sorting5 <- NA
  }
  return(lv.df)
}



### part 2: one single model run----------------------------------------------------------
lv.df <- BEF_beta(alpha = 5, beta = 0.5,  patches = 10,
                  env_niche = 0.2, heterogeneity = 0.5, show_plot = T, dispersal.amount = 0.01,
                  max_alpha = 0.125, threshold=0.05, stochasticity = 0.01)
lv.df


### part 3: Simulation to explore interaction effect of B diversity and heterogeneity --
heterogeneity_V <-seq(0, 1, 0.1) 
repetitions <- 1000
show_plot <- F
dispersal.amount_V <- c(0, 0.001, 0.01, 0.1)

set.seed(5)
results.df<-data.frame()
pb <- txtProgressBar(min = 0, max = length(dispersal.amount_V), style = 3)
for(dispersal.amount in dispersal.amount_V){
  setTxtProgressBar(pb, which(dispersal.amount_V == dispersal.amount))
    for(rep in 1:repetitions){
      for(hetero in heterogeneity_V){
          lv.df <- BEF_beta(heterogeneity = hetero, dispersal.amount = dispersal.amount,
          show_plot = show_plot)
          lv.df$rep <- rep
          results.df <- rbind(results.df,lv.df)
      }
    }
}  


### part 4: Data cleansing -----------------------------------------------------------------

#write.csv (results.df, file="Beta_div_data_uncleaned.csv")
#results.df <- read.csv (file="Beta_div_data_uncleaned.csv")
# replace NA with 0
results.df$beta.BC.realized.1[which(is.na(results.df$beta.BC.realized.1))] <- 0
results.df$beta.J.realized.1[which(is.na(results.df$beta.J.realized.1))] <- 0
results.df$alpha.realized.1[which(is.na(results.df$alpha.realized.1))] <- 0

#remove biomass of 0
results.df <- results.df[-which(results.df$total.biomass==0),]

#create columns for factorial heterogeneity and dispersal
results.df$heterogeneity.factor <- as.factor(results.df$heterogeneity)
results.df$disp.factor <- as.factor(results.df$disp)
#write.csv (results.df, file="Beta_div_data_cleaned.csv")

### part 5: Plots -------------------------------------------------------------------------

# FIGURE 2
plot1summary.het <- summarySE(results.df, measurevar = "alpha.realized.1", groupvars = c("heterogeneity","disp.factor"), na.rm=T)
plot1error.het <- ggplot(plot1summary.het, aes(x=heterogeneity, y=alpha.realized.1, color=disp.factor)) +
  geom_point(aes(fill=disp.factor), size=5, shape=21,color="black", position=position_dodge(width=0.07)) +
  geom_errorbar(aes(ymin=alpha.realized.1-sd, ymax=alpha.realized.1+sd), width=0.04, size=1, position=position_dodge(width=0.07))+
  scale_color_viridis(name="Dispersal", option="B", end=0.9, discrete=T)+
  scale_fill_viridis(name="Dispersal", option="B",end=0.9, discrete=T)+
  labs(x="Heterogeneity", y=expression("Realized "*alpha*"-diversity"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1),
        strip.text.x = element_text(size = 32),legend.text=element_blank(),legend.position="none", 
        legend.title=element_blank(), axis.title.y = element_text(size=32), axis.title.x = element_text(size=32),
        axis.text.y=element_text(color="black", size=32, angle=90, hjust=0.6), axis.text.x=element_text(size=32, color="black"), axis.ticks.length = unit(0.25, "cm"))
plot1error.het
#ggsave("global_alpha~heterogeneity.pdf", width= 10, height=8.5,dpi=600)

plot2summary.het <- summarySE(results.df, measurevar = "beta.BC.realized.1", groupvars = c("heterogeneity","disp.factor"), na.rm=T)
plot2error.het <- ggplot(plot2summary.het, aes(x=heterogeneity, y=beta.BC.realized.1, color=disp.factor)) +
  geom_point(aes(fill=disp.factor), size=5, shape=21,color="black", position=position_dodge(width=0.07)) +
  geom_errorbar(aes(ymin=beta.BC.realized.1-sd, ymax=beta.BC.realized.1+sd), width=0.04, size=1, position=position_dodge(width=0.07))+
  scale_color_viridis(name="Dispersal", option="B", end=0.9, discrete=T)+
  scale_fill_viridis(name="Dispersal", option="B",end=0.9, discrete=T)+
  labs(x="Heterogeneity", y=expression("Realized "*beta*"-diversity (Bray-Curtis)"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1),
        strip.text.x = element_text(size = 32), legend.text =element_text(size=32), legend.position = "top",
        legend.title=element_text(size=32), axis.title.y = element_text(size=32), axis.title.x = element_text(size=32), legend.key = element_rect(fill = NA),
        axis.text.y=element_text(color="black", size=32, angle=90, hjust=0.5), axis.text.x=element_text(size=32, color="black"), axis.ticks.length = unit(0.25, "cm"))
plot2error.het
#ggsave("global_beta~heterogeneity.png", width= 10, height=8.5,dpi=600)

plot3summary.het <- summarySE(results.df, measurevar = "env.fit", groupvars = c("heterogeneity","disp.factor"), na.rm=T)
plot3error.het <- ggplot(plot3summary.het, aes(x=heterogeneity, y=env.fit, color=disp.factor)) +
  geom_point(aes(fill=disp.factor), size=5, shape=21,color="black", position=position_dodge(width=0.07)) +
  geom_errorbar(aes(ymin=env.fit-sd, ymax=env.fit+sd), width=0.04, size=1, position=position_dodge(width=0.07))+
  scale_color_viridis(name="Dispersal", option="B", end=0.9, discrete=T)+
  scale_fill_viridis(name="Dispersal", option="B",end=0.9, discrete=T)+
  labs(x="Heterogeneity", y="Environmental matching")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1),
        strip.text.x = element_text(size = 32),legend.position="none", 
        legend.title=element_blank(), axis.title.y = element_text(size=32), axis.title.x = element_text(size=32),
        axis.text.y=element_text(color="black", size=32, angle=90, hjust=0.6), axis.text.x=element_text(size=32, color="black"), axis.ticks.length = unit(0.25, "cm"))
plot3error.het
#ggsave("global_env_fit~heterogeneity.png", width= 10, height=8.5,dpi=600)


plot1error.het+ plot2error.het+  plot3error.het+ plot_layout(ncol=3)
#ggsave("figure2_12012020.png", width= 26, height=9, dpi=300)
#ggsave("figure2_12012020.svg", width= 26, height=9, dpi=600)


# FIGURE 3
#3A: alpha diversity
results.df$alpha.cut <- cut(results.df$alpha.realized.1, breaks=c(seq(0,7,0.5)), 
                            labels = c("0.25","0.75","1.25","1.75","2.25","2.75","3.25","3.75","4.25","4.75","5.25","5.75","6.25","6.75"))
mean.biomass <- c()
biomass.sd <- c()
lower.5 <- c()
higher.5 <- c()
alpha.cat <- c()
for(i in sort(unique(results.df$alpha.cut))){
  subset.data <- results.df$total.biomass[which(results.df$alpha.cut==i)]
  mean.biomass <- append(mean.biomass, mean(subset.data, na.rm=T))
  alpha.cat <- append(alpha.cat, i)
  lower.5 <- append(lower.5, quantile(subset.data, 0.05))
  higher.5 <- append(higher.5, quantile(subset.data, 0.95))
  biomass.sd <- append(biomass.sd, sd(subset.data))
  output.data.alpha <- data.frame(alpha.cat,mean.biomass, biomass.sd, lower.5, higher.5)
}

output.data.alpha$alpha.cat <- as.numeric(output.data.alpha$alpha.cat)

plot4<- ggplot(output.data.alpha, aes(x=alpha.cat, y=mean.biomass)) +
  geom_point(size=5, shape=21,color="black", fill="#39568CFF") +
  geom_errorbar(aes(ymin=mean.biomass-biomass.sd, ymax=mean.biomass+biomass.sd), color="#39568CFF", width=0.2, size=1,)+
  ylim(-0.8, 125)+
  scale_x_continuous(breaks=c(0,2,4,6))+
  labs(x=expression("Realized "*alpha*"-diversity"), y="Metacommunity biomass")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1),
        strip.text.x = element_text(size = 32),legend.position="none", 
        legend.title=element_blank(), axis.title.y = element_text(size=32), axis.title.x = element_text(size=32),
        axis.text.y=element_text(color="black", size=32), axis.text.x=element_text(size=32, color="black"), 
        axis.ticks.length = unit(0.25, "cm"))
plot4


#3B: beta diversity
results.df$beta.cut <- cut(results.df$beta.BC.realized.1, breaks=seq(0,1,0.1),
                           labels = c("0.05","0.15","0.25","0.35","0.45","0.55","0.65","0.75","0.85","0.95"))

mean.biomass <- c()
biomass.sd <- c()
lower.5 <- c()
higher.5 <- c()
beta.cat <- c()
for(i in sort(unique(results.df$beta.cut))){
  subset.data <- results.df$total.biomass[which(results.df$beta.cut==i)]
  mean.biomass <- append(mean.biomass, mean(subset.data, na.rm=T))
  beta.cat <- append(beta.cat, i)
  lower.5 <- append(lower.5, quantile(subset.data, 0.05))
  higher.5 <- append(higher.5, quantile(subset.data, 0.95))
  biomass.sd <- append(biomass.sd, sd(subset.data))
  output.data.beta <- data.frame(beta.cat,mean.biomass, biomass.sd, lower.5, higher.5)
}

output.data.beta$beta.cat <- as.numeric (output.data.beta$beta.cat)

plot5<- ggplot(output.data.beta, aes(x=beta.cat, y=mean.biomass)) +
  geom_point(size=5, shape=21,color="black", fill="#39568CFF") +
  geom_errorbar(aes(ymin=mean.biomass-biomass.sd, ymax=mean.biomass+biomass.sd), color="#39568CFF", width=0.03, size=1,)+
  xlim(0.0,1)+
  ylim(-0.8, 125)+
  labs(x=expression("Realized "*beta*"-diversity (Bray-Curtis)"), y="Metacommunity biomass")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1),
        strip.text.x = element_text(size = 32),legend.position="none", 
        legend.title=element_blank(), axis.title.y = element_blank(), axis.title.x = element_text(size=32),
        axis.text.y=element_blank(), axis.text.x=element_text(size=32, color="black"), axis.ticks.length = unit(0.25, "cm"))
plot5

#env.matching
results.df$matching.cut <- cut(results.df$env.fit, breaks=seq(0.675,1,0.025), 
                               labels = c("0.6875","0.7125","0.7375","0.7875","0.8125","0.8375","0.8625","0.8875","0.9125","0.9375","0.9625","0.9875","1"))

mean.biomass <- c()
biomass.sd <- c()
lower.5 <- c()
higher.5 <- c()
matching.cat <- c()
for(i in sort(unique(results.df$matching.cut))){
  subset.data <- results.df$total.biomass[which(results.df$matching.cut==i)]
  mean.biomass <- append(mean.biomass, mean(subset.data, na.rm=T))
  matching.cat <- append(matching.cat, i)
  lower.5 <- append(lower.5, quantile(subset.data, 0.05))
  higher.5 <- append(higher.5, quantile(subset.data, 0.95))
  biomass.sd <- append(biomass.sd, sd(subset.data))
  output.data.env <- data.frame(matching.cat,mean.biomass, biomass.sd, lower.5, higher.5)
}
output.data.env$matching.cat <- as.numeric (output.data.env$matching.cat)

output.data.env <- output.data.env [-which(is.na(output.data.env$biomass.sd)),]

plot6<- ggplot(output.data.env, aes(x=matching.cat, y=mean.biomass)) +
  geom_point(size=5, shape=21,color="black", fill="#39568CFF") +
  geom_errorbar(aes(ymin=mean.biomass-biomass.sd, ymax=mean.biomass+biomass.sd), color="#39568CFF", width=0.01, size=1,)+
  ylim(-0.8, 125)+
  labs(x=expression("Environmental matching"), y="Metacommunity biomass")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1),
        strip.text.x = element_text(size = 32),legend.position="none", 
        legend.title=element_blank(), axis.title.y = element_blank(), axis.title.x = element_text(size=32),
        axis.text.y=element_blank(), axis.text.x=element_text(size=32, color="black"), axis.ticks.length = unit(0.25, "cm"))
plot6

plot4 + plot5 + plot6 + plot_layout(ncol=3)
#ggsave("figure3_12012020.png", width= 26, height=9, dpi=300)
#ggsave("figure3_12012020.svg", width= 26, height=9, dpi=600)

# FIGURE 4

# without env.matching
hetero.c <- c()
disp.c <- c()
beta.effect.c <- c()
alpha.effect.c <- c()
alpha.error.c <- c()
beta.error.c <- c()
interaction.r2 <- c()
additive.r2 <- c()
interaction.AIC <- c()
additive.AIC <- c()

for(i in 1:length(unique(results.df$disp))){
  for(j in 1:length(unique(results.df$heterogeneity))){
    data.subset <- results.df[which(results.df$disp==unique(results.df$disp)[i] &
                                      results.df$heterogeneity==unique(results.df$heterogeneity)[j]),]
    m1 <- lm(scale(data.subset$total.biomass) ~ scale(data.subset$alpha.realized.1) + scale(data.subset$beta.BC.realized.1))
    m1.interactions <- lm(scale(data.subset$total.biomass) ~ scale(data.subset$alpha.realized.1) * scale(data.subset$beta.BC.realized.1))
    beta.effect <- summary(m1)$coefficient[3,1]
    beta.error <- summary(m1)$coefficient[3,2]
    alpha.effect <- summary(m1)$coefficient[2,1]
    alpha.error <- summary(m1)$coefficient[2,2]
    hetero.c <- append(hetero.c, unique(results.df$heterogeneity)[j])
    disp.c <- append (disp.c, unique(results.df$disp)[i])
    beta.effect.c <- append(beta.effect.c, beta.effect)
    alpha.effect.c <- append(alpha.effect.c, alpha.effect)
    alpha.error.c <- append(alpha.error.c, alpha.error)
    beta.error.c <- append(beta.error.c, beta.error)
    interaction.r2 <- append(interaction.r2, summary(m1.interactions)$adj.r.squared)
    additive.r2 <- append(additive.r2, summary(m1)$adj.r.squared)
    interaction.AIC <- append (interaction.AIC, AIC(m1.interactions))
    additive.AIC <- append (additive.AIC, AIC(m1))
    df.lm <- data.frame(hetero.c, disp.c, alpha.effect.c, beta.effect.c, alpha.error.c, beta.error.c)
  }
}

interaction.r2 ; mean(interaction.r2)
additive.r2 ; mean(additive.r2)
interaction.AIC; mean(interaction.AIC)
additive.AIC; mean(additive.AIC)

df.lm.long  <- melt(df.lm, id.vars=c("hetero.c", "disp.c","beta.error.c","alpha.error.c"))
df.lm.long$disp.c <- factor(df.lm.long$disp.c, labels = c("Dispersal = 0", "Dispersal = 0.001" , "Dispersal = 0.01", "Dispersal = 0.1"))
plot7 <- ggplot(df.lm.long, aes(x=value, y=hetero.c, fill=variable)) +
  geom_point(shape=21, size=7)+
  scale_fill_manual (values=c("#433C84FF","#29AF7FFF"), labels=c(expression("Realized "*alpha*"-diversity"),expression("Realized "*beta*"-diversity (Bray-Curtis)")))+
  geom_vline(xintercept= 0, linetype="dotted")+
  xlim(-1,1)+
  facet_wrap(~disp.c,ncol=4)+
  labs(x="Effect on metacommunity biomass", y="Heterogeneity")+
  theme(panel.spacing = unit (1.5, "lines"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1), legend.spacing.x = unit(1, 'cm'),
        strip.text.x = element_text(size = 32),legend.text=element_text(size=32),legend.position="none", legend.key = element_rect(fill = NA),
        legend.title=element_blank(),axis.title.x = element_blank(), axis.title.y = element_text(size=36), axis.text.x=element_blank(),
        axis.text.y=element_text(color="black", size=28, angle=90, hjust=0.5), 
        axis.ticks.length = unit(0.25, "cm"))
plot7

# with env.matching
hetero.c <- c()
disp.c <- c()
beta.effect.c <- c()
alpha.effect.c <- c()
alpha.error.c <- c()
beta.error.c <- c()
envfit.effect.c <- c()
envfit.error.c <- c()
interaction.r2 <- c()
additive.r2 <- c()
interaction.AIC <- c()
additive.AIC <- c()

for(i in 1:length(unique(results.df$disp))){
  for(j in 1:length(unique(results.df$heterogeneity))){
    data.subset <- results.df[which(results.df$disp==unique(results.df$disp)[i] &
                                      results.df$heterogeneity==unique(results.df$heterogeneity)[j]),]
    m1 <- lm(scale(data.subset$total.biomass) ~ scale(data.subset$alpha.realized.1) + scale(data.subset$beta.BC.realized.1)+
               scale(data.subset$env.fit))
    m1.interactions <- lm(scale(data.subset$total.biomass) ~ scale(data.subset$alpha.realized.1) * scale(data.subset$beta.BC.realized.1) 
                          * scale(data.subset$env.fit))
    beta.effect <- summary(m1)$coefficient[3,1]
    beta.error <- summary(m1)$coefficient[3,2]
    alpha.effect <- summary(m1)$coefficient[2,1]
    alpha.error <- summary(m1)$coefficient[2,2]
    envfit.effect <- summary(m1)$coefficient[4,1]
    envfit.error <- summary(m1)$coefficient[4,2]
    hetero.c <- append(hetero.c, unique(results.df$heterogeneity)[j])
    disp.c <- append (disp.c, unique(results.df$disp)[i])
    beta.effect.c <- append(beta.effect.c, beta.effect)
    alpha.effect.c <- append(alpha.effect.c, alpha.effect)
    alpha.error.c <- append(alpha.error.c, alpha.error)
    beta.error.c <- append(beta.error.c, beta.error)
    envfit.effect.c <- append(envfit.effect.c, envfit.effect)
    envfit.error.c <- append(envfit.error.c, envfit.error)
    interaction.r2 <- append(interaction.r2, summary(m1.interactions)$adj.r.squared)
    additive.r2 <- append(additive.r2, summary(m1)$adj.r.squared)
    interaction.AIC <- append (interaction.AIC, AIC(m1.interactions))
    additive.AIC <- append (additive.AIC, AIC(m1))
    df.lm <- data.frame(hetero.c, disp.c, alpha.effect.c, beta.effect.c, alpha.error.c, beta.error.c, envfit.effect.c, envfit.error.c)
  }
}

interaction.r2 ; mean(interaction.r2)
additive.r2 ; mean(additive.r2)
interaction.AIC; mean(interaction.AIC)
additive.AIC; mean(additive.AIC)

df.lm.long  <- melt(df.lm, id.vars=c("hetero.c", "disp.c","beta.error.c","alpha.error.c", "envfit.error.c"))
df.lm.long$disp.c <- factor(df.lm.long$disp.c, labels = c("Dispersal = 0", "Dispersal = 0.001" , "Dispersal = 0.01", "Dispersal = 0.1"))

plot8 <- ggplot(df.lm.long, aes(x=value, y=hetero.c, fill=variable)) +
  geom_point(shape=21, size=7)+
  #scale_fill_viridis(discrete=T, labels=c(expression("Realized "*alpha*"-diversity"),expression("Realized "*beta*"-diversity (Bray-Curtis)", "Environmental matching")))+
  scale_fill_manual (values=c("#433C84FF","#29AF7FFF","#f6d645"),labels=c(expression("Realized "*alpha*"-diversity"),expression("Realized "*beta*"-diversity (Bray-Curtis)"),"Environmental matching"))+
  xlim(-1,1)+
  geom_vline(xintercept= 0, linetype="dotted")+
  facet_wrap(~disp.c,ncol=4)+
  labs(x="Effect on metacommunity biomass", y="Heterogeneity")+
  theme(panel.spacing = unit (1.5, "lines"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1), legend.spacing.x = unit(1, 'cm'),
        strip.text.x = element_blank(),legend.text=element_text(size=32), legend.position="top", legend.key = element_rect(fill = NA),
        legend.title=element_blank(), axis.title = element_text(size=36), axis.text.x=element_text(color="black", size=26),
        axis.text.y=element_text(color="black", size=28, angle=90, hjust=0.5), axis.ticks.length = unit(0.25, "cm"))
plot8

plot7 + plot8 + plot_layout(ncol=1)
#ggsave("figure4_12012020.png", width= 21, height=17,dpi=300)
#ggsave("figure4_12012020.svg", width= 21, height=17,dpi=600)

