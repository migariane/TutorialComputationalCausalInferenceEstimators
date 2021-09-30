###################################################################################################

# Tutorial: causal inference methods made easy for applied resarchers/epidemiologists/statisticians 

# ICON-LSHTM, LONDON, 30th June 2020

# Miguel Angel Luque Fernandez, PhD
# Assistant Professor of Epidemiology and Biostatistics
# Matthew Smith, PhD
# Research Fellow

# Inequalities in Cancer Outcomes Network, LSHTM, London, UK

# Copyright (c) 2020 Permission is hereby granted, free of charge, to any person obtaining a copy 
# of this software and associated documentation files (the "Software"), to deal in the Software 
# without restriction, including without limitation the rights to use, copy, modify, merge, 
# publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to 
# whom the Software is furnished to do so, subject to the following conditions: The above 
# copyright notice and this permission notice shall be included in all copies or substantial 
# portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING 
# BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NON 
# INFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES 
# OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR 
# IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# Bug reports: miguel-angel.luque at lshtm.ac.uk	

# The rhc dataset can be dowloaded at http://biostat.mc.vanderbilt.edu/wiki/Main/DataSets

###################################################################################################

# Preliminaries
  rm(list=ls())
  set.seed(1)  

  ### Box 1: Setting the data
  setwd("your path")
  #setwd("~/Dropbox/ESTIMATORSCIproject/R_Stata_master_files/Data")
  library(haven)
  #data <- read_dta("~/Dropbox/ESTIMATORSCIproject/R_Stata_master_files/Data/rhc.dta")
  data <- read_dta("rhc.dta")
  # Define the outcome (Y), exposure (A), confounder (C), and confounders (W)
  data$Y  <- data$death_d30; data$Y <- as.numeric(data$Y); Y <- data$Y 
  data$A  <- data$rhc; data$A <- as.numeric(data$A); A <- data$A
  data$C  <- data$sex; data$C <- as.numeric(data$C); C <- data$C
  data$w1 <- data$age; data$w1 <- as.numeric(data$w1); w1 <- data$w1
  data$w2 <- data$edu; data$w2 <- as.numeric(data$w2); w2 <- data$w2
  data$w3 <- data$race; data$w3 <- as.numeric(data$w3); w3 <- data$w3
  data$w4 <- data$carcinoma; data$w4 <- as.numeric(data$w4); w4 <- data$w4
  data2   <- as.data.frame(Y); data2$A <- A; data2$C <- C; data2$w1 <- w1; data2$w2 <- w2; data2$w3 <- w3; data2$w4 <- w4
    
  
  ### Box 2: Naive estimate of the ATE
  naive <- lm(Y ~ A + C, data=data); naive  # Naive estimate of the ATE is 0.07352

# 3. G-Formula
  
  ## 3.1 Non-parametric G-formula
  
    ### Box 3: Non-parametric G-formula for the ATE 
    mean(data$A[data$C==1], na.rm=TRUE)       # 
    mean(data$A[data$C==0], na.rm=TRUE)       # 
    mean(data$Y[data$A==1], na.rm=TRUE) - mean(data$Y[data$A==0],na.rm=TRUE)    # Unadjusted Estimate
    reg <- lm(Y ~ A, data=data); reg          # Unadjusted Estimate Regression 
    pr.l <- prop.table(table(data$C)); pr.l   # Marginal probability of C
    tab.out <- aggregate(Y ~ A + C, data, mean); tab.out    # Table of Means in 
    ATE <- ((mean(data$Y[data$A==1 & data$C==1]) - mean(data$Y[data$A==0 & data$C==1]))*pr.l[2]) + 
     (mean(data$Y[data$A==1 & data$C==0]) - mean(data$Y[data$A==0 & data$C==0]))*pr.l[1]   # G-formula Non-parametric ATE 
    ATE; rm(ATE)    # The ATE from the non-parametric estimator is 0.073692
    
    ### Box 4: Bootstrap the 95% confidence intervals (CI) for the ATE/ATT estimated using the non-parametric G-Formula
    # ATE
    library(boot)
    g.comp = function(data,indices)     # Define the function to estimate the ATE
    {
      dat=data[indices,]
      pr.l <- prop.table(table(dat$C))
      
      ATE = ((mean(dat$Y[dat$A==1 & dat$C==1]) - mean(dat$Y[dat$A==0 & dat$C==1]))*pr.l[2]) + 
       (mean(dat$Y[dat$A==1 & dat$C==0]) - mean(dat$Y[dat$A==0 & dat$C==0]))*pr.l[1] ; ATE
    }
    g.comp(data,indices=1:nrow(data))     # Can get original estimate, by plugging in indices 1:n
    boot.out=boot(data,g.comp,200)      # Draw 200 bootstrap sample estimates 
    boot.ci(boot.out,type="perc",conf=0.95)     # compute confidence intervals using percentile method
    boot.ci(boot.out,type="norm",conf=0.95)
    
    ### Box 5: Non-parametric G-Formula using a fully saturated regression model in Stata (A)
        # Method 1: conditional probabilities
    data$A1 <- ifelse(data$A == 1, 1, 0)
    data$A0 <- ifelse(data$A == 0, 1, 0)
    data$C1 <- ifelse(data$C == 1, 1, 0)
    data$C0 <- ifelse(data$C == 0, 1, 0)
    reg <- glm(Y ~ -1 + (A1 + A0) + A1:(C1) + A0:(C1), data=data); summary(reg)
    ATE <- mean((reg$coefficients[1] + reg$coefficients[3]*C) - (reg$coefficients[2] + reg$coefficients[4]*C)); ATE
    rm(ATE)
    
    ### Box 6: Non-parametric G-Formula using a fully saturated regression model in Stata (B)
        # Method 2: Marginal probabilities
    install.packages("margins")
    library(margins)
    reg <- glm(Y ~ -1 + (A1 + A0) + A1:(C1) + A0:(C1), data=data); summary(reg)
    Y1 <- margins(reg, variables="A1"); Y1
    Y0 <- margins(reg, variables="A0"); Y0
    ATE <- mean(Y1$fitted[A==1]) - mean(Y0$fitted[A==0]); ATE
    rm(ATE)
    
  ## 3.2 Parametric G-formula
    ### Box 7: Parametric G-formula by hand
    mod1  <- glm(Y ~  C, family="binomial", data=data[data$A==1,])    # Expected probability amongst those with RHC
    mod0  <- glm(Y ~  C, family="binomial", data=data[data$A==0,])    # Expected probability amongst those without RHC
    GcompRA  <- cbind(Y1 = predict(mod1, newdata=data.frame(A = 1, C), type="response"),
                   Y0 = predict(mod0, newdata=data.frame(A = 0, C), type="response"))
    GcompRA <- as.data.frame(GcompRA)
    Y.1 <- GcompRA$Y1
    Y.0 <- GcompRA$Y0
    ATE <- mean((Y.1) - (Y.0), na.rm=TRUE); ATE     # Difference between expected probabilities (ATE) 
    rm(ATE)
  
    
    ### Box 8: Parametric regression adjustment (one confounder) using stdReg R-package
    install.packages("stdReg")
    library(stdReg)
    reg <- glm(Y ~ A + C, data = data, family = poisson(link="log")); summary(reg)
    reg.std <- stdGlm(fit=reg, data = data, X = "A", x=seq(0,1))
    print(summary(reg.std, contrast = "difference", reference=0))
    plot(reg.std)

    ### Box 9: Bootstrap for the parametric regression adjustment one confounder
    library(boot)           # Install the Bootstrap package
    attach(data)          
    g.comp=function(data,indices)       # Define the function to estimate the ATE
    {
      dat=data[indices,]
      glm1  <- glm(Y ~ C, family="binomial", dat=dat[dat$A==1,])
      glm2  <- glm(Y ~ C, family="binomial", dat=dat[dat$A==0,])
      Y.1 = predict(glm1, newdata=data.frame(A = 1, C), type="response")
      Y.0 = predict(glm2, newdata=data.frame(A = 0, C), type="response")
      ATE <- mean((Y.1) - mean(Y.0)); ATE
    }
    g.comp(data,indices=1:nrow(data))     # Can get original estimate, by plugging in indices 1:n
    boot.out=boot(data,g.comp,200)        # Draw 1000 bootstrap sample estimates of RD
    boot.ci(boot.out,type="norm",conf=0.95)     # Bootstrapped 95% CI based on normal approximation
    boot.ci(boot.out,type="perc",conf=0.95)     # Bootstrapped 95% CI based on percentiles of the bootstrap replicates
    
    # Now with more than one confounder
    
    ### Box 10: Parametric multivariable regression adjustment implementation of the G-Formula
    mod1  <- glm(Y ~  C + w1 + w2 + as.factor(w3) + as.factor(w4), family="binomial", data=data[data$A==1,])    # Expected probability amongst those with RHC
    mod0  <- glm(Y ~  C + w1 + w2 + as.factor(w3) + as.factor(w4), family="binomial", data=data[data$A==0,])    # Expected probability amongst those without RHC
    GcompRA  <- cbind(Y1 = predict(mod1, newdata=data.frame(A = 1, C, w1, w2, as.factor(w3), as.factor(w4)), type="response"),
                   Y0 = predict(mod0, newdata=data.frame(A = 0, C, w1, w2, as.factor(w3), as.factor(w4)), type="response"))
    GcompRA <- as.data.frame(GcompRA)
    Y.1 <- GcompRA$Y1
    Y.0 <- GcompRA$Y0
    ATE <- mean((Y.1) - (Y.0), na.rm=TRUE); ATE     # ATE 
    rm(ATE)
    
    ### Box 11: Parametric multivariable regression adjustment using "stdReg" R-package
    install.packages("stdReg")
    library(stdReg)
    reg <- glm(Y ~ A + C + w1 + w2 + as.factor(w3) + as.factor(w4), data = data, family = poisson(link="log")); summary(reg)
    reg.std <- stdGlm(fit=reg, data=data, X="A", x=seq(0,1))
    print(summary(reg.std, contrast="difference", reference=0))
    plot(reg.std)
    
    ### Box 12: Parametric multivariable regression adjustment using "margins" R-package
    reg1 <- glm(Y ~ -1 + (A1 + A0) + A1:(C1 + w1 + w2 + as.factor(w3) + as.factor(w4)) + A0:(C0 + w1 + w2 + as.factor(w3) + as.factor(w4)) , data=data); summary(reg1)
    poY1m <- margins(reg1, variables="A1"); poY1m
    poY0m <- margins(reg1, variables="A0"); poY0m
    ATE2 <- poY1m$fitted[A==1] - poY0m$fitted[A==0]; mean(ATE2)
    rm(ATE2)    

    ### Box 13 Bootstrap for the multivariable parametric regression adjustment
    library(boot)           # Install the Bootstrap package
    attach(data)          
    g.comp=function(data,indices)       # Define the function to estimate the ATE
    {
      dat=data[indices,]
      glm1  <- glm(Y ~ C + w1 + w2 + as.factor(w3) + as.factor(w4), family="binomial", dat=dat[dat$A==1,])
      glm2  <- glm(Y ~ C + w1 + w2 + as.factor(w3) + as.factor(w4), family="binomial", dat=dat[dat$A==0,])
      Y.1 = predict(glm1, newdata=data.frame(A = 1, C, w1, w2, as.factor(w3), as.factor(w4)), type="response")
      Y.0 = predict(glm2, newdata=data.frame(A = 0, C, w1, w2, as.factor(w3), as.factor(w4)), type="response")
      mean((Y.1) - mean(Y.0))
    }
    g.comp(data,indices=1:nrow(data))     # Can get original estimate, by plugging in indices 1:n
    boot.out=boot(data,g.comp,200)        # Draw 1000 bootstrap sample estimates of RD
    boot.ci(boot.out,type="norm",conf=0.95)     # Bootstrapped 95% CI based on normal approximation
    boot.ci(boot.out,type="perc",conf=0.95)     # Bootstrapped 95% CI based on percentiles of the bootstrap replicates
  
# 4. Inverse Probability of Treatment Weighting
  
  ## 4.1 Inverse probability of treatment weighting based on the propensity score plus regression adjustment
    
    # Box 14 (IPTW by hand)
    p.s <- glm(A ~ as.factor(C) + w1 + w2 + as.factor(w3) + as.factor(w4), data=data, family=binomial)      # Propensity score mmodel for the exposure
    p.score <- ifelse(data$A == 0, 1 - predict(p.s, type = "response"), predict(p.s, type = "response"))  # Assign Propensity score weights
    #table(p.score)      # Table of Propensity Scores
    data$w <- 1/p.score       # Generate IP Weights 
    data2$w <- 1/p.score
    #table(data$w); summary(data$w); sd(data$w)
    
    ATE <- mean(data$w*as.numeric(data$A==1)*data$Y) - mean(data$w*as.numeric(data$A==0)*data$Y);ATE  # Estimate ATE
    rm(ATE)
  
    # Box 15 Bootstrap computation for the IPTW estimator
    library(boot)
    iptw.w = function(data,indices)     # Define the function to estimate the ATE
    {
      dat=data[indices,]
      mean(dat$w*as.numeric(dat$A==1)*dat$Y) - mean(dat$w*as.numeric(dat$A==0)*dat$Y)
    }
    iptw.w(data,indices=1:nrow(data))     # Can get original estimate, by plugging in indices 1:n
    boot.out=boot(data,iptw.w,100)      # Draw 200 bootstrap sample estimates 
    boot.ci(boot.out,type="perc",conf=0.95)     # compute confidence intervals using percentile method
    boot.ci(boot.out,type="norm",conf=0.95)

    ### Box 16: Computation of the IPTW estimator for the ATE using IPW R-package
    install.packages("ipw", "survey")
    library(ipw)
    library(survey)
    
    # Univariable
    ipw.ATE <- ipwpoint(exposure = A, family = "binomial", link = "logit", 
                        numerator = ~ 1, 
                        denominator = ~ C, 
                        data = data2)
    summary(ipw.ATE$ipw.weights)
    ipwplot(weights = ipw.ATE$ipw.weights, logscale = FALSE, main = "Unstabilized weights", xlim = c(0.5, 2))
    summary(ipw.ATE$num.mod)
    summary(ipw.ATE$den.mod)
    data2$usw <- ipw.ATE$ipw.weights
    msm <- (svyglm(Y ~ A, design = svydesign(~ 1, weights = ~ usw, data = data2)))
    coef(msm);  confint(msm)
    
    # Multivariable
    ipw.ATE <- ipwpoint(exposure = A, family = "binomial", link = "logit", 
                        numerator = ~ 1, 
                        denominator = ~ C + w1 + w2 + as.factor(w3) + as.factor(w4), 
                        data = data2)
    summary(ipw.ATE$ipw.weights)
    ipwplot(weights = ipw.ATE$ipw.weights, logscale = FALSE, main = "Unstabilized weights", xlim = c(0.5, 2))
    summary(ipw.ATE$num.mod)
    summary(ipw.ATE$den.mod)
    data2$usw <- ipw.ATE$ipw.weights
    msm <- (svyglm(Y ~ A, design = svydesign(~ 1, weights = ~ usw, data = data2)))
    coef(msm);  confint(msm)
  
    ### Box 17: Assessing IPTW balance
    install.packages("twang")
    library(twang)
    ps.balance <- ps(A ~ C + w1 + w2 + as.factor(w3) + as.factor(w4), data = data2, 
                     n.trees=1000, interaction.depth=2, shrinkage=0.01, perm.test.iters=0,
                     stop.method=c("es.mean","ks.max"), estimand = "ATE", verbose=FALSE)
    plot(ps.balance)
    summary(ps.balance$gbm.obj, n.trees=ps.balance$desc$ks.max.ATE$n.trees, plot=FALSE)
    data2.balance <- bal.table(ps.balance); data2.balance  
    
    ### Box 18: Assessing IPTW overlap by hand
    install.packages("xtable")
    library(xtable)
    pretty.tab <- data2.balance$ks.max.ATE[,c("tx.mn","ct.mn","ks")]
    pretty.tab <- cbind(pretty.tab, data2.balance$unw[,"ct.mn"])
    names(pretty.tab) <- c("E(Y1|t=1)","E(Y0|t=1)","KS","E(Y0|t=0)")
    xtable(pretty.tab, caption = "Balance of the treatment and comparison groups",
           label = "tab:balance", digits = c(0, 2, 2, 2, 2), align=c("l","r","r","r","r"))
    plot(ps.balance, plots = 6)
    
    ### Box 19: Assessing overlap using plots
    # Fit a propensity score model
    m_PS<-glm(A ~ C + w1 + w2 + as.factor(w3) + as.factor(w4), data = data2, family=binomial(link="logit"))
    summary(m_PS)
    
    # Estimate the propensity score
    data$PS<-fitted.values(m_PS)
    
    # Histogram of the PS
    hist(data$PS[data$rhc==0])
    hist(data$PS[data$rhc==1])
    plot(density(data$PS[data$rhc==0]),col="red",lwd=2, xlab="PS")
    lines(density(data$PS[data$rhc==1]),col="blue",lwd=2)
    legend("topright", legend=c("No RHC", "RHC"), pch="--", col=c("red","blue"), bty="n", lwd=2)
    
    # Look at minimum and maximum PS in each exposure group
    min(data$PS[data$rhc==0])
    min(data$PS[data$rhc==1])
    max(data$PS[data$rhc==0])
    max(data$PS[data$rhc==1])
    
    # Investigate overlap (i.e. positivity) 
    data$overlap <- ifelse(data$PS>=min(data$PS[data$rhc==1]) & data$PS<=max(data$PS[data$rhc==0]),1,0); table(data$overlap,data$rhc)
    
  ## 4.2 Marginal structural model with stabilised weights
    ### Box 20: Computation of the IPTW estimator for the ATE using a MSM
    
    # Unstabilized weights 
    msm <- lm(Y  ~ A + C +  w1 + w2 + as.factor(w3) + as.factor(w4), data = data, weights = data$w)     # MSM
    library(sandwich)
    SE <-sqrt(diag(vcovHC(msm, type="HC0")))       # robust standard errors
    beta <- coef(msm)
    lcl <- beta-1.96*SE 
    ucl <- beta+1.96*SE
    cbind(beta, lcl, ucl)[2,]
  
    # Stabilized weights
    denom.fit <- glm(A ~ as.factor(C) + w1 + w2 + as.factor(w3) + as.factor(w4), 
                     family = binomial(), data = data)
    denom.p <- predict(denom.fit, type = "response")    # Stablized Weights  
  
    numer.fit <- glm(A ~ 1, family = binomial(), data = data)
    summary(numer.fit)
    numer.p <- predict(numer.fit, type = "response")    # estimation of numerator of ip weights
  
    data$sw <- ifelse(data$A == 0, ((1-numer.p)/(1-denom.p)), (numer.p/denom.p))
  
    msm <- lm(Y  ~ A, data = data, weights = sw)
  
    SE <-sqrt(diag(vcovHC(msm, type="HC0"))) # robust standard errors
    beta <- coef(msm)
    lcl <- beta-1.96*SE 
    ucl <- beta+1.96*SE
    cbind(beta, lcl, ucl)[2,]  
    
## Double-robust methods
  ## 5.1 IPTW with regression adjustment  
    ### Box 21: Computation of the IPTW-RA estimator for the ATE and bootstrap for statistical inference
    glm1  <- glm(Y ~ C + w1 + w2 + as.factor(w3) + as.factor(w4),  weights = data$w[data$A==1], data=data[data$A==1,])
    Y.1 = predict(glm1,  newdata=data.frame(A = 1, C, w1, w2, as.factor(w3), as.factor(w4)), type="response")
    glm2  <- glm(Y ~ C + w1 + w2 + as.factor(w3) + as.factor(w4),  weights = data$w[data$A==0], data=data[data$A==0,])
    Y.0 = predict(glm2,  newdata=data.frame(A = 0, C, w1, w2, as.factor(w3), as.factor(w4)), type="response")
    ATE <- mean(Y.1 - Y.0); ATE
    ATE2 <- mean(data$w*as.numeric(data$A==1)*Y.1)/mean( data$w*as.numeric(data$A==1)) - mean(data$w*as.numeric(data$A==0)*Y.0)/mean(data$w*as.numeric(data$A==0));ATE2
    rm(ATE, ATE2)
  
    ### Box 22: Computation of the IPTW-RA estimator for the ATE using the ipw R-package
    library(ipw)
    ipw.ATE <- ipwpoint(exposure = A, family = "binomial", link = "logit", 
                        numerator = ~ C, 
                        denominator = ~ C + w1 + w2 + as.factor(w3) + as.factor(w4), 
                        data = data2)
    summary(ipw.ATE$ipw.weights)
    ipwplot(weights = ipw.ATE$ipw.weights, logscale = FALSE, main = "Stabilized weights", xlim = c(0.5, 2))
    summary(ipw.ATE$num.mod)
    summary(ipw.ATE$den.mod)
    
    data2$sw <- ipw.ATE$ipw.weights
    msm <- (svyglm(Y ~ A, design = svydesign(~ 1, weights = ~ sw, data = data2)))
    coef(msm);  confint(msm)
  
  ## 5.2 Augmented inverse probability weighting
    
    ### Box 23: Computation of the AIPTW estimator for the ATE and bootstrap for statistical inference
    mod  <- glm(Y ~ A + C + w1 + w2 + as.factor(w3) + as.factor(w4), family="binomial", data=data)
    PO   <- cbind(Yhat = predict(mod),
                  Y1 = predict(mod, newdata=data.frame(A = 1, C, w1, w2, as.factor(w3), as.factor(w4)), type="response"),
                  Y0 = predict(mod, newdata=data.frame(A = 0, C, w1, w2, as.factor(w3), as.factor(w4)), type="response"))
    RA <- as.data.frame(PO)   # Potential Outcomes
    Yhat <- RA$Yhat
    Y.1a <- RA$Y1
    Y.0a <- RA$Y0
  
    g <- glm(A ~ C + w1 + w2 + as.factor(w3) + as.factor(w4), family = binomial(), data = data)
    gw <- predict(g, type = "response")
    gws <- ifelse(data$A == 0, (-(1 - data$A)/(1 - gw)),(data$A/gw)); sum(gws)    # estimation of weights
    AIPTW <- mean(gws*(data$Y - plogis(RA$Yhat)) + ((Y.1a) - (Y.0a))); AIPTW  # ATE 
    RR <- mean(Y.1a/Y.0a); RR # RR
    
    IC <- (gws*(data$Y - plogis(RA$Yhat)) + ((Y.1a) - (Y.0a)))-AIPTW    # Estimate the influence function (functional Delta method)
    n  <- nrow(data)
    varHat.IC <- var(IC)/n; varHat.IC
    lci <- AIPTW-1.96*sqrt(varHat.IC)
    uci <- AIPTW+1.96*sqrt(varHat.IC)
    cat(AIPTW,lci,uci)      # Inference Influence function
  
    AIPTW.b = function(data,indices)        # Inference using Bootstrap
    {
      dat=data[indices,]
      mod  <- glm(Y ~ A + C + w1 + w2 + as.factor(w3) + as.factor(w4), family="binomial", data=data)
      Yhat = predict(mod)
      Y1 = predict(mod, newdata=data.frame(A = 1, C, w1, w2, as.factor(w3), as.factor(w4)))
      Y0 = predict(mod, newdata=data.frame(A = 0, C, w1, w2, as.factor(w3), as.factor(w4)))
      g <- glm(A ~ C + w1 + w2 + as.factor(w3) + as.factor(w4), family="binomial",  data = data)
      gw <- predict(g,type="response")
      gws <- ifelse(A == 0, (-(1 - A)/(1 - gw)),(A/gw))
      mean(gws*(Y - plogis(Yhat)) + (plogis(Y1) - plogis(Y0)))
    }
    AIPTW.b(data,indices=1:nrow(data))     # Can get original estimate, by plugging in indices 1:n
    boot.out=boot(data,AIPTW.b,200)        # Draw 200 bootstrap sample estimates 
    boot.ci(boot.out,type="perc",conf=0.95)     # compute confidence intervals using percentile method
    boot.ci(boot.out,type="norm",conf=0.95)
    
    ### Box 23 bis: Computation of the AIPTW estimator using drtmle package for the ATE
    w <- subset(data, select=c(C, w1, w2, as.factor(w3), as.factor(w4)))
    fit1 <- drtmle(W = w, A = A, Y = Y, # input data
                   a_0 = c(0, 1), # return estimates for A = 0 and A = 1
                   SL_Q = "SL.npreg", # use kernel regression for E(Y | A = a, W)
                   glm_g = "C + w1 + w2 + w3 + w4", # use misspecified main terms glm for E(A | W)
                   SL_Qr = "SL.npreg", # use kernel regression to guard against
                   # misspecification of outcome regression
                   #SL_gr = "SL.npreg", # use kernel regression to guard against
                   # misspecification of propensity score
                   returnModels = TRUE # for visualizing fits later
    )
    ATE <- ci(fit1, contrast = c(-1,1)); ATE
    rm(ATE)
    
# 6. DATA-ADAPTIVE ESTIMATION: ENSEMBLE LEARNING TARGETED MAXIMUMLIKELIHOOD ESTIMATION
    
    ### Box 25: Computational implementation of TMLE by hand
    # Step 1
    Gcomp <- glm(Y ~ A + C + w1 + w2 + as.factor(w3) + as.factor(w4), family="binomial", data=data2)
    # Prediction for A, A=1 and, A=0
    QAW <- predict(Gcomp)
    Q1W = predict(Gcomp, newdata=data.frame(A = 1, data2[,c("C", "w1","w2","w3","w4")])) 
    Q0W = predict(Gcomp, newdata=data.frame(A = 0, data2[,c("C", "w1","w2","w3","w4")]))
    # Step 2 estimation of the propensity score (ps)
    psm <- glm(A ~ C + w1 + w2 + as.factor(w3) + as.factor(w4), family = binomial, data=data2) 
    gW = predict(psm, type = "response")
    g1W = (1 / gW)
    g0W = (-1 / (1-gW))
    # Step 3 computation of H and estimation of epsilon
    HAW <- (data2$A / gW -(1-data2$A) / (1 - gW))
    H1W = (1/gW)
    H0W = (-1 / (1 - gW))
    epsilon <- coef(glm(data2$Y ~ -1 + HAW + offset(QAW), family = "binomial"))
    # Step 4 update from Q0 to Q1 ATE
    Q1W_1 <- plogis(Q1W + epsilon * H1W) 
    Q0W_1 <- plogis(Q0W + epsilon * H0W)
    # Step 5 targeted estimate of the ATE
    ATE <- mean(Q1W_1 - Q0W_1); ATE
    # Step 6 statistical inference
    d1 <- ((data2$A * (Y - Q1W_1)/gW)) + Q1W_1 - mean(Q1W_1)
    d0 <- ((1 - data2$A) * (Y - Q0W_1)/(1 - gW)) + Q0W_1 - mean(Q0W_1)
    IF <- d1 - d0
    n <- nrow(data2)
    varIF <- var(IF)/n
    LCI <- ATE - 1.96*sqrt(varIF)
    UCI <- ATE + 1.96*sqrt(varIF)
    cbind(ATE, LCI, UCI)
    rm(ATE)
    
    ### Box 26: TMLE with data-adaptive estimation using the R package
    set.seed(777)
    library(tmle)
    w <- subset(data, select=c(C, w1, w2, w3, w4))
    fittmle <- tmle(data$Y, data$A, W=w, family="binomial", 
                    Q.SL.library = c("SL.glm","SL.glm.interaction","SL.step.interaction","SL.gam","SL.randomForest"),
                    g.SL.library = c("SL.glm","SL.glm.interaction","SL.step.interaction","SL.gam","SL.randomForest"))
    fittmle
  
#  7. Simulation
    ### Box 27: Data generation for the Monte Carlo experiment
    
    rm(list=ls())
    
    # Super Learner libraries
    SL.library <- c("SL.glm","SL.step","SL.step.interaction","SL.glm.interaction","SL.gam") #"SL.randomForest","SL.glmnet" 
                  
    # Data generation A: dual misspecification for the model of the outcome and treatment 
    set.seed(7777)
    generateData <- function(n){
      w1 <- round(runif(n, min=1, max=5), digits=0) 
      w2 <- rbinom(n, size=1, prob=0.45)
      w3 <- round(runif(n, min=0, max=1), digits=0 + 0.75*w2 + 0.8*w1)
      w4 <- round(runif(n, min=0, max=1), digits=0 + 0.75*w2 + 0.2*w1)
      A  <- rbinom(n, size=1, prob= plogis(-1 -  0.15*w4 + 1.5*w2 + 0.75*w3 + 0.25*w1 + 0.8*w2*w4))
      # Counterfactuals
      Y.1 <- rbinom(n, size=1, prob = plogis(-3 + 1 + 0.25*w4 + 0.75*w3 + 0.8*w2*w4 + 0.05*w1))
      Y.0 <- rbinom(n, size=1, prob = plogis(-3 + 0 + 0.25*w4 + 0.75*w3 + 0.8*w2*w4 + 0.05*w1)) 
      # Observed outcome
      Y <- Y.1*A + Y.0*(1 - A)
      # return data.frame
      data.frame(w1, w2, w3, w4, A, Y, Y.1, Y.0)
      }
    
    # True ATE
    ObsDataTrueATE <- generateData(n=5000000)
    True_ATE <- mean(ObsDataTrueATE$Y.1 - ObsDataTrueATE$Y.0);True_ATE 
    True_EY.1 <- mean(ObsDataTrueATE$Y.1)
    True_EY.0 <- mean(ObsDataTrueATE$Y.0)
    True_RR <- (True_EY.1 / True_EY.0);True_RR
    
    #Simulations 
    library(tmle) 
    library(SuperLearner) 
    #install.packages("dbarts")
    R <- 1000
    #Empty vectors
    naive_RR <- rep(NA,R) 
    ATEtmle1 <- rep(NA,R) 
    RRtmle1 <- rep(NA,R) 
    ATE_AIPTW   <- rep(NA,R) 
    RR_AIPTW <- rep(NA,R) 
    ATEtmle2 <- rep(NA,R) 
    RRtmle2 <- rep(NA,R) 
    ATEtmle3 <- rep(NA,R) 
    RRtmle3 <- rep(NA,R)
    for(r in 1:R){
      print(paste("This is simulation run number",r)) 
      CancerData <- generateData(n=1000)
      # ATE naive approach
      naive_RR[r] <- exp(glm(data = CancerData, Y ~ A + w1 + w2 + w3 + w4,  family = poisson(link="log"))$coef[2])
      # TMLE implementation by hand
      # Step 1
      gm <- glm(Y ~ A + w1 + w2 + w3 + w4, family="binomial", data=CancerData)
      # Prediction for A, A=1 and, A=0
      QAW <- predict(gm)
      Q1W = predict(gm, newdata=data.frame(A = 1, CancerData[,c("w1","w2","w3","w4")])) 
      Q0W = predict(gm, newdata=data.frame(A = 0, CancerData[,c("w1","w2","w3","w4")]))
      # Step 2 estimation of the propensity score (ps)
      psm <- glm(A ~ w1 + w2 + w3 + w4, family = binomial, data=CancerData) 
      gW = predict(psm, type = "response")
      g1W = (1 / gW)
      g0W = (-1 / (1-gW))
      # Step 3 computation of H and estimation of epsilon
      HAW <- (CancerData$A / gW -(1-CancerData$A) / (1 - gW))
      H1W = (1/gW)
      H0W = (-1 / (1 - gW))
      epsilon <- coef(glm(CancerData$Y ~ -1 + HAW + offset(QAW), family = "binomial"))
      # Step 4 updated ATE
      ATEtmle1[r] <- mean(plogis(Q1W + epsilon * H1W) - plogis(Q0W + epsilon * H0W))
      # Step 5 updated MOR
      T1.EY1 <- mean(plogis(Q1W + epsilon * H1W)) 
      T1.EY0 <- mean(plogis(Q0W + epsilon * H0W))
      RRtmle1[r] <- (T1.EY1 / T1.EY0)
      
      # Augmented inverse probability treatment weight (AIPTW) estimator
      ATE_AIPTW[r] <- mean((HAW*(CancerData$Y - plogis(QAW)) + (plogis(Q1W)-plogis(Q0W))))
      AIPTW1 <- mean(CancerData$A * (CancerData$Y - plogis(Q1W)) / gW + plogis(Q1W) )
      AIPTW0 <- mean((1- CancerData$A) * (CancerData$Y - plogis(Q0W)) / (1-gW) + plogis(Q0W)) 
      RR_AIPTW[r] <- mean( AIPTW1 / AIPTW0)
      
      # R-package tmle (base implementation includes SL.step, SL.glm and SL.glm.interaction) 
      ATE2 <- tmle(Y=CancerData$Y, A=CancerData$A, W=CancerData[,c("w1","w2","w3","w4")], family="binomial")
      ATEtmle2[r] <- ATE2$estimates$ATE$psi
      RRtmle2[r] <- ATE2$estimates$RR$psi
      
      # Improved Super learner
      ATE3 <- tmle(Y = CancerData$Y, A=CancerData$A, W=CancerData[,c("w1","w2","w3","w4")], family="binomial", Q.SL.library=SL.library, g.SL.library=SL.library)
      ATEtmle3[r] <- ATE3$estimates$ATE$psi
      RRtmle3[r] <- ATE3$estimates$RR$psi
    }
    # Mean naive
    mean(naive_RR)
    # Mean AIPTW
    mean(ATE_AIPTW)
    mean(RR_AIPTW)
    # Estimate of TMLE by hand
    mean(ATEtmle1)
    mean(RRtmle1)
    # Estimate of TMLE + SL default implementation
    mean(ATEtmle2)
    mean(RRtmle2)
    # Estimate of TMLE + SL2 default plus more algorithms
    mean(ATEtmle3)
    mean(RRtmle3)
    save.image("your path\results.RData")
    
# Relative Bias ATE
abs(mean((True_ATE - ATE_AIPTW) / True_ATE)*100)
abs(mean((True_ATE - ATEtmle1) / True_ATE)*100)
abs(mean((True_ATE - ATEtmle2) / True_ATE)*100)
abs(mean((True_ATE - ATEtmle3) / True_ATE)*100)

# Relative Bias RR 
abs(mean((True_RR - naive_RR) / True_RR)*100)
abs(mean((True_RR - RR_AIPTW) / True_RR)*100)
abs(mean((True_RR - RRtmle1) / True_RR)*100)
abs(mean((True_RR - RRtmle2) / True_RR)*100)
abs(mean((True_RR - RRtmle3) / True_RR)*100)    
    
    
    
    
    
    
    
    
    
    
    
    
    



