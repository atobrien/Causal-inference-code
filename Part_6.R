# Import data
library(dplyr)
library(multcomp)
library(BSagri)
require(sandwich)
df.original<-read.csv("C:/Users/Me/Desktop/hivpoint.csv")
# Scan names of columns
names(df.original)

# Consider the following list of variables measured at baseline: age, sex, 
# geographic origin, calendar year, mode of HIV transmission, CD4 cell count 
# (a measure of immunosuppression), and viral load. 
# For all questions, assume this full set of confounders is 
# sufficient and necessary to achieve conditional exchangeability 
# for treatment and for censoring. 
# 
# At baseline (confounders/censoring):
#   
# 1.age, 
# 2.sex, 
# 3.geographic origin, 
# 4.calendar year, 
# 5.mode of HIV transmission, 
# 6.CD4 cell count (a measure of immunosuppression), 
# 7.and viral load.
# 
# Predictor:
# 8.treatment
#
# Outcome: 
# 9.logrna

# Review variables of interest
sum(is.na(df.original$treatment)) # 0 NAs categorical 

sum(is.na(df.original$age_0_cat)) # 0 NAs categorical 
sum(is.na(df.original$SEX)) # 0 NAs categorical
sum(is.na(df.original$origin)) # 0 NAs categorical 
sum(is.na(df.original$year_0_cat)) # 0 NAs categorical
sum(is.na(df.original$mode)) # 0 NAs categorical
sum(is.na(df.original$cd4_0_cat)) # 0 NAs categorical
sum(is.na(df.original$rna_0_cat)) # 0 NAs categorical

sum(is.na(df.original$logrna)) # 452 missing/censored & continuous


df.original<-df.original %>% mutate(age_0_cat= as.factor(age_0_cat),
                       origin= as.factor(origin),
                       year_0_cat= as.factor(year_0_cat),
                       mode= as.factor(mode),
                       cd4_0_cat= as.factor(cd4_0_cat),
                       rna_0_cat= as.factor(rna_0_cat))




#______________________________________________________________________________#
#______________________________________________________________________________#
#______________________________________________________________________________#
#______________________________________________________________________________#
#______________________________________________________________________________#
# PART 1
# Question 2

#
###############################################################################
# CENSORING WEIGHT ESTIMATION #################################################
###############################################################################
#
# Estimation of denominator of ip weights for not being censored
denom.cens <- glm(censoring ~                       
                    as.factor(treatment)+
                    as.factor(age_0_cat) + 
                    as.factor(SEX) + 
                    as.factor(origin) + 
                    as.factor(year_0_cat) +
                    as.factor(mode) +
                    as.factor(cd4_0_cat) +
                    as.factor(rna_0_cat),
                  family = binomial(link="logit"), 
                  data = df.original)
summary(denom.cens)
denom.p.cens <- predict(denom.cens, df.original, type = "response")

# Compute non-stabilized weight for not being ccccensored
df.original$w.c <- ifelse(df.original$cens == 0, 
                          ((1)/(1-denom.p.cens)),
                          1)  
#
################################################################################
# NON-STABILIZED IP WEIGHTS ####################################################
################################################################################
#
#
# Estimation of denominator of stabilized ip weights
denom.fit <- glm(treatment ~ 
                   as.factor(age_0_cat) + 
                   as.factor(SEX) + 
                   as.factor(origin) + 
                   as.factor(year_0_cat) +
                   as.factor(mode) +
                   as.factor(cd4_0_cat) +
                   as.factor(rna_0_cat),
                 family = binomial(link="logit"), 
                 data = df.original)
denom.p <- predict(denom.fit, df.original,type = "response")

# Compute non-stabilized weightcs
df.original$obs.ip <- ifelse(df.original$treatment == 0, 
                         1-denom.p,
                         denom.p)
df.original$w.ip <- 1/df.original$obs.ip
summary(df.original$w.ip)

# Compute the combined weights for censoring and ip weights
df.original$w.cip <- df.original$w.ip * df.original$w.c

################################################################################
# Obtaining final estimates ####################################################
################################################################################
glm.obj.q2 <- glm(logrna~
                 as.factor(treatment) + 
                 cluster(id), 
               data = df.original, 
               weights = w.cip)
summary(glm.obj.q2)
comp<-glht(glm.obj.q2)
CIadj<-CIGLM(comp,method="Adj")
CIadj
vcov(glm.obj.q2)

beta <- coef(glm.obj.q2)
SE <-sqrt(diag(vcovHC(glm.obj.q2, type="HC0"))) # robust standard errors
lcl <- beta-1.96*SE 
ucl <- beta+1.96*SE
p.value <- 2*(1-pnorm(abs(beta/SE)))
round(cbind(beta, SE, lcl, ucl, p.value),3)[2,]

#______________________________________________________________________________#
#______________________________________________________________________________#
#______________________________________________________________________________#
#______________________________________________________________________________#
#______________________________________________________________________________#
# Question 3

require(geepack)

df.0<-subset(df.original, censoring == 0)
data <- df.0
grid <- seq(from = -2,to = 1, by = 0.01) # set by = 0.001 for finer estimate

j = 0
store.Hpsi.coefs <- double(length(grid))
for (i in grid){
  psi = i
  j = j+1
  data$Hpsi <- data$logrna - psi * data$treatment 
  
gee.obj.q3 <- geeglm(treatment ~
                      as.factor(age_0_cat) + 
                      as.factor(SEX) + 
                      as.factor(origin) + 
                      as.factor(year_0_cat) +
                      as.factor(mode) +
                      as.factor(cd4_0_cat) +
                      as.factor(rna_0_cat)+
                      Hpsi, 
                    data = data,
                    weight = w.c, 
                    id=id, 
                    corstr="independence", 
                    family = binomial(logit))
  store.Hpsi.coefs[j] <- coef(gee.obj.q3)["Hpsi"]
  cat("Iteration", j, "completed\n")
}

store.results <- as.data.frame(cbind(grid, abs(store.Hpsi.coefs)))
names(store.results)
names(store.results) <- c("grid", "Hpsi.est")  
store.results[store.results$Hpsi.est == min(store.results$Hpsi.est),]

# Using dplyr
# calipered range
store.results %>% filter((floor(((store.results$Hpsi.est*1000)))/1000) %in% c(seq(0.01, 0.06, 0.0000001)))
# upper boundries (UCI 95%)
head(store.results %>% filter((floor(((store.results$Hpsi.est*1000)))/1000) %in% c(seq(0.01, 0.06, 0.0000001))), 1)
# lower boundries (LCI 95%)
tail(store.results %>% filter((floor(((store.results$Hpsi.est*1000)))/1000) %in% c(seq(0.01, 0.06, 0.0000001))), 1)

#______________________________________________________________________________#
#______________________________________________________________________________#
#______________________________________________________________________________#
#______________________________________________________________________________#
#______________________________________________________________________________#
# Question 4

lm.q4 <- lm(logrna ~
           as.factor(treatment)+
           as.factor(age_0_cat) + 
           as.factor(SEX) +
           as.factor(origin) + 
           as.factor(year_0_cat) +
           as.factor(mode) +
           as.factor(cd4_0_cat) +
           as.factor(rna_0_cat),
         data = df.original,
         weights = w.c)
summary(lm.q4)
confint(lm.q4)

comp<-glht(lm.q4)
CIadj<-CIGLM(comp,method="Adj")
CIadj
vcov(lm.q4)

beta <- coef(lm.q4)
SE <-sqrt(diag(vcovHC(lm.q4, type="HC0"))) # robust standard errors
lcl <- beta-1.96*SE 
ucl <- beta+1.96*SE
p.value <- 2*(1-pnorm(abs(beta/SE)))
round(cbind(beta, SE, lcl, ucl, p.value),3)[2,]

#______________________________________________________________________________#
#______________________________________________________________________________#
#______________________________________________________________________________#
#______________________________________________________________________________#
#______________________________________________________________________________#
# Question 5

# Propensity score
glm.ps <- glm(treatment ~ 
                as.factor(age_0_cat) + 
                as.factor(SEX) + 
                as.factor(origin) + 
                as.factor(year_0_cat) +
                as.factor(mode) +
                as.factor(cd4_0_cat) +
                as.factor(rna_0_cat),
              family = binomial(link="logit"), 
              data = df.original)

## Extract propensity score 
df.original$pscores <- fitted(glm.ps)

# Fit linear outcome model with propensity score
lm.q5 <- lm(logrna ~
              as.factor(treatment)+
              pscores,
            data = df.original,
            weights = w.c)
summary(lm.q5)
confint(lm.q5)

comp<-glht(lm.q5)
CIadj<-CIGLM(comp,method="Adj")
CIadj
vcov(lm.q5)

beta <- coef(lm.q5)
SE <-sqrt(diag(vcovHC(lm.q5, type="HC0"))) # robust standard errors
lcl <- beta-1.96*SE 
ucl <- beta+1.96*SE
p.value <- 2*(1-pnorm(abs(beta/SE)))
round(cbind(beta, SE, lcl, ucl, p.value),3)[2,]

#______________________________________________________________________________#
#______________________________________________________________________________#
#______________________________________________________________________________#
#______________________________________________________________________________#
#______________________________________________________________________________#
# Question 6

# create table 1 with indicator
df.original$interv <- rep(-1, nrow(df.original))

#create table 2 for untreated and label outcome as NA
df.untreat <- df.original
df.untreat$interv <- 0
df.untreat$treatment<- 0
df.untreat$logrna <- rep(NA, nrow(df.untreat))

#create table 3 for untreated and label outcome as NA
df.treat <- df.original
df.treat$interv <- 1
df.treat$treatment<- 1
df.treat$logrna <- rep(NA, nrow(df.treat))

# create a dataset with 3 copies of each subject
df.q6<- as.data.frame(rbind(df.original, df.untreat, df.treat))

# estimate non-parametric model
fit.q6<-glm(logrna~
              as.factor(treatment) +
              as.factor(age_0_cat) + 
              as.factor(SEX) + 
              as.factor(origin) + 
              as.factor(year_0_cat) +
              as.factor(mode) +
              as.factor(cd4_0_cat) +
              as.factor(rna_0_cat),
           data=df.q6)
summary(fit.q6)

# predict outcome Y in tables 2 and 3 for outcome
df.q6$meanY <- predict(fit.q6, df.q6)

# get the means of the tables
with(df.q6, tapply(meanY, list(interv), mean))

#############################95% CI BOOTSTRAP #################################
B <- 100 # 10,000 bootstraps
N <- 5542 #replicate the size of a table
start.time <- Sys.time() 
mean_list<-replicate(B,{
  tab2<-sample_n(df.untreat,N,replace=TRUE) #create table 2
  tab3<-sample_n(df.treat,N,replace=TRUE) #create table 3
  tab1<-as.data.frame(rbind(df.original, tab2, tab3)) #fuse tables
  fit.q6BS<-glm(logrna~
                  as.factor(treatment) +
                  as.factor(age_0_cat) + 
                  as.factor(SEX) + 
                  as.factor(origin) + 
                  as.factor(year_0_cat) +
                  as.factor(mode) +
                  as.factor(cd4_0_cat) +
                  as.factor(rna_0_cat),
                data=df.q6) #fit model
  tab1$meanY <- predict(fit.q6BS, tab1) #predict death
  results<-with(tab1, tapply(meanY, list(interv), mean))# get mean of predictions
  return(results) #store results in matrix and repeat
})
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
mean_df<-as.data.frame(mean_list)
ceff<-mean_df[3,]-mean_df[2,]
tceff<-t(ceff)
t.test(tceff)

#______________________________________________________________________________#
#______________________________________________________________________________#
#______________________________________________________________________________#
#______________________________________________________________________________#
#______________________________________________________________________________#
# Question 7

# create table 1 with indicator
df.original$interv <- rep(-1, nrow(df.original))

#create table 2 for untreated and label outcome as NA
df.untreat <- df.original
df.untreat$interv <- 0
df.untreat$treatment<- 0
df.untreat$logrna <- rep(NA, nrow(df.untreat))

#create table 3 for untreated and label outcome as NA
df.treat <- df.original
df.treat$interv <- 1
df.treat$treatment<- 1
df.treat$logrna <- rep(NA, nrow(df.treat))

# create a dataset with 3 copies of each subject
df.q7<- as.data.frame(rbind(df.original, df.untreat, df.treat))

# estimate non-parametric model
fit.q7<-glm(logrna~
              as.factor(treatment) +
              pscores,
            data=df.q7)
summary(fit.q7)

# predict outcome Y in tables 2 and 3 for outcome
df.q7$meanY <- predict(fit.q7, df.q7)

# get the means of the tables
with(df.q7, tapply(meanY, list(interv), mean))

#############################95% CI BOOTSTRAP #################################
B <- 100 # 10,000 bootstraps
N <- 5542 #replicate the size of a table
start.time <- Sys.time() 
mean_list<-replicate(B,{
  tab2<-sample_n(df.untreat,N,replace=TRUE) #create table 2
  tab3<-sample_n(df.treat,N,replace=TRUE) #create table 3
  tab1<-as.data.frame(rbind(df.original, tab2, tab3)) #fuse tables
  fit.q7BS<-glm(logrna~
                  as.factor(treatment) +
                  pscores,
                data=df.q7) #fit model
  tab1$meanY <- predict(fit.q7BS, tab1) #predict death
  results<-with(tab1, tapply(meanY, list(interv), mean))# get mean of predictions
  return(results) #store results in matrix and repeat
})
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
mean_df<-as.data.frame(mean_list)
ceff<-mean_df[3,]-mean_df[2,]
tceff<-t(ceff)
t.test(tceff)


#______________________________________________________________________________#
#______________________________________________________________________________#
#______________________________________________________________________________#
#______________________________________________________________________________#
#______________________________________________________________________________#
# PART 2
# Question 3

df.part2<-read.csv("C:/Users/Me/Desktop/ohie.csv")

library(sem) 
model1 <- tsls(bmi ~ medicaid, 
               ~ lottery, 
               data = df.part2)
summary(model1)
confint(model1) 

comp<-glht(model1)
CIadj<-CIGLM(comp,method="Adj")
CIadj
vcov(model1)

beta <- coef(model1)
SE <-sqrt(diag(vcovHC(model1, type="HC0"))) # robust standard errors
lcl <- beta-1.96*SE 
ucl <- beta+1.96*SE
p.value <- 2*(1-pnorm(abs(beta/SE)))
round(cbind(beta, SE, lcl, ucl, p.value),3)[2,]


################################################################################
# first try it "manually"
a <- mean(df.part2$bmi[df.part2$lottery == 1], na.rm = TRUE) # E[Y|Z=1]
a # print
b <- mean(df.part2$bmi[df.part2$lottery == 0], na.rm = TRUE) # E[Y|Z=0]
b # print
c <- mean(df.part2$medicaid[df.part2$lottery == 1], na.rm = TRUE)    # Pr(A=1|Z=1)
c # print
d <- mean(df.part2$medicaid[df.part2$lottery == 0], na.rm = TRUE)    # Pr(A=1|Z=0)
d # print
(a-b)/(c-d)                                                  # the IV estimate

# then wrap it in a function
wald <- function(outcome, treatment, iv, data) {
  # the numerator, i.e., E[Y|Z=1] - E[Y|Z=0]
  num <- coef(lm(data[,outcome] ~ data[,iv]))[2]
  # the denominator, i.e., Pr[A=1|Z=1] - Pr[A=1|Z=0]
  denom <- coef(lm(data[,treatment] ~ data[,iv]))[2]
  return(as.numeric(num/denom))
}

wald('bmi', 'medicaid', 'lottery', df.part2)

###############################################################################
#P2-Q3
table(df.part2$lottery, df.part2$medicaid)
summary(table(df.part2$lottery, df.part2$medicaid))

