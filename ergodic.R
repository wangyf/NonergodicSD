#plot
library(ggplot2)
library(patchwork)
library(INLA)
library(dplyr)
library(lme4)
library("ggpubr")

rm(list = ls())

setwd("/home/yow004/Projects/NonergodicSD")

#flatfile_fname <- 'data/new_SD_PGA.csv'

flatfile_fname<- "/home/yow004/Projects/NonergodicSD/data/bayarea.groundmotion.fitted.csv"
source('R_lib/regression/inla/regression_inla_model1_unbounded_hyp.R')

#load flatfile
utmzone = 10
df_flatfile <- read.csv(flatfile_fname)
names(df_flatfile) <- c('eqid','date','eqlat','eqlon','eqZ','mag','SD','SDD','Site','R','Vs30','Vs30class','PGA','pPGA','qPGA','sPGA','stalat','stalon')
df_flatfile$ssn <- df_flatfile$Site
df_flatfile$UTMzone <- utmzone #north california

#df_flatfile$PGA = df_flatfile$PGA/log(10)

df_flatfile<- subset(df_flatfile,mag<6)

eq<-LongLatToUTM(df_flatfile$eqlat,df_flatfile$eqlon,utmzone)
df_flatfile[,c('eqX','eqY')] <- eq[,c('X','Y')]/1000

sta<-LongLatToUTM(df_flatfile$stalat,df_flatfile$stalon,utmzone)
df_flatfile[,c('staX','staY')] <- sta[,c('X','Y')]/1000

# Preprocess Input Data
# ---------------------------
n_data <- nrow(df_flatfile)
#earthquake data
data_eq_all <- df_flatfile[,c('eqid','mag','eqX', 'eqY')]
out_unq  <- UniqueIdxInv(df_flatfile[,'eqid'])
eq_idx   <- out_unq$idx
eq_inv   <- out_unq$inv
data_eq  <- data_eq_all[eq_idx,]
X_eq     <- data_eq[,c(3,4)] #earthquake coordinates
X_eq_all <- data_eq_all[,c(3,4)]
#create earthquake ids for all records (1 to n_eq)
eq_id <- eq_inv
n_eq  <- nrow(data_eq)



#station data
data_sta_all <- df_flatfile[,c('ssn','Vs30','staX','staY')]
out_unq   <- UniqueIdxInv(df_flatfile[,'ssn'])
sta_idx   <- out_unq$idx
sta_inv   <- out_unq$inv
data_sta  <- data_sta_all[sta_idx,]
X_sta     <- data_sta[,c(3,4)] #station coordinates
X_sta_all <- data_sta_all[,c(3,4)]
#create station indices for all records (1 to n_sta)
sta_id <- sta_inv
n_sta  <- nrow(data_sta)

##
df_flatfile[,c('eq_id','sta_id')] = c(eq_id,sta_id)

#df <- na.omit(df_flatfile)
df <- df_flatfile
Vref = 760
y     <- df[,'PGA']
M     <- df[,'mag']
M2    <- (8.5-M)**2
R     <- df[,'R']
lnRef <- log((R**2+4.5**2)**0.5)
lnvs  <- log(df[,'Vs30']/Vref)
eqid  <- df[,'eq_id']
stid  <- df[,'sta_id']

inladata <- data.frame(M,M2,lnRef,R,lnvs,y,eqid,stid)

fit_inla2 <-inla(y ~ 1 + M + M2 + lnRef + R + lnvs +
                   f(eqid, model="iid") + f(stid, model="iid"), data=inladata,
                 num.threads = 8,quantiles = c(0.05,0.5,0.95),
                 control.predictor=list(compute=TRUE),verbose = TRUE)

#fit_inla2 <-inla(y ~ 1 + M + log(R) +
#                  f(eqid, model="iid"), data=inladata,
#                num.threads = 4,quantiles = c(0.05,0.5,0.95),
#                control.predictor=list(compute=TRUE),verbose = TRUE)

summary(fit_inla2)



event<-inladata %>% group_by(eqid) %>% filter(row_number()==1)
dataM=event[,'M']
deltaB<-fit_inla2$summary.random$eqid

event2<-df_flatfile %>% group_by(eqid) %>% filter(row_number()==1)
sig<-event2$SD
dE_trugman<-event2$qPGA


station<-inladata %>% group_by(stid) %>% filter(row_number()==1)
dataS=station[,'stid']
deltaS<-fit_inla2$summary.random$stid

dE <- deltaB$`0.5quant`
df_E<-data.frame(dataM,dE)
p0<-ggplot(df_E,aes(x=dataM,y=dE))+geom_point()+
  stat_summary_bin(bins=10,color='red',width=1,geom='errorbar',fun.min=~quantile(.x,probs = .25),fun.max=~quantile(.x,probs = .75))+
  stat_summary_bin(bins=10,color='red',geom='point',fun=median,size=4)



#create list




dS <- deltaS$`0.5quant`
SS <- dataS$stid


df_S<-data.frame(SS,dS)
p1<-ggplot(df_S,aes(x=SS,y=dS))+geom_point()


y_pred <-fit_inla2$summary.fitted.values
inladata[,'resid']<-inladata$y - y_pred$mean



p2<-ggplot(inladata,aes(x=M,y=resid))+geom_point()+
  stat_summary_bin(bins=10,color='red',width=1,geom='errorbar',fun.min=~quantile(.x,probs = .25),fun.max=~quantile(.x,probs = .75))+
  stat_summary_bin(bins=10,color='red',geom='point',fun=median,size=4)

p3<-ggplot(inladata,aes(x=R,y=resid))+geom_point()+
  stat_summary_bin(bins=10,color='red',width=1,geom='errorbar',fun.min=~quantile(.x,probs = .25),fun.max=~quantile(.x,probs = .75))+
  stat_summary_bin(bins=10,color='red',geom='point',fun=median,size=4)

p4<-ggplot(inladata,aes(x=exp(lnvs)*Vref,y=resid))+geom_point()+
  stat_summary_bin(bins=10,color='red',width=1,geom='errorbar',fun.min=~quantile(.x,probs = .25),fun.max=~quantile(.x,probs = .75))+
  stat_summary_bin(bins=10,color='red',geom='point',fun=median,size=4)


p0+p1+p2+p3+p4




cofit<-fit_inla2$summary.fixed$mean
pred<-matrix(nrow=nrow(inladata),ncol=1)
event_term<-matrix(nrow=nrow(inladata),ncol=1)
site_term<-matrix(nrow=nrow(inladata),ncol=1)
for (i in 1:nrow(inladata)){
  pred[i,]<-cofit %*% c(1, inladata[i,'M'],inladata[i,'M2'],inladata[i,'lnRef'],inladata[i,'R'],inladata[i,"lnvs"])
  event_term[i,]<-dE[eq_id[i]]
  site_term[i,]<-dS[sta_id[i]]
}

pred_all<-pred+event_term+site_term
pred_fixed<-y_pred$mean-event_term-site_term
tmp<-data.frame(inladata$y,pred,event_term,site_term,y_pred$mean,pred_all,pred_fixed)

total_resid = inladata$y - pred_fixed

#ggplot(as.data.frame(M,total_resid),aes(x=M,y=total_resid))+geom_point()

df_plot<-data.frame(dE,sig)
ggplot(df_plot,aes(x=dE,y=sig))+geom_point()


res <- cor.test(df_plot$dE, df_plot$sig, 
                method = "pearson")
res

p1<-ggscatter(df_plot, x = "dE", y = "sig", 
          add = "reg.line", conf.int = FALSE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Event Term", ylab = "Stress Drop")



df_plot<-data.frame(dE_trugman,sig)
res <- cor.test(df_plot$dE_trugman, df_plot$sig, 
                method = "pearson")
res
p2<-ggscatter(df_plot, x = "dE_trugman", y = "sig", 
              add = "reg.line", conf.int = FALSE, 
              cor.coef = TRUE, cor.method = "pearson",
              xlab = "Event Term", ylab = "Stress Drop")

p1+p2
