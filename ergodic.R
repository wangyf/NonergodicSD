#plot
library(ggplot2)
library(patchwork)
library(INLA)
library(dplyr)
library(lme4)

rm(list = ls())

setwd("~/Projects/ngmm_tools/nonergodicSD/")

flatfile_fname <- 'data/new_SD_PGA.csv'
source('R_lib/regression/inla/regression_inla_model1_unbounded_hyp.R')

#load flatfile
utmzone = 10
df_flatfile <- read.csv(flatfile_fname)
names(df_flatfile) <- c('eqid','date','eqlat','eqlon','eqZ','mag','SD','SDD','Site','R','Vs30','Vs30class','PGA','stalat','stalon')
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





Y <- df_flatfile[,'PGA']
M <- df_flatfile[,'mag']
R <- df_flatfile[,'R']




# linear regression using lm package
#lmod <- lm(PGA ~ mag + log(R), data=df_flatfile)
#summary(lmod)
#df_flatfile$predicted<-predict(lmod)
#df_flatfile$residuals<-residuals(lmod)
#ggplot()+
#  geom_point(df_flatfile,mapping=aes(x=mag,y=PGA),color='black') +
#  geom_point(df_flatfile,mapping=aes(x=mag,y=predicted),color='red')
#  stat_smooth(method = "lm", color = 'red')




#fm1 <- lmer(data=df_flatfile,formula='PGA ~ 1 + mag + log(R) + (1|eqid) + (1|Site)')
#summary(fm1)
#total residual
#df_flatfile$residuals<-residuals(fm1)
#ggplot()+
#    geom_point(df_flatfile,mapping=aes(x=mag,y=residuals),color='black') 



fit_inla <-inla(PGA ~ 1 + mag + log(R) + 
                  f(eqid, model="iid") + f(Site,model='iid'), data=df_flatfile,
                num.threads = 4,quantiles = c(0.05,0.5,0.95),
                control.predictor=list(compute=TRUE))
summary(fit_inla)

event<-df_flatfile %>% group_by(eq_id) %>% filter(row_number()==1)
dataM=event[,'mag']
deltaB<-fit_inla$summary.random$eqid

station<-df_flatfile %>% group_by(sta_id) %>% filter(row_number()==1)
dataS=station[,'sta_id']
deltaS<-fit_inla$summary.random$Site


df_plot <-cbind(dataM,deltaB)
names(df_plot)[c(1,2)] <- c('M','DeltaB')
names(df_plot)[c(6,7,8)] <- c('q05','q50','q95')
p1<-ggplot(df_plot,aes(x=M,y=q50))+geom_point()

dS <- deltaS$`0.5quant`
SS <- dataS$sta_id
df_plot<-data.frame(SS,dS)
p2<-ggplot(df_plot,aes(x=SS,y=dS))+geom_point()

y_pred <-fit_inla$summary.fitted.values
df_flatfile[,'resid']<-df_flatfile$PGA - y_pred$mean
p3<-ggplot(df_flatfile,aes(x=mag,y=resid))+geom_point()

p1+p2+p3
