
library(rstan)
library(data.table)
# setwd("~/Downloads/HFMD")

case1 <- readRDS("case_data_preprocessed.rds")
casepo <- readRDS("hochiminhcity_population.rds")

case1 <- as.data.table(case1)
casepo <- as.data.table(casepo)

lambda_index_group <- readRDS("lambda_index_group.rds")
lambda_index_group <- as.data.table(lambda_index_group)

#### model LinPw ####
model_fit <- readRDS("model_linpw_fit.rds")
po <- rstan::extract(model_fit)
##### convergence evaluation #####
model_su <- summary(model_fit)$summary
min(model_su[,"n_eff"]) # minimum effective sample size
max(model_su[,"Rhat"]) # maximum effective sample size

trace_name <- rownames(model_su[order(model_su[,"n_eff"]),][1:4,])
# Trace plot
po_trace <- rstan:::extract(model_fit, inc_warmup = TRUE, permuted = FALSE,
                            pars = trace_name)
bayesplot:::color_scheme_set("mix-blue-pink")
trace <- bayesplot:::mcmc_trace(po_trace, n_warmup = 1e3,
                                facet_args = list(nrow=1,ncol=4,labeller = label_parsed))+
  theme_bw()+
  theme(legend.position='bottom')

##### force of infection #####
log_beta0 <- po$log_beta0
log_beta1 <- po$log_beta1
log_beta2 <- po$log_beta2

lambda_index_group_apre <-lambda_index_group[,{
  lambda1 <- exp(log_beta0 + log_beta1[,time_virus_index] + log_beta2[,age_index] )
  est = quantile(lambda1,probs = c(0.5, 0.025, 0.975))
  list(mean=mean(lambda1),median=est[1],CI_L=est[2],CI_U=est[3])
},by=c("age_group","year_group","virus","index")]

lambda_index_group_apre[,virus1 := factor(virus,levels = c("CV-A6","CV-A10","CV-A16","EV-A71"),labels = c("CV-A6","CV-A10","CV-A16","EV-A71"))]

saveRDS(lambda_index_group_apre,"est_foi_model_linpw.rds")

##### estimated cases #####
xi0 <- po$xi0
log_one_over_phi1 <- po$log_one_over_phi

case1[,age:=floor(`Age in years`)]
case1_0 <- case1[which(case1[,age==0]),c(2,5,6,7,10,13)]
case1_0[,month:=0]
set(case1_0,which(case1_0[,(`Age in years`>=1/12)&(`Age in years`<2/12)]),"month",1)
set(case1_0,which(case1_0[,(`Age in years`>=2/12)&(`Age in years`<3/12)]),"month",2)
set(case1_0,which(case1_0[,(`Age in years`>=3/12)&(`Age in years`<4/12)]),"month",3)
set(case1_0,which(case1_0[,(`Age in years`>=4/12)&(`Age in years`<5/12)]),"month",4)
set(case1_0,which(case1_0[,(`Age in years`>=5/12)&(`Age in years`<6/12)]),"month",5)
set(case1_0,which(case1_0[,(`Age in years`>=6/12)&(`Age in years`<7/12)]),"month",6)
set(case1_0,which(case1_0[,(`Age in years`>=7/12)&(`Age in years`<8/12)]),"month",7)
set(case1_0,which(case1_0[,(`Age in years`>=8/12)&(`Age in years`<9/12)]),"month",8)
set(case1_0,which(case1_0[,(`Age in years`>=9/12)&(`Age in years`<10/12)]),"month",9)
set(case1_0,which(case1_0[,(`Age in years`>=10/12)&(`Age in years`<11/12)]),"month",10)
set(case1_0,which(case1_0[,(`Age in years`>=11/12)&(`Age in years`<1)]),"month",11)

case2_0 <- case1_0[,list(count=.N),by=c("virus","year_sampling","age","month")]
colnames(case2_0)[2]="year_group"
colnames(case2_0)[3]="age_group"

tmp <- expand.grid(virus=c("CV-A6","CV-A10","CV-A16","EV-A71"),year_group=2013:2018,age_group=0,month=1:11)
tmp$self_index <- 1:dim(tmp)[1]
tmp <- as.data.table(tmp)

case2_0 <- merge(case2_0,tmp,by=c("virus","year_group","age_group","month"))

zero_index1 <- tmp[which(tmp[,!self_index%in%case2_0$self_index]),]
zero_index1[,count:=0]

case3_0 <- rbind(case2_0,zero_index1)
case3_0 <- case3_0[order(case3_0[,self_index]),]

case4_0 <- merge(case3_0,lambda_index_group,by=c("virus","year_group","age_group"))
case4_0 <- case4_0[order(case4_0[,self_index]),]

colnames(casepo)[1]="age_group"
colnames(casepo)[2]="year_group"

case5_0 <- merge(case4_0,casepo[,1:3],by=c("age_group","year_group"))

virus_index <- data.table(virus=c("CV-A6","CV-A10","CV-A16","EV-A71"),virus_index=1:4)
case6_0 <- merge(case5_0,virus_index,by="virus")
case6_0 <- case6_0[order(case6_0[,self_index]),]

case_0_pre <- case6_0[,{
  xi_this = xi0[,virus_index]+1
  phi_this = exp(-log_one_over_phi1[,virus_index])
  lambda1 = exp(log_beta0 + log_beta1[,time_virus_index] + log_beta2[,age_index])
  case_pre1 = pop*xi_this/lambda1*exp(-month*lambda1/12)*(exp(min(month/12,1/xi_this)*lambda1)-1)*(1-exp(-lambda1/12))*phi_this
  est = quantile(case_pre1,probs = c(0.5, 0.025, 0.975))
  list(mean=mean(case_pre1),median=est[1],CI_L=est[2],CI_U=est[3])
},by=c("age_group","year_group","virus","month","count","self_index")]


case7_0 <- case6_0[,list(count1=sum(count)),by=c("age_group","year_group","virus","index","time_virus_index","age_index","virus_index","pop")]
i=1
xi_this = xi0[,case7_0[i,]$virus_index]+1
phi_this = exp(-log_one_over_phi1[,case7_0[i,]$virus_index])
lambda1 = exp(log_beta0 + log_beta1[,case7_0[i,]$time_virus_index] + log_beta2[,case7_0[i,]$age_index])
case_pre1 = 0
for (x in 1:11){
  case_pre2 = case7_0[i,]$pop*xi_this/lambda1*exp(-x*lambda1/12)*(exp(min(x/12,1/xi_this)*lambda1)-1)*(1-exp(-lambda1/12))*phi_this
  case_pre1 = case_pre1+case_pre2
}
est = quantile(case_pre1,probs = c(0.5, 0.025, 0.975))
case_0_pre1 = data.table(mean=mean(case_pre1),median=est[1],CI_L=est[2],CI_U=est[3],age_group=0,year_group=case7_0[i,]$year_group,virus=case7_0[i,]$virus,count=case7_0[i,]$count1)
for (i in 2:dim(case7_0)[1]){
  xi_this = xi0[,case7_0[i,]$virus_index]+1
  phi_this = exp(-log_one_over_phi1[,case7_0[i,]$virus_index])
  lambda1 = exp(log_beta0 + log_beta1[,case7_0[i,]$time_virus_index] + log_beta2[,case7_0[i,]$age_index])
  case_pre1 = 0
  for (x in 1:11){
    case_pre2 = case7_0[i,]$pop*xi_this/lambda1*exp(-x*lambda1/12)*(exp(min(x/12,1/xi_this)*lambda1)-1)*(1-exp(-lambda1/12))*phi_this
    case_pre1 = case_pre1+case_pre2
  }
  est = quantile(case_pre1,probs = c(0.5, 0.025, 0.975))
  case_0_pre2 = data.table(mean=mean(case_pre1),median=est[1],CI_L=est[2],CI_U=est[3],age_group=0,year_group=case7_0[i,]$year_group,virus=case7_0[i,]$virus,count=case7_0[i,]$count1)
  case_0_pre1 = rbind(case_0_pre1,case_0_pre2)
}



case1_1 <- case1[which(case1[,age>0]),c(2,5,6,7,10,13)]
case1_1[,age_group:=age]
set(case1_1,which(case1_1[,age_group>=7]),"age_group",7)
case2_1 <- case1_1[,list(count=.N),by=c("virus","year_sampling","age_group")]
colnames(case2_1)[2]="year_group"



tmp2 <- expand.grid(virus=c("CV-A6","CV-A10","CV-A16","EV-A71"),year_group=2013:2018,age_group=1:7)
tmp2$self_index <- 1:dim(tmp2)[1]
tmp2 <- as.data.table(tmp2)

case2_1 <- merge(case2_1, tmp2, by=c("virus","year_group","age_group"))
zero_index2 <- tmp2[which(tmp2[,!self_index%in%case2_1$self_index]),]
zero_index2[,count:=0]

case3_1 <- rbind(case2_1,zero_index2)
case3_1 <- case3_1[order(case3_1[,self_index]),]

case4_1 <- merge(case3_1,lambda_index_group,by=c("virus","year_group","age_group"))
case4_1 <- case4_1[order(case4_1[,self_index]),]

virus_index
case5_1 <- merge(case4_1,virus_index,by="virus")
case5_1 <- case5_1[order(case5_1[,self_index]),]

pop_group <- expand.grid(year_group=2013:2018,age_group=1:7)
pop_group$pop <- NA
for (i in 1:dim(pop_group)[1]){
  age <- pop_group[i,]$age_group
  year1 <- pop_group[i,]$year_group
  if (age<7){
    tmp = data.table(age_group=age,year_group=year1)
    tmp1 = merge(tmp,casepo,by=c("age_group","year_group"))
    pop_group[i,]$pop = tmp1$pop
  }else{
    tmp = expand.grid(age_group=7:11,year_group=year1)
    tmp1 = merge(tmp,casepo,by=c("age_group","year_group"))
    pop_group[i,]$pop = sum(tmp1$pop)
  }
}

case6_1 <- merge(case5_1,pop_group,by=c("year_group","age_group"))
case6_1 <- case6_1[order(case6_1[,self_index]),]


i=1
case6_1[i,]
virus1 <- case6_1[i,]$virus
age1 <- case6_1[i,]$age_group
year1 <- case6_1[i,]$year_group

if (age1==7) age1=9

data1 <- data.table(year_group = (year1-age1):(year1-1),age_group = 0:(age1-1),virus=virus1)
set(data1,which(data1[,age_group>7]),"age_group",7)
set(data1,which(data1[,year_group<2011]),"year_group",2011)
data1[,order1 := 1:dim(data1)[1]]
data1 = merge(data1,lambda_index_group,by=c("virus","age_group","year_group"))
data1 <- data1[order(data1[,order1]),]

sum_lambda=0
for (j in 1:dim(data1)[1]){
  sum_lambda1 <- exp(log_beta0 + log_beta1[,data1[j,]$time_virus_index] + log_beta2[,data1[j,]$age_index])
  sum_lambda <- sum_lambda + sum_lambda1
}

xi_this <- xi0[,case6_1[i,]$virus_index]+1
phi_this <- exp(-log_one_over_phi1[,case6_1[i,]$virus_index])
stopifnot(data1[1,]$age_group==0)
lambda0_this <- exp(log_beta0 + log_beta1[,data1[1,]$time_virus_index] + log_beta2[,data1[1,]$age_index])
lambda_this <- exp(log_beta0 + log_beta1[,case6_1[i,]$time_virus_index] + log_beta2[,case6_1[i,]$age_index])
case_pre1 <- case6_1[i,]$pop*xi_this/lambda0_this*(exp(lambda0_this/xi_this)-1)*exp(-sum_lambda)*(1-exp(-lambda_this))*phi_this
est = quantile(case_pre1,probs = c(0.5, 0.025, 0.975))
case_1_pre <- data.table(mean=mean(case_pre1),median=est[1],CI_L=est[2],CI_U=est[3],age_group=age1,year_group=year1,virus=virus1,count=case6_1[i,]$count)


for (i in 2:dim(case6_1)[1]){
  virus1 <- case6_1[i,]$virus
  age1 <- case6_1[i,]$age_group
  year1 <- case6_1[i,]$year_group
  
  if (age1==7) age1=9
  
  data1 <- data.table(year_group = (year1-age1):(year1-1),age_group = 0:(age1-1),virus=virus1)
  set(data1,which(data1[,age_group>7]),"age_group",7)
  set(data1,which(data1[,year_group<2011]),"year_group",2011)
  data1[,order1 := 1:dim(data1)[1]]
  data1 = merge(data1,lambda_index_group,by=c("virus","age_group","year_group"))
  data1 <- data1[order(data1[,order1]),]
  
  sum_lambda=0
  for (j in 1:dim(data1)[1]){
    sum_lambda1 <- exp(log_beta0 + log_beta1[,data1[j,]$time_virus_index] +log_beta2[,data1[j,]$age_index])
    sum_lambda <- sum_lambda + sum_lambda1
  }
  
  xi_this <- xi0[,case6_1[i,]$virus_index]+1
  phi_this <- exp(-log_one_over_phi1[,case6_1[i,]$virus_index])
  stopifnot(data1[1,]$age_group==0)
  lambda0_this <- exp(log_beta0 + log_beta1[,data1[1,]$time_virus_index] + log_beta2[,data1[1,]$age_index])
  lambda_this <- exp(log_beta0 + log_beta1[,case6_1[i,]$time_virus_index] + log_beta2[,case6_1[i,]$age_index])
  case_pre1 <- case6_1[i,]$pop*xi_this/lambda0_this*(exp(lambda0_this/xi_this)-1)*exp(-sum_lambda)*(1-exp(-lambda_this))*phi_this
  est = quantile(case_pre1,probs = c(0.5, 0.025, 0.975))
  case_1_pre1 <- data.table(mean=mean(case_pre1),median=est[1],CI_L=est[2],CI_U=est[3],age_group=age1,year_group=year1,virus=virus1,count=case6_1[i,]$count)
  case_1_pre <- rbind(case_1_pre,case_1_pre1)
}
set(case_1_pre,which(case_1_pre[,age_group>7]),"age_group",7)

exactPoiCI <- function (X, conf.level=0.95) {
  alpha = 1 - conf.level
  upper <- 0.5 * qchisq(1-alpha/2, 2*X+2)
  lower <- 0.5 * qchisq(alpha/2, 2*X)
  return(c(lower, upper))
}
exactPoiCI(30,0.95)
exactPoiCI(0,0.95)

case_1_pre[,count_L := 0]
case_1_pre[,count_U := 0]

for (i in 1:dim(case_1_pre)[1]){
  x=exactPoiCI(case_1_pre[i,count],0.95)
  set(case_1_pre,i,"count_L",x[1])
  set(case_1_pre,i,"count_U",x[2])
}

case_0_pre[,count_L := 0]
case_0_pre[,count_U := 0]

for (i in 1:dim(case_0_pre)[1]){
  x=exactPoiCI(case_0_pre[i,count],0.95)
  set(case_0_pre,i,"count_L",x[1])
  set(case_0_pre,i,"count_U",x[2])
}

case_0_pre1[,count_L := 0]
case_0_pre1[,count_U := 0]

for (i in 1:dim(case_0_pre1)[1]){
  x=exactPoiCI(case_0_pre1[i,count],0.95)
  set(case_0_pre1,i,"count_L",x[1])
  set(case_0_pre1,i,"count_U",x[2])
}


case_1_pre[,virus1 := factor(virus,levels = c("CV-A6","CV-A10","CV-A16","EV-A71"),labels = c("CV-A6","CV-A10","CV-A16","EV-A71"))]
case_0_pre[,virus1 := factor(virus,levels = c("CV-A6","CV-A10","CV-A16","EV-A71"),labels = c("CV-A6","CV-A10","CV-A16","EV-A71"))]
case_0_pre1[,virus1 := factor(virus,levels = c("CV-A6","CV-A10","CV-A16","EV-A71"),labels = c("CV-A6","CV-A10","CV-A16","EV-A71"))]
case_1_pre1 = rbind(case_1_pre,case_0_pre1)


# notes: estimated values: mean, median, CI_L, CI_U; data: count, count_L, count_U
saveRDS(case_1_pre1,"est_case_model_linpw.rds")

##### estimated seroprevalence #####
log_beta0 <- po$log_beta0
log_beta1 <- po$log_beta1
log_beta2 <- po$log_beta2

xi0 <- po$xi0
log_one_over_phi1 <- po$log_one_over_phi

sero_pre_0 <- expand.grid(year_end = 2013:2018,age=0,virus=c("CV-A6","CV-A10","CV-A16","EV-A71"))
sero_pre_0 <- as.data.table(sero_pre_0)
sero_pre_0[,birth_year:=year_end]
colnames(sero_pre_0)[1]="year_group"
colnames(sero_pre_0)[2]="age_group"
sero_pre_0 <- merge(sero_pre_0,lambda_index_group,by=c("year_group","age_group","virus"))
virus_index <- data.table(virus=c("CV-A6","CV-A10","CV-A16","EV-A71"),virus_index=1:4)
sero_pre_0 <- merge(sero_pre_0,virus_index,by=c("virus"))

k=1
virus_index <- sero_pre_0[k,]$virus_index
time_virus_index <- sero_pre_0[k,]$time_virus_index
age_index <- sero_pre_0[k,]$age_index
xi <- 1+xi0[,virus_index]
lambda1 <- exp(log_beta0 + log_beta1[,time_virus_index] + log_beta2[,age_index])

i <- 1
x <- i/12
pr_sus <- xi/lambda1*exp(-x*lambda1)*(exp(min(x,1/xi)*lambda1)-1)
quantile(pr_sus,probs=c(0.5,0.025,0.975))
for (i in 2:11){
  x = i/12
  pr_sus1 <- xi/lambda1*exp(-x*lambda1)*(exp(min(x,1/xi)*lambda1)-1)
  quantile(pr_sus1,probs=c(0.5,0.025,0.975))
  pr_sus <- pr_sus+pr_sus1
}
pr_sero <- 1-pr_sus/12
est <- quantile(pr_sero,probs=c(0.5,0.025,0.975))
sero_pre_0_1 <- data.table(mean=mean(pr_sero),median=est[1],CI_L=est[2],CI_U=est[3],
                           virus=sero_pre_0[k,]$virus,year_group=sero_pre_0[k,]$year_group,age_group=sero_pre_0[k,]$age_group)


for (k in 2:dim(sero_pre_0)[1]){
  virus_index <- sero_pre_0[k,]$virus_index
  time_virus_index <- sero_pre_0[k,]$time_virus_index
  age_index <- sero_pre_0[k,]$age_index
  xi <- 1+xi0[,virus_index]
  lambda1 <- exp(log_beta0 + log_beta1[,time_virus_index] + log_beta2[,age_index])
  i <- 1
  x <- i/12
  pr_sus <- xi/lambda1*exp(-x*lambda1)*(exp(min(x,1/xi)*lambda1)-1)
  quantile(pr_sus,probs=c(0.5,0.025,0.975))
  for (i in 2:11){
    x = i/12
    pr_sus1 <- xi/lambda1*exp(-x*lambda1)*(exp(min(x,1/xi)*lambda1)-1)
    quantile(pr_sus1,probs=c(0.5,0.025,0.975))
    pr_sus <- pr_sus+pr_sus1
  }
  pr_sero <- 1-pr_sus/12
  est <- quantile(pr_sero,probs=c(0.5,0.025,0.975))
  sero_pre_0_2 <- data.table(mean=mean(pr_sero),median=est[1],CI_L=est[2],CI_U=est[3],
                             virus=sero_pre_0[k,]$virus,year_group=sero_pre_0[k,]$year_group,age_group=sero_pre_0[k,]$age_group)
  sero_pre_0_1 <- rbind(sero_pre_0_1,sero_pre_0_2)
}


sero_pre_1 <- expand.grid(year_end = 2013:2018,age=1:6,virus=c("CV-A6","CV-A10","CV-A16","EV-A71"))
sero_pre_1 <- as.data.table(sero_pre_1)
sero_pre_1[,birth_year:=0]
set(sero_pre_1,NULL,"birth_year",sero_pre_1[,year_end-age])
colnames(sero_pre_1)[1]="year_group"
colnames(sero_pre_1)[2]="age_group"

i=1
sero_pre_1[i,]
tmp1 = data.table(year_group=sero_pre_1[i,]$birth_year:sero_pre_1[i,]$year_group,age_group=0:sero_pre_1[i,]$age_group,virus=sero_pre_1[i,]$virus)
set(tmp1,which(tmp1[,year_group<2011]),"year_group",2011)
tmp1[,order1:=1:dim(tmp1)[1]]
tmp1 <- merge(tmp1,lambda_index_group,by=c("age_group","year_group","virus"))
virus_index <- data.table(virus_index=1:4,virus=c("CV-A6","CV-A10","CV-A16","EV-A71"))
tmp1 <- merge(tmp1,virus_index,by="virus")
tmp1 <- tmp1[order(tmp1[,order1]),]

stopifnot(tmp1[1,]$age_group==0)
lambda0_sero <- exp(log_beta0 + log_beta1[,tmp1[1,]$time_virus_index] + log_beta2[,tmp1[1,]$age_index])
xi_this_sero <- 1+xi0[,tmp1[1,]$virus_index]

sum_lambda=0
for (j in 1:dim(tmp1)[1]){
  sum_lambda1 <- exp(log_beta0 + log_beta1[,tmp1[j,]$time_virus_index] + log_beta2[,tmp1[j,]$age_index])
  sum_lambda <- sum_lambda+sum_lambda1
}

sero_pre <- 1-xi_this_sero/lambda0_sero*(exp(lambda0_sero/xi_this_sero)-1)*exp(-sum_lambda)
est <- quantile(sero_pre,probs=c(0.5,0.025,0.975))
sero_pre_1_1 <- data.table(mean=mean(sero_pre),median=est[1],CI_L=est[2],CI_U=est[3],virus=sero_pre_1[i,]$virus,year_group=sero_pre_1[i,]$year_group,age_group=sero_pre_1[i,]$age_group)

for (i in 2:dim(sero_pre_1)[1]){
  # print(i)
  tmp1 = data.table(year_group=sero_pre_1[i,]$birth_year:sero_pre_1[i,]$year_group,age_group=0:sero_pre_1[i,]$age_group,virus=sero_pre_1[i,]$virus)
  set(tmp1,which(tmp1[,year_group<2011]),"year_group",2011)
  tmp1[,order1:=1:dim(tmp1)[1]]
  tmp1 <- merge(tmp1,lambda_index_group,by=c("age_group","year_group","virus"))
  tmp1 <- merge(tmp1,virus_index,by="virus")
  tmp1 <- tmp1[order(tmp1[,order1]),]
  
  stopifnot(tmp1[1,]$age_group==0)
  lambda0_sero <- exp(log_beta0 + log_beta1[,tmp1[1,]$time_virus_index] + log_beta2[,tmp1[1,]$age_index])
  xi_this_sero <- 1+xi0[,tmp1[1,]$virus_index]
  
  sum_lambda=0
  for (j in 1:dim(tmp1)[1]){
    sum_lambda1 <- exp(log_beta0 + log_beta1[,tmp1[j,]$time_virus_index] + log_beta2[,tmp1[j,]$age_index])
    sum_lambda <- sum_lambda+sum_lambda1
  }
  
  sero_pre <- 1-xi_this_sero/lambda0_sero*(exp(lambda0_sero/xi_this_sero)-1)*exp(-sum_lambda)
  est <- quantile(sero_pre,probs=c(0.5,0.025,0.975))
  sero_pre_1_2 <- data.table(mean=mean(sero_pre),median=est[1],CI_L=est[2],CI_U=est[3],virus=sero_pre_1[i,]$virus,year_group=sero_pre_1[i,]$year_group,age_group=sero_pre_1[i,]$age_group)
  sero_pre_1_1 <- rbind(sero_pre_1_1,sero_pre_1_2)
}

sero_pre <- rbind(sero_pre_0_1,sero_pre_1_1)

sero_pre[,virus1 := factor(virus,levels = c("CV-A6","CV-A10","CV-A16","EV-A71"),labels = c("CV-A6","CV-A10","CV-A16","EV-A71"))]
saveRDS(sero_pre,"est_seroprevalence_model_linpw.rds")

#### model LinFt ####
model_fit <- readRDS("model_linft_fit.rdss")
po <- rstan::extract(model_fit)

##### convergence evaluation #####
model_su <- summary(model_fit)$summary
min(model_su[,"n_eff"])
max(model_su[,"Rhat"])
trace_name <- rownames(model_su[order(model_su[,"n_eff"]),][1:4,])
# Trace plot
po_trace <- rstan:::extract(model_fit, inc_warmup = TRUE, permuted = FALSE,
                            pars = trace_name)
bayesplot:::color_scheme_set("mix-blue-pink")
trace <- bayesplot:::mcmc_trace(po_trace, n_warmup = 1e3,
                                facet_args = list(nrow=1,ncol=4,labeller = label_parsed))+
  theme_bw()+
  theme(legend.position='bottom')

##### force of infection #####
log_beta0 <- po$log_beta0
log_beta1 <- po$log_beta1
beta2 <- po$beta2

lambda_index_group_apre <-lambda_index_group[,{
  lambda1 <- exp(log_beta0 + log_beta1[,time_virus_index] - beta2*age_group-log(beta2))*((age_group+1/beta2)*(1-exp(-beta2)) - exp(-beta2))
  #lambda1 <- exp(log_beta0 + log_beta1[,time_virus_index] - beta2*age_group-log(beta2))*age_group
  est = quantile(lambda1,probs = c(0.5, 0.025, 0.975))
  list(mean=mean(lambda1),median=est[1],CI_L=est[2],CI_U=est[3])
},by=c("age_group","year_group","virus","index")]

lambda_index_group_apre[,virus1 := factor(virus,levels = c("CV-A6","CV-A10","CV-A16","EV-A71"),labels = c("CV-A6","CV-A10","CV-A16","EV-A71"))]
saveRDS(lambda_index_group_apre,"est_foi_model_linft.rds")

##### estimated cases #####
xi0 <- po$xi0
log_one_over_phi1 <- po$log_one_over_phi

case1[,age:=floor(`Age in years`)]
case1_0 <- case1[which(case1[,age==0]),c(2,5,6,7,10,13)]
case1_0[,month:=0]
set(case1_0,which(case1_0[,(`Age in years`>=1/12)&(`Age in years`<2/12)]),"month",1)
set(case1_0,which(case1_0[,(`Age in years`>=2/12)&(`Age in years`<3/12)]),"month",2)
set(case1_0,which(case1_0[,(`Age in years`>=3/12)&(`Age in years`<4/12)]),"month",3)
set(case1_0,which(case1_0[,(`Age in years`>=4/12)&(`Age in years`<5/12)]),"month",4)
set(case1_0,which(case1_0[,(`Age in years`>=5/12)&(`Age in years`<6/12)]),"month",5)
set(case1_0,which(case1_0[,(`Age in years`>=6/12)&(`Age in years`<7/12)]),"month",6)
set(case1_0,which(case1_0[,(`Age in years`>=7/12)&(`Age in years`<8/12)]),"month",7)
set(case1_0,which(case1_0[,(`Age in years`>=8/12)&(`Age in years`<9/12)]),"month",8)
set(case1_0,which(case1_0[,(`Age in years`>=9/12)&(`Age in years`<10/12)]),"month",9)
set(case1_0,which(case1_0[,(`Age in years`>=10/12)&(`Age in years`<11/12)]),"month",10)
set(case1_0,which(case1_0[,(`Age in years`>=11/12)&(`Age in years`<1)]),"month",11)

case2_0 <- case1_0[,list(count=.N),by=c("virus","year_sampling","age","month")]
colnames(case2_0)[2]="year_group"
colnames(case2_0)[3]="age_group"

tmp <- expand.grid(virus=c("CV-A6","CV-A10","CV-A16","EV-A71"),year_group=2013:2018,age_group=0,month=1:11)
tmp$self_index <- 1:dim(tmp)[1]
tmp <- as.data.table(tmp)

case2_0 <- merge(case2_0,tmp,by=c("virus","year_group","age_group","month"))

zero_index1 <- tmp[which(tmp[,!self_index%in%case2_0$self_index]),]
zero_index1[,count:=0]

case3_0 <- rbind(case2_0,zero_index1)
case3_0 <- case3_0[order(case3_0[,self_index]),]

case4_0 <- merge(case3_0,lambda_index_group,by=c("virus","year_group","age_group"))
case4_0 <- case4_0[order(case4_0[,self_index]),]

casepo
colnames(casepo)[1]="age_group"
colnames(casepo)[2]="year_group"

case5_0 <- merge(case4_0,casepo[,1:3],by=c("age_group","year_group"))

virus_index <- data.table(virus=c("CV-A6","CV-A10","CV-A16","EV-A71"),virus_index=1:4)
case6_0 <- merge(case5_0,virus_index,by="virus")
case6_0 <- case6_0[order(case6_0[,self_index]),]

case_0_pre <- case6_0[,{
  xi_this = xi0[,virus_index]+1
  phi_this = exp(-log_one_over_phi1[,virus_index])
  lambda1 = exp(log_beta0 + log_beta1[,time_virus_index]-log(beta2))*(1/beta2*(1-exp(-beta2)) - exp(-beta2))
  case_pre1 = pop*xi_this/lambda1*exp(-month*lambda1/12)*(exp(min(month/12,1/xi_this)*lambda1)-1)*(1-exp(-lambda1/12))*phi_this
  est = quantile(case_pre1,probs = c(0.5, 0.025, 0.975))
  list(mean=mean(case_pre1),median=est[1],CI_L=est[2],CI_U=est[3])
},by=c("age_group","year_group","virus","month","count","self_index")]


case7_0 <- case6_0[,list(count1=sum(count)),by=c("age_group","year_group","virus","index","time_virus_index","age_index","virus_index","pop")]
i=1
xi_this = xi0[,case7_0[i,]$virus_index]+1
phi_this = exp(-log_one_over_phi1[,case7_0[i,]$virus_index])
lambda1 = exp(log_beta0 + log_beta1[,case7_0[i,]$time_virus_index]-log(beta2))*(1/beta2*(1-exp(-beta2)) - exp(-beta2))
case_pre1 = 0
for (x in 1:11){
  case_pre2 = case7_0[i,]$pop*xi_this/lambda1*exp(-x*lambda1/12)*(exp(min(x/12,1/xi_this)*lambda1)-1)*(1-exp(-lambda1/12))*phi_this
  case_pre1 = case_pre1+case_pre2
}
est = quantile(case_pre1,probs = c(0.5, 0.025, 0.975))
case_0_pre1 = data.table(mean=mean(case_pre1),median=est[1],CI_L=est[2],CI_U=est[3],age_group=0,year_group=case7_0[i,]$year_group,virus=case7_0[i,]$virus,count=case7_0[i,]$count1)
for (i in 2:dim(case7_0)[1]){
  xi_this = xi0[,case7_0[i,]$virus_index]+1
  phi_this = exp(-log_one_over_phi1[,case7_0[i,]$virus_index])
  lambda1 = exp(log_beta0 + log_beta1[,case7_0[i,]$time_virus_index]-log(beta2))*(1/beta2*(1-exp(-beta2)) - exp(-beta2))
  case_pre1 = 0
  for (x in 1:11){
    case_pre2 = case7_0[i,]$pop*xi_this/lambda1*exp(-x*lambda1/12)*(exp(min(x/12,1/xi_this)*lambda1)-1)*(1-exp(-lambda1/12))*phi_this
    case_pre1 = case_pre1+case_pre2
  }
  est = quantile(case_pre1,probs = c(0.5, 0.025, 0.975))
  case_0_pre2 = data.table(mean=mean(case_pre1),median=est[1],CI_L=est[2],CI_U=est[3],age_group=0,year_group=case7_0[i,]$year_group,virus=case7_0[i,]$virus,count=case7_0[i,]$count1)
  case_0_pre1 = rbind(case_0_pre1,case_0_pre2)
}



case1_1 <- case1[which(case1[,age>0]),c(2,5,6,7,10,13)]
case1_1[,age_group:=age]
set(case1_1,which(case1_1[,age_group>=7]),"age_group",7)
case2_1 <- case1_1[,list(count=.N),by=c("virus","year_sampling","age_group")]
colnames(case2_1)[2]="year_group"



tmp2 <- expand.grid(virus=c("CV-A6","CV-A10","CV-A16","EV-A71"),year_group=2013:2018,age_group=1:7)
tmp2$self_index <- 1:dim(tmp2)[1]
tmp2 <- as.data.table(tmp2)

case2_1 <- merge(case2_1, tmp2, by=c("virus","year_group","age_group"))
zero_index2 <- tmp2[which(tmp2[,!self_index%in%case2_1$self_index]),]
zero_index2[,count:=0]

case3_1 <- rbind(case2_1,zero_index2)
case3_1 <- case3_1[order(case3_1[,self_index]),]

case4_1 <- merge(case3_1,lambda_index_group,by=c("virus","year_group","age_group"))
case4_1 <- case4_1[order(case4_1[,self_index]),]

virus_index
case5_1 <- merge(case4_1,virus_index,by="virus")
case5_1 <- case5_1[order(case5_1[,self_index]),]

pop_group <- expand.grid(year_group=2013:2018,age_group=1:7)
pop_group$pop <- NA
for (i in 1:dim(pop_group)[1]){
  age <- pop_group[i,]$age_group
  year1 <- pop_group[i,]$year_group
  if (age<7){
    tmp = data.table(age_group=age,year_group=year1)
    tmp1 = merge(tmp,casepo,by=c("age_group","year_group"))
    pop_group[i,]$pop = tmp1$pop
  }else{
    tmp = expand.grid(age_group=7:11,year_group=year1)
    tmp1 = merge(tmp,casepo,by=c("age_group","year_group"))
    pop_group[i,]$pop = sum(tmp1$pop)
  }
}

case6_1 <- merge(case5_1,pop_group,by=c("year_group","age_group"))
case6_1 <- case6_1[order(case6_1[,self_index]),]


i=1
case6_1[i,]
virus1 <- case6_1[i,]$virus
age1 <- case6_1[i,]$age_group
year1 <- case6_1[i,]$year_group

if (age1==7) age1=9

data1 <- data.table(year_group = (year1-age1):(year1-1),age_group = 0:(age1-1),virus=virus1)
set(data1,which(data1[,age_group>7]),"age_group",7)
set(data1,which(data1[,year_group<2011]),"year_group",2011)
data1[,order1 := 1:dim(data1)[1]]
data1 = merge(data1,lambda_index_group,by=c("virus","age_group","year_group"))
data1 <- data1[order(data1[,order1]),]

sum_lambda=0
for (j in 1:dim(data1)[1]){
  sum_lambda1 <- exp(log_beta0 + log_beta1[,data1[j,]$time_virus_index] - beta2*data1[j,]$age_group-log(beta2))*((data1[j,]$age_group+1/beta2)*(1-exp(-beta2)) - exp(-beta2))
  sum_lambda <- sum_lambda + sum_lambda1
}

xi_this <- xi0[,case6_1[i,]$virus_index]+1
phi_this <- exp(-log_one_over_phi1[,case6_1[i,]$virus_index])
stopifnot(data1[1,]$age_group==0)
lambda0_this <- exp(log_beta0 + log_beta1[,data1[1,]$time_virus_index]-log(beta2))*(1/beta2*(1-exp(-beta2)) - exp(-beta2))
lambda_this <- exp(log_beta0 + log_beta1[,case6_1[i,]$time_virus_index] - beta2*age1-log(beta2))*((age1+1/beta2)*(1-exp(-beta2)) - exp(-beta2))
case_pre1 <- case6_1[i,]$pop*xi_this/lambda0_this*(exp(lambda0_this/xi_this)-1)*exp(-sum_lambda)*(1-exp(-lambda_this))*phi_this
est = quantile(case_pre1,probs = c(0.5, 0.025, 0.975))
case_1_pre <- data.table(mean=mean(case_pre1),median=est[1],CI_L=est[2],CI_U=est[3],age_group=age1,year_group=year1,virus=virus1,count=case6_1[i,]$count)


for (i in 2:dim(case6_1)[1]){
  virus1 <- case6_1[i,]$virus
  age1 <- case6_1[i,]$age_group
  year1 <- case6_1[i,]$year_group
  
  if (age1==7) age1=9
  
  data1 <- data.table(year_group = (year1-age1):(year1-1),age_group = 0:(age1-1),virus=virus1)
  set(data1,which(data1[,age_group>7]),"age_group",7)
  set(data1,which(data1[,year_group<2011]),"year_group",2011)
  data1[,order1 := 1:dim(data1)[1]]
  data1 = merge(data1,lambda_index_group,by=c("virus","age_group","year_group"))
  data1 <- data1[order(data1[,order1]),]
  
  sum_lambda=0
  for (j in 1:dim(data1)[1]){
    sum_lambda1 <- exp(log_beta0 + log_beta1[,data1[j,]$time_virus_index] - beta2*data1[j,]$age_group-log(beta2))*((data1[j,]$age_group+1/beta2)*(1-exp(-beta2)) - exp(-beta2))
    sum_lambda <- sum_lambda + sum_lambda1
  }
  
  xi_this <- xi0[,case6_1[i,]$virus_index]+1
  phi_this <- exp(-log_one_over_phi1[,case6_1[i,]$virus_index])
  stopifnot(data1[1,]$age_group==0)
  lambda0_this <- exp(log_beta0 + log_beta1[,data1[1,]$time_virus_index]-log(beta2))*(1/beta2*(1-exp(-beta2)) - exp(-beta2))
  lambda_this <- exp(log_beta0 + log_beta1[,case6_1[i,]$time_virus_index] - beta2*age1-log(beta2))*((age1+1/beta2)*(1-exp(-beta2)) - exp(-beta2))
  case_pre1 <- case6_1[i,]$pop*xi_this/lambda0_this*(exp(lambda0_this/xi_this)-1)*exp(-sum_lambda)*(1-exp(-lambda_this))*phi_this
  est = quantile(case_pre1,probs = c(0.5, 0.025, 0.975))
  case_1_pre1 <- data.table(mean=mean(case_pre1),median=est[1],CI_L=est[2],CI_U=est[3],age_group=age1,year_group=year1,virus=virus1,count=case6_1[i,]$count)
  case_1_pre <- rbind(case_1_pre,case_1_pre1)
}
set(case_1_pre,which(case_1_pre[,age_group>7]),"age_group",7)

exactPoiCI <- function (X, conf.level=0.95) {
  alpha = 1 - conf.level
  upper <- 0.5 * qchisq(1-alpha/2, 2*X+2)
  lower <- 0.5 * qchisq(alpha/2, 2*X)
  return(c(lower, upper))
}
exactPoiCI(30,0.95)
exactPoiCI(0,0.95)

case_1_pre[,count_L := 0]
case_1_pre[,count_U := 0]

for (i in 1:dim(case_1_pre)[1]){
  x=exactPoiCI(case_1_pre[i,count],0.95)
  set(case_1_pre,i,"count_L",x[1])
  set(case_1_pre,i,"count_U",x[2])
}

case_0_pre[,count_L := 0]
case_0_pre[,count_U := 0]

for (i in 1:dim(case_0_pre)[1]){
  x=exactPoiCI(case_0_pre[i,count],0.95)
  set(case_0_pre,i,"count_L",x[1])
  set(case_0_pre,i,"count_U",x[2])
}

case_0_pre1[,count_L := 0]
case_0_pre1[,count_U := 0]

for (i in 1:dim(case_0_pre1)[1]){
  x=exactPoiCI(case_0_pre1[i,count],0.95)
  set(case_0_pre1,i,"count_L",x[1])
  set(case_0_pre1,i,"count_U",x[2])
}

case_1_pre[,virus1 := factor(virus,levels = c("CV-A6","CV-A10","CV-A16","EV-A71"),labels = c("CV-A6","CV-A10","CV-A16","EV-A71"))]
case_0_pre[,virus1 := factor(virus,levels = c("CV-A6","CV-A10","CV-A16","EV-A71"),labels = c("CV-A6","CV-A10","CV-A16","EV-A71"))]
case_0_pre1[,virus1 := factor(virus,levels = c("CV-A6","CV-A10","CV-A16","EV-A71"),labels = c("CV-A6","CV-A10","CV-A16","EV-A71"))]
case_1_pre1 = rbind(case_1_pre,case_0_pre1)

saveRDS(case_1_pre1,"est_case_model_linft.rds")

##### estimated seroprevalence #####
log_beta0 <- po$log_beta0
log_beta1 <- po$log_beta1
beta2 <- po$beta2

xi0 <- po$xi0
log_one_over_phi1 <- po$log_one_over_phi

sero_pre_0 <- expand.grid(year_end = 2013:2018,age=0,virus=c("CV-A6","CV-A10","CV-A16","EV-A71"))
sero_pre_0 <- as.data.table(sero_pre_0)
sero_pre_0[,birth_year:=year_end]
colnames(sero_pre_0)[1]="year_group"
colnames(sero_pre_0)[2]="age_group"
sero_pre_0 <- merge(sero_pre_0,lambda_index_group,by=c("year_group","age_group","virus"))
virus_index <- data.table(virus=c("CV-A6","CV-A10","CV-A16","EV-A71"),virus_index=1:4)
sero_pre_0 <- merge(sero_pre_0,virus_index,by=c("virus"))

k=1
virus_index <- sero_pre_0[k,]$virus_index
time_virus_index <- sero_pre_0[k,]$time_virus_index
age_group <- sero_pre_0[k,]$age_group
xi <- 1+xi0[,virus_index]
lambda1 <- exp(log_beta0 + log_beta1[,time_virus_index] - beta2*age_group-log(beta2))*((age_group+1/beta2)*(1-exp(-beta2)) - exp(-beta2))

i <- 1
x <- i/12
pr_sus <- xi/lambda1*exp(-x*lambda1)*(exp(min(x,1/xi)*lambda1)-1)
quantile(pr_sus,probs=c(0.5,0.025,0.975))
for (i in 2:11){
  x = i/12
  pr_sus1 <- xi/lambda1*exp(-x*lambda1)*(exp(min(x,1/xi)*lambda1)-1)
  quantile(pr_sus1,probs=c(0.5,0.025,0.975))
  pr_sus <- pr_sus+pr_sus1
}
pr_sero <- 1-pr_sus/12
est <- quantile(pr_sero,probs=c(0.5,0.025,0.975))
sero_pre_0_1 <- data.table(mean=mean(pr_sero),median=est[1],CI_L=est[2],CI_U=est[3],
                           virus=sero_pre_0[k,]$virus,year_group=sero_pre_0[k,]$year_group,age_group=sero_pre_0[k,]$age_group)


for (k in 2:dim(sero_pre_0)[1]){
  virus_index <- sero_pre_0[k,]$virus_index
  time_virus_index <- sero_pre_0[k,]$time_virus_index
  age_group <- sero_pre_0[k,]$age_group
  xi <- 1+xi0[,virus_index]
  #lambda1 <- exp(log_beta0 + log_beta1[,time_virus_index] + log_beta2[,age_index])
  lambda1 <- exp(log_beta0 + log_beta1[,time_virus_index] - beta2*age_group-log(beta2))*((age_group+1/beta2)*(1-exp(-beta2)) - exp(-beta2))
  
  i <- 1
  x <- i/12
  pr_sus <- xi/lambda1*exp(-x*lambda1)*(exp(min(x,1/xi)*lambda1)-1)
  quantile(pr_sus,probs=c(0.5,0.025,0.975))
  for (i in 2:11){
    x = i/12
    pr_sus1 <- xi/lambda1*exp(-x*lambda1)*(exp(min(x,1/xi)*lambda1)-1)
    quantile(pr_sus1,probs=c(0.5,0.025,0.975))
    pr_sus <- pr_sus+pr_sus1
  }
  pr_sero <- 1-pr_sus/12
  est <- quantile(pr_sero,probs=c(0.5,0.025,0.975))
  sero_pre_0_2 <- data.table(mean=mean(pr_sero),median=est[1],CI_L=est[2],CI_U=est[3],
                             virus=sero_pre_0[k,]$virus,year_group=sero_pre_0[k,]$year_group,age_group=sero_pre_0[k,]$age_group)
  sero_pre_0_1 <- rbind(sero_pre_0_1,sero_pre_0_2)
}


sero_pre_1 <- expand.grid(year_end = 2013:2018,age=1:6,virus=c("CV-A6","CV-A10","CV-A16","EV-A71"))
sero_pre_1 <- as.data.table(sero_pre_1)
sero_pre_1[,birth_year:=0]
set(sero_pre_1,NULL,"birth_year",sero_pre_1[,year_end-age])
colnames(sero_pre_1)[1]="year_group"
colnames(sero_pre_1)[2]="age_group"

i=1
sero_pre_1[i,]
tmp1 = data.table(year_group=sero_pre_1[i,]$birth_year:sero_pre_1[i,]$year_group,age_group=0:sero_pre_1[i,]$age_group,virus=sero_pre_1[i,]$virus)
set(tmp1,which(tmp1[,year_group<2011]),"year_group",2011)
tmp1[,order1:=1:dim(tmp1)[1]]
tmp1 <- merge(tmp1,lambda_index_group,by=c("age_group","year_group","virus"))
virus_index <- data.table(virus_index=1:4,virus=c("CV-A6","CV-A10","CV-A16","EV-A71"))
tmp1 <- merge(tmp1,virus_index,by="virus")
tmp1 <- tmp1[order(tmp1[,order1]),]

stopifnot(tmp1[1,]$age_group==0)
#lambda0_sero <- exp(log_beta0 + log_beta1[,tmp1[1,]$time_virus_index] + log_beta2[,tmp1[1,]$age_index])
lambda0_sero <- exp(log_beta0 + log_beta1[,tmp1[1,]$time_virus_index] - beta2*tmp1[1,]$age_group-log(beta2))*((tmp1[1,]$age_group+1/beta2)*(1-exp(-beta2)) - exp(-beta2))

xi_this_sero <- 1+xi0[,tmp1[1,]$virus_index]

sum_lambda=0
for (j in 1:dim(tmp1)[1]){
  sum_lambda1 <- exp(log_beta0 + log_beta1[,tmp1[j,]$time_virus_index] - beta2*tmp1[j,]$age_group-log(beta2))*((tmp1[j,]$age_group+1/beta2)*(1-exp(-beta2)) - exp(-beta2))
  sum_lambda <- sum_lambda+sum_lambda1
}

sero_pre <- 1-xi_this_sero/lambda0_sero*(exp(lambda0_sero/xi_this_sero)-1)*exp(-sum_lambda)
est <- quantile(sero_pre,probs=c(0.5,0.025,0.975))
sero_pre_1_1 <- data.table(mean=mean(sero_pre),median=est[1],CI_L=est[2],CI_U=est[3],virus=sero_pre_1[i,]$virus,year_group=sero_pre_1[i,]$year_group,age_group=sero_pre_1[i,]$age_group)

for (i in 2:dim(sero_pre_1)[1]){
  # print(i)
  tmp1 = data.table(year_group=sero_pre_1[i,]$birth_year:sero_pre_1[i,]$year_group,age_group=0:sero_pre_1[i,]$age_group,virus=sero_pre_1[i,]$virus)
  set(tmp1,which(tmp1[,year_group<2011]),"year_group",2011)
  tmp1[,order1:=1:dim(tmp1)[1]]
  tmp1 <- merge(tmp1,lambda_index_group,by=c("age_group","year_group","virus"))
  tmp1 <- merge(tmp1,virus_index,by="virus")
  tmp1 <- tmp1[order(tmp1[,order1]),]
  
  stopifnot(tmp1[1,]$age_group==0)
  lambda0_sero <- exp(log_beta0 + log_beta1[,tmp1[1,]$time_virus_index] - beta2*tmp1[1,]$age_group-log(beta2))*((tmp1[1,]$age_group+1/beta2)*(1-exp(-beta2)) - exp(-beta2))
  xi_this_sero <- 1+xi0[,tmp1[1,]$virus_index]
  
  sum_lambda=0
  for (j in 1:dim(tmp1)[1]){
    sum_lambda1 <- exp(log_beta0 + log_beta1[,tmp1[j,]$time_virus_index] - beta2*tmp1[j,]$age_group-log(beta2))*((tmp1[j,]$age_group+1/beta2)*(1-exp(-beta2)) - exp(-beta2))
    sum_lambda <- sum_lambda+sum_lambda1
  }
  
  sero_pre <- 1-xi_this_sero/lambda0_sero*(exp(lambda0_sero/xi_this_sero)-1)*exp(-sum_lambda)
  est <- quantile(sero_pre,probs=c(0.5,0.025,0.975))
  sero_pre_1_2 <- data.table(mean=mean(sero_pre),median=est[1],CI_L=est[2],CI_U=est[3],virus=sero_pre_1[i,]$virus,year_group=sero_pre_1[i,]$year_group,age_group=sero_pre_1[i,]$age_group)
  sero_pre_1_1 <- rbind(sero_pre_1_1,sero_pre_1_2)
}

sero_pre <- rbind(sero_pre_0_1,sero_pre_1_1)

sero_pre[,virus1 := factor(virus,levels = c("CV-A6","CV-A10","CV-A16","EV-A71"),labels = c("CV-A6","CV-A10","CV-A16","EV-A71"))]
saveRDS(sero_pre,"est_seroprevalence_model_linft.rds")


#### model ExpPw ####
model_fit <- readRDS("model_exppw_fit.rds")
po <- rstan::extract(model_fit)

##### convergence evaluation #####
model_su <- summary(model_fit)$summary
min(model_su[,"n_eff"])
max(model_su[,"Rhat"])
trace_name <- rownames(model_su[order(model_su[,"n_eff"]),][1:4,])
# Trace plot
po_trace <- rstan:::extract(model_fit, inc_warmup = TRUE, permuted = FALSE,
                            pars = trace_name)
bayesplot:::color_scheme_set("mix-blue-pink")
trace <- bayesplot:::mcmc_trace(po_trace, n_warmup = 1e3,
                                facet_args = list(nrow=1,ncol=4,labeller = label_parsed))+
  theme_bw()+
  theme(legend.position='bottom')

##### force of infection #####
log_beta0 <- po$log_beta0
log_beta1 <- po$log_beta1
log_beta2 <- po$log_beta2

lambda_index_group_apre <-lambda_index_group[,{
  lambda1 <- exp(log_beta0 + log_beta1[,time_virus_index] + log_beta2[,age_index] )
  est = quantile(lambda1,probs = c(0.5, 0.025, 0.975))
  list(mean=mean(lambda1),median=est[1],CI_L=est[2],CI_U=est[3])
},by=c("age_group","year_group","virus","index")]

lambda_index_group_apre[,virus1 := factor(virus,levels = c("CV-A6","CV-A10","CV-A16","EV-A71"),labels = c("CV-A6","CV-A10","CV-A16","EV-A71"))]
saveRDS(lambda_index_group_apre,"est_foi_model_exppw.rds")

##### estimated cases #####
xi <- po$xi
log_one_over_phi1 <- po$log_one_over_phi

case1[,age:=floor(`Age in years`)]
case1_0 <- case1[which(case1[,age==0]),c(2,5,6,7,10,13)]
case1_0[,month:=0]
set(case1_0,which(case1_0[,(`Age in years`>=1/12)&(`Age in years`<2/12)]),"month",1)
set(case1_0,which(case1_0[,(`Age in years`>=2/12)&(`Age in years`<3/12)]),"month",2)
set(case1_0,which(case1_0[,(`Age in years`>=3/12)&(`Age in years`<4/12)]),"month",3)
set(case1_0,which(case1_0[,(`Age in years`>=4/12)&(`Age in years`<5/12)]),"month",4)
set(case1_0,which(case1_0[,(`Age in years`>=5/12)&(`Age in years`<6/12)]),"month",5)
set(case1_0,which(case1_0[,(`Age in years`>=6/12)&(`Age in years`<7/12)]),"month",6)
set(case1_0,which(case1_0[,(`Age in years`>=7/12)&(`Age in years`<8/12)]),"month",7)
set(case1_0,which(case1_0[,(`Age in years`>=8/12)&(`Age in years`<9/12)]),"month",8)
set(case1_0,which(case1_0[,(`Age in years`>=9/12)&(`Age in years`<10/12)]),"month",9)
set(case1_0,which(case1_0[,(`Age in years`>=10/12)&(`Age in years`<11/12)]),"month",10)
set(case1_0,which(case1_0[,(`Age in years`>=11/12)&(`Age in years`<1)]),"month",11)

summary(case1_0)

case2_0 <- case1_0[,list(count=.N),by=c("virus","year_sampling","age","month")]
colnames(case2_0)[2]="year_group"
colnames(case2_0)[3]="age_group"

tmp <- expand.grid(virus=c("CV-A6","CV-A10","CV-A16","EV-A71"),year_group=2013:2018,age_group=0,month=1:11)
tmp$self_index <- 1:dim(tmp)[1]
tmp <- as.data.table(tmp)

case2_0 <- merge(case2_0,tmp,by=c("virus","year_group","age_group","month"))

zero_index1 <- tmp[which(tmp[,!self_index%in%case2_0$self_index]),]
zero_index1[,count:=0]

case3_0 <- rbind(case2_0,zero_index1)
case3_0 <- case3_0[order(case3_0[,self_index]),]

case4_0 <- merge(case3_0,lambda_index_group,by=c("virus","year_group","age_group"))
summary(case4_0)
case4_0 <- case4_0[order(case4_0[,self_index]),]

casepo
colnames(casepo)[1]="age_group"
colnames(casepo)[2]="year_group"

case5_0 <- merge(case4_0,casepo[,1:3],by=c("age_group","year_group"))

virus_index <- data.table(virus=c("CV-A6","CV-A10","CV-A16","EV-A71"),virus_index=1:4)
case6_0 <- merge(case5_0,virus_index,by="virus")
case6_0 <- case6_0[order(case6_0[,self_index]),]

case_0_pre <- case6_0[,{
  xi_this = xi[,virus_index]
  phi_this = exp(-log_one_over_phi1[,virus_index])
  lambda1 = exp(log_beta0 + log_beta1[,time_virus_index] + log_beta2[,age_index])
  case_pre1 = pop*xi_this/(xi_this-lambda1)*(exp(-month*lambda1/12)-exp(-xi_this*month/12))*(1-exp(-lambda1/12))*phi_this
  est = quantile(case_pre1,probs = c(0.5, 0.025, 0.975))
  list(mean=mean(case_pre1),median=est[1],CI_L=est[2],CI_U=est[3])
},by=c("age_group","year_group","virus","month","count","self_index")]


case7_0 <- case6_0[,list(count1=sum(count)),by=c("age_group","year_group","virus","index","time_virus_index","age_index","virus_index","pop")]
i=1
xi_this = xi[,case7_0[i,]$virus_index]
phi_this = exp(-log_one_over_phi1[,case7_0[i,]$virus_index])
lambda1 = exp(log_beta0 + log_beta1[,case7_0[i,]$time_virus_index] + log_beta2[,case7_0[i,]$age_index])
case_pre1 = 0
for (x in 1:11){
  case_pre2 = case7_0[i,]$pop*xi_this/(xi_this-lambda1)*(exp(-x*lambda1/12)-exp(-xi_this*x/12))*(1-exp(-lambda1/12))*phi_this
  case_pre1 = case_pre1+case_pre2
}
est = quantile(case_pre1,probs = c(0.5, 0.025, 0.975))
case_0_pre1 = data.table(mean=mean(case_pre1),median=est[1],CI_L=est[2],CI_U=est[3],age_group=0,year_group=case7_0[i,]$year_group,virus=case7_0[i,]$virus,count=case7_0[i,]$count1)
for (i in 2:dim(case7_0)[1]){
  xi_this = xi[,case7_0[i,]$virus_index]
  phi_this = exp(-log_one_over_phi1[,case7_0[i,]$virus_index])
  lambda1 = exp(log_beta0 + log_beta1[,case7_0[i,]$time_virus_index] + log_beta2[,case7_0[i,]$age_index])
  case_pre1 = 0
  for (x in 1:11){
    case_pre2 = case7_0[i,]$pop*xi_this/(xi_this-lambda1)*(exp(-x*lambda1/12)-exp(-xi_this*x/12))*(1-exp(-lambda1/12))*phi_this
    case_pre1 = case_pre1+case_pre2
  }
  est = quantile(case_pre1,probs = c(0.5, 0.025, 0.975))
  case_0_pre2 = data.table(mean=mean(case_pre1),median=est[1],CI_L=est[2],CI_U=est[3],age_group=0,year_group=case7_0[i,]$year_group,virus=case7_0[i,]$virus,count=case7_0[i,]$count1)
  case_0_pre1 = rbind(case_0_pre1,case_0_pre2)
}



case1_1 <- case1[which(case1[,age>0]),c(2,5,6,7,10,13)]
case1_1[,age_group:=age]
set(case1_1,which(case1_1[,age_group>=7]),"age_group",7)
case2_1 <- case1_1[,list(count=.N),by=c("virus","year_sampling","age_group")]
colnames(case2_1)[2]="year_group"



tmp2 <- expand.grid(virus=c("CV-A6","CV-A10","CV-A16","EV-A71"),year_group=2013:2018,age_group=1:7)
tmp2$self_index <- 1:dim(tmp2)[1]
tmp2 <- as.data.table(tmp2)

case2_1 <- merge(case2_1, tmp2, by=c("virus","year_group","age_group"))
zero_index2 <- tmp2[which(tmp2[,!self_index%in%case2_1$self_index]),]
zero_index2[,count:=0]

case3_1 <- rbind(case2_1,zero_index2)
case3_1 <- case3_1[order(case3_1[,self_index]),]

case4_1 <- merge(case3_1,lambda_index_group,by=c("virus","year_group","age_group"))
case4_1 <- case4_1[order(case4_1[,self_index]),]

virus_index
case5_1 <- merge(case4_1,virus_index,by="virus")
case5_1 <- case5_1[order(case5_1[,self_index]),]

pop_group <- expand.grid(year_group=2013:2018,age_group=1:7)
pop_group$pop <- NA
for (i in 1:dim(pop_group)[1]){
  age <- pop_group[i,]$age_group
  year1 <- pop_group[i,]$year_group
  if (age<7){
    tmp = data.table(age_group=age,year_group=year1)
    tmp1 = merge(tmp,casepo,by=c("age_group","year_group"))
    pop_group[i,]$pop = tmp1$pop
  }else{
    tmp = expand.grid(age_group=7:11,year_group=year1)
    tmp1 = merge(tmp,casepo,by=c("age_group","year_group"))
    pop_group[i,]$pop = sum(tmp1$pop)
  }
}

case6_1 <- merge(case5_1,pop_group,by=c("year_group","age_group"))
case6_1 <- case6_1[order(case6_1[,self_index]),]


i=1
case6_1[i,]
virus1 <- case6_1[i,]$virus
age1 <- case6_1[i,]$age_group
year1 <- case6_1[i,]$year_group

if (age1==7) age1=9

data1 <- data.table(year_group = (year1-age1):(year1-1),age_group = 0:(age1-1),virus=virus1)
set(data1,which(data1[,age_group>7]),"age_group",7)
set(data1,which(data1[,year_group<2011]),"year_group",2011)
data1[,order1 := 1:dim(data1)[1]]
data1 = merge(data1,lambda_index_group,by=c("virus","age_group","year_group"))
data1 <- data1[order(data1[,order1]),]

sum_lambda=0
for (j in 1:dim(data1)[1]){
  sum_lambda1 <- exp(log_beta0 + log_beta1[,data1[j,]$time_virus_index] + log_beta2[,data1[j,]$age_index])
  sum_lambda <- sum_lambda + sum_lambda1
}

quantile(sum_lambda,probs = c(0.5, 0.025, 0.975))
quantile(lambda_this,probs = c(0.5, 0.025, 0.975))

xi_this <- xi[,case6_1[i,]$virus_index]
phi_this <- exp(-log_one_over_phi1[,case6_1[i,]$virus_index])
stopifnot(data1[1,]$age_group==0)
lambda0_this <- exp(log_beta0 + log_beta1[,data1[1,]$time_virus_index] + log_beta2[,data1[1,]$age_index])
lambda_this <- exp(log_beta0 + log_beta1[,case6_1[i,]$time_virus_index] + log_beta2[,case6_1[i,]$age_index])
case_pre1 <- case6_1[i,]$pop*xi_this/(xi_this-lambda0_this)*exp(-sum_lambda)*(1-exp(-lambda_this))*phi_this
est = quantile(case_pre1,probs = c(0.5, 0.025, 0.975))
case_1_pre <- data.table(mean=mean(case_pre1),median=est[1],CI_L=est[2],CI_U=est[3],age_group=age1,year_group=year1,virus=virus1,count=case6_1[i,]$count)


for (i in 2:dim(case6_1)[1]){
  virus1 <- case6_1[i,]$virus
  age1 <- case6_1[i,]$age_group
  year1 <- case6_1[i,]$year_group
  
  if (age1==7) age1=9
  
  data1 <- data.table(year_group = (year1-age1):(year1-1),age_group = 0:(age1-1),virus=virus1)
  set(data1,which(data1[,age_group>7]),"age_group",7)
  set(data1,which(data1[,year_group<2011]),"year_group",2011)
  data1[,order1 := 1:dim(data1)[1]]
  data1 = merge(data1,lambda_index_group,by=c("virus","age_group","year_group"))
  data1 <- data1[order(data1[,order1]),]
  
  sum_lambda=0
  for (j in 1:dim(data1)[1]){
    sum_lambda1 <- exp(log_beta0 + log_beta1[,data1[j,]$time_virus_index] +log_beta2[,data1[j,]$age_index])
    sum_lambda <- sum_lambda + sum_lambda1
  }
  
  xi_this <- xi[,case6_1[i,]$virus_index]
  phi_this <- exp(-log_one_over_phi1[,case6_1[i,]$virus_index])
  stopifnot(data1[1,]$age_group==0)
  lambda0_this <- exp(log_beta0 + log_beta1[,data1[1,]$time_virus_index] + log_beta2[,data1[1,]$age_index])
  lambda_this <- exp(log_beta0 + log_beta1[,case6_1[i,]$time_virus_index] + log_beta2[,case6_1[i,]$age_index])
  case_pre1 <- case6_1[i,]$pop*xi_this/(xi_this-lambda0_this)*exp(-sum_lambda)*(1-exp(-lambda_this))*phi_this
  est = quantile(case_pre1,probs = c(0.5, 0.025, 0.975))
  case_1_pre1 <- data.table(mean=mean(case_pre1),median=est[1],CI_L=est[2],CI_U=est[3],age_group=age1,year_group=year1,virus=virus1,count=case6_1[i,]$count)
  case_1_pre <- rbind(case_1_pre,case_1_pre1)
}
set(case_1_pre,which(case_1_pre[,age_group>7]),"age_group",7)

exactPoiCI <- function (X, conf.level=0.95) {
  alpha = 1 - conf.level
  upper <- 0.5 * qchisq(1-alpha/2, 2*X+2)
  lower <- 0.5 * qchisq(alpha/2, 2*X)
  return(c(lower, upper))
}
exactPoiCI(30,0.95)
exactPoiCI(0,0.95)

case_1_pre[,count_L := 0]
case_1_pre[,count_U := 0]

for (i in 1:dim(case_1_pre)[1]){
  x=exactPoiCI(case_1_pre[i,count],0.95)
  set(case_1_pre,i,"count_L",x[1])
  set(case_1_pre,i,"count_U",x[2])
}

case_0_pre[,count_L := 0]
case_0_pre[,count_U := 0]

for (i in 1:dim(case_0_pre)[1]){
  x=exactPoiCI(case_0_pre[i,count],0.95)
  set(case_0_pre,i,"count_L",x[1])
  set(case_0_pre,i,"count_U",x[2])
}

case_0_pre1[,count_L := 0]
case_0_pre1[,count_U := 0]

for (i in 1:dim(case_0_pre1)[1]){
  x=exactPoiCI(case_0_pre1[i,count],0.95)
  set(case_0_pre1,i,"count_L",x[1])
  set(case_0_pre1,i,"count_U",x[2])
}

case_1_pre[,virus1 := factor(virus,levels = c("CV-A6","CV-A10","CV-A16","EV-A71"),labels = c("CV-A6","CV-A10","CV-A16","EV-A71"))]
case_0_pre[,virus1 := factor(virus,levels = c("CV-A6","CV-A10","CV-A16","EV-A71"),labels = c("CV-A6","CV-A10","CV-A16","EV-A71"))]
case_0_pre1[,virus1 := factor(virus,levels = c("CV-A6","CV-A10","CV-A16","EV-A71"),labels = c("CV-A6","CV-A10","CV-A16","EV-A71"))]
case_1_pre1 = rbind(case_1_pre,case_0_pre1)

saveRDS(case_1_pre1,"est_case_model_exppw.rds")

##### estimated seroprevalence #####
log_beta0 <- po$log_beta0
log_beta1 <- po$log_beta1
log_beta2 <- po$log_beta2

xi <- po$xi
log_one_over_phi1 <- po$log_one_over_phi

sero_pre_0 <- expand.grid(year_end = 2013:2018,age=0,virus=c("CV-A6","CV-A10","CV-A16","EV-A71"))
sero_pre_0 <- as.data.table(sero_pre_0)
sero_pre_0[,birth_year:=year_end]
colnames(sero_pre_0)[1]="year_group"
colnames(sero_pre_0)[2]="age_group"
sero_pre_0 <- merge(sero_pre_0,lambda_index_group,by=c("year_group","age_group","virus"))
virus_index <- data.table(virus=c("CV-A6","CV-A10","CV-A16","EV-A71"),virus_index=1:4)
sero_pre_0 <- merge(sero_pre_0,virus_index,by=c("virus"))

k=1
virus_index <- sero_pre_0[k,]$virus_index
time_virus_index <- sero_pre_0[k,]$time_virus_index
age_index <- sero_pre_0[k,]$age_index
xi_this <- xi[,virus_index]
lambda1 <- exp(log_beta0 + log_beta1[,time_virus_index] + log_beta2[,age_index])

i <- 1
x <- i/12
pr_sus <- xi_this/(xi_this-lambda1)*(exp(-x*lambda1)-exp(-xi_this*x))
quantile(pr_sus,probs=c(0.5,0.025,0.975))
for (i in 2:11){
  x = i/12
  pr_sus1 <- xi_this/(xi_this-lambda1)*(exp(-x*lambda1)-exp(-xi_this*x))
  quantile(pr_sus1,probs=c(0.5,0.025,0.975))
  pr_sus <- pr_sus+pr_sus1
}
pr_sero <- 1-pr_sus/12
est <- quantile(pr_sero,probs=c(0.5,0.025,0.975))
sero_pre_0_1 <- data.table(mean=mean(pr_sero),median=est[1],CI_L=est[2],CI_U=est[3],
                           virus=sero_pre_0[k,]$virus,year_group=sero_pre_0[k,]$year_group,age_group=sero_pre_0[k,]$age_group)


for (k in 2:dim(sero_pre_0)[1]){
  virus_index <- sero_pre_0[k,]$virus_index
  time_virus_index <- sero_pre_0[k,]$time_virus_index
  age_index <- sero_pre_0[k,]$age_index
  xi_this <- xi[,virus_index]
  lambda1 <- exp(log_beta0 + log_beta1[,time_virus_index] + log_beta2[,age_index])
  i <- 1
  x <- i/12
  pr_sus <- xi_this/(xi_this-lambda1)*(exp(-x*lambda1)-exp(-x*xi_this))
  quantile(pr_sus,probs=c(0.5,0.025,0.975))
  for (i in 2:11){
    x = i/12
    pr_sus1 <- xi_this/(xi_this-lambda1)*(exp(-x*lambda1)-exp(-x*xi_this))
    quantile(pr_sus1,probs=c(0.5,0.025,0.975))
    pr_sus <- pr_sus+pr_sus1
  }
  pr_sero <- 1-pr_sus/12
  est <- quantile(pr_sero,probs=c(0.5,0.025,0.975))
  sero_pre_0_2 <- data.table(mean=mean(pr_sero),median=est[1],CI_L=est[2],CI_U=est[3],
                             virus=sero_pre_0[k,]$virus,year_group=sero_pre_0[k,]$year_group,age_group=sero_pre_0[k,]$age_group)
  sero_pre_0_1 <- rbind(sero_pre_0_1,sero_pre_0_2)
}


sero_pre_1 <- expand.grid(year_end = 2013:2018,age=1:6,virus=c("CV-A6","CV-A10","CV-A16","EV-A71"))
sero_pre_1 <- as.data.table(sero_pre_1)
sero_pre_1[,birth_year:=0]
set(sero_pre_1,NULL,"birth_year",sero_pre_1[,year_end-age])
colnames(sero_pre_1)[1]="year_group"
colnames(sero_pre_1)[2]="age_group"

i=1
sero_pre_1[i,]
tmp1 = data.table(year_group=sero_pre_1[i,]$birth_year:sero_pre_1[i,]$year_group,age_group=0:sero_pre_1[i,]$age_group,virus=sero_pre_1[i,]$virus)
set(tmp1,which(tmp1[,year_group<2011]),"year_group",2011)
tmp1[,order1:=1:dim(tmp1)[1]]
tmp1 <- merge(tmp1,lambda_index_group,by=c("age_group","year_group","virus"))
virus_index <- data.table(virus_index=1:4,virus=c("CV-A6","CV-A10","CV-A16","EV-A71"))
tmp1 <- merge(tmp1,virus_index,by="virus")
tmp1 <- tmp1[order(tmp1[,order1]),]

stopifnot(tmp1[1,]$age_group==0)
lambda0_sero <- exp(log_beta0 + log_beta1[,tmp1[1,]$time_virus_index] + log_beta2[,tmp1[1,]$age_index])
xi_this_sero <- xi[,tmp1[1,]$virus_index]


sum_lambda=0
for (j in 1:dim(tmp1)[1]){
  sum_lambda1 <- exp(log_beta0 + log_beta1[,tmp1[j,]$time_virus_index] + log_beta2[,tmp1[j,]$age_index])
  sum_lambda <- sum_lambda+sum_lambda1
}

sero_pre <- 1-xi_this_sero/(xi_this_sero-lambda0_sero)*exp(-sum_lambda)
est <- quantile(sero_pre,probs=c(0.5,0.025,0.975))
sero_pre_1_1 <- data.table(mean=mean(sero_pre),median=est[1],CI_L=est[2],CI_U=est[3],virus=sero_pre_1[i,]$virus,year_group=sero_pre_1[i,]$year_group,age_group=sero_pre_1[i,]$age_group)

for (i in 2:dim(sero_pre_1)[1]){
  print(i)
  tmp1 = data.table(year_group=sero_pre_1[i,]$birth_year:sero_pre_1[i,]$year_group,age_group=0:sero_pre_1[i,]$age_group,virus=sero_pre_1[i,]$virus)
  set(tmp1,which(tmp1[,year_group<2011]),"year_group",2011)
  tmp1[,order1:=1:dim(tmp1)[1]]
  tmp1 <- merge(tmp1,lambda_index_group,by=c("age_group","year_group","virus"))
  tmp1 <- merge(tmp1,virus_index,by="virus")
  tmp1 <- tmp1[order(tmp1[,order1]),]
  
  stopifnot(tmp1[1,]$age_group==0)
  lambda0_sero <- exp(log_beta0 + log_beta1[,tmp1[1,]$time_virus_index] + log_beta2[,tmp1[1,]$age_index])
  xi_this_sero <- xi[,tmp1[1,]$virus_index]
  
  sum_lambda=0
  for (j in 1:dim(tmp1)[1]){
    sum_lambda1 <- exp(log_beta0 + log_beta1[,tmp1[j,]$time_virus_index] + log_beta2[,tmp1[j,]$age_index])
    sum_lambda <- sum_lambda+sum_lambda1
  }
  
  sero_pre <- 1-xi_this_sero/(xi_this_sero-lambda0_sero)*exp(-sum_lambda)
  est <- quantile(sero_pre,probs=c(0.5,0.025,0.975))
  sero_pre_1_2 <- data.table(mean=mean(sero_pre),median=est[1],CI_L=est[2],CI_U=est[3],virus=sero_pre_1[i,]$virus,year_group=sero_pre_1[i,]$year_group,age_group=sero_pre_1[i,]$age_group)
  sero_pre_1_1 <- rbind(sero_pre_1_1,sero_pre_1_2)
}

sero_pre <- rbind(sero_pre_0_1,sero_pre_1_1)

sero_pre[,virus1 := factor(virus,levels = c("CV-A6","CV-A10","CV-A16","EV-A71"),labels = c("CV-A6","CV-A10","CV-A16","EV-A71"))]
saveRDS(sero_pre,"est_seroprevalence_model_exppw.rds")

#### model ExpFt ####
model_fit <- readRDS("model_expft_fit.rds")
po <- rstan::extract(model_fit)

##### convergence evaluation #####
model_su <- summary(model_fit)$summary
min(model_su[,"n_eff"])
max(model_su[,"Rhat"])
trace_name <- rownames(model_su[order(model_su[,"n_eff"]),][1:4,])
# Trace plot
po_trace <- rstan:::extract(model_fit, inc_warmup = TRUE, permuted = FALSE,
                            pars = trace_name)
bayesplot:::color_scheme_set("mix-blue-pink")
trace <- bayesplot:::mcmc_trace(po_trace, n_warmup = 1e3,
                                facet_args = list(nrow=1,ncol=4,labeller = label_parsed))+
  theme_bw()+
  theme(legend.position='bottom')

##### force of infection #####
log_beta0 <- po$log_beta0
log_beta1 <- po$log_beta1
beta2 <- po$beta2

lambda_index_group_apre <-lambda_index_group[,{
  lambda1 <- exp(log_beta0 + log_beta1[,time_virus_index] - beta2*age_group-log(beta2))*((age_group+1/beta2)*(1-exp(-beta2)) - exp(-beta2))
  est = quantile(lambda1,probs = c(0.5, 0.025, 0.975))
  list(mean=mean(lambda1),median=est[1],CI_L=est[2],CI_U=est[3])
},by=c("age_group","year_group","virus","index")]

lambda_index_group_apre[,virus1 := factor(virus,levels = c("CV-A6","CV-A10","CV-A16","EV-A71"),labels = c("CV-A6","CV-A10","CV-A16","EV-A71"))]
saveRDS(lambda_index_group_apre,"est_foi_model_expft.rds")

##### estimated cases #####
xi <- po$xi
log_one_over_phi1 <- po$log_one_over_phi

case1[,age:=floor(`Age in years`)]
case1_0 <- case1[which(case1[,age==0]),c(2,5,6,7,10,13)]
case1_0[,month:=0]
set(case1_0,which(case1_0[,(`Age in years`>=1/12)&(`Age in years`<2/12)]),"month",1)
set(case1_0,which(case1_0[,(`Age in years`>=2/12)&(`Age in years`<3/12)]),"month",2)
set(case1_0,which(case1_0[,(`Age in years`>=3/12)&(`Age in years`<4/12)]),"month",3)
set(case1_0,which(case1_0[,(`Age in years`>=4/12)&(`Age in years`<5/12)]),"month",4)
set(case1_0,which(case1_0[,(`Age in years`>=5/12)&(`Age in years`<6/12)]),"month",5)
set(case1_0,which(case1_0[,(`Age in years`>=6/12)&(`Age in years`<7/12)]),"month",6)
set(case1_0,which(case1_0[,(`Age in years`>=7/12)&(`Age in years`<8/12)]),"month",7)
set(case1_0,which(case1_0[,(`Age in years`>=8/12)&(`Age in years`<9/12)]),"month",8)
set(case1_0,which(case1_0[,(`Age in years`>=9/12)&(`Age in years`<10/12)]),"month",9)
set(case1_0,which(case1_0[,(`Age in years`>=10/12)&(`Age in years`<11/12)]),"month",10)
set(case1_0,which(case1_0[,(`Age in years`>=11/12)&(`Age in years`<1)]),"month",11)

summary(case1_0)

case2_0 <- case1_0[,list(count=.N),by=c("virus","year_sampling","age","month")]
colnames(case2_0)[2]="year_group"
colnames(case2_0)[3]="age_group"

tmp <- expand.grid(virus=c("CV-A6","CV-A10","CV-A16","EV-A71"),year_group=2013:2018,age_group=0,month=1:11)
tmp$self_index <- 1:dim(tmp)[1]
tmp <- as.data.table(tmp)

case2_0 <- merge(case2_0,tmp,by=c("virus","year_group","age_group","month"))

zero_index1 <- tmp[which(tmp[,!self_index%in%case2_0$self_index]),]
zero_index1[,count:=0]

case3_0 <- rbind(case2_0,zero_index1)
case3_0 <- case3_0[order(case3_0[,self_index]),]

case4_0 <- merge(case3_0,lambda_index_group,by=c("virus","year_group","age_group"))
summary(case4_0)
case4_0 <- case4_0[order(case4_0[,self_index]),]

casepo
colnames(casepo)[1]="age_group"
colnames(casepo)[2]="year_group"

case5_0 <- merge(case4_0,casepo[,1:3],by=c("age_group","year_group"))

virus_index <- data.table(virus=c("CV-A6","CV-A10","CV-A16","EV-A71"),virus_index=1:4)
case6_0 <- merge(case5_0,virus_index,by="virus")
case6_0 <- case6_0[order(case6_0[,self_index]),]

case_0_pre <- case6_0[,{
  xi_this = xi[,virus_index]
  phi_this = exp(-log_one_over_phi1[,virus_index])
  lambda1 = exp(log_beta0 + log_beta1[,time_virus_index]-log(beta2))*(1/beta2*(1-exp(-beta2)) - exp(-beta2))
  case_pre1 = pop*xi_this/(xi_this-lambda1)*(exp(-month*lambda1/12)-exp(-xi_this*month/12))*(1-exp(-lambda1/12))*phi_this
  est = quantile(case_pre1,probs = c(0.5, 0.025, 0.975))
  list(mean=mean(case_pre1),median=est[1],CI_L=est[2],CI_U=est[3])
},by=c("age_group","year_group","virus","month","count","self_index")]


case7_0 <- case6_0[,list(count1=sum(count)),by=c("age_group","year_group","virus","index","time_virus_index","age_index","virus_index","pop")]
i=1
xi_this = xi[,case7_0[i,]$virus_index]
phi_this = exp(-log_one_over_phi1[,case7_0[i,]$virus_index])
lambda1 = exp(log_beta0 + log_beta1[,case7_0[i,]$time_virus_index]-log(beta2))*(1/beta2*(1-exp(-beta2)) - exp(-beta2))
case_pre1 = 0
for (x in 1:11){
  case_pre2 = case7_0[i,]$pop*xi_this/(xi_this-lambda1)*(exp(-x*lambda1/12)-exp(-xi_this*x/12))*(1-exp(-lambda1/12))*phi_this
  case_pre1 = case_pre1+case_pre2
}
est = quantile(case_pre1,probs = c(0.5, 0.025, 0.975))
case_0_pre1 = data.table(mean=mean(case_pre1),median=est[1],CI_L=est[2],CI_U=est[3],age_group=0,year_group=case7_0[i,]$year_group,virus=case7_0[i,]$virus,count=case7_0[i,]$count1)
for (i in 2:dim(case7_0)[1]){
  xi_this = xi[,case7_0[i,]$virus_index]
  phi_this = exp(-log_one_over_phi1[,case7_0[i,]$virus_index])
  lambda1 = exp(log_beta0 + log_beta1[,case7_0[i,]$time_virus_index]-log(beta2))*(1/beta2*(1-exp(-beta2)) - exp(-beta2))
  case_pre1 = 0
  for (x in 1:11){
    case_pre2 = case7_0[i,]$pop*xi_this/(xi_this-lambda1)*(exp(-x*lambda1/12)-exp(-xi_this*x/12))*(1-exp(-lambda1/12))*phi_this
    case_pre1 = case_pre1+case_pre2
  }
  est = quantile(case_pre1,probs = c(0.5, 0.025, 0.975))
  case_0_pre2 = data.table(mean=mean(case_pre1),median=est[1],CI_L=est[2],CI_U=est[3],age_group=0,year_group=case7_0[i,]$year_group,virus=case7_0[i,]$virus,count=case7_0[i,]$count1)
  case_0_pre1 = rbind(case_0_pre1,case_0_pre2)
}



case1_1 <- case1[which(case1[,age>0]),c(2,5,6,7,10,13)]
case1_1[,age_group:=age]
set(case1_1,which(case1_1[,age_group>=7]),"age_group",7)
case2_1 <- case1_1[,list(count=.N),by=c("virus","year_sampling","age_group")]
colnames(case2_1)[2]="year_group"



tmp2 <- expand.grid(virus=c("CV-A6","CV-A10","CV-A16","EV-A71"),year_group=2013:2018,age_group=1:7)
tmp2$self_index <- 1:dim(tmp2)[1]
tmp2 <- as.data.table(tmp2)

case2_1 <- merge(case2_1, tmp2, by=c("virus","year_group","age_group"))
zero_index2 <- tmp2[which(tmp2[,!self_index%in%case2_1$self_index]),]
zero_index2[,count:=0]

case3_1 <- rbind(case2_1,zero_index2)
case3_1 <- case3_1[order(case3_1[,self_index]),]

case4_1 <- merge(case3_1,lambda_index_group,by=c("virus","year_group","age_group"))
case4_1 <- case4_1[order(case4_1[,self_index]),]

virus_index
case5_1 <- merge(case4_1,virus_index,by="virus")
case5_1 <- case5_1[order(case5_1[,self_index]),]

pop_group <- expand.grid(year_group=2013:2018,age_group=1:7)
pop_group$pop <- NA
for (i in 1:dim(pop_group)[1]){
  age <- pop_group[i,]$age_group
  year1 <- pop_group[i,]$year_group
  if (age<7){
    tmp = data.table(age_group=age,year_group=year1)
    tmp1 = merge(tmp,casepo,by=c("age_group","year_group"))
    pop_group[i,]$pop = tmp1$pop
  }else{
    tmp = expand.grid(age_group=7:11,year_group=year1)
    tmp1 = merge(tmp,casepo,by=c("age_group","year_group"))
    pop_group[i,]$pop = sum(tmp1$pop)
  }
}

case6_1 <- merge(case5_1,pop_group,by=c("year_group","age_group"))
case6_1 <- case6_1[order(case6_1[,self_index]),]


i=1
case6_1[i,]
virus1 <- case6_1[i,]$virus
age1 <- case6_1[i,]$age_group
year1 <- case6_1[i,]$year_group

if (age1==7) age1=9

data1 <- data.table(year_group = (year1-age1):(year1-1),age_group = 0:(age1-1),virus=virus1)
set(data1,which(data1[,age_group>7]),"age_group",7)
set(data1,which(data1[,year_group<2011]),"year_group",2011)
data1[,order1 := 1:dim(data1)[1]]
data1 = merge(data1,lambda_index_group,by=c("virus","age_group","year_group"))
data1 <- data1[order(data1[,order1]),]

sum_lambda=0
for (j in 1:dim(data1)[1]){
  sum_lambda1 <- exp(log_beta0 + log_beta1[,data1[j,]$time_virus_index] - beta2*data1[j,]$age_group-log(beta2))*((data1[j,]$age_group+1/beta2)*(1-exp(-beta2)) - exp(-beta2))
  sum_lambda <- sum_lambda + sum_lambda1
}

quantile(sum_lambda,probs = c(0.5, 0.025, 0.975))
quantile(lambda_this,probs = c(0.5, 0.025, 0.975))

xi_this <- xi[,case6_1[i,]$virus_index]
phi_this <- exp(-log_one_over_phi1[,case6_1[i,]$virus_index])
stopifnot(data1[1,]$age_group==0)
lambda0_this <- exp(log_beta0 + log_beta1[,data1[1,]$time_virus_index]-log(beta2))*(1/beta2*(1-exp(-beta2)) - exp(-beta2))
lambda_this <- exp(log_beta0 + log_beta1[,case6_1[i,]$time_virus_index] - beta2*age1-log(beta2))*((age1+1/beta2)*(1-exp(-beta2)) - exp(-beta2))
case_pre1 <- case6_1[i,]$pop*xi_this/(xi_this-lambda0_this)*exp(-sum_lambda)*(1-exp(-lambda_this))*phi_this
est = quantile(case_pre1,probs = c(0.5, 0.025, 0.975))
case_1_pre <- data.table(mean=mean(case_pre1),median=est[1],CI_L=est[2],CI_U=est[3],age_group=age1,year_group=year1,virus=virus1,count=case6_1[i,]$count)


for (i in 2:dim(case6_1)[1]){
  virus1 <- case6_1[i,]$virus
  age1 <- case6_1[i,]$age_group
  year1 <- case6_1[i,]$year_group
  
  if (age1==7) age1=9
  
  data1 <- data.table(year_group = (year1-age1):(year1-1),age_group = 0:(age1-1),virus=virus1)
  set(data1,which(data1[,age_group>7]),"age_group",7)
  set(data1,which(data1[,year_group<2011]),"year_group",2011)
  data1[,order1 := 1:dim(data1)[1]]
  data1 = merge(data1,lambda_index_group,by=c("virus","age_group","year_group"))
  data1 <- data1[order(data1[,order1]),]
  
  sum_lambda=0
  for (j in 1:dim(data1)[1]){
    sum_lambda1 <- exp(log_beta0 + log_beta1[,data1[j,]$time_virus_index] - beta2*data1[j,]$age_group-log(beta2))*((data1[j,]$age_group+1/beta2)*(1-exp(-beta2)) - exp(-beta2))
    sum_lambda <- sum_lambda + sum_lambda1
  }
  
  xi_this <- xi[,case6_1[i,]$virus_index]
  phi_this <- exp(-log_one_over_phi1[,case6_1[i,]$virus_index])
  stopifnot(data1[1,]$age_group==0)
  lambda0_this <- exp(log_beta0 + log_beta1[,data1[1,]$time_virus_index]-log(beta2))*(1/beta2*(1-exp(-beta2)) - exp(-beta2))
  lambda_this <- exp(log_beta0 + log_beta1[,case6_1[i,]$time_virus_index] - beta2*age1-log(beta2))*((age1+1/beta2)*(1-exp(-beta2)) - exp(-beta2))
  case_pre1 <- case6_1[i,]$pop*xi_this/(xi_this-lambda0_this)*exp(-sum_lambda)*(1-exp(-lambda_this))*phi_this
  est = quantile(case_pre1,probs = c(0.5, 0.025, 0.975))
  case_1_pre1 <- data.table(mean=mean(case_pre1),median=est[1],CI_L=est[2],CI_U=est[3],age_group=age1,year_group=year1,virus=virus1,count=case6_1[i,]$count)
  case_1_pre <- rbind(case_1_pre,case_1_pre1)
}
set(case_1_pre,which(case_1_pre[,age_group>7]),"age_group",7)

exactPoiCI <- function (X, conf.level=0.95) {
  alpha = 1 - conf.level
  upper <- 0.5 * qchisq(1-alpha/2, 2*X+2)
  lower <- 0.5 * qchisq(alpha/2, 2*X)
  return(c(lower, upper))
}
exactPoiCI(30,0.95)
exactPoiCI(0,0.95)

case_1_pre[,count_L := 0]
case_1_pre[,count_U := 0]

for (i in 1:dim(case_1_pre)[1]){
  x=exactPoiCI(case_1_pre[i,count],0.95)
  set(case_1_pre,i,"count_L",x[1])
  set(case_1_pre,i,"count_U",x[2])
}

case_0_pre[,count_L := 0]
case_0_pre[,count_U := 0]

for (i in 1:dim(case_0_pre)[1]){
  x=exactPoiCI(case_0_pre[i,count],0.95)
  set(case_0_pre,i,"count_L",x[1])
  set(case_0_pre,i,"count_U",x[2])
}

case_0_pre1[,count_L := 0]
case_0_pre1[,count_U := 0]

for (i in 1:dim(case_0_pre1)[1]){
  x=exactPoiCI(case_0_pre1[i,count],0.95)
  set(case_0_pre1,i,"count_L",x[1])
  set(case_0_pre1,i,"count_U",x[2])
}

summary(case_0_pre)
case_1_pre[,virus1 := factor(virus,levels = c("CV-A6","CV-A10","CV-A16","EV-A71"),labels = c("CV-A6","CV-A10","CV-A16","EV-A71"))]
case_0_pre[,virus1 := factor(virus,levels = c("CV-A6","CV-A10","CV-A16","EV-A71"),labels = c("CV-A6","CV-A10","CV-A16","EV-A71"))]
case_0_pre1[,virus1 := factor(virus,levels = c("CV-A6","CV-A10","CV-A16","EV-A71"),labels = c("CV-A6","CV-A10","CV-A16","EV-A71"))]
case_1_pre1 = rbind(case_1_pre,case_0_pre1)

saveRDS(case_1_pre1,"est_case_model_expft.rds")

##### estimated seroprevalence #####
log_beta0 <- po$log_beta0
log_beta1 <- po$log_beta1
beta2 <- po$beta2

xi <- po$xi
log_one_over_phi1 <- po$log_one_over_phi

sero_pre_0 <- expand.grid(year_end = 2013:2018,age=0,virus=c("CV-A6","CV-A10","CV-A16","EV-A71"))
sero_pre_0 <- as.data.table(sero_pre_0)
sero_pre_0[,birth_year:=year_end]
colnames(sero_pre_0)[1]="year_group"
colnames(sero_pre_0)[2]="age_group"
sero_pre_0 <- merge(sero_pre_0,lambda_index_group,by=c("year_group","age_group","virus"))
virus_index <- data.table(virus=c("CV-A6","CV-A10","CV-A16","EV-A71"),virus_index=1:4)
sero_pre_0 <- merge(sero_pre_0,virus_index,by=c("virus"))

k=1
virus_index <- sero_pre_0[k,]$virus_index
time_virus_index <- sero_pre_0[k,]$time_virus_index
age_group <- sero_pre_0[k,]$age_group
xi_this <- xi[,virus_index]
lambda1 <- exp(log_beta0 + log_beta1[,time_virus_index] - beta2*age_group-log(beta2))*((age_group+1/beta2)*(1-exp(-beta2)) - exp(-beta2))

i <- 1
x <- i/12
pr_sus <- xi_this/(xi_this-lambda1)*(exp(-x*lambda1)-exp(-xi_this*x))
quantile(pr_sus,probs=c(0.5,0.025,0.975))
for (i in 2:11){
  x = i/12
  pr_sus1 <- xi_this/(xi_this-lambda1)*(exp(-x*lambda1)-exp(-xi_this*x))
  quantile(pr_sus1,probs=c(0.5,0.025,0.975))
  pr_sus <- pr_sus+pr_sus1
}
pr_sero <- 1-pr_sus/12
est <- quantile(pr_sero,probs=c(0.5,0.025,0.975))
sero_pre_0_1 <- data.table(mean=mean(pr_sero),median=est[1],CI_L=est[2],CI_U=est[3],
                           virus=sero_pre_0[k,]$virus,year_group=sero_pre_0[k,]$year_group,age_group=sero_pre_0[k,]$age_group)


for (k in 2:dim(sero_pre_0)[1]){
  virus_index <- sero_pre_0[k,]$virus_index
  time_virus_index <- sero_pre_0[k,]$time_virus_index
  age_group <- sero_pre_0[k,]$age_group
  xi_this <- xi[,virus_index]
  lambda1 <- exp(log_beta0 + log_beta1[,time_virus_index] - beta2*age_group-log(beta2))*((age_group+1/beta2)*(1-exp(-beta2)) - exp(-beta2))
  i <- 1
  x <- i/12
  pr_sus <- xi_this/(xi_this-lambda1)*(exp(-x*lambda1)-exp(-x*xi_this))
  quantile(pr_sus,probs=c(0.5,0.025,0.975))
  for (i in 2:11){
    x = i/12
    pr_sus1 <- xi_this/(xi_this-lambda1)*(exp(-x*lambda1)-exp(-x*xi_this))
    quantile(pr_sus1,probs=c(0.5,0.025,0.975))
    pr_sus <- pr_sus+pr_sus1
  }
  pr_sero <- 1-pr_sus/12
  est <- quantile(pr_sero,probs=c(0.5,0.025,0.975))
  sero_pre_0_2 <- data.table(mean=mean(pr_sero),median=est[1],CI_L=est[2],CI_U=est[3],
                             virus=sero_pre_0[k,]$virus,year_group=sero_pre_0[k,]$year_group,age_group=sero_pre_0[k,]$age_group)
  sero_pre_0_1 <- rbind(sero_pre_0_1,sero_pre_0_2)
}


sero_pre_1 <- expand.grid(year_end = 2013:2018,age=1:6,virus=c("CV-A6","CV-A10","CV-A16","EV-A71"))
sero_pre_1 <- as.data.table(sero_pre_1)
sero_pre_1[,birth_year:=0]
set(sero_pre_1,NULL,"birth_year",sero_pre_1[,year_end-age])
colnames(sero_pre_1)[1]="year_group"
colnames(sero_pre_1)[2]="age_group"

i=1
sero_pre_1[i,]
tmp1 = data.table(year_group=sero_pre_1[i,]$birth_year:sero_pre_1[i,]$year_group,age_group=0:sero_pre_1[i,]$age_group,virus=sero_pre_1[i,]$virus)
set(tmp1,which(tmp1[,year_group<2011]),"year_group",2011)
tmp1[,order1:=1:dim(tmp1)[1]]
tmp1 <- merge(tmp1,lambda_index_group,by=c("age_group","year_group","virus"))
virus_index <- data.table(virus_index=1:4,virus=c("CV-A6","CV-A10","CV-A16","EV-A71"))
tmp1 <- merge(tmp1,virus_index,by="virus")
tmp1 <- tmp1[order(tmp1[,order1]),]

stopifnot(tmp1[1,]$age_group==0)
lambda0_sero <- exp(log_beta0 + log_beta1[,tmp1[1,]$time_virus_index] - beta2*tmp1[1,]$age_group-log(beta2))*((tmp1[1,]$age_group+1/beta2)*(1-exp(-beta2)) - exp(-beta2))
xi_this_sero <- xi[,tmp1[1,]$virus_index]


sum_lambda=0
for (j in 1:dim(tmp1)[1]){
  sum_lambda1 <- exp(log_beta0 + log_beta1[,tmp1[j,]$time_virus_index] - beta2*tmp1[j,]$age_group-log(beta2))*((tmp1[j,]$age_group+1/beta2)*(1-exp(-beta2)) - exp(-beta2))
  sum_lambda <- sum_lambda+sum_lambda1
}

sero_pre <- 1-xi_this_sero/(xi_this_sero-lambda0_sero)*exp(-sum_lambda)
est <- quantile(sero_pre,probs=c(0.5,0.025,0.975))
sero_pre_1_1 <- data.table(mean=mean(sero_pre),median=est[1],CI_L=est[2],CI_U=est[3],virus=sero_pre_1[i,]$virus,year_group=sero_pre_1[i,]$year_group,age_group=sero_pre_1[i,]$age_group)

for (i in 2:dim(sero_pre_1)[1]){
  print(i)
  tmp1 = data.table(year_group=sero_pre_1[i,]$birth_year:sero_pre_1[i,]$year_group,age_group=0:sero_pre_1[i,]$age_group,virus=sero_pre_1[i,]$virus)
  set(tmp1,which(tmp1[,year_group<2011]),"year_group",2011)
  tmp1[,order1:=1:dim(tmp1)[1]]
  tmp1 <- merge(tmp1,lambda_index_group,by=c("age_group","year_group","virus"))
  tmp1 <- merge(tmp1,virus_index,by="virus")
  tmp1 <- tmp1[order(tmp1[,order1]),]
  
  stopifnot(tmp1[1,]$age_group==0)
  lambda0_sero <- exp(log_beta0 + log_beta1[,tmp1[1,]$time_virus_index] - beta2*tmp1[1,]$age_group-log(beta2))*((tmp1[1,]$age_group+1/beta2)*(1-exp(-beta2)) - exp(-beta2))
  xi_this_sero <- xi[,tmp1[1,]$virus_index]
  
  sum_lambda=0
  for (j in 1:dim(tmp1)[1]){
    sum_lambda1 <- exp(log_beta0 + log_beta1[,tmp1[j,]$time_virus_index] - beta2*tmp1[j,]$age_group-log(beta2))*((tmp1[j,]$age_group+1/beta2)*(1-exp(-beta2)) - exp(-beta2))
    sum_lambda <- sum_lambda+sum_lambda1
  }
  
  sero_pre <- 1-xi_this_sero/(xi_this_sero-lambda0_sero)*exp(-sum_lambda)
  est <- quantile(sero_pre,probs=c(0.5,0.025,0.975))
  sero_pre_1_2 <- data.table(mean=mean(sero_pre),median=est[1],CI_L=est[2],CI_U=est[3],virus=sero_pre_1[i,]$virus,year_group=sero_pre_1[i,]$year_group,age_group=sero_pre_1[i,]$age_group)
  sero_pre_1_1 <- rbind(sero_pre_1_1,sero_pre_1_2)
}

sero_pre <- rbind(sero_pre_0_1,sero_pre_1_1)

sero_pre[,virus1 := factor(virus,levels = c("CV-A6","CV-A10","CV-A16","EV-A71"),labels = c("CV-A6","CV-A10","CV-A16","EV-A71"))]
saveRDS(sero_pre,"est_seroprevalence_model_expft.rds")

