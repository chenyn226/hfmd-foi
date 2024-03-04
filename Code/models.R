
library(rstan)
library(readxl)
library(data.table)

#### load data ####
# case data
case1 <- readRDS("case_data_preprocessed.rds")

# population size data
casepo <- readRDS("hochiminhcity_population.rds")

# serological survey data
sero1 <- read.csv("SampleHFMD_sero_2015.csv")
sero2 <- read_xls("SampleHFMD_sero_2017.xls")

# lambda index
lambda_index_group <- readRDS("lambda_index_group.rds")

case1 <- as.data.table(case1)
casepo <- as.data.table(casepo)
lambda_index_group <- as.data.table(lambda_index_group)
sero1 <- as.data.table(sero1)
sero2 <- as.data.table(sero2)

#### model construction ####
## notes about the components in "data" of the code below
# K: the total number of age-time-serotype specific force of infection (lambda)
# K_time_virus: the total number of parameters varied by time and serotypes
# K_age: the total number of parameters varied by age
# time_virus_index: a vector recording the corresponding time-virus index of each force of infection (lambda)
# age_index: a vector recording the corresponding age index of each force of infection (lambda)
# K_virus: the total number of parameters varied by serotypes

# N_group1: the total number of cases varied by age, time and serotypes for those with age > 1 year old
# lambda_index_group1: a vector recording the corresponding force of infection (lambda) index for cases with age > 1 year old
# log_pop_group1: the age-time specific log(population size) for cases with age > 1 year old
# virus_index1: a vector recording the corresponding virus index for cases with age > 1 year old
# case_group11: the number of age-time-serotype specific number of cases for those with age > 1 year old
# matrix_length_group1: the number of columns of matrix "lambda_sum_index_group1"
# lambda_sum_index_group1: the calculation Pr(susceptible) includes exp(-sum(lambda)), this matrix records the lambda to be summed for each cases aged > 1 year old
# lambda_length_group1: for different case, the number of lambdas in exp(-sum(lambda)) is different, this vector records the corresponding number of lambda needed to be summed for each case aged > 1 year old
# lambda_0_index_group1: the calculation Pr(susceptible) includes lambda_{(0,1),t}^S, this vector records the index corresponding lambda_{(0,1),t}^S for each case aged > 1 year old

# N_group0: the total number of cases varied by age, time and serotypes for those with age < 1 year old
# log_pop_group0: the age-time specific log(population size) for cases with age < 1 year old
# age_x: for age < 1 year old, the model consider month unit for the age to estimate maternal immunity, this is a vector recording month
# virus_index0: a vector recording the corresponding virus index for cases with age < 1 year old
# lambda_index0: the calculation Pr(susceptible) and Pr(infection) includes lambda_{(0,1),t}^S, this vector records the index corresponding lambda_{(0,1),t}^S for each case aged < 1 year old
# case_group00: the number of age-time-serotype specific number of cases for those with age < 1 year old

# Nsero0neg: in serological data, the total number of subjects varied by age, time and serotypes for those with age < 1 year old
# age_x_sero0neg: a vector recording month for subjects with age < 1 year old in serological data
# virus_index_sero0neg: a vector recording the corresponding virus index for subjects with age < 1 year old in serological data
# lambda_index_sero0neg: a vector recording the corresponding lambda index for subjects with age < 1 year old in serological data

# Nsero0pos, age_x_sero0pos, virus_index_sero0pos, lambda_index_sero0pos: 
# similar to Nsero0neg, age_x_sero0neg, virus_index_sero0neg, lambda_index_sero0neg, but for positive results

# Nsero1neg, lambda_index_group_sero1neg, virus_index_sero1neg:
# similar to Nsero0neg, age_x_sero0neg, virus_index_sero0neg, but for age > 1 year old
# Csero1neg: number of age-time-serotype specific sujects aged > 1 year old with negative antibody testing result
# matrix_length_group_sero1neg, lambda_sum_index_group_sero1neg, lambda_length_group_sero1neg, lambda_0_index_group_sero1neg:
# similar to matrix_length_group1, lambda_sum_index_group1, lambda_length_group1, lambda_0_index_group1, but for serological data with negative antibody

# Nsero1pos, lambda_index_group_sero1pos, virus_index_sero1pos, Csero1pos, matrix_length_group_sero1pos, lambda_sum_index_group_sero1pos, lambda_length_group_sero1pos, lambda_0_index_group_sero1pos
# similar to Nsero1neg, lambda_index_group_sero1neg, virus_index_sero1neg, ..., but for seropositive results

##### model LinPw #####
# tricks in the code to improve efficiency:
# log1m_exp() is more efficient than log(1-exp()), so use lambda/xi + log1m_exp(-lambda*xi) to replace log(exp(lambda/xi)-1)
# -lambda[lambda_index0[i]]/20*log_sum_exp(-1.667*age_x[i],-20/xi[virus_index0[i]]) is the smooth version of lambda*min{a,1/xi}

model_linpw_txt <- "
data{
int<lower=1> K;
int<lower=1> K_time_virus;
int<lower=1> K_age;
int<lower=1,upper=K_time_virus> time_virus_index[K];
int<lower=1,upper=K_age> age_index[K];
int<lower=1> K_virus;


int<lower=1> N_group1;
int<lower=1,upper=K> lambda_index_group1[N_group1];
real<lower=0> log_pop_group1[N_group1];
int<lower=1,upper=K_virus> virus_index1[N_group1];
int<lower=0> case_group11[N_group1];

int<lower=1> matrix_length_group1;
int<lower=0,upper=K> lambda_sum_index_group1[N_group1,matrix_length_group1];
int<lower=1,upper=matrix_length_group1> lambda_length_group1[N_group1];
int<lower=0,upper=K> lambda_0_index_group1[N_group1];


int<lower=1> N_group0;
vector<lower=0>[N_group0] log_pop_group0;
vector<lower=0.1>[N_group0] age_x;
int<lower=1,upper=K_virus> virus_index0[N_group0];
int<lower=1,upper=K> lambda_index0[N_group0];
int<lower=0> case_group00[N_group0];


int<lower=1> Nsero0neg;
vector<lower=0.1>[Nsero0neg] age_x_sero0neg;
int<lower=1,upper=K_virus> virus_index_sero0neg[Nsero0neg];
int<lower=1,upper=K> lambda_index_sero0neg[Nsero0neg];

int<lower=1> Nsero0pos;
vector<lower=0.1>[Nsero0pos] age_x_sero0pos;
int<lower=1,upper=K_virus> virus_index_sero0pos[Nsero0pos];
int<lower=1,upper=K> lambda_index_sero0pos[Nsero0pos];

int<lower=1> Nsero1neg;
int<lower=1,upper=K> lambda_index_group_sero1neg[Nsero1neg];
int<lower=1,upper=K_virus> virus_index_sero1neg[Nsero1neg];
int<lower=1> Csero1neg[Nsero1neg];

int<lower=1> matrix_length_group_sero1neg;
int<lower=0,upper=K> lambda_sum_index_group_sero1neg[Nsero1neg,matrix_length_group_sero1neg];
int<lower=1,upper=matrix_length_group_sero1neg> lambda_length_group_sero1neg[Nsero1neg];
int<lower=0,upper=K> lambda_0_index_group_sero1neg[Nsero1neg];

int<lower=1> Nsero1pos;
int<lower=1,upper=K> lambda_index_group_sero1pos[Nsero1pos];
int<lower=1,upper=K_virus> virus_index_sero1pos[Nsero1pos];
int<lower=1> Csero1pos[Nsero1pos];

int<lower=1> matrix_length_group_sero1pos;
int<lower=0,upper=K> lambda_sum_index_group_sero1pos[Nsero1pos,matrix_length_group_sero1pos];
int<lower=1,upper=matrix_length_group_sero1pos> lambda_length_group_sero1pos[Nsero1pos];
int<lower=0,upper=K> lambda_0_index_group_sero1pos[Nsero1pos];

}
parameters{
real log_beta0;
vector[K_time_virus] log_beta1;
vector[K_age] log_beta2;

real<lower=0> sigma_beta1;
real<lower=0> sigma_beta2;
vector<lower=0>[K_virus] log_one_over_phi; 
vector<lower=0>[K_virus] xi0;

}
transformed parameters{
vector<lower=0>[K] lambda; 
vector[K_virus] xi=xi0+1;
vector[N_group1] log_mean1;
vector[N_group0] log_mean0;

lambda = exp(log_beta0 + log_beta1[time_virus_index] + log_beta2[age_index]);

for (j in 1:N_group1){
log_mean1[j] =  - sum(lambda[lambda_sum_index_group1[j,1:lambda_length_group1[j]]]) + log_pop_group1[j] + log(xi[virus_index1[j]]) - log(lambda[lambda_0_index_group1[j]]) + lambda[lambda_0_index_group1[j]]/xi[virus_index1[j]] + log1m_exp(-lambda[lambda_0_index_group1[j]]/xi[virus_index1[j]]) + log1m_exp(-lambda[lambda_index_group1[j]]) - log_one_over_phi[virus_index1[j]];
}

for (i in 1:N_group0){
log_mean0[i] = log_pop_group0[i] + log(xi[virus_index0[i]]) - log(lambda[lambda_index0[i]]) - age_x[i]/12*lambda[lambda_index0[i]] + log(exp(-lambda[lambda_index0[i]]/20*log_sum_exp(-1.667*age_x[i],-20/xi[virus_index0[i]]))-1) + log1m_exp(-lambda[lambda_index0[i]]/12) - log_one_over_phi[virus_index0[i]];
}

}
model{
case_group11 ~ poisson_log(log_mean1);
case_group00 ~ poisson_log(log_mean0);

for (i in 1:Nsero0neg){
target += log(xi[virus_index_sero0neg[i]]) - log(lambda[lambda_index_sero0neg[i]]) - age_x_sero0neg[i]*lambda[lambda_index_sero0neg[i]] + log(exp(-lambda[lambda_index_sero0neg[i]]/20*log_sum_exp(-20*age_x_sero0neg[i],-20/xi[virus_index_sero0neg[i]]))-1);
}
for (i in 1:Nsero0pos){
target += log1m(xi[virus_index_sero0pos[i]]/lambda[lambda_index_sero0pos[i]]*exp(-age_x_sero0pos[i]*lambda[lambda_index_sero0pos[i]])*(exp(-lambda[lambda_index_sero0pos[i]]/20*log_sum_exp(-20*age_x_sero0pos[i],-20/xi[virus_index_sero0pos[i]]))-1));
}
for (j in 1:Nsero1neg){
target += Csero1neg[j]*(- sum(lambda[lambda_sum_index_group_sero1neg[j,1:lambda_length_group_sero1neg[j]]]) + log(xi[virus_index_sero1neg[j]]) - log(lambda[lambda_0_index_group_sero1neg[j]]) + lambda[lambda_0_index_group_sero1neg[j]]/xi[virus_index_sero1neg[j]] + log1m_exp(-lambda[lambda_0_index_group_sero1neg[j]]/xi[virus_index_sero1neg[j]]));
}
for (j in 1:Nsero1pos){
target += Csero1pos[j]*log1m(exp(- sum(lambda[lambda_sum_index_group_sero1pos[j,1:lambda_length_group_sero1pos[j]]]))*xi[virus_index_sero1pos[j]]/lambda[lambda_0_index_group_sero1pos[j]]*(exp(lambda[lambda_0_index_group_sero1pos[j]]/xi[virus_index_sero1pos[j]])-1));
}
log_beta0 ~ normal(0,5);
log_beta1 ~ normal(0,sigma_beta1);
log_beta2 ~ normal(0,sigma_beta2);
sigma_beta1 ~ cauchy(0,1);
sigma_beta2 ~ cauchy(0,1);

log_one_over_phi ~ exponential(0.1);
xi0 ~ exponential(1);
}
"
model_linpw_compiled <- rstan::stan_model(
  model_name = 'model_linpw',
  model_code = gsub('\t','',model_linpw_txt)
)

##### model LinFt #####
model_linft_txt <- "
data{
int<lower=1> K;
vector<lower=0>[K] lambda_age;
int<lower=1> K_time_virus;
int<lower=1,upper=K_time_virus> time_virus_index[K];
int<lower=1> K_virus;


int<lower=1> N_group1;
int<lower=1,upper=K> lambda_index_group1[N_group1];
real<lower=0> log_pop_group1[N_group1];
int<lower=1,upper=K_virus> virus_index1[N_group1];
int<lower=0> case_group11[N_group1];

int<lower=1> matrix_length_group1;
int<lower=0,upper=K> lambda_sum_index_group1[N_group1,matrix_length_group1];
int<lower=1,upper=matrix_length_group1> lambda_length_group1[N_group1];
int<lower=0,upper=K> lambda_0_index_group1[N_group1];


int<lower=1> N_group0;
vector<lower=0>[N_group0] log_pop_group0;
vector<lower=0.1>[N_group0] age_x;
int<lower=1,upper=K_virus> virus_index0[N_group0];
int<lower=1,upper=K> lambda_index0[N_group0];
int<lower=0> case_group00[N_group0];


int<lower=1> Nsero0neg;
vector<lower=0.1>[Nsero0neg] age_x_sero0neg;
int<lower=1,upper=K_virus> virus_index_sero0neg[Nsero0neg];
int<lower=1,upper=K> lambda_index_sero0neg[Nsero0neg];

int<lower=1> Nsero0pos;
vector<lower=0.1>[Nsero0pos] age_x_sero0pos;
int<lower=1,upper=K_virus> virus_index_sero0pos[Nsero0pos];
int<lower=1,upper=K> lambda_index_sero0pos[Nsero0pos];

int<lower=1> Nsero1neg;
int<lower=1,upper=K> lambda_index_group_sero1neg[Nsero1neg];
int<lower=1,upper=K_virus> virus_index_sero1neg[Nsero1neg];
int<lower=1> Csero1neg[Nsero1neg];

int<lower=1> matrix_length_group_sero1neg;
int<lower=0,upper=K> lambda_sum_index_group_sero1neg[Nsero1neg,matrix_length_group_sero1neg];
int<lower=1,upper=matrix_length_group_sero1neg> lambda_length_group_sero1neg[Nsero1neg];
int<lower=0,upper=K> lambda_0_index_group_sero1neg[Nsero1neg];

int<lower=1> Nsero1pos;
int<lower=1,upper=K> lambda_index_group_sero1pos[Nsero1pos];
int<lower=1,upper=K_virus> virus_index_sero1pos[Nsero1pos];
int<lower=1> Csero1pos[Nsero1pos];

int<lower=1> matrix_length_group_sero1pos;
int<lower=0,upper=K> lambda_sum_index_group_sero1pos[Nsero1pos,matrix_length_group_sero1pos];
int<lower=1,upper=matrix_length_group_sero1pos> lambda_length_group_sero1pos[Nsero1pos];
int<lower=0,upper=K> lambda_0_index_group_sero1pos[Nsero1pos];

}
parameters{
real log_beta0;
vector[K_time_virus] log_beta1;
real<lower=0> beta2;

real<lower=0> sigma_beta;
vector<lower=0>[K_virus] log_one_over_phi; 
vector<lower=0>[K_virus] xi0;

}
transformed parameters{
vector<lower=0>[K] lambda; 
vector[K_virus] xi=xi0+1;
vector[N_group1] log_mean1;
vector[N_group0] log_mean0;

lambda = exp(log_beta0 + log_beta1[time_virus_index] - beta2*lambda_age - log(beta2)) .* ((lambda_age+1/beta2)*(1-exp(-beta2)) - exp(-beta2));

for (j in 1:N_group1){
log_mean1[j] =  - sum(lambda[lambda_sum_index_group1[j,1:lambda_length_group1[j]]]) + log_pop_group1[j] + log(xi[virus_index1[j]]) - log(lambda[lambda_0_index_group1[j]]) + lambda[lambda_0_index_group1[j]]/xi[virus_index1[j]] + log1m_exp(-lambda[lambda_0_index_group1[j]]/xi[virus_index1[j]]) + log1m_exp(-lambda[lambda_index_group1[j]]) - log_one_over_phi[virus_index1[j]];
}

for (i in 1:N_group0){
log_mean0[i] = log_pop_group0[i] + log(xi[virus_index0[i]]) - log(lambda[lambda_index0[i]]) - age_x[i]/12*lambda[lambda_index0[i]] + log(exp(-lambda[lambda_index0[i]]/20*log_sum_exp(-1.667*age_x[i],-20/xi[virus_index0[i]]))-1) + log1m_exp(-lambda[lambda_index0[i]]/12) - log_one_over_phi[virus_index0[i]];
}

}
model{
case_group11 ~ poisson_log(log_mean1);
case_group00 ~ poisson_log(log_mean0);

for (i in 1:Nsero0neg){
target += log(xi[virus_index_sero0neg[i]]) - log(lambda[lambda_index_sero0neg[i]]) - age_x_sero0neg[i]*lambda[lambda_index_sero0neg[i]] + log(exp(-lambda[lambda_index_sero0neg[i]]/20*log_sum_exp(-20*age_x_sero0neg[i],-20/xi[virus_index_sero0neg[i]]))-1);
}
for (i in 1:Nsero0pos){
target += log1m(xi[virus_index_sero0pos[i]]/lambda[lambda_index_sero0pos[i]]*exp(-age_x_sero0pos[i]*lambda[lambda_index_sero0pos[i]])*(exp(-lambda[lambda_index_sero0pos[i]]/20*log_sum_exp(-20*age_x_sero0pos[i],-20/xi[virus_index_sero0pos[i]]))-1));
}
for (j in 1:Nsero1neg){
target += Csero1neg[j]*(- sum(lambda[lambda_sum_index_group_sero1neg[j,1:lambda_length_group_sero1neg[j]]]) + log(xi[virus_index_sero1neg[j]]) - log(lambda[lambda_0_index_group_sero1neg[j]]) + lambda[lambda_0_index_group_sero1neg[j]]/xi[virus_index_sero1neg[j]] + log1m_exp(-lambda[lambda_0_index_group_sero1neg[j]]/xi[virus_index_sero1neg[j]]));
}
for (j in 1:Nsero1pos){
target += Csero1pos[j]*log1m(exp(- sum(lambda[lambda_sum_index_group_sero1pos[j,1:lambda_length_group_sero1pos[j]]]))*xi[virus_index_sero1pos[j]]/lambda[lambda_0_index_group_sero1pos[j]]*(exp(lambda[lambda_0_index_group_sero1pos[j]]/xi[virus_index_sero1pos[j]])-1));
}
log_beta0 ~ normal(0,5);
log_beta1 ~ normal(0,sigma_beta);
sigma_beta ~ cauchy(0,1);

beta2 ~ exponential(1);
log_one_over_phi ~ exponential(0.1);
xi0 ~ exponential(1);
}
"
model_linft_compiled <- rstan::stan_model(
  model_name = 'model_linft',
  model_code = gsub('\t','',model_linft_txt)
)

##### model ExpPw #####
model_exppw_txt <- "
data{
int<lower=1> K;
int<lower=1> K_time_virus;
int<lower=1> K_age;
int<lower=1,upper=K_time_virus> time_virus_index[K];
int<lower=1,upper=K_age> age_index[K];
int<lower=1> K_virus;


int<lower=1> N_group1;
int<lower=1,upper=K> lambda_index_group1[N_group1];
real<lower=0> log_pop_group1[N_group1];
int<lower=1,upper=K_virus> virus_index1[N_group1];
int<lower=0> case_group11[N_group1];

int<lower=1> matrix_length_group1;
int<lower=0,upper=K> lambda_sum_index_group1[N_group1,matrix_length_group1];
int<lower=1,upper=matrix_length_group1> lambda_length_group1[N_group1];
int<lower=0,upper=K> lambda_0_index_group1[N_group1];


int<lower=1> N_group0;
vector<lower=0>[N_group0] log_pop_group0;
vector<lower=0.1>[N_group0] age_x;
int<lower=1,upper=K_virus> virus_index0[N_group0];
int<lower=1,upper=K> lambda_index0[N_group0];
int<lower=0> case_group00[N_group0];


int<lower=1> Nsero0neg;
vector<lower=0.1>[Nsero0neg] age_x_sero0neg;
int<lower=1,upper=K_virus> virus_index_sero0neg[Nsero0neg];
int<lower=1,upper=K> lambda_index_sero0neg[Nsero0neg];

int<lower=1> Nsero0pos;
vector<lower=0.1>[Nsero0pos] age_x_sero0pos;
int<lower=1,upper=K_virus> virus_index_sero0pos[Nsero0pos];
int<lower=1,upper=K> lambda_index_sero0pos[Nsero0pos];

int<lower=1> Nsero1neg;
int<lower=1,upper=K> lambda_index_group_sero1neg[Nsero1neg];
int<lower=1,upper=K_virus> virus_index_sero1neg[Nsero1neg];
int<lower=1> Csero1neg[Nsero1neg];

int<lower=1> matrix_length_group_sero1neg;
int<lower=0,upper=K> lambda_sum_index_group_sero1neg[Nsero1neg,matrix_length_group_sero1neg];
int<lower=1,upper=matrix_length_group_sero1neg> lambda_length_group_sero1neg[Nsero1neg];
int<lower=0,upper=K> lambda_0_index_group_sero1neg[Nsero1neg];

int<lower=1> Nsero1pos;
int<lower=1,upper=K> lambda_index_group_sero1pos[Nsero1pos];
int<lower=1,upper=K_virus> virus_index_sero1pos[Nsero1pos];
int<lower=1> Csero1pos[Nsero1pos];

int<lower=1> matrix_length_group_sero1pos;
int<lower=0,upper=K> lambda_sum_index_group_sero1pos[Nsero1pos,matrix_length_group_sero1pos];
int<lower=1,upper=matrix_length_group_sero1pos> lambda_length_group_sero1pos[Nsero1pos];
int<lower=0,upper=K> lambda_0_index_group_sero1pos[Nsero1pos];

}
parameters{
real log_beta0;
vector[K_time_virus] log_beta1;
vector[K_age] log_beta2;

real<lower=0> sigma_beta1;
real<lower=0> sigma_beta2;
vector<lower=0>[K_virus] log_one_over_phi; 
vector<lower=0>[K_virus] xi;

}
transformed parameters{
vector<lower=0>[K] lambda; 
vector[N_group1] log_mean1;
vector[N_group0] log_mean0;

lambda = exp(log_beta0 + log_beta1[time_virus_index] + log_beta2[age_index]);

for (j in 1:N_group1){
log_mean1[j] =  - sum(lambda[lambda_sum_index_group1[j,1:lambda_length_group1[j]]]) + log_pop_group1[j] + log(xi[virus_index1[j]]) - log(xi[virus_index1[j]]-lambda[lambda_0_index_group1[j]]) + log1m_exp(-lambda[lambda_index_group1[j]]) - log_one_over_phi[virus_index1[j]];
}

for (i in 1:N_group0){
log_mean0[i] = log_pop_group0[i] + log(xi[virus_index0[i]]) + log(exp(-lambda[lambda_index0[i]]*age_x[i]/12)-exp(-xi[virus_index0[i]]*age_x[i]/12)) - log(xi[virus_index0[i]]-lambda[lambda_index0[i]]) + log1m_exp(-lambda[lambda_index0[i]]/12) - log_one_over_phi[virus_index0[i]];
}

}
model{
case_group11 ~ poisson_log(log_mean1);
case_group00 ~ poisson_log(log_mean0);

for (i in 1:Nsero0neg){
target += log(xi[virus_index_sero0neg[i]]) + log(exp(-lambda[lambda_index_sero0neg[i]]*age_x_sero0neg[i]) - exp(-xi[virus_index_sero0neg[i]]*age_x_sero0neg[i])) -  log( xi[virus_index_sero0neg[i]] - lambda[lambda_index_sero0neg[i]] );


}
for (i in 1:Nsero0pos){
target += log1m( xi[virus_index_sero0pos[i]]/(xi[virus_index_sero0pos[i]]-lambda[lambda_index_sero0pos[i]])*(exp(-lambda[lambda_index_sero0pos[i]]*age_x_sero0pos[i]) - exp(-xi[virus_index_sero0pos[i]]*age_x_sero0pos[i]))  );
}
for (j in 1:Nsero1neg){
target += Csero1neg[j]*(- sum(lambda[lambda_sum_index_group_sero1neg[j,1:lambda_length_group_sero1neg[j]]]) + log(xi[virus_index_sero1neg[j]]) - log(xi[virus_index_sero1neg[j]]-lambda[lambda_0_index_group_sero1neg[j]]) );
}
for (j in 1:Nsero1pos){
target += Csero1pos[j]*log1m(exp(- sum(lambda[lambda_sum_index_group_sero1pos[j,1:lambda_length_group_sero1pos[j]]]))*xi[virus_index_sero1pos[j]]/(xi[virus_index_sero1pos[j]]-lambda[lambda_0_index_group_sero1pos[j]]) );
}
log_beta0 ~ normal(0,5);
log_beta1 ~ normal(0,sigma_beta1);
log_beta2 ~ normal(0,sigma_beta2);
sigma_beta1 ~ cauchy(0,1);
sigma_beta2 ~ cauchy(0,1);

log_one_over_phi ~ exponential(0.1);
xi ~ gamma(4,2);
}
"

model_exppw_compiled <- rstan::stan_model(
  model_name = 'model_exppw',
  model_code = gsub('\t','',model_exppw_txt)
)

##### model ExpFt #####
model_expft_txt <- "
data{
int<lower=1> K;
vector<lower=0>[K] lambda_age;
int<lower=1> K_time_virus;
int<lower=1,upper=K_time_virus> time_virus_index[K];
int<lower=1> K_virus;


int<lower=1> N_group1;
int<lower=1,upper=K> lambda_index_group1[N_group1];
real<lower=0> log_pop_group1[N_group1];
int<lower=1,upper=K_virus> virus_index1[N_group1];
int<lower=0> case_group11[N_group1];

int<lower=1> matrix_length_group1;
int<lower=0,upper=K> lambda_sum_index_group1[N_group1,matrix_length_group1];
int<lower=1,upper=matrix_length_group1> lambda_length_group1[N_group1];
int<lower=0,upper=K> lambda_0_index_group1[N_group1];


int<lower=1> N_group0;
vector<lower=0>[N_group0] log_pop_group0;
vector<lower=0.1>[N_group0] age_x;
int<lower=1,upper=K_virus> virus_index0[N_group0];
int<lower=1,upper=K> lambda_index0[N_group0];
int<lower=0> case_group00[N_group0];


int<lower=1> Nsero0neg;
vector<lower=0.1>[Nsero0neg] age_x_sero0neg;
int<lower=1,upper=K_virus> virus_index_sero0neg[Nsero0neg];
int<lower=1,upper=K> lambda_index_sero0neg[Nsero0neg];

int<lower=1> Nsero0pos;
vector<lower=0.1>[Nsero0pos] age_x_sero0pos;
int<lower=1,upper=K_virus> virus_index_sero0pos[Nsero0pos];
int<lower=1,upper=K> lambda_index_sero0pos[Nsero0pos];

int<lower=1> Nsero1neg;
int<lower=1,upper=K> lambda_index_group_sero1neg[Nsero1neg];
int<lower=1,upper=K_virus> virus_index_sero1neg[Nsero1neg];
int<lower=1> Csero1neg[Nsero1neg];

int<lower=1> matrix_length_group_sero1neg;
int<lower=0,upper=K> lambda_sum_index_group_sero1neg[Nsero1neg,matrix_length_group_sero1neg];
int<lower=1,upper=matrix_length_group_sero1neg> lambda_length_group_sero1neg[Nsero1neg];
int<lower=0,upper=K> lambda_0_index_group_sero1neg[Nsero1neg];

int<lower=1> Nsero1pos;
int<lower=1,upper=K> lambda_index_group_sero1pos[Nsero1pos];
int<lower=1,upper=K_virus> virus_index_sero1pos[Nsero1pos];
int<lower=1> Csero1pos[Nsero1pos];

int<lower=1> matrix_length_group_sero1pos;
int<lower=0,upper=K> lambda_sum_index_group_sero1pos[Nsero1pos,matrix_length_group_sero1pos];
int<lower=1,upper=matrix_length_group_sero1pos> lambda_length_group_sero1pos[Nsero1pos];
int<lower=0,upper=K> lambda_0_index_group_sero1pos[Nsero1pos];

}
parameters{
real log_beta0;
vector[K_time_virus] log_beta1;
real<lower=0> beta2;

real<lower=0> sigma_beta;
vector<lower=0>[K_virus] log_one_over_phi; 
vector<lower=0>[K_virus] xi;

}
transformed parameters{
vector<lower=0>[K] lambda; 
vector[N_group1] log_mean1;
vector[N_group0] log_mean0;

lambda = exp(log_beta0 + log_beta1[time_virus_index] - beta2*lambda_age - log(beta2)) .* ((lambda_age+1/beta2)*(1-exp(-beta2)) - exp(-beta2));

for (j in 1:N_group1){
log_mean1[j] =  - sum(lambda[lambda_sum_index_group1[j,1:lambda_length_group1[j]]]) + log_pop_group1[j] + log(xi[virus_index1[j]]) - log(xi[virus_index1[j]]-lambda[lambda_0_index_group1[j]]) + log1m_exp(-lambda[lambda_index_group1[j]]) - log_one_over_phi[virus_index1[j]];
}

for (i in 1:N_group0){
log_mean0[i] = log_pop_group0[i] + log(xi[virus_index0[i]]) + log(exp(-lambda[lambda_index0[i]]*age_x[i]/12)-exp(-xi[virus_index0[i]]*age_x[i]/12)) - log(xi[virus_index0[i]]-lambda[lambda_index0[i]]) + log1m_exp(-lambda[lambda_index0[i]]/12) - log_one_over_phi[virus_index0[i]];
}

}
model{
case_group11 ~ poisson_log(log_mean1);
case_group00 ~ poisson_log(log_mean0);

for (i in 1:Nsero0neg){
target += log(xi[virus_index_sero0neg[i]]) + log(exp(-lambda[lambda_index_sero0neg[i]]*age_x_sero0neg[i]) - exp(-xi[virus_index_sero0neg[i]]*age_x_sero0neg[i])) -  log( xi[virus_index_sero0neg[i]] - lambda[lambda_index_sero0neg[i]] );
}
for (i in 1:Nsero0pos){
target += log1m( xi[virus_index_sero0pos[i]]/(xi[virus_index_sero0pos[i]]-lambda[lambda_index_sero0pos[i]])*(exp(-lambda[lambda_index_sero0pos[i]]*age_x_sero0pos[i]) - exp(-xi[virus_index_sero0pos[i]]*age_x_sero0pos[i]))  );
}
for (j in 1:Nsero1neg){
target += Csero1neg[j]*(- sum(lambda[lambda_sum_index_group_sero1neg[j,1:lambda_length_group_sero1neg[j]]]) + log(xi[virus_index_sero1neg[j]]) - log(xi[virus_index_sero1neg[j]]-lambda[lambda_0_index_group_sero1neg[j]]) );
}
for (j in 1:Nsero1pos){
target += Csero1pos[j]*log1m(exp(- sum(lambda[lambda_sum_index_group_sero1pos[j,1:lambda_length_group_sero1pos[j]]]))*xi[virus_index_sero1pos[j]]/(xi[virus_index_sero1pos[j]]-lambda[lambda_0_index_group_sero1pos[j]]) );
}
log_beta0 ~ normal(0,5);
log_beta1 ~ normal(0,sigma_beta);
sigma_beta ~ cauchy(0,1);

beta2 ~ exponential(1);

log_one_over_phi ~ exponential(0.1);
xi ~ gamma(4,2);
}
"

model_expft_compiled <- rstan::stan_model(
  model_name = 'model_expft',
  model_code = gsub('\t','',model_expft_txt)
)

#### model input: stan data ####
##### case data #####
## when age=0
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

# for those with no recording of cases, set the case number = 0
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

stan_data <- list()
stan_data$N_group0 <- dim(case6_0)[1]
stan_data$log_pop_group0 <- case6_0[,log(pop)]
stan_data$age_x <- case6_0$month
stan_data$virus_index0 <- case6_0$virus_index
stan_data$lambda_index0 <- case6_0$index
stan_data$case_group00 <- case6_0$count

## when age>0
case1_1 <- case1[which(case1[,age>0]),c(2,5,6,7,10,13)]
case1_1[,age_group:=age]
set(case1_1,which(case1_1[,age_group>=7]),"age_group",7)
case2_1 <- case1_1[,list(count=.N),by=c("virus","year_sampling","age_group")]
colnames(case2_1)[2]="year_group"

# for those with no recording of cases, set the case number = 0
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


matrix_length_group1 = dim(data1)[1]

lambda_sum_index_group1 = matrix(0,nrow = dim(case6_1)[1],ncol=30)
lambda_length_group1 <- rep(0,dim(case6_1)[1])


lambda_length_group1[1] <- dim(data1)[1]
lambda_sum_index_group1[i,1:dim(data1)[1]] = data1$index


for (i in 2:dim(case6_1)[1]){
  virus1 <- case6_1[i,]$virus
  age1 <- case6_1[i,]$age_group
  year1 <- case6_1[i,]$year_group
  
  if (age1==7) age1=9
  
  data1 <- data.table(year_group = (year1-age1):(year1-1),age_group = 0:(age1-1),virus=virus1)
  
  stopifnot(max(data1$age_group)==age1-1)
  set(data1,which(data1[,age_group>7]),"age_group",7)
  set(data1,which(data1[,year_group<2011]),"year_group",2011)
  
  data1[,order1 := 1:dim(data1)[1]]
  data1 = merge(data1,lambda_index_group,by=c("virus","age_group","year_group"))
  data1 <- data1[order(data1[,order1]),]
  
  matrix_length_group2 = dim(data1)[1]
  if (matrix_length_group2>matrix_length_group1){
    matrix_length_group1 <- matrix_length_group2
  }
  
  
  lambda_length_group1[i] <- dim(data1)[1]
  lambda_sum_index_group1[i,1:dim(data1)[1]] = data1$index
  
  
}

lambda_sum_index_group1 <- lambda_sum_index_group1[,1:9]

stan_data$matrix_length_group1 <- matrix_length_group1
stan_data$lambda_sum_index_group1 <- lambda_sum_index_group1
stan_data$lambda_length_group1 <- lambda_length_group1


stan_data$N_group1 <- dim(case6_1)[1]
stan_data$lambda_index_group1 <- case6_1$index
stan_data$log_pop_group1 <- case6_1[,log(pop)]
stan_data$virus_index1 <- case6_1$virus_index
stan_data$case_group11 <- case6_1$count


zero_tmp <- case6_1[,c(1:3,5)]
set(zero_tmp,which(zero_tmp[,age_group==7]),"age_group",9)
zero_tmp[,year:=year_group-age_group]
zero_tmp <- zero_tmp[,c(3,4,5)]
colnames(zero_tmp)[3]="year_group"
zero_tmp[,age_group:=0]

set(zero_tmp,which(zero_tmp[,year_group<2011]),"year_group",2011)
zero_tmp <- merge(zero_tmp,lambda_index_group,by=c("year_group","age_group","virus"))
zero_tmp <- zero_tmp[order(zero_tmp[,self_index]),]

stan_data$lambda_0_index_group1 <- zero_tmp$index

stan_data$K <- max(lambda_index_group$index)
stan_data$lambda_age <- lambda_index_group$age_group
stan_data$K_time_virus <- max(lambda_index_group$time_virus_index)
stan_data$time_virus_index <- lambda_index_group$time_virus_index
stan_data$K_virus <- max(virus_index$virus_index)

##### serological data #####
sero1[,age:=floor(AGE_MAX)]
sero1[,birth_year:=YEAR-age]
sero1[which(sero1[,age==0]),]
sero1 <- sero1[,c(7,9,10,11,12,13,14,1,3)]
colnames(sero1)[1]="year_end"
colnames(sero1)[2]="CV-A6"
colnames(sero1)[3]="CV-A10"
colnames(sero1)[4]="CV-A16"
colnames(sero1)[5]="EV-A71"

sero2[,age:=floor(AGE_MAX)]
sero2 <- sero2[,c(3,9,11:15,2,5)]
colnames(sero2)[1]="birth_year"
colnames(sero2)[2]="year_end"

sero3 <- rbind(sero1,sero2)
sero4 <- melt(sero3,measure.vars = c("CV-A6","CV-A10","CV-A16","EV-A71"))
colnames(sero4)[6]="virus"

sero4[,pos1:=0]
set(sero4,which(sero4[,value>=8]),"pos1",1)
sero5 <- sero4[,list(count=.N),by=c("year_end","age","birth_year","virus","pos1")]

sero5_0_neg <- sero4[which(sero4[,(age==0)&(pos1==0)]),]
colnames(sero5_0_neg)[1] = "year_group"
colnames(sero5_0_neg)[2] = "age_group"
sero5_0_neg <- merge(sero5_0_neg,lambda_index_group,by=c("year_group","age_group","virus"))
sero5_0_neg <- merge(sero5_0_neg,virus_index,by="virus")

stan_data$Nsero0neg <- dim(sero5_0_neg)[1]
stan_data$age_x_sero0neg <- sero5_0_neg$AGE_MAX
stan_data$virus_index_sero0neg <- sero5_0_neg$virus_index
stan_data$lambda_index_sero0neg <- sero5_0_neg$index

sero5_0_pos <- sero4[which(sero4[,(age==0)&(pos1==1)]),]
colnames(sero5_0_pos)[1] = "year_group"
colnames(sero5_0_pos)[2] = "age_group"
sero5_0_pos <- merge(sero5_0_pos,lambda_index_group,by=c("year_group","age_group","virus"))

virus_index
sero5_0_pos <- merge(sero5_0_pos,virus_index,by="virus")

stan_data$Nsero0pos <- dim(sero5_0_pos)[1]
stan_data$age_x_sero0pos <- sero5_0_pos$AGE_MAX
stan_data$virus_index_sero0pos <- sero5_0_pos$virus_index
stan_data$lambda_index_sero0pos <- sero5_0_pos$index



sero5_1_neg <- sero5[which(sero5[,(age>0)&(pos1==0)]),]
colnames(sero5_1_neg)[1] = "year_group"
colnames(sero5_1_neg)[2] = "age_group"
sero5_1_neg <- merge(sero5_1_neg,lambda_index_group,by=c("year_group","age_group","virus"))

virus_index
sero5_1_neg <- merge(sero5_1_neg,virus_index,by="virus")
sero5_1_neg[,order2:=1:dim(sero5_1_neg)[1]]

lambda_sum_index_group_sero1neg <- matrix(0,nrow = dim(sero5_1_neg)[1],ncol=20)
lambda_length_group_sero1neg <- rep(0, dim(sero5_1_neg)[1])

i=1
virus1 <- sero5_1_neg[i,]$virus
year_end1 <- sero5_1_neg[i,]$year_group
year_start1 <- sero5_1_neg[i,]$birth_year
age1 <- sero5_1_neg[i,]$age_group
tmp1 <- data.table(age_group=0:age1,year_group=year_start1:year_end1)
set(tmp1,which(tmp1[,year_group<2011]),"year_group",2011)
tmp1$virus = virus1
tmp1 <- merge(tmp1,lambda_index_group,by=c("year_group","age_group","virus"))
matrix_length_group_sero1neg1 <- dim(tmp1)[1]
lambda_sum_index_group_sero1neg[i,1:dim(tmp1)[1]]=tmp1$index
lambda_length_group_sero1neg[i]=dim(tmp1)[1]

for (i in 2:dim(sero5_1_neg)[1]){
  virus1 <- sero5_1_neg[i,]$virus
  year_end1 <- sero5_1_neg[i,]$year_group
  year_start1 <- sero5_1_neg[i,]$birth_year
  age1 <- sero5_1_neg[i,]$age_group
  tmp1 <- data.table(age_group=0:age1,year_group=year_start1:year_end1)
  set(tmp1,which(tmp1[,year_group<2011]),"year_group",2011)
  tmp1$virus = virus1
  tmp1 <- merge(tmp1,lambda_index_group,by=c("year_group","age_group","virus"))
  
  if (dim(tmp1)[1] > matrix_length_group_sero1neg1){
    matrix_length_group_sero1neg1 = dim(tmp1)[1]
  }
  
  lambda_sum_index_group_sero1neg[i,1:dim(tmp1)[1]]=tmp1$index
  lambda_length_group_sero1neg[i]=dim(tmp1)[1]
}

lambda_sum_index_group_sero1neg <- lambda_sum_index_group_sero1neg[,1:matrix_length_group_sero1neg1]

tmp0 <- sero5_1_neg[,c(1,4,11)]
colnames(tmp0)[2] = "year_group"
tmp0[,age_group:=0]
set(tmp0,which(tmp0[,year_group<2011]),"year_group",2011)
tmp0 <- merge(tmp0,lambda_index_group,by=c("year_group","age_group","virus"))
tmp0 <- tmp0[order(tmp0[,order2]),]

stan_data$Nsero1neg <- dim(sero5_1_neg)[1]
stan_data$lambda_index_group_sero1neg <-sero5_1_neg$index
stan_data$virus_index_sero1neg <- sero5_1_neg$virus_index
stan_data$Csero1neg <- sero5_1_neg$count

stan_data$matrix_length_group_sero1neg <- matrix_length_group_sero1neg1
stan_data$lambda_sum_index_group_sero1neg <- lambda_sum_index_group_sero1neg
stan_data$lambda_length_group_sero1neg <- lambda_length_group_sero1neg
stan_data$lambda_0_index_group_sero1neg <- tmp0$index


sero5_1_pos <- sero5[which(sero5[,(age>0)&(pos1==1)]),]
colnames(sero5_1_pos)[1] = "year_group"
colnames(sero5_1_pos)[2] = "age_group"
sero5_1_pos <- merge(sero5_1_pos,lambda_index_group,by=c("year_group","age_group","virus"))

virus_index
sero5_1_pos <- merge(sero5_1_pos,virus_index,by="virus")
sero5_1_pos[,order2:=1:dim(sero5_1_pos)[1]]

lambda_sum_index_group_sero1pos <- matrix(0,nrow = dim(sero5_1_pos)[1],ncol=20)
lambda_length_group_sero1pos <- rep(0, dim(sero5_1_pos)[1])

i=1
virus1 <- sero5_1_pos[i,]$virus
year_end1 <- sero5_1_pos[i,]$year_group
year_start1 <- sero5_1_pos[i,]$birth_year
age1 <- sero5_1_pos[i,]$age_group
tmp1 <- data.table(age_group=0:age1,year_group=year_start1:year_end1)
set(tmp1,which(tmp1[,year_group<2011]),"year_group",2011)
tmp1$virus = virus1
tmp1 <- merge(tmp1,lambda_index_group,by=c("year_group","age_group","virus"))
matrix_length_group_sero1pos1 <- dim(tmp1)[1]
lambda_sum_index_group_sero1pos[i,1:dim(tmp1)[1]]=tmp1$index
lambda_length_group_sero1pos[i]=dim(tmp1)[1]

for (i in 2:dim(sero5_1_pos)[1]){
  virus1 <- sero5_1_pos[i,]$virus
  year_end1 <- sero5_1_pos[i,]$year_group
  year_start1 <- sero5_1_pos[i,]$birth_year
  age1 <- sero5_1_pos[i,]$age_group
  tmp1 <- data.table(age_group=0:age1,year_group=year_start1:year_end1)
  set(tmp1,which(tmp1[,year_group<2011]),"year_group",2011)
  tmp1$virus = virus1
  tmp1 <- merge(tmp1,lambda_index_group,by=c("year_group","age_group","virus"))
  
  if (dim(tmp1)[1] > matrix_length_group_sero1pos1){
    matrix_length_group_sero1pos1 = dim(tmp1)[1]
  }
  
  lambda_sum_index_group_sero1pos[i,1:dim(tmp1)[1]]=tmp1$index
  lambda_length_group_sero1pos[i]=dim(tmp1)[1]
}

lambda_sum_index_group_sero1pos <- lambda_sum_index_group_sero1pos[,1:matrix_length_group_sero1pos1]

tmp0 <- sero5_1_pos[,c(1,4,11)]
colnames(tmp0)[2] = "year_group"
tmp0[,age_group:=0]
set(tmp0,which(tmp0[,year_group<2011]),"year_group",2011)
tmp0 <- merge(tmp0,lambda_index_group,by=c("year_group","age_group","virus"))
tmp0 <- tmp0[order(tmp0[,order2]),]

stan_data$Nsero1pos <- dim(sero5_1_pos)[1]
stan_data$lambda_index_group_sero1pos <- sero5_1_pos$index
stan_data$virus_index_sero1pos <- sero5_1_pos$virus_index
stan_data$Csero1pos <- sero5_1_pos$count

stan_data$matrix_length_group_sero1pos <- matrix_length_group_sero1pos1
stan_data$lambda_sum_index_group_sero1pos <- lambda_sum_index_group_sero1pos
stan_data$lambda_length_group_sero1pos <- lambda_length_group_sero1pos
stan_data$lambda_0_index_group_sero1pos <- tmp0$index



stan_data1 <- stan_data
stan_data1 <- stan_data1[-17]
stan_data1$K_age <- max(lambda_index_group$age_index)
stan_data1$age_index <- lambda_index_group$age_index

saveRDS(stan_data1,"stan_data_pw.rds")
saveRDS(stan_data,"stan_data_ft.rds")

#### model fitting ####
stan_data_pw <- readRDS("stan_data_pw.rds")
stan_data_ft <- readRDS("stan_data_ft.rds")

model_linpw_fit <- rstan::sampling(model_linpw_compiled,
                                   data=stan_data_pw,
                                   chains=4,
                                   init=list(list(log_beta0=-5,log_beta1=rep(-1,stan_data_pw$K_time_virus),log_beta2=rep(0,stan_data_pw$K_age),xi0=rep(0.01,stan_data_pw$K_virus)),
                                             list(log_beta0=-3,log_beta1=rep(-2,stan_data_pw$K_time_virus),log_beta2=rep(-0.5,stan_data_pw$K_age),xi0=rep(0.1,stan_data_pw$K_virus)),
                                             list(log_beta0=-1,log_beta1=rep(-3,stan_data_pw$K_time_virus),log_beta2=rep(0.5,stan_data_pw$K_age),xi0=rep(0.2,stan_data_pw$K_virus)),
                                             list(log_beta0=-2,log_beta1=rep(-2.5,stan_data_pw$K_time_virus),log_beta2=rep(-1,stan_data_pw$K_age),xi0=rep(0.3,stan_data_pw$K_virus))),
                                   warmup=1e3,
                                   iter=1e4,
                                   cores=4)

model_linft_fit <- rstan::sampling(model_linft_compiled,
                                   data=stan_data_ft,
                                   chains=4,
                                   warmup=1e3,
                                   iter=1e4,
                                   cores=4)

model_exppw_fit <- rstan::sampling(model_exppw_compiled,
                                   data=stan_data_pw,
                                   chains=4,
                                   init=list(list(log_beta0=-5,log_beta1=rep(-1,stan_data_pw$K_time_virus),log_beta2=rep(0,stan_data_pw$K_age)),
                                             list(log_beta0=-3,log_beta1=rep(-2,stan_data_pw$K_time_virus),log_beta2=rep(-0.5,stan_data_pw$K_age)),
                                             list(log_beta0=-1,log_beta1=rep(-3,stan_data_pw$K_time_virus),log_beta2=rep(0.5,stan_data_pw$K_age)),
                                             list(log_beta0=-2,log_beta1=rep(-2.5,stan_data_pw$K_time_virus),log_beta2=rep(-1,stan_data_pw$K_age))),
                                   warmup=1e3,
                                   iter=1e4,
                                   cores=4)


model_expft_fit <- rstan::sampling(model_expft_compiled,
                                   data=stan_data_ft,
                                   chains=4,
                                   warmup=1e3,
                                   iter=1e4,
                                   cores=4)


saveRDS(model_linpw_fit,"model_linpw_fit.rds")
saveRDS(model_linft_fit,"model_linft_fit.rds")
saveRDS(model_exppw_fit,"model_exppw_fit.rds")
saveRDS(model_expft_fit,"model_expft_fit.rds")

