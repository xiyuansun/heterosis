library(stringr)


#' A function to create a STAN model for either
#' differential expression or heterosis with appropriate priors on sparse signals. 
#'
#' @param type a string that determines whether a differential expression or heterosis model is constructed
#' @param average a string that determines the prior on parental average
#' @param difference a string that determines the prior on parental difference
#' @param hybrid a string that determines the prior on hybrid difference from the parental average (heterosis model only)
#' @return a string that provides a STAN model 

construct_model = function(type=c("differential expression","heterosis"),
                           average=c("normal","t"),
                           difference=c("normal","t"),
                           hybrid=c("normal","t"))
{
  isHeterosis   = pmatch(type,       c("differential expression","heterosis"))==2
  isAverageT    = pmatch(average,    c("normal","t"))==2
  isDifferenceT = pmatch(difference, c("normal","t"))==2
  isHybridT     = pmatch(hybrid,     c("normal","t"))==2

  data = 
  str_join("data {\n",
           "  int<lower=1> N;\n",
           "  int<lower=1> G;\n",
           "  int<lower=1> S;\n",
           "  int<lower=0> y[N];\n",
           "  int<lower=1,upper=G> gene[N];\n",
           "  int<lower=1,upper=S> sample[N];\n",
           "  int<lower=-1,upper=1> parent[N];\n",
           ifelse(isHeterosis  ,"  int<lower=0,upper=1> hybrid[N];\n",""),
           ifelse(isAverageT   ,"  real<lower=0> v_mean;\n",""),
           ifelse(isDifferenceT,"  real<lower=0> v_diff;\n",""),
           ifelse(isHybridT    ,"  real<lower=0> v_diffh;\n",""),
           "}\n")

  parameters = 
  str_join("parameters {\n",
           "  real e[N];\n",
           "  real mean[G];\n",
           "  real diff[G];\n",
           ifelse(isHeterosis,"  real diffh[G];\n",""),
           "  real<lower=0> sigma_od[G];\n", 
           "  real norm[S];\n",
           "  real theta_mean;\n",
           "  real theta_diff;\n",
           ifelse(isHeterosis,"  real theta_diff;\n",""),
           "  real<lower=0> sigma_mean;\n",
           "  real<lower=0> sigma_diff;\n",
           ifelse(isHeterosis,"  real<lower=0> sigma_diffh;\n",""),
           "  real<lower=0> a_od;\n",
           "  real<lower=0> b_od;\n",
           "}\n")

  transformed_parameters = 
  str_join("transformed parameters {\n",
           "  for (i in 1:N){\n",
           "    lambda[i] <- norm[sample[i]] + mean[gene[i]] + ",

           switch(pmatch(type, c("differential expression","heterosis")),
             "parent[i]*diff[gene[i]]",
             "parent[i]*diff[gene[i]] + hybrid[i]*diffh[gene[i]]"), 

           " + e[i]\n}\n")

  model = 
  str_join("model {\n  for (i in 1:N) {\n    y[i] ~ poison(exp(lambda[i]))\n",
           "    e[i] ~ normal(0, sigma_od[gene[i]])\n  }\n  for (g in 1:G) {\n",

           switch(pmatch(average, c("normal","t")),
             "    mean[g] ~ normal(theta_mean, sigma_mean)\n",
             "    mean[g] ~ student_t(v_mean, theta_mean, sigma_mean)\n"),

           switch(pmatch(difference, c("normal","t")),
             "    diff[g] ~ normal(theta_diff, sigma_diff)\n",
             "    diff[g] ~ student_t(v_diff, theta_diff, sigma_diff)\n"),

           ifelse(isHeterosis, 
             switch(pmatch(hybrid, c("normal","t")),
               "    diffh[g] ~ normal(theta_diffh, sigma_diffh)\n",
               "    diffh[g] ~ student_t(v_diffh, theta_diffh, sigma_diffh)\n"), ""),

           "    sigma_od[g] ~ gamma(a_od, b_od)\n  }\n  for(s in 1:S) {\n    norm[s] ~ normal(theta_norm, sigma_norm)\n  }\n",
           "  theta_mean ~ normal(0,1)\n  theta_diff ~ normal(0,1)\n",

           ifelse(isHeterosis, "  theta_diffh ~ normal(0,1)\n",""),

           "  sigma_mean ~ uniform(0,1)\n  sigma_diff ~ uniform(0,1)\n",

           ifelse(isHeterosis, "  sigma_diffh ~ uniform(0,1)\n",""),

           "}\n"
          )

  return(str_join(data,parameters,transformed_parameters,model))
}


#' Creates a simulated data set for either differential expression or heterosis data.
#'
#' @param type a string that determines whether a differential expression or heterosis model is constructed
#' @param average a string that determines the prior on parental average
#' @param difference a string that determines the prior on parental difference
#' @param hybrid a string that determines the prior on hybrid difference from the parental average (heterosis model only)
#' @param params a list with all the parameter values with required elements 
#' @return a string that provides a STAN model 

simulate_data = function(type=c("differential expression","heterosis"),
                           average=c("normal","t"),
                           difference=c("normal","t"),
                           hybrid=c("normal","t"),
                           params=NULL)
{
  attach(params)
  
}

