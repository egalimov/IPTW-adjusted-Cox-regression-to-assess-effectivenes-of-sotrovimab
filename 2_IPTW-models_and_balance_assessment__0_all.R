library(twang)
library(cobalt)
library(xtable)
library(survey)
library(stringr)
library(WeightIt)
library(stats)

#load data and set var types 
t1 = read.csv('S://Evgeny/1_Sotrovimab/221209_target3_CE.csv')
t1$Gender = as.factor(t1$Gender)
t1$Time_period = as.factor(t1$Time_period)
ordered(t1$Time_period, levels = c('1','2','3') )
t1$RenalDisease = as.factor(t1$RenalDisease)
t1$Mult_high_risk_bin = as.factor(t1$Mult_high_risk_bin)
t1$Moderate_risk_bin = as.factor(t1$Moderate_risk_bin)
t1$SolidOrganTransplants = as.factor(t1$SolidOrganTransplants)
t1$Vacc_status_ordered = as.factor(t1$Vacc_status_ordered)
ordered(t1$Vacc_status_ordered, levels = c('0','1','2') )
t1$Vacc_stcatus_bin = as.factor(t1$Vacc_stcatus_bin)
t1$EthnicCategory = as.factor(t1$EthnicCategory)





#### function to get weights from propensity score models
  ### the function output an updated df with added column with weights
  # data - dataset name
  # treat_col_name - col name for treatment: ie 'Treatment'
  # out_weights_col_name - name of column with output weights: ie 'w_GBM_stab'
  # model_m: specify the model - 1) GBM, 2) logr
  # formula_m - formula with outcome and factors, ie y~x1+x2
  # GBM_stop.method: 
          # if model == GBM, the options are: 1) 'es.stat.mean', 2) 'ks.stat.max', 3) 'es.stat.max', 4) 'ks.stat.mean'
          # if model != GBM. set 'none'
  # estimand: 1) 'ATE', 2) 'ATT'
  # weights_stab: 1) raw, 2) stabilized
  # trimming: TRUE or FALSE 
  # trim_values = c(0.01, 0.99)  #set lower and upper quintiles to trim: ie c(0.01, 0.99)
get_iptw_weights = function(data, out_weights_col_name, treat_col_name, model_m, formula_m, GBM_stop.method,estimand_m , weights_stab, trimming, trim_values){
  if (model_m=='GBM'){
    ps_w = ps( formula_m,
               data = data,
               estimand = estimand_m,
               stop.method = stop.methods[c(GBM_stop.method)],
               n.trees = 10000,
               interactio.depth = 2,
               shrinkage = 0.01,
               perm.test.iters = 0,
               verbose = FALSE)
    data$ps = ps_w$ps[[1]]
        if (weights_stab == 'raw'){ data[out_weights_col_name] = get_w_from_ps(data$ps, data[,c(treat_col_name)], estimand = estimand_m, stabilize = FALSE)     }
        if (weights_stab == 'stabilized'){ data[out_weights_col_name] = get_w_from_ps(data$ps, data[,c(treat_col_name)], estimand = estimand_m, stabilize = TRUE)       }
  }
  if (model_m=='logr'){
    ps_log = glm(  formula_m  ,  data = data, family = binomial) 
    data$ps = predict(ps_log, type = 'response')
    
        if (weights_stab == 'raw'){ data[out_weights_col_name] = get_w_from_ps(data$ps, data[,c(treat_col_name)], estimand = estimand_m, stabilize = FALSE)     }
        if (weights_stab == 'stabilized'){ data[out_weights_col_name] = get_w_from_ps(data$ps, data[,c(treat_col_name)], estimand = estimand_m, stabilize = TRUE)       }
  }
  
  # trimming weight
  if (trimming) { 
     q_low = quantile(  unlist( as.vector( data[out_weights_col_name] ) ), trim_values[1] )[1]
     q_high = quantile(  unlist( as.vector( data[out_weights_col_name] ) ), trim_values[2] )[1]
     #data2 = data
     data = data[ data[c(out_weights_col_name)]>q_low,] 
     data = data[ data[c(out_weights_col_name)]<q_high,] 
  }
  return (data)
}



##################### function removes p-values repeated for other categories in case of categorical variables
remove_wrong_p = function(a){
  prev_str = '1'
  for (i in 1:nrow(a))
  {
    #i=3
    cur_str = row.names(a)[i]
    
    if (unlist(gregexpr(':',cur_str)) == -1){next}
    
    if (unlist(gregexpr(':',cur_str)) != -1){
      if (str_sub(cur_str, 1, unlist(gregexpr(':',cur_str))-1) == prev_str){
        a$ks.pval[i] = NA
        prev_str = str_sub(cur_str, 1, unlist(gregexpr(':',cur_str))-1)
      }
      prev_str = str_sub(cur_str, 1, unlist(gregexpr(':',cur_str))-1)
      
      
    }
  }     
  return(a)
}   



############################### Function to get balanced tables: saves to types of tables:
    # obtained through dx.wts method (twang)
    # obtained through cobalt package
# confounders = c("Age","Gender","EthnicCategory","Time_period","RenalDisease","Mult_high_risk_bin","Moderate_risk_bin","SolidOrganTransplants","Vacc_stcatus_bin","Days_SinceLastVaccination")
# # weights_col: 'w'
# weights_col = 'w_GBM_stab' 
# path = 'S://Evgeny/1_Sotrovimab/PS/1116/'
# name_m = '1_'
# estimand_m = 'ATE'
# treat_col_name = 'Treatment'
# #interactions_m: do we include interactions between confounders - False or True
# interactions_m = FALSE

####
get_balance_tables = function(data, weights_col,interactions_m, confounders, estimand_m,treat_col_name, path, name_m ){
  # create dir
  dir.create(path)
  
  # dx.wts method
  bal_obj = dx.wts(data[,c(weights_col)], data=data,estimand = estimand_m, vars = confounders, 
                   treat.var = treat_col_name, perm.test.iters = 0)
  
  bal_obj$desc[[1]]$bal.tab$results = remove_wrong_p(bal_obj$desc[[1]]$bal.tab$results)
  bal_obj$desc[[2]]$bal.tab$results = remove_wrong_p(bal_obj$desc[[2]]$bal.tab$results)
  
  write.csv(bal_obj$desc[[1]]$bal.tab, paste(path,name_m,'unw_bal.tab.csv', sep = '')   )
  write.csv(bal_obj$desc[[2]]$bal.tab, paste(path,name_m,'weighted_bal.tab.csv', sep = '')   )
  
  # cobalt method
  bal_cob_obj = bal.tab(data[,confounders], treat = data[,c(treat_col_name)], weights=data[,c(weights_col)], thresholds = c(m=.1, v=2), binary='std', int=TRUE)
  sink(file = paste(path,name_m,'cobalt_bal.tab.txt', sep = ''))
  print( bal_cob_obj )
  sink()
}

###############################   Function to plot diagnostic plots from twang
    # 1) boxplots of ps
    # 2) standardized differences betwen treated and unterated
    # 3) p values for t and chi-square tests
    # 4) p values for KS tests
# path = 'S://Evgeny/1_Sotrovimab/PS/1116/'
# name_m = '1_'
# # weights_col: 'w'
# weights_col = 'w_GBM_stab' 
# path = 'S://Evgeny/1_Sotrovimab/PS/1116/'
# name_m = '1_'
# estimand_m = 'ATE'
# treat_col_name = 'Treatment'

plot_diagnostics = function(data, weights_col, estimand_m, confounders, treat_col_name, path, name_m){
  # get the dx.wts object
    bal_obj = dx.wts(data[,c(weights_col)], data=data,estimand = estimand_m, vars = confounders, 
                   treat.var = treat_col_name, perm.test.iters = 0)
    # remove repeated p-values
    bal_obj$desc[[1]]$bal.tab$results = remove_wrong_p(bal_obj$desc[[1]]$bal.tab$results)
    bal_obj$desc[[2]]$bal.tab$results = remove_wrong_p(bal_obj$desc[[2]]$bal.tab$results)
    # get formula for boxplots
    formula_boxt_plot = formula(paste('ps~', treat_col_name, sep='') )
  # boxplots for ps
    boxplot(formula_boxt_plot, data, horizontal = TRUE)
    p2 = recordPlot()
    png(paste(path,name_m,'_1_ps_boxplot.png', sep = ''), width = 12, height = 8, units = 'cm', res = 200  )
    print(p2)
    dev.off() 
  # Plots of standardized effect size of difference between means for confounding variables before and after weighting
    class(bal_obj) = 'ps'
    p3 = plot(bal_obj, plots = 3)
    png(paste(path,name_m,'_2_Plots of standardized diff.png', sep = ''), width = 12, height = 8, units = 'cm', res = 200  )
    print(p3)
    dev.off() 
  # plots of p-values from t test comparing means for confounding variables before and after weighting
    p4 = plot(bal_obj, plots = 4) 
    png(paste(path,name_m,'_3_p-values from t or chi-square tests.png', sep = ''), width = 12, height = 8, units = 'cm', res = 200  )
    print(p4)
    dev.off() 
  # plots of p-values from KS test comparing means for confounding variables before and after weighting
    p5 = plot(bal_obj, plots = 5) 
    png(paste(path,name_m,'_4_p-values from KS tests.png', sep = ''), width = 12, height = 8, units = 'cm', res = 200  )
    print(p5)
    dev.off() 
}




###################################################################
########Running the analysis

data = t1
path = 'S://Evgeny/1_Sotrovimab/PS/1209_0_all/'
formula_m = formula('Treatment ~ Age+Gender+EthnicCategory+Time_period+RenalDisease+Mult_high_risk_bin+Moderate_risk_bin+SolidOrganTransplants+Vacc_stcatus_bin+Days_SinceLastVaccination')
treat_col_name = 'Treatment'
confounders = c("Age","Gender","EthnicCategory","Time_period","RenalDisease","Mult_high_risk_bin","Moderate_risk_bin","SolidOrganTransplants","Vacc_stcatus_bin","Days_SinceLastVaccination")
estimand_m = 'ATE'

dir.create(path)
# get dataset with weights
data_GBM_es.mean_raw = get_iptw_weights(data, 'w_GBM_es.mean_raw', 'Treatment', 'GBM',formula_m,'es.stat.mean','ATE','raw',FALSE, c(0.01,0.99)  )
data_GBM_es.mean_raw_trim = get_iptw_weights(data, 'w_GBM_es.mean_raw', 'Treatment', 'GBM',formula_m,'es.stat.mean','ATE','raw',TRUE, c(0.01,0.99)  )
data_GBM_es.mean_stab = get_iptw_weights(data, 'w_GBM_es.mean_stab', 'Treatment', 'GBM',formula_m,'es.stat.mean','ATE','stabilized',FALSE,c(0.01,0.99)  )
data_GBM_es.mean_stab_trim = get_iptw_weights(data, 'w_GBM_es.mean_stab', 'Treatment', 'GBM',formula_m,'es.stat.mean','ATE','stabilized',TRUE,c(0.01,0.99)  )
data_logr_raw = get_iptw_weights(data, 'w_logr_raw', 'Treatment', 'logr',formula_m,'es.stat.mean','ATE','raw',FALSE,c(0.01,0.99)  )
data_logr_raw_trim = get_iptw_weights(data, 'w_logr_raw', 'Treatment', 'logr',formula_m,'es.stat.mean','ATE','raw',TRUE,c(0.01,0.99)  )
data_logr_stab = get_iptw_weights(data, 'w_logr_stab', 'Treatment', 'logr',formula_m,'es.stat.mean','ATE','stabilized',FALSE,c(0.01,0.99)  )
data_logr_stab_trim = get_iptw_weights(data, 'w_logr_stab', 'Treatment', 'logr',formula_m,'es.stat.mean','ATE','stabilized',TRUE,c(0.01,0.99)  )

write.csv(data_GBM_es.mean_raw, paste(path, '1_GBM_raw_', 'data.csv', sep=''))
write.csv(data_GBM_es.mean_raw_trim, paste(path, '2_GBM_raw_trim_', 'data.csv', sep=''))
write.csv(data_GBM_es.mean_stab, paste(path, '3_GBM_stab_', 'data.csv', sep=''))
write.csv(data_GBM_es.mean_stab_trim, paste(path, '4_GBM_stab_trim_', 'data.csv', sep=''))
write.csv(data_logr_raw, paste(path, '5_logr_raw_', 'data.csv', sep=''))
write.csv(data_logr_raw_trim, paste(path, '6_logr_raw_trim_', 'data.csv', sep=''))
write.csv(data_logr_stab, paste(path, '7_logr_stab_', 'data.csv', sep=''))
write.csv(data_logr_stab_trim, paste(path, '8_logr_stab_trim_', 'data.csv', sep=''))



# get balance tables
get_balance_tables(data_GBM_es.mean_raw, 'w_GBM_es.mean_raw', FALSE, confounders, estimand_m,treat_col_name, path, '1_GBM_raw_' )
get_balance_tables(data_GBM_es.mean_raw_trim, 'w_GBM_es.mean_raw', FALSE, confounders, estimand_m,treat_col_name, path, '2_GBM_raw_trim_' )
get_balance_tables(data_GBM_es.mean_stab, 'w_GBM_es.mean_stab', FALSE, confounders, estimand_m,treat_col_name, path, '3_GBM_stab_' )
get_balance_tables(data_GBM_es.mean_stab_trim, 'w_GBM_es.mean_stab', FALSE, confounders, estimand_m,treat_col_name, path, '4_GBM_stab_trim_' )
get_balance_tables(data_logr_raw, 'w_logr_raw', FALSE, confounders, estimand_m,treat_col_name, path, '5_logr_raw_' )
get_balance_tables(data_logr_raw_trim, 'w_logr_raw', FALSE, confounders, estimand_m,treat_col_name, path, '6_logr_raw_trim_' )
get_balance_tables(data_logr_stab, 'w_logr_stab', FALSE, confounders, estimand_m,treat_col_name, path, '7_logr_stab_' )
get_balance_tables(data_logr_stab_trim, 'w_logr_stab', FALSE, confounders, estimand_m,treat_col_name, path, '8_logr_stab_trim_' )

# get twang plots
plot_diagnostics(data_GBM_es.mean_raw, 'w_GBM_es.mean_raw', estimand_m, confounders, treat_col_name, path, '1_GBM_raw_')
plot_diagnostics(data_GBM_es.mean_raw_trim, 'w_GBM_es.mean_raw', estimand_m, confounders, treat_col_name, path, '2_GBM_raw_trim_')
plot_diagnostics(data_GBM_es.mean_stab, 'w_GBM_es.mean_stab', estimand_m, confounders, treat_col_name, path, '3_GBM_stab_')
plot_diagnostics(data_GBM_es.mean_stab_trim, 'w_GBM_es.mean_stab', estimand_m, confounders, treat_col_name, path, '4_GBM_stab_trim_')
plot_diagnostics(data_logr_raw, 'w_logr_raw', estimand_m, confounders, treat_col_name, path, '5_logr_raw_')
plot_diagnostics(data_logr_raw_trim, 'w_logr_raw', estimand_m, confounders, treat_col_name, path, '6_logr_raw_trim_')
plot_diagnostics(data_logr_stab, 'w_logr_stab', estimand_m, confounders, treat_col_name, path, '7_logr_stab_')
plot_diagnostics(data_logr_stab_trim, 'w_logr_stab', estimand_m, confounders, treat_col_name, path, '8_logr_stab_trim_')



# get summary statistics for weights
path = 'S://Evgeny/1_Sotrovimab/PS/1209_0_all/'
gbm_raw = read.csv(paste(path, '1_GBM_raw_', 'data.csv', sep=''))
gbm_stab = read.csv(paste(path, '3_GBM_stab_', 'data.csv', sep=''))
logr_raw = read.csv(paste(path, '5_logr_raw_', 'data.csv', sep=''))
logr_stab = read.csv(paste(path, '7_logr_stab_', 'data.csv', sep=''))


data_names = c('gbm_raw','gbm_stab','logr_raw','logr_stab')
weight_names = c('w_GBM_es.mean_raw','w_GBM_es.mean_stab','w_logr_raw','w_logr_stab')

df2 = setNames( data.frame(matrix(ncol = 6, nrow = 0)), c('weight_name','mean','sd','median', 'min','max')   )
for (i in 1:length(data_names)){

  d = eval(parse(text = data_names[i]))
  d$w_logr_stab
    
  df2[nrow(df2) + 1,] = c( weight_names[i], round( mean(  d[,c(weight_names[i])]),2 ), 
                          round( sd(  d[,c(weight_names[i])] ),2  ), round( median(  d[,c(weight_names[i])] ),2  ),
                          round( min(  d[,c(weight_names[i])]),2   ),round( max(  d[,c(weight_names[i])]),2   )            )

  }


write.csv(df2, paste(path,'0','_summary_stats_for_weights.csv', sep = '') )












































































  



