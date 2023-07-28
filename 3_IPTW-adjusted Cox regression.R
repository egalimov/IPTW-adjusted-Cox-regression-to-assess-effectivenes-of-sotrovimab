library(survey)
library(survminer)
library(survival)
library(powerSurvEpi)
library(stringr)
      
# return the model or summary of the model, specified in the output ('model' or 'summary')

# function to run weighted cox regression
    # data1 = d_raw   # dataset
    # form_str = 'Treatment'    # variable for which we assess the effect
    # weights_m ='w_GBM_es.mean_raw'     # file with weights
    # output = 'summary'     # output format ('model' or 'summary')
surv_function_weighted = function(form_str, data1, weights_m, output){
  design.obj = svydesign(ids=~1, weights = ~data1[,c(weights_m)], data = data1)

  a = svycoxph(form_str, data = data1, design = design.obj)
  
  if (output == 'model') {out = a}
  if (output == 'summary') {out = summary(a)}
  return(out)
}

# function to run cox regression without weights
    # data1 = d_raw   # dataset
    # form_str = 'Treatment'    # variable for which we assess the effect
    # output = 'summary'     # output format ('model' or 'summary')
surv_function_unweighted = function(form_str, data1, output){
  a = coxph(form_str, data = data1)
  if (output == 'model') {out = a}
  if (output == 'summary') {out = summary(a)}
  return(out)
}



# function to test the proportional hazards assumption for a Cox regression
    # data1 = d_raw   # dataset
    # form_str = 'Treatment'    # variable for which we assess the effect
    # weights_m ='w_GBM_es.mean_raw'     # file with weights
    surv_coxzph_function = function(form_str, data1, weights_m){
  design.obj = svydesign(ids=~1, weights = ~data1[,c(weights_m)], data = data1)
  t = cox.zph(svycoxph(form_str, data = data1, design = design.obj))
  return(t)
}




# function to run weighted cox regression and get all the output tables
    # data -  dataset
    # form_str = formula in string format
    # weights_m ='w_GBM_es.mean_raw'     # file with weights
    # target = name of file for saving
univar_weighted = function(form_str, target, data, weights_m){
  # get summary of weighted regression   
  a = surv_function_weighted(form_str, data,weights_m,  "summary")
  
  # make df from the summary and populate it
  c = data.frame(stringsAsFactors = FALSE)
  c = data.frame(matrix(NA, nrow = 3, ncol = 14))
  colnames(c) = c("target","formula", "n", "events", "concordance", "se_concordance", 
                  "test", "test_value", "test_df", "test_p") #, 
  c$target[1] = target
  c$formula[1] = Reduce( paste, deparse(form_str) )
  c$n[1] = a$n
  c$events[1] = a$nevent
  c$concordance[1] = a$concordance[1]
  c$se_concordance[1] = a$concordance[2]
  c$test[1] = 'Likelihood ratio test'
  c$test_value[1] = a$logtest[1]
  c$test_df[1] = a$logtest[2]    
  c$test_p[1] = a$logtest[3] 
  
  c$test[2] = 'Wald test'
  c$test_value[2] = a$waldtest[1]
  c$test_df[2] = a$waldtest[2]
  c$test_p[2] = a$waldtest[3]   
  
  c$test[3] = 'Score (logrank) test'
  c$test_value[3] = a$sctest[1]
  c$test_df[3] = a$sctest[2]  
  c$test_p[3] = a$sctest[3]  
  
  
  # get summary of the test for the proportional hazards assumption 
  t = surv_coxzph_function(form_str, data, weights_m)
  # make df from the summary and populate it
  t2 = as.data.frame(t$table)
  t2 = cbind(feature = rownames(t2), t2)
  rownames(t2) = 1:nrow(t2)
  colnames(t2) = c('GOF_PH_assumption_test','GOF_chisq','GOF_DF','GOF_p')
  
  if (nrow(c)>=nrow(t2)){
    row_d = nrow(c) - nrow(t2)
    if (row_d != 0) {for (i in 1:row_d) {t2[nrow(t2)+1,] = NA}}
    c3 = cbind(c, t2)
  }
  if (nrow(c)<nrow(t2)){
    row_d = nrow(t2) - nrow(c)
    for (i in 1:row_d) {c[nrow(c)+1,] = NA}
    c3 = cbind(c, t2)
  }
  

  # make df with coefficients/CI/power/ sample size to reach alpha0.05 
  c2 = cbind(as.data.frame(a$coefficients), as.data.frame(a$conf.int))
  
      # get power of the analysis
                # only consider treatment
            if  (str_split_fixed(as.character(form_str)[3], " \\+ ", 10)[2]==''){
            data$group = ifelse(data$Treatment==1,'E','C')
            form_str2 = as.formula(paste(as.character(form_str)[2],' ~ group',sep=''))
            p_test = powerCT(form_str2, dat=data, nE=nrow(data[data$group=='E',]), nC = nrow(data[data$group=='C',]),RR=c2[1,7],alpha=0.05 )
            colnames(c3)[11] = 'Power of the model'
            c3[1,11] = p_test[5]
            }
                # consider treatment and one covariate (bin or cont)
            if  (str_split_fixed(as.character(form_str)[3], " \\+ ", 10)[2]!=''){
            treat_var = str_split_fixed(as.character(form_str)[3], " \\+ ", 10)[1]
            covar_var = str_split_fixed(as.character(form_str)[3], " \\+ ", 10)[2]
            
            p_test_epi = powerEpi(X1=data[,c(treat_var)],
                                  X2=data[,c(covar_var)],
                                  failureFlag=data[,c(status)],
                                  n=nrow(data), 
                                  theta=c2[1,7],
                                  alpha=0.05)
            colnames(c3)[11] = 'Power of the model'
            c3[1,11] = p_test_epi[1]
            }
       # get sample size to reach alpha0.05 power of the analysis
            if  (str_split_fixed(as.character(form_str)[3], " \\+ ", 10)[2]==''){
              data$group = ifelse(data$Treatment==1,'E','C')
              form_str2 = as.formula(paste(as.character(form_str)[2],' ~ group',sep=''))
              p_size = ssizeCT(form_str2, dat=data, power=0.8,
                               k=nrow(data[data$group=='E',])/nrow(data[data$group=='C',]),
                               RR=c2[1,7],alpha=0.05 )
              colnames(c3)[12] = 'Sample size to detect the effect, pow=0.8,alpha=0.05'
              c3[1,12] = 'nE'
              c3[2,12] = p_size[[5]][1]
              c3[1,13] = 'nC'
              c3[2,13] = p_size[[5]][2]
            }
            # consider treatment and one covariate (bin or cont)
            if  (str_split_fixed(as.character(form_str)[3], " \\+ ", 10)[2]!=''){
              treat_var = str_split_fixed(as.character(form_str)[3], " \\+ ", 10)[1]
              covar_var = str_split_fixed(as.character(form_str)[3], " \\+ ", 10)[2]
              
              s_test_epi = ssizeEpi(X1=data[,c(treat_var)],
                                    X2=data[,c(covar_var)],
                                    failureFlag=data[,c(status)],
                                    power=0.8, 
                                    theta=c2[1,7],
                                    alpha=0.05)
              colnames(c3)[12] = 'Sample size to detect the effect, pow=0.8,alpha=0.05'
              c3[1,12] = s_test_epi[1]
            }
    
  c2 = c2[ ,c("coef", "exp(coef)", "se(coef)", "lower .95", "upper .95", "z", "Pr(>|z|)") ]
  c2 = cbind(feature = rownames(c2), c2)
  rownames(c2) = 1:nrow(c2)
  
  if (nrow(c3)>=nrow(c2)){
    row_d = nrow(c3) - nrow(c2)
    if (row_d != 0) {for (i in 1:row_d) {c2[nrow(c2)+1,] = NA}}
    c4 = cbind(c3, c2)
  }
  if (nrow(c3)<nrow(c2)){
    row_d = nrow(c2) - nrow(c3)
    for (i in 1:row_d) {c3[nrow(c3)+1,] = NA}
    c4 = cbind(c3, c2)
  }
  return(c4)
}



# function to run cox regression without weights and get all the output tables
    # data -  dataset
    # form_str = formula in string format
    # weights_m ='w_GBM_es.mean_raw'     # file with weights
    # target = name of file for saving
univar_unweighted = function(form_str, target, data){
  # get summary of cox regression   
  a = surv_function_unweighted(form_str, data,  "summary")

  # make df from the summary and populate it
  c = data.frame(stringsAsFactors = FALSE)
  c = data.frame(matrix(NA, nrow = 3, ncol = 14))
  colnames(c) = c("target","formula", "n", "events", "concordance", "se_concordance", 
                  "test", "test_value", "test_df", "test_p") #, 
  # "GOF_PH_assumption_test", "GOF_chisq","GOF_DF", "GOF_p")
  c$target[1] = target
  c$formula[1] = Reduce( paste, deparse(form_str) )
  c$n[1] = a$n
  c$events[1] = a$nevent
  c$concordance[1] = a$concordance[1]
  c$se_concordance[1] = a$concordance[2]
  c$test[1] = 'Likelihood ratio test'
  c$test_value[1] = a$logtest[1]
  c$test_df[1] = a$logtest[2]    
  c$test_p[1] = a$logtest[3] 
  
  c$test[2] = 'Wald test'
  c$test_value[2] = a$waldtest[1]
  c$test_df[2] = a$waldtest[2]
  c$test_p[2] = a$waldtest[3]   
  
  c$test[3] = 'Score (logrank) test'
  c$test_value[3] = a$sctest[1]
  c$test_df[3] = a$sctest[2]  
  c$test_p[3] = a$sctest[3]  
  
  # get summary of the test for the proportional hazards assumption 
  t = surv_coxzph_function(form_str, data, weights_m)
  # make df from the summary and populate it
  t2 = as.data.frame(t$table)
  t2 = cbind(feature = rownames(t2), t2)
  rownames(t2) = 1:nrow(t2)
  colnames(t2) = c('GOF_PH_assumption_test','GOF_chisq','GOF_DF','GOF_p')
  
  if (nrow(c)>=nrow(t2)){
    row_d = nrow(c) - nrow(t2)
    if (row_d != 0) {for (i in 1:row_d) {t2[nrow(t2)+1,] = NA}}
    c3 = cbind(c, t2)
  }
  if (nrow(c)<nrow(t2)){
    row_d = nrow(t2) - nrow(c)
    for (i in 1:row_d) {c[nrow(c)+1,] = NA}
    c3 = cbind(c, t2)
  }
  

  # make df with coefficients/CI/power/ sample size to reach alpha0.05 
  c2 = cbind(as.data.frame(a$coefficients), as.data.frame(a$conf.int))
            # get power of the analysis
            # only consider treatment
            if  (str_split_fixed(as.character(form_str)[3], " \\+ ", 10)[2]==''){
              data$group = ifelse(data$Treatment==1,'E','C')
              form_str2 = as.formula(paste(as.character(form_str)[2],' ~ group',sep=''))
              p_test = powerCT(form_str2, dat=data, nE=nrow(data[data$group=='E',]), nC = nrow(data[data$group=='C',]),RR=c2[1,7],alpha=0.05 )
              colnames(c3)[11] = 'Power of the model'
              c3[1,11] = p_test[5]
            }
            # consider treatment and one covariate (bin or cont)
            if  (str_split_fixed(as.character(form_str)[3], " \\+ ", 10)[2]!=''){
              treat_var = str_split_fixed(as.character(form_str)[3], " \\+ ", 10)[1]
              covar_var = str_split_fixed(as.character(form_str)[3], " \\+ ", 10)[2]
              
              p_test_epi = powerEpi(X1=data[,c(treat_var)],
                                    X2=data[,c(covar_var)],
                                    failureFlag=data[,c(status)],
                                    n=nrow(data), 
                                    theta=c2[1,7],
                                    alpha=0.05)
              colnames(c3)[11] = 'Power of the model'
              c3[1,11] = p_test_epi[1]
            }
            # get sample size to reach alpha0.05 power of the analysis
            if  (str_split_fixed(as.character(form_str)[3], " \\+ ", 10)[2]==''){
              data$group = ifelse(data$Treatment==1,'E','C')
              form_str2 = as.formula(paste(as.character(form_str)[2],' ~ group',sep=''))
              p_size = ssizeCT(form_str2, dat=data, power=0.8,
                               k=nrow(data[data$group=='E',])/nrow(data[data$group=='C',]),
                               RR=c2[1,7],alpha=0.05 )
              colnames(c3)[12] = 'Sample size to detect the effect, pow=0.8,alpha=0.05'
              c3[1,12] = 'nE'
              c3[2,12] = p_size[[5]][1]
              c3[1,13] = 'nC'
              c3[2,13] = p_size[[5]][2]
            }
            # consider treatment and one covariate (bin or cont)
            if  (str_split_fixed(as.character(form_str)[3], " \\+ ", 10)[2]!=''){
              treat_var = str_split_fixed(as.character(form_str)[3], " \\+ ", 10)[1]
              covar_var = str_split_fixed(as.character(form_str)[3], " \\+ ", 10)[2]
              
              s_test_epi = ssizeEpi(X1=data[,c(treat_var)],
                                    X2=data[,c(covar_var)],
                                    failureFlag=data[,c(status)],
                                    power=0.8, 
                                    theta=c2[1,7],
                                    alpha=0.05)
              colnames(c3)[12] = 'Sample size to detect the effect, pow=0.8,alpha=0.05'
              c3[1,12] = s_test_epi[1]
            }

  c2 = c2[ ,c("coef", "exp(coef)", "se(coef)", "lower .95", "upper .95", "z", "Pr(>|z|)") ]
  c2 = cbind(feature = rownames(c2), c2)
  rownames(c2) = 1:nrow(c2)
  
  if (nrow(c3)>=nrow(c2)){
    row_d = nrow(c3) - nrow(c2)
    if (row_d != 0) {for (i in 1:row_d) {c2[nrow(c2)+1,] = NA}}
    c4 = cbind(c3, c2)
  }
  if (nrow(c3)<nrow(c2)){
    row_d = nrow(c2) - nrow(c3)
    for (i in 1:row_d) {c3[nrow(c3)+1,] = NA}
    c4 = cbind(c3, c2)
  }
  return(c4)
}




#################################################################################
### Running the analysis

# lists to iterate over the main cohort and 7 subgroups
paths = c('S://Evgeny/1_Sotrovimab/PS/1209_0_all/',
          'S://Evgeny/1_Sotrovimab/PS/1209_1_age_before65/',
          'S://Evgeny/1_Sotrovimab/PS/1209_2_age_65andabove/',
          'S://Evgeny/1_Sotrovimab/PS/1209_3_time_period_1/',
          'S://Evgeny/1_Sotrovimab/PS/1209_3_time_period_2/',
          'S://Evgeny/1_Sotrovimab/PS/1209_3_time_period_3/',
          'S://Evgeny/1_Sotrovimab/PS/1209_6_renal_disease_0/',
          'S://Evgeny/1_Sotrovimab/PS/1209_6_renal_disease_1/'
          )

datasets = c('3_GBM_stab_', '7_logr_stab_','7_logr_stab_','7_logr_stab_','7_logr_stab_','7_logr_stab_','7_logr_stab_','7_logr_stab_')

weights_list = c('w_GBM_es.mean_stab', 'w_logr_stab','w_logr_stab','w_logr_stab','w_logr_stab','w_logr_stab','w_logr_stab','w_logr_stab')

vars_list = c('Treatment + Time_period', 
              'Treatment + Mult_high_risk_bin', 
              'Treatment',
              'Treatment + Mult_high_risk_bin + Age + EthnicCategory',
              'Treatment',
              'Treatment + Days_SinceLastVaccination',
              'Treatment + Mult_high_risk_bin + Days_SinceLastVaccination',
              'Treatment'
              )




# lists to iterate over 5 outcomes: Death, All types of hospitalisation, Covid-related hospitalisation, All types of hospitalisation+Death,  Covid-related hospitalisation + Death
surv_names = c('0_1_WholeData__AllHosp_death','0_2_WholeData__CovidHosp_death',
               '0_3_WholeData__AllHosp','0_4_WholeData__CovidHosp', '0_5_WholeData__death')
statuses_list = c("status_time_All_hosp_and_death","status_Covid_hosp_and_death",
             "status_All_hosp","status_Covid_hosp","status_Death")
events_list = c("survival_time_All_hosp_and_death","survival_time_Covid_hosp_and_death",
           "survival_time_All_hosp","survival_time_Covid_hosp","survival_time_Death")



###  iterating over the cohorts and 5 outcomes for each cohort
# iterating over cohorts
for c in seq(1:8){

path = paths[i]
d_stab = read.csv(paste(path,datasets[c], 'data.csv', sep=''))   # loading the dataset
weights_m = weights_list[c]   # file with weights
vars = vars_list[c]      # variables to include in the formula

# create dir for each cohort
dir.create(path)


# iterating over 5 outcomes
for (i in seq(1:5)){
  
events = events_list[i]
status = statuses_list[i]
surv_name = surv_names[i]

# create formula for cox analysis
form_str = as.formula(paste('Surv(', events,',', status,')~',vars, sep=''))

# create results for weighted cox regression
surv_out_weighted = univar_weighted(form_str, surv_name, d_stab, weights_m)
write.csv(surv_out_weighted, paste(path, surv_name,'_1_weighted.csv', sep=''))

# create results for unweighted cox regression
surv_out_unweighted = univar_unweighted(form_str, surv_name, d_stab)
write.csv(surv_out_unweighted, paste(path, surv_name,'_2_unweighted.csv', sep=''))




# get median time
survfit_temp = survfit(form_str, data = d_stab)
sink(file = paste(path, surv_name,'_surv_summary.txt', sep=''))
print( survfit_temp )
sink()

# get the table with events
sink(file = paste(path, surv_name,'_surv_events.txt', sep=''))
print( summary(survfit_temp) )
sink()



# plot survival curves
survfit_temp$call$formula = form_str
ggsurvplot(survfit_temp)
p2 = ggsurvplot(survfit_temp, pval=TRUE, conf.int=TRUE, risk.table="abs_pct",risk.table.col="strata",
                linetype="strata",surv.median.line="hv",ggtheme=theme_bw(),ncensor.plot=TRUE,
                # palette=c("#81B622","#FA26A0","#29A0B1", "#FF4500"),
                font.title =  c(16, 'bold'),
                font.x =  c(16, 'bold'),
                font.y =  c(16, 'bold'),
                #fontsize = 12
                font.tickslab = c(16, 'bold')           )

# add risk table
p2$table = ggrisktable(survfit_temp, #tables_theme = theme_cleantable(),
                       color = 'strata',
                       risk.table.type = 'percentage',
                       ggtheme = theme(axis.title = element_text(size = 16, face = 'bold'),
                                       axis.text = element_text(size = 10),
                                       panel.background = element_rect(colour = 'white',
                                                                       fill = 'white'),
                                       #panel.grid.major = element(color = 'black',
                                       #                           size = 0.5),
                                       panel.border = element_rect(colour='black', size=0.5,
                                                                   fill = NA)
                       ),
                       fontsize = 3)

# add censor plot
p2$ncensor.plot = p2$ncensor.plot + theme(axis.text = element_text(size = 10),
                                          axis.title = element_text(size = 16, face = 'bold'))

png(paste(path, surv_name,'_surv_plot.png', sep=''), width = 20, height = 15, units = 'cm', res = 600  )
print(p2)
dev.off() 



}






}



































































