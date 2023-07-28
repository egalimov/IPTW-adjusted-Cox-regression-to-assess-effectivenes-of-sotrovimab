# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 15:03:57 2022

@author: galimove
"""

import pandas as pd
import numpy as np
import copy


# load the data
data = pd.read_csv('S://Evgeny/1_Sotrovimab/221017_without_excluded_final.csv')

# get df with missing data
miss_null = data.isnull()
miss_null_sum = miss_null.sum(axis=0)
miss_null_sum = pd.DataFrame({'col':miss_null_sum.index, 'null':miss_null_sum.values})

miss_na = data.isna()
miss_na_sum = miss_na.sum(axis=0)
miss_na_sum = pd.DataFrame({'col':miss_na_sum.index, 'na':miss_na_sum.values})

miss = pd.merge(miss_null_sum,miss_na_sum, on = ['col'])
miss_val = miss[ (miss['na']>0) | (miss['null']>0) ]



# list of columns of interest
col_list_of_interest = ['Age','Gender','EthnicCategory','EthnicityDescription','COVIDDate','SotrovimabDate','DaysToTreatment','IndexDate','VaccinationStatus','SolidOrganTransplants',
  'RenalDisease', 'HighRisk', 'ModerateRisk']

# dictionary with value counts for neach column of interest
val_counts = {}
for i in col_list_of_interest:
    val_counts[i] = a5[i].value_counts()




################################################################################
# data preproessing
data3 = copy.deepcopy(data)

# dervive time period - we set up 3 time periods for the study: 
    # 1) before 13/02/2022
    # 2) 13/02/2022 - 31/05/2022
    # 3) 01/06/2022 and after

def time_period_derivation(index_date):  
#for i in range(0,100):
    #cdform = pd.to_datetime(data3['COVIDDate'][i], format = '%d/%m/%Y %H:%M')
    cdform = pd.to_datetime(index_date, format = '%d/%m/%Y %H:%M')
    #print(i)
    #print(cdform)
    time_period = np.nan
    if (type(cdform) != pd._libs.tslibs.timestamps.Timestamp):
       # print('nan')
       return np.nan
    if ( ( cdform < pd.to_datetime('13/02/2022', format = '%d/%m/%Y') ) ):
        time_period = 1
    if ( ( cdform >= pd.to_datetime('13/02/2022', format = '%d/%m/%Y') ) & ( cdform <= pd.to_datetime('31/05/2022', format = '%d/%m/%Y' )) ):
        time_period = 2
    if (  ( cdform >= pd.to_datetime('01/06/2022', format = '%d/%m/%Y' )) ):
        time_period = 3
    #print(time_period)
    #print('---------------')
    return time_period

data3['Time_period'] = list(  map(time_period_derivation, data3['IndexDate'])   )


# saving table with value counts for time periods
a7 = data3[['PatientKey', 'IndexDate3','Time_period']]
a7['Time_period'].value_counts().to_csv('S://Evgeny/1_Sotrovimab/1_Time_period_value_counts.csv')
a7_2 = a7.sort_values(by = 'IndexDate3', ascending = [True])
# getting earliest and latest dates in the study
time_periods_df = pd.DataFrame( {'earliest - period 0': [ a7_2.iloc[0,1] ] , 'latest - period 3': [a7_2.iloc[a7_2.shape[0]-1,1] ]} )
time_periods_df.to_csv('S://Evgeny/1_Sotrovimab/1_Time_period_dates_updated.csv')
 

###### deriving variables

# derive mult_high_risk_bin
data3['Mult_high_risk_bin'] = list(  map(lambda x: 1 if x>=2 else 0, data3['HighRisk'] )   )

# derive moderate_risk_bin
data3['Moderate_risk_bin'] = list(  map(lambda x: 1 if x>=1 else 0, data3['ModerateRisk'] )   )

# derive vacc_status_bin
data3['Vacc_stcatus_bin'] = list(  map(lambda x: 1 if x=='Fully vaccinated' else 0, data3['VaccinationStatus'] )   )

# derive vacc_status_ordered
def vaccinated_status(x):
    out = 0
    if x == 'Unvaccinated':
        out = 0
    if x == 'Partially vaccinated':
        out = 1
    if x == 'Fully vaccinated':
        out = 2
    return out        
data3['Vacc_status_ordered'] = list(  map(vaccinated_status, data3['VaccinationStatus'] )   )


# derive target 1 - all_hosp_and_death
def taregt1(APC_StartDate, DateofDeath ):
     
    APC_StartDate_f = pd.to_datetime(APC_StartDate, format = '%d/%m/%Y %H:%M')
    DateofDeath_f = pd.to_datetime(DateofDeath, format = '%d/%m/%Y %H:%M')
    
    target1 = np.nan
    if ( (type(APC_StartDate_f) == pd._libs.tslibs.timestamps.Timestamp) & (type(DateofDeath_f) == pd._libs.tslibs.timestamps.Timestamp) ):
        if ( APC_StartDate_f < DateofDeath_f):
            target1 = APC_StartDate_f
        if ( DateofDeath_f <= APC_StartDate_f):
            target1 = DateofDeath_f

    if ( (type(APC_StartDate_f) == pd._libs.tslibs.timestamps.Timestamp) & (type(DateofDeath_f) != pd._libs.tslibs.timestamps.Timestamp) ):
        target1 = APC_StartDate_f

    if ( (type(APC_StartDate_f) != pd._libs.tslibs.timestamps.Timestamp) & (type(DateofDeath_f) == pd._libs.tslibs.timestamps.Timestamp) ):
        target1 = DateofDeath_f
    return target1
    

data3['All_hosp_and_death'] = list(  map(taregt1, data3['APC_StartDate'], data3['DateofDeath'])   )


# derive target 2 - covid_hosp_and_death
def taregt2(APC_COVIDPrimaryStartDate, DateofDeath ):

    APC_COVIDPrimaryStartDate_f = pd.to_datetime(APC_COVIDPrimaryStartDate, format = '%d/%m/%Y %H:%M')
    DateofDeath_f = pd.to_datetime(DateofDeath, format = '%d/%m/%Y %H:%M')
    
    target2 = np.nan
    if ( (type(APC_COVIDPrimaryStartDate_f) == pd._libs.tslibs.timestamps.Timestamp) & (type(DateofDeath_f) == pd._libs.tslibs.timestamps.Timestamp) ):
        if ( APC_COVIDPrimaryStartDate_f < DateofDeath_f):
            target2 = APC_COVIDPrimaryStartDate_f
        if ( DateofDeath_f <= APC_COVIDPrimaryStartDate_f):
            target2 = DateofDeath_f

    if ( (type(APC_COVIDPrimaryStartDate_f) == pd._libs.tslibs.timestamps.Timestamp) & (type(DateofDeath_f) != pd._libs.tslibs.timestamps.Timestamp) ):
        target2 = APC_COVIDPrimaryStartDate_f

    if ( (type(APC_COVIDPrimaryStartDate_f) != pd._libs.tslibs.timestamps.Timestamp) & (type(DateofDeath_f) == pd._libs.tslibs.timestamps.Timestamp) ):
        target2 = DateofDeath_f
    return target2
    

data3['Covid_hosp_and_death'] = list(  map(taregt2, data3['APC_COVIDPrimaryStartDate'], data3['DateofDeath'])   )


# derive treatment group 
data3.columns
data3['Treatment'] = data3['SotrovimabDate'].apply(lambda x: 1 if type(x)==str else 0)


# define survival time and status for survival analysis | All_hosp_and_death
# acute period = index_date+1 - index_date+28
def get_survival_time_All_hosp_and_death(index_date, All_hosp_and_death):
    #index_date = data3['IndexDate'][12]
    #All_hosp_and_death = data3['All_hosp_and_death'][12]
    
    index_date_f = pd.to_datetime(index_date, format = '%d/%m/%Y %H:%M')
    All_hosp_and_death_f = pd.to_datetime(All_hosp_and_death, format = '%d/%m/%Y %H:%M')
    
    if ( type(All_hosp_and_death_f) == pd._libs.tslibs.timestamps.Timestamp ):
        out = (All_hosp_and_death_f - index_date_f).days

    if ( type(All_hosp_and_death_f) != pd._libs.tslibs.timestamps.Timestamp ):
        out = np.nan

    return out



data3['survival_time_All_hosp_and_death'] = list(  map(get_survival_time_All_hosp_and_death, data3['IndexDate'], data3['All_hosp_and_death'])   )

# get status
data3['status_time_All_hosp_and_death'] = data3['survival_time_All_hosp_and_death'].apply(lambda x: 0 if np.isnan(x) else 1)

# get survival time for censored patients
data3['survival_time_All_hosp_and_death'] = data3['survival_time_All_hosp_and_death'].apply(lambda x: 28 if np.isnan(x) else x)
    


# define survival time and status for survival analysis | Covid_hosp_and_death
# acute period = index_date+1 - index_date+28
def get_survival_time_Covid_hosp_and_death(index_date, Covid_hosp_and_death):
    #index_date = data3['IndexDate'][12]
    #All_hosp_and_death = data3['All_hosp_and_death'][12]
    
    index_date_f = pd.to_datetime(index_date, format = '%d/%m/%Y %H:%M')
    Covid_hosp_and_death_f = pd.to_datetime(Covid_hosp_and_death, format = '%d/%m/%Y %H:%M')
    
    if ( type(Covid_hosp_and_death_f) == pd._libs.tslibs.timestamps.Timestamp ):
        out = (Covid_hosp_and_death_f - index_date_f).days

    if ( type(Covid_hosp_and_death_f) != pd._libs.tslibs.timestamps.Timestamp ):
        out = np.nan

    return out
    
data3['survival_time_Covid_hosp_and_death'] = list(  map(get_survival_time_All_hosp_and_death, data3['IndexDate'], data3['Covid_hosp_and_death'])   )

# get status
data3['status_Covid_hosp_and_death'] = data3['survival_time_Covid_hosp_and_death'].apply(lambda x: 0 if np.isnan(x) else 1)

# get survival time for censored patients
data3['survival_time_Covid_hosp_and_death'] = data3['survival_time_Covid_hosp_and_death'].apply(lambda x: 28 if np.isnan(x) else x)


######
# get status and survival time 

###     APC_StartDate
# get survival time
data3['survival_time_All_hosp'] = list(  map(get_survival_time_All_hosp_and_death, data3['IndexDate'], data3['APC_StartDate'])   )

# get status      
data3['status_All_hosp'] = data3['survival_time_All_hosp'].apply(lambda x: 0 if np.isnan(x) else 1)

# get survival time for censored patients
data3['survival_time_All_hosp'] = data3['survival_time_All_hosp'].apply(lambda x: 28 if np.isnan(x) else x)



###     APC_COVIDPrimaryStartDate
# get survival time
data3['survival_time_Covid_hosp'] = list(  map(get_survival_time_All_hosp_and_death, data3['IndexDate'], data3['APC_COVIDPrimaryStartDate'])   )


# get status      
data3['status_Covid_hosp'] = data3['survival_time_Covid_hosp'].apply(lambda x: 0 if np.isnan(x) else 1)

# get survival time for censored patients
data3['survival_time_Covid_hosp'] = data3['survival_time_Covid_hosp'].apply(lambda x: 28 if np.isnan(x) else x)





###     DateofDeath
# get survival time
data3['survival_time_Death'] = list(  map(get_survival_time_All_hosp_and_death, data3['IndexDate'], data3['DateofDeath'])   )

# get status      
data3['status_Death'] = data3['survival_time_Death'].apply(lambda x: 0 if np.isnan(x) else 1)

# get survival time for censored patients
data3['survival_time_Death'] = data3['survival_time_Death'].apply(lambda x: 28 if np.isnan(x) else x)




# for DaysSinceLastVaccination convert nans to zeros
data3['Days_SinceLastVaccination'] = data3['DaysSinceLastVaccination'].apply(lambda x: 0 if np.isnan(x) else x)




    
# check data3 missing data
miss_null = data3.isnull()
miss_null_sum = miss_null.sum(axis=0)
miss_null_sum = pd.DataFrame({'col':miss_null_sum.index, 'null':miss_null_sum.values})

miss_na = data3.isna()
miss_na_sum = miss_na.sum(axis=0)
miss_na_sum = pd.DataFrame({'col':miss_na_sum.index, 'na':miss_na_sum.values})

miss = pd.merge(miss_null_sum,miss_na_sum, on = ['col'])
miss_val = miss[ (miss['na']>0) | (miss['null']>0) ]







#####################
# get the final datset with the outcome:  all-cause-related hospitalisation and death
final_target1 = data3[['Treatment','Age', 'Gender','EthnicCategory','Time_period','RenalDisease','Mult_high_risk_bin',
               'Moderate_risk_bin','SolidOrganTransplants', 'Vacc_stcatus_bin','Vacc_status_ordered','Days_SinceLastVaccination',
               'status_time_All_hosp_and_death', 'survival_time_All_hosp_and_death']]

miss_null = final_target1.isnull()
miss_null_sum = miss_null.sum(axis=0)
miss_null_sum = pd.DataFrame({'col':miss_null_sum.index, 'null':miss_null_sum.values})

miss_na = final_target1.isna()
miss_na_sum = miss_na.sum(axis=0)
miss_na_sum = pd.DataFrame({'col':miss_na_sum.index, 'na':miss_na_sum.values})

miss = pd.merge(miss_null_sum,miss_na_sum, on = ['col'])
miss_val = miss[ (miss['na']>0) | (miss['null']>0) ]

final_target1.to_csv('S://Evgeny/1_Sotrovimab/221209_target1_CE.csv')

### get the final datset with the outcome:  COVID-19-related hospitalisation and death
final_target2 = data3[['Treatment','Age', 'Gender','EthnicCategory','Time_period','RenalDisease','Mult_high_risk_bin',
               'Moderate_risk_bin','SolidOrganTransplants', 'Vacc_stcatus_bin','Vacc_status_ordered','Days_SinceLastVaccination',
               'status_Covid_hosp_and_death', 'survival_time_Covid_hosp_and_death']]

miss_null = final_target2.isnull()
miss_null_sum = miss_null.sum(axis=0)
miss_null_sum = pd.DataFrame({'col':miss_null_sum.index, 'null':miss_null_sum.values})

miss_na = final_target2.isna()
miss_na_sum = miss_na.sum(axis=0)
miss_na_sum = pd.DataFrame({'col':miss_na_sum.index, 'na':miss_na_sum.values})

miss = pd.merge(miss_null_sum,miss_na_sum, on = ['col'])
miss_val = miss[ (miss['na']>0) | (miss['null']>0) ]

final_target2.to_csv('S://Evgeny/1_Sotrovimab/221209_target2_CE.csv')



###  get the final datset with the outcomes:
                                            # COVID-19-related hospitalisation   
                                            # all-cause-related hospitalisation   
                                            # death
                                            # COVID-19-related hospitalisation and death
                                            # all-cause-related hospitalisation and death

final_target3 = data3[['Treatment','Age', 'Gender','EthnicCategory','Time_period','RenalDisease','Mult_high_risk_bin',
               'Moderate_risk_bin','SolidOrganTransplants', 'Vacc_stcatus_bin','Vacc_status_ordered','Days_SinceLastVaccination',
               'status_time_All_hosp_and_death', 'survival_time_All_hosp_and_death','status_Covid_hosp_and_death', 'survival_time_Covid_hosp_and_death',
               'status_All_hosp','survival_time_All_hosp', 'status_Covid_hosp','survival_time_Covid_hosp','status_Death','survival_time_Death']]

miss_null = final_target3.isnull()
miss_null_sum = miss_null.sum(axis=0)
miss_null_sum = pd.DataFrame({'col':miss_null_sum.index, 'null':miss_null_sum.values})

miss_na = final_target3.isna()
miss_na_sum = miss_na.sum(axis=0)
miss_na_sum = pd.DataFrame({'col':miss_na_sum.index, 'na':miss_na_sum.values})

miss = pd.merge(miss_null_sum,miss_na_sum, on = ['col'])
miss_val = miss[ (miss['na']>0) | (miss['null']>0) ]

final_target3.to_csv('S://Evgeny/1_Sotrovimab/221209_target3_CE.csv')


################################################################################
################################################################################
################################################################################
################################################################################
### Getting and saving datasets for subgroup analysis

data3_age_before65 = final_target3[final_target3['Age'] < 65]
data3_age_65andabove = final_target3[final_target3['Age'] >= 65]

data3_time_period_1 = final_target3[final_target3['Time_period'] == 1]
data3_time_period_2 = final_target3[final_target3['Time_period'] == 2]
data3_time_period_3 = final_target3[final_target3['Time_period'] == 3]

data3_renal_disease_0 = final_target3[final_target3['RenalDisease'] == 0]
data3_renal_disease_1 = final_target3[final_target3['RenalDisease'] == 1]



data3_age_before65.to_csv('S://Evgeny/1_Sotrovimab/221209_1_age_before65.csv')
data3_age_65andabove.to_csv('S://Evgeny/1_Sotrovimab/221209_2_age_65andabove.csv')

data3_time_period_1.to_csv('S://Evgeny/1_Sotrovimab/221209_3_time_period_1.csv')
data3_time_period_2.to_csv('S://Evgeny/1_Sotrovimab/221209_4_time_period_2.csv')
data3_time_period_3.to_csv('S://Evgeny/1_Sotrovimab/221209_5_time_period_3.csv')

data3_renal_disease_0.to_csv('S://Evgeny/1_Sotrovimab/221209_renal_disease_0.csv')
data3_renal_disease_1.to_csv('S://Evgeny/1_Sotrovimab/221209_renal_disease_1.csv')


















