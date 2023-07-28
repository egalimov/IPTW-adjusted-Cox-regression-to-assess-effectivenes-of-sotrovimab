# IPTW-adjusted Cox regression to assess effectivenes of Sotrovimab

## Project Description
This repository contains the code for a retrospective cohort study to assess the effectivenes of Sotrovimab in the highest-risk patients

The main project aim was to assess the risk of COVID-19-related hospitalisation and/or COVID-19-related death within 28 days of the observed/imputed treatment date between highest-risk patients treated and not treated with sotrovimab.

Inverse probability of treatment weighting (IPTW) was used to balance baseline patient characteristics in the treated and untreated cohorts. IPTW based on propensity scores was used to adjust for measured confounders between the treated and untreated cohorts. Propensity scores (probability of treatment based on baseline covariates) were obtained using logistic regression or gradient boosting machine models. Propensity score models were used to predict the probability of treatment based on the following covariates: age, gender, time period of COVID-19 diagnosis (i.e., Omicron BA.1, BA.2 or BA.5, as defined above), presence of renal disease (binary), presence of multiple highest-risk conditions (≥2, binary), presence of high-risk conditions (binary), solid organ transplant (binary), COVID vaccination status (binary), time since vaccination, and ethnicity (see full list of variables and models in the publication). To obtain an appropriate estimation of the variance of the treatment effect and better control the type I error rate, inverse probability of treatment weights were stabilised. The balance in baseline characteristics between weighted treated and untreated groups was assessed using standardised differences.

Cox proportional hazards models with stabilised weights were performed to assess the hazard ratio (HR) of COVID-19-related hospitalisation and/or COVID-19-related death. Covariates not balanced after weighting (standardised differences >0.1) were included in the Cox proportional hazards model. IPTWs and accordingly doubly robust estimation was performed separately for each Cox model.



## About the repository
The repository contains scripts for the following steps:
1) Data preprocessing, checks for missing data. The output dataset is used in the next step 2. Script name: 1_data_preprocessing.py
2) Running models (logistic regression and GBM) to get propensity scores and inverse probabiility of treatment weights (various sets of weights were created: stabilised/not stabilised and trimmid/not trimmed). Assessing balance between control and treated groups after weighting using standardised differences and p-value distributions. Data with weights are saved in files endinbg data.csv and used in the step 3. 
Script name: 2_IPTW-models_and_balance_assessment__0_all.R 
	For subgroups analyses: 
	* 2_IPTW-models_and_balance_assessment__1_age_below65
	* 2_IPTW-models_and_balance_assessment__2_age_65andabove  
	* 2_IPTW-models_and_balance_assessment__3_timeperiod1
	* 2_IPTW-models_and_balance_assessment__3_timeperiod2
	* 2_IPTW-models_and_balance_assessment__3_timeperiod3
	* 2_IPTW-models_and_balance_assessment__6_renal_disease0
	* 2_IPTW-models_and_balance_assessment__6_renal_disease1

3) Performing IPTW Cox regression for the main cohort and 7 subgroups for each of 5 outcomes: 
	* COVID-19-related hospitalisation   
    * all-cause-related hospitalisation   
    * death
    * COVID-19-related hospitalisation and death
    * all-cause-related hospitalisation and death
Script name: 3_IPTW-adjusted Cox regression.R 



## Data sources
This analysis utilised the Discover dataset (https://www.discover-now.co.uk/), which is accessible via the Discover-NOW Health Data Research Hub for Real World Evidence through their data scientist specialists and IG committee-approved analysts, hosted by Imperial College Health Partners.


## Cohort derivation
This retrospective cohort study was based on data from the Discover dataset.

In the sotrovimab-treated cohort, the index date was defined as the date of sotrovimab prescription. Patients in the treated cohort must have had a recorded prescription for sotrovimab within 28 days of their COVID-19 diagnosis. In the untreated comparator cohort, the index date was the imputed treatment date. Dates were imputed based on the distribution of time to treatment in the treated cohort. The baseline period was defined as the 365 days immediately prior to index. Patients were followed up for 28 days from the index date (acute period), during which time patient outcomes were evaluated.

Patients in both cohorts were eligible for inclusion if they were aged ≥12 years on the index date and met at least one of the NHS highest-risk criteria for receiving early treatment with sotrovimab. At the time of study, these criteria included Down’s syndrome, solid cancer, haematological diseases (including cancers), renal disease, liver disease, immune-mediated inflammatory disorders, immune deficiencies, HIV/AIDS, solid organ and stem cell transplant recipients and rare neurological conditions.2,3 Patients meeting the NHS highest-risk criteria were identified via the presence of International Classification of Disease version 10 (ICD-10) and Systematized Nomenclature of Medicine (SNOMED) codes appearing at any time in the patient’s records since first registration in North West London.

Patients were also required to be non-hospitalised; to be considered non-hospitalised at the time of an event, patients must not have had an inpatient hospital visit (event from admission to discharge) starting on or before the date of the event, unless: the visit was a day case (in the NHS, a day case is a planned elective admission without a planned overnight stay, used to administer treatments under medical supervision or to conduct minor procedures), or the visit did not incur an overnight stay (the start and end time of the visit do not cross midnight).

Patients were excluded if they received more than one COVID-19 treatment (sotrovimab, nirmatrelvir/ritonavir, molnupiravir or remdesivir) in an outpatient setting before the index date; or were diagnosed with COVID-19 while hospitalised.


## Requirements
The work was done using:
Python 3.9.5	
	Packages:
	* Numpy	1.21
	* Pandas 1.3.4

R version 4.2.1	
	Packages:
	* Cobalt 4.4.1
	* PowerSurvEpi 0.1.3
	* Stats 4.2.1
	* Stringr 1.4.1
	* Survey 4.1-1
	* Survival 3.3-1
	* Survminer 0.4.9
	* Twang	2.5
	* Xtable 1.8-4
	* WeightIt 0.13.1


## Publication 
The current version of the manuscript in pre-print form is availbale at 

 
## Authors
- Evgeniy Galimov |   e.r.galimov@gmail.com   |  https://www.linkedin.com/in/evgeniygalimov


## License
This project is licensed under the [MIT License](https://opensource.org/licenses/MIT).




