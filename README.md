# *BMT Risk Prediction*

## Objective

This project evaluates the performance of ris triage tools for predicting clinical deterioration among hospitalized hematopoietic stem cell transplantation recipients. The research was performed on data from patients hospitalized at OHSU 2019-2022. 

## Cohort identification
All hospitalizations of hematopoietic stem cell transplants are included in the cohort for the study period of 2019-2022. Hospitalizations must include admissions to the inpatient wards or ICU. 

## Scripts

*01_cohort* applies inclusion/exclusion criteria to the BMT patient cohort

*02_scores* cleans EDI data extracted from the EHR

*03_data_prep* utilizes vitals, labs, respiratory criteria, and patient assessments from the EHR to calculate SIRS, MEWS, NEWS, and qSOFA continuously throughout cohort hospitalizations. Using admission-discharge-transfer data, determines the patient outcomes of ICU transfer, and/or death on wards or transfer to hospice. 

*04_analysis* performs hospitalization level and time series logistic regression for each triage score, and creates final tables and figures. 

