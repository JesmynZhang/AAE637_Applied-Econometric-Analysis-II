* Load the dataset
**#
use "/Users/jesmyn/Downloads/schoolexp.dta", clear

**Section 1: Propensity Score Matching and Inverse Probability Weights

*1.1

* Calculate the 25th percentile of the expenditure per pupil
summarize exppp, detail
local p25 = r(p25)
* Generate the binary variable lowexpend: 1 if exppp is below the 25th percentile, 0 otherwise
gen lowexpend = exppp < `p25'
* To view the first few rows and check the new variable
list exppp lowexpend in 1/5

*1.2

* Run the OLS regression with specified controls and cluster standard errors at the district level
regress math4 lowexpend lunch enroll, vce(cluster dcode)
* Store the regression results for later use
estimates store Model1

*1.4
* Estimate a logit model to calculate propensity scores
logit lowexpend lunch enroll
* Generate propensity scores based on the logit model
predict p_scores, pr
* View the generated propensity scores
list p_scores in 1/5
* Summarize the propensity scores to get mean, SD, min, and max
summarize p_scores

*1.6
* Create a histogram of the propensity scores for treatment and control groups
histogram p_scores, by(lowexpend) title("Common Support for Propensity Scores") xlabel(0(0.1)1) name(CommonSupport, replace)

*1.7
* Estimate the propensity score matching model with a single nearest neighbor
 psmatch2 lexppp, outcome(math4) pscore(p_scores) neighbor(1)
* Save the dataset with the matching results
save matched_data.dta, replace
* Store the results for later use
estimates store PSM1

*1.8
* Estimate the propensity score matching model using kmatch with a single nearest neighbor
kmatch ps lowexpend, pscore(p_score)
* Store the results
estimates store PSM2

*1.9
* Re-estimate the propensity score matching model using kmatch with common support adjustment
kmatch ps lowexpend, pscore(p_score) comsup
* Store the results for later use
estimates store PSM3

*1.10
* Determine the common support region
sum p_score if lowexpend == 1
local p_min_treated = r(min)
local p_max_treated = r(max)

sum p_score if lowexpend == 0
local p_min_control = r(min)
local p_max_control = r(max)

local p_common_min = max(`p_min_treated', `p_min_control')
local p_common_max = min(`p_max_treated', `p_max_control')

* Exclude units outside of the common support
gen common_support = p_score >= `p_common_min' & p_score <= `p_common_max'

* Re-estimate the propensity score matching model using kmatch with units within the common support
kmatch ps lowexpend if common_support, pscore(p_score)

* Store the results for later
estimates store common_support_model

*1.12

* Use esttab to create a regression table for all stored models
esttab Model1 PSM1 PSM2 common_support_model , ///
    b(3) se star(* 0.1 ** 0.05 *** 0.01) ///
    nobaselevels label booktabs nonumber ///

*1.13
* Calculate the mean outcome for the treated group
gen treated_outcome = math4 if lowexpend == 1 & common_support == 1
summarize treated_outcome, meanonly
scalar mean_treated = r(mean)

* Calculate the mean outcome for the control group within the matched pairs
gen control_outcome = math4 if lowexpend == 0 & common_support == 1
summarize control_outcome, meanonly
scalar mean_control = r(mean)

* Calculate the ATT
scalar att = mean_treated - mean_control

* Display the ATT
disp "The Average Treatment Effect on the Treated (ATT) is " att

*1.14
*Single Nearest Neighbor
kmatch ps lowexpend, pscore(p_scores) nn(1) vce(cluster dcode)
estimates store single_nn
*Five Nearest Neighbors
kmatch ps lowexpend, pscore(p_scores) nn(5) vce(cluster dcode)
estimates store five_nn
*Five Nearest Neighbors with a Caliper of 0.05
kmatch ps lowexpend, pscore(p_scores) nn(5) caliper(0.05) vce(cluster dcode)
estimates store five_nn_caliper
*Epanechnikov Kernel Function
kmatch ps lowexpend, pscore(p_scores) kernel(epan) vce(cluster dcode)
estimates store epan_kernel

*1.15
*PSM
logit lowexpend enroll
predict propensity_score, pr

*IPW
gen ipw_weights = 1/(propensity_score*(lowexpend==1) + (1-propensity_score)*(lowexpend==0))
regress math4 lowexpend [pw=ipw_weights], vce(cluster dcode)

*DR
regress math4 i.enroll if lowexpend==1
predict yhat_treated, xb
regress math4 i.enroll if lowexpend==0
predict yhat_control, xb
gen dr_estimate = . 
replace dr_estimate = yhat_control + (math4 - yhat_control)*ipw_weights if lowexpend==0
replace dr_estimate = yhat_treated + (math4 - yhat_treated)*ipw_weights if lowexpend==1
summarize dr_estimate, meanonly
display   dr_estimate


**Section 2: Constructing Simulated Data
*2.1
clear
set obs 10000
set seed 123456
gen id = _n
foreach decade in 1980 1990 2000 {
    gen measure_`decade' = runiform()
}


*2.2
set seed 123456
* 𝜆 = 5
gen city = floor(-ln(runiform()) / (1/5) + 1)
tab city, m

*2.4
set seed 123456
gen neighborhood = ceil(_n/100)
tab neighborhood

*2.5 
set seed 123456
egen unique_neighborhood = group(neighborhood)
bysort unique_neighborhood: gen temp_treated = runiform() if _n == 1
by unique_neighborhood: egen neighborhood_treated = max(temp_treated)
replace neighborhood_treated = neighborhood_treated < 0.5
drop temp_treated unique_neighborhood
tab neighborhood_treated



*2.6
egen total_neighborhoods = count(neighborhood), by(neighborhood)
egen treated_neighborhoods = total(neighborhood_treated), by(neighborhood)
gen percent_neighborhoods_treated = treated_neighborhoods / total_neighborhoods * 100
summarize percent_neighborhoods_treated, meanonly
display "Percentage of neighborhoods treated: " r(mean)

gen is_treated_individual = neighborhood_treated
summarize is_treated_individual, mean
display "Percentage of the population treated: " r(mean)*100

*2.7
gen Age = runiformint(25, 55)
gen Male = runiform(0, 1) > 0.5

*2.8
gen NeighborInc = rnormal(10000, 2500)
gen CityInc = rnormal(5000, 1000)
gen AgeInc = rgamma(Age, 1/100)
gen Income = 65000 + (5000 * Male) + NeighborInc + CityInc + AgeInc + rnormal(10000, 5000)
summarize Age Male NeighborInc CityInc AgeInc Income

*2.10
* Set the seed for reproducibility
set seed 123456

* Assume city and neighborhood identifiers are already created in the data
* Create city-specific growth rates
gen GrowCity = rnormal(0.01, 0.0025)
bysort city: egen MeanGrowCity = mean(GrowCity)

* Create neighborhood-specific growth rates
gen GrowNeigh = rnormal(0.025, 0.01)
bysort neighborhood: egen MeanGrowNeigh = mean(GrowNeigh)

* Generate growth rates for treated and untreated neighborhoods
gen GrowTreat = .
gen GrowUntreat = .
replace GrowTreat = rnormal(1.05 + MeanGrowCity + MeanGrowNeigh - 0.05, 0.05) if neighborhood_treated == 1
replace GrowUntreat = rnormal(1 + MeanGrowCity + MeanGrowNeigh, 0.01) if neighborhood_treated == 0

* Combine the growth rates into one variable for simplicity
gen GrowthRate = GrowTreat
replace GrowthRate = GrowUntreat if missing(GrowTreat)

* Generate a second measure of income adjusted for growth
gen IncomeSecondPeriod = Income * GrowthRate

* Check the results
summarize GrowCity GrowNeigh GrowTreat GrowUntreat GrowthRate Income IncomeSecondPeriod

**#
*2.11
* Duplicate the dataset for two time periods
expand 2

* Generate the year variable
gen year = 2021
by id (year), sort: replace year = 2023 if _n == 2

* Adjust age for the second period
replace Age = Age + 2 if year == 2023

* Adjust treatment to only occur in 2023
replace is_treated_individual = 0 if year == 2021
by id (year), sort: gen income_adj = Income if year == 2021
by id (year), sort: replace income_adj = IncomeSecondPeriod if year == 2023

* Check the data structure
browse

*Section 3: Panel Models, Identifying Variation, and Inference
*3.1
gen ln_income = ln(IncomeSecondPeriod)

preserve
keep if year == 2023

areg ln_income i.city i.Age i.Male, absorb(city)

restore

*3.7
* Set the panel data structure
xtset id year

* Estimate the panel regression with individual and time fixed effects
* and cluster standard errors at the individual level
xtreg ln_income is_treated_individual Male i.year, fe vce(cluster id)

*3.9
*Robust Standard Errors
xtreg ln_income i.year, fe robust
estimates store robust
*Cluster by Neighborhood
xtreg ln_income i.year, fe vce(cluster neighborhood)
estimates store cluster_neighborhood
*Two-Way Cluster by Neighborhood and Year
reghdfe ln_income, absorb(id year) vce(cluster neighborhood)
estimates store twoway_cluster
*Bootstrap by Neighborhood
bootstrap, reps(500): xtreg ln_income i.year, fe
xtreg ln_income i.year, fe vce(cluster neighborhood)

*3.10
* Fixed-effects regression with clustering by neighborhood
xtreg ln_income i.year, fe vce(cluster neighborhood)
estimates store Model_FE_ClusterN

* Compile results into a table
esttab Model_FE_ClusterN, cell("b(star) t(par)") noobs compress nogap wide label














