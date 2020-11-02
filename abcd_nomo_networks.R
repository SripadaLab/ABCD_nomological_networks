library(lme4)
library(lmerTest)
library(lavaan)
library(lavaan.survey)
library(survey)
library(psych)
library(stringr)

# helper function to extract standardized coefficient estimates
# from lmer mixed-effects models
stdCoef.merMod <- function(object) {
    sdy <- sd(getME(object,"y"))
    sdx <- apply(getME(object,"X"), 2, sd)
    sc <- fixef(object)*sdx/sdy
    se.fixef <- coef(summary(object))[,"Std. Error"]
    se <- se.fixef*sdx/sdy
    return(data.frame(stdcoef=sc, stdse=se))
}

# helper function to read abcd phenotype files
read.abcd = function(file,sep="\t",skip=1,cols=NULL,descriptions=FALSE,descriptions.only=FALSE) {
    headers = names(read.table(file,sep=sep,header=T,stringsAsFactors=F)[-1,])
    if (descriptions) {
        descrip = names(read.table(file,sep=sep,header=T,stringsAsFactors=T,skip=1, na.strings = c("NA", "", " ")))
    }
    data = read.table(file,sep=sep,header=T,stringsAsFactors=T,skip=skip, na.strings = c("NA", "", " "))
    names(data) = headers
    if (!is.null(cols)) {
        data = subset(data,select=cols)
    }
    if (descriptions) {
        if (descriptions.only) {
            temp = data[1,]
            temp[1,] = descrip
            temp
        } else {
            list(data=data,descrip=descrip)
        }
    } else {
        data
    }
}

# load demographics and phenotype data
# two data files below can be requested from the authors (permitting DUA)
dat = read.csv('./ABCD_phenotypic_data_BPM_withPscores_7.15.20.csv', na.strings = c("NA", "", " "))
cbcl_dat= read.abcd('./abcd_cbcls01.txt')

# merge phenotype/demographics
cbcl_cols = c("cbcl_scr_syn_anxdep_r","cbcl_scr_syn_withdep_r","cbcl_scr_syn_somatic_r",
              "cbcl_scr_syn_social_r","cbcl_scr_syn_thought_r","cbcl_scr_syn_attention_r",
              "cbcl_scr_syn_rulebreak_r","cbcl_scr_syn_aggressive_r")
dat = merge(dat, cbcl_dat[, c('subjectkey', cbcl_cols)], by='subjectkey')
for (covariate in cbcl_cols) {
    dat = dat[!is.na(dat[, covariate]), ]
}
print(dim(dat))


# scale sample where to assist lmer convergence (where applicable)
dat['acs_raked_propensity_score_scaled'] = dat$acs_raked_propensity_score/1000
# set reference levels for categorical variables
dat = within(dat, gender <- relevel(gender, ref = 'M'))
dat = within(dat, RaceEthnicity <- relevel(RaceEthnicity, ref = 'White'))
dat = within(dat, pcgeduc_recode <- relevel(pcgeduc_recode, ref = 'zLess than High School'))
dat = within(dat, HouseholdMaritalStatus <- relevel(HouseholdMaritalStatus, ref = 'yes'))
dat = within(dat, HouseholdIncome <- relevel(HouseholdIncome, ref = '[>=100K]'))

# controls and nesting structure for mixed-effects models
controls = c('C(gender)', 'C(RaceEthnicity)', 'C(pcgeduc_recode)',
             'C(HouseholdMaritalStatus)', 'C(HouseholdIncome)', 'Age')
nested_str = ' + (1| site_id_l/rel_family_id)'


# lavaan modelling needs categorical variables to be coded as dummy variables 
race_dummy = dummy.code(dat$RaceEthnicity, na.rm=TRUE)[, c("Hispanic", "Black", "Other", "Asian")]
covariates = cbind('gender_dummy' = dummy.code(dat$gender, na.rm=TRUE)[, c('F')],
                   'race_dummy_1' = race_dummy[, 1], 
                   'race_dummy_2' = race_dummy[, 2],
                   'race_dummy_3' = race_dummy[, 3],
                   'race_dummy_4' = race_dummy[, 4])
educ_dummy = dummy.code(dat$pcgeduc_recode, na.rm=TRUE)[, c('College Degree', 'Masters or Professional Degree', 'Some College', 'Associates or Occupational Degree', 'High School Degree or Equivalent')]
covariates = cbind(covariates, 
                   'educ_dummy_1' = educ_dummy[, 1],
                   'educ_dummy_2' = educ_dummy[, 2],
                   'educ_dummy_3' = educ_dummy[, 3],
                   'educ_dummy_4' = educ_dummy[, 4],
                   'educ_dummy_5' = educ_dummy[, 5])
covariates = cbind(covariates, 
                   'marital_dummy' = dummy.code(dat$HouseholdMaritalStatus, na.rm=TRUE)[, c('no')])
income_dummy = dummy.code(dat$HouseholdIncome, na.rm=TRUE)[, c('[<50K]', '[>=50K & <100K]')]
covariates = cbind(covariates, 
                   'income_dummy_1' = income_dummy[, 1],
                   'income_dummy_2' = income_dummy[, 2])
dat = cbind(dat, covariates)


# drop rows with missing values or duplicated subject-key's
dat = dat[!is.na(dat$rel_family_id), ]
for (covariate in colnames(covariates)) {
    dat = dat[!is.na(dat[, covariate]), ]
}
dat = dat[!is.na(dat[, 'acs_raked_propensity_score_scaled']), ]
dat = dat[!duplicated(dat$subjectkey), ]
print(dim(dat))


# declare lavvan formula for modelling P/INT/EXT-Factors
bifactor_formula_raw = '
P =~ cbcl_scr_syn_anxdep_r + cbcl_scr_syn_withdep_r + cbcl_scr_syn_somatic_r + cbcl_scr_syn_social_r + cbcl_scr_syn_thought_r + cbcl_scr_syn_attention_r + cbcl_scr_syn_rulebreak_r + cbcl_scr_syn_aggressive_r
int =~ cbcl_scr_syn_withdep_r + cbcl_scr_syn_somatic_r + cbcl_scr_syn_anxdep_r
ext =~ 1*cbcl_scr_syn_rulebreak_r + 1*cbcl_scr_syn_aggressive_r
P ~~ 0*int
P ~~ 0*ext
int ~~ 0*ext

DEPENDENT ~ PHENOTYPE + gender_dummy + race_dummy_1 + race_dummy_2 + race_dummy_3 + race_dummy_4 + educ_dummy_1 + educ_dummy_2 + educ_dummy_3 + educ_dummy_4 + educ_dummy_5 + marital_dummy + income_dummy_1 + income_dummy_2
'


# declare list of phenotypes we're interested in modelling (as seen in tables)
phenotypes_nonksads = c('pps_y_ss_severity_score', 'bpm_t_scr_attention_t', 'bpm_t_scr_internal_t', 'bpm_t_scr_external_t', 'bis_y_ss_bis_sum', 'bis_y_ss_bas_rr', 'bis_y_ss_bas_drive', 'bis_y_ss_bas_fs', 'upps_y_ss_negative_urgency', 'upps_y_ss_positive_urgency', 'upps_y_ss_lack_of_planning', 'upps_y_ss_lack_of_perseverance', 'upps_y_ss_sensation_seeking', 'G_Dave', 'srpf_y_ss_ses', 'srpf_y_ss_iiss', 'srpf_y_ss_dfs', 'psb_y_ss_mean', 'psb_p_ss_mean', 'fes_p_ss_fc_pr', 'fes_y_ss_fc_pr', 'pmq_y_ss_mean')
phenotypes_ksads = c('Any_DepDx', 'ANY_SUISH', 'ANY_GAD', 
                     'Any_SocAnx', 'Any_SepAnx', 'Any_Phobia',
                     'Any_PanicAgorDx', 'Any_BipolarDx', 'Any_ADHD',
                     'Any_CD', 'Any_ODD')
phenotypes_demo = c('fh_parent_alcdrug_binary', 'FHtotal', 'adi_weightavg_kind')


# prepare a new dataframe for cbcl-internalizing/externalizing modelling 
dat2 = data.frame(dat)
dat2['Zcbcl_scr_syn_internal_t'] = scale(dat$cbcl_scr_syn_internal_t)
dat2['Zcbcl_scr_syn_external_t'] = scale(dat$cbcl_scr_syn_external_t)
dat2['Any_DepDx'] = as.factor(dat[['Any_DepDx']])
dat2['ANY_SUISH'] = as.factor(dat[['ANY_SUISH']])
dat2['ANY_GAD'] = as.factor(dat[['ANY_GAD']])
dat2['Any_SocAnx'] = as.factor(dat[['Any_SocAnx']])
dat2['Any_SepAnx'] = as.factor(dat[['Any_SepAnx']])
dat2['Any_Phobia'] = as.factor(dat[['Any_Phobia']])
dat2['Any_PanicAgorDx'] = as.factor(dat[['Any_PanicAgorDx']])
dat2['Any_BipolarDx'] = as.factor(dat[['Any_BipolarDx']])
dat2['Any_ADHD'] = as.factor(dat[['Any_ADHD']])
dat2['Any_CD'] = as.factor(dat[['Any_CD']])
dat2['Any_ODD'] = as.factor(dat[['Any_ODD']])
dat2 = within(dat2, Any_DepDx <- relevel(Any_DepDx, ref = '0'))
dat2 = within(dat2, ANY_SUISH <- relevel(ANY_SUISH, ref = '0'))
dat2 = within(dat2, ANY_GAD <- relevel(ANY_GAD, ref = '0'))
dat2 = within(dat2, Any_SocAnx <- relevel(Any_SocAnx, ref = '0'))
dat2 = within(dat2, Any_SepAnx <- relevel(Any_SepAnx, ref = '0'))
dat2 = within(dat2, Any_Phobia <- relevel(Any_Phobia, ref = '0'))
dat2 = within(dat2, Any_PanicAgorDx <- relevel(Any_PanicAgorDx, ref = '0'))
dat2 = within(dat2, Any_BipolarDx <- relevel(Any_BipolarDx, ref = '0'))
dat2 = within(dat2, Any_ADHD <- relevel(Any_ADHD, ref = '0'))
dat2 = within(dat2, Any_CD <- relevel(Any_CD, ref = '0'))
dat2 = within(dat2, Any_ODD <- relevel(Any_ODD, ref = '0'))


# declare the phenotype set we're interested in modelling
# 'phenotypes' below can be any one of the three list of phenotypes defined above
phenotypes = phenotypes_nonksads  #phenotypes_ksads  #phenotypes_demo
# dependents to consider (P, int, ext come from lavaan models)
dependents = c('P', 'int', 'ext', 
               'cbcl_scr_syn_internal_t', 'cbcl_scr_syn_external_t')
# column names of output (should not need to be changed)
default_res_names = c('dependent', 'phenotype', 'raw_beta', 'raw_pvalue', 'raw_lower_ci', 'raw_upper_ci', 'std_beta', 'std_pvalue', 'std_lower_ci', 'std_upper_ci')
n_done = 0
full_res = c()
fdr_raw_pvals = c()
fdr_std_pvals = c()
for (dependent in dependents) {
    dep_raw_pvals = c()
    dep_std_pvals = c()
    for (phenotype in phenotypes) {
        # if the dependent is cbcl-int/ext, then lavaan modelling is not necessary
        # directly construct mixed-effects models using predefined factor scores
        if ((dependent == 'cbcl_scr_syn_internal_t') | (dependent == 'cbcl_scr_syn_external_t')) {
            fmla = paste(dependent, ' ~ ', phenotype, ' + ', 
                         paste(controls[controls!='C(gender)'], collapse=' + '),
                         nested_str, sep='')
            model = lmer(formula=fmla, data=dat2, 
                         na.action='na.exclude', 
                         weights = acs_raked_propensity_score_scaled)
            std_model = stdCoef.merMod(model)
            raw_se = coef(summary(model))[phenotype, 'Std. Error']
            std_se = std_model[phenotype, 'stdse']
            
            lavaan_raw_beta = fixef(model)[[phenotype]]
            lavaan_std_beta = std_model[phenotype, 'stdcoef']
            lavaan_pval_raw = coef(summary(model))[phenotype, 'Pr(>|t|)']
            lavaan_pval_std = coef(summary(model))[phenotype, 'Pr(>|t|)']
            lavaan_lowerci_raw = lavaan_raw_beta - (qnorm(0.975)*raw_se)
            lavaan_lowerci_std = lavaan_std_beta - (qnorm(0.975)*std_se)
            lavaan_upperci_raw = lavaan_raw_beta + (qnorm(0.975)*raw_se)
            lavaan_upperci_std = lavaan_std_beta + (qnorm(0.975)*std_se)
            
        } else {
            # use lavaan to construct P, INT, EXT factors and compute regression weights
            bifactor_formula = str_replace(bifactor_formula_raw, 'PHENOTYPE', phenotype)
            bifactor_formula = str_replace(bifactor_formula, 'DEPENDENT', dependent)
            model = cfa(bifactor_formula, 
                        data=dat, std.lv = TRUE) 
            design = svydesign(ids=~rel_family_id,
                               strata=~site_id_l, 
                               nest=TRUE, 
                               weights=~acs_raked_propensity_score_scaled,
                               data=dat) 
            survey.model = lavaan.survey(lavaan.fit = model, survey.design = design)
            std_solution = standardizedSolution(survey.model)
            result = (std_solution[which(std_solution$lhs==dependent & std_solution$op == '~' & std_solution$rhs == phenotype), c('est.std', 'pvalue', 'ci.lower', 'ci.upper')])
            unstd_solution = parameterEstimates(survey.model)
            result_2 = (unstd_solution[which(unstd_solution$lhs==dependent & unstd_solution$op == '~' & unstd_solution$rhs == phenotype), c('est', 'pvalue', 'ci.lower', 'ci.upper')])
            
            # result_2 is unstandardized model, and result is standardized
            lavaan_raw_beta = result_2['est']
            lavaan_std_beta = result['est.std']
            lavaan_pval_raw = result_2['pvalue']
            lavaan_pval_std = result['pvalue']
            lavaan_lowerci_raw = result_2['ci.lower']
            lavaan_lowerci_std = result['ci.lower']
            lavaan_upperci_raw = result_2['ci.upper']
            lavaan_upperci_std = result['ci.upper']
        }
        # gather results to match default_res_names above
        res = cbind(dependent,
                    phenotype,
                    lavaan_raw_beta,
                    lavaan_pval_raw,
                    lavaan_lowerci_raw,
                    lavaan_upperci_raw,
                    lavaan_std_beta,
                    lavaan_pval_std,
                    lavaan_lowerci_std,
                    lavaan_upperci_std)
        colnames(res) = default_res_names
        full_res = rbind(full_res, res)
        
        dep_raw_pvals = c(dep_raw_pvals, lavaan_pval_raw)
        dep_std_pvals = c(dep_std_pvals, lavaan_pval_std)
        n_done = n_done + 1
        print(paste('(', n_done, '/', length(phenotypes)*length(dependents), ')   Done ', phenotype, ' - ', dependent, '...', sep=''))
    }
    # adjust pvalues using FDR
    fdr_raw_pvals = c(fdr_raw_pvals, p.adjust(dep_raw_pvals, method='fdr'))
    fdr_std_pvals = c(fdr_std_pvals, p.adjust(dep_std_pvals, method='fdr'))
}
colnames(full_res) = default_res_names
# include FDR adjust pvalues
full_res[, 'FDR_raw_pvalue'] = fdr_raw_pvals
full_res[, 'FDR_std_pvalue'] = fdr_std_pvals
# save results as csv
write.csv(full_res, file='./nonksads_lavaan_sample_weights_FDR.csv', quote=FALSE, row.names=FALSE)









# similar modeling for demographics table
race_dummy = dummy.code(dat$RaceEthnicity, na.rm=TRUE)[, c("Hispanic", "Black", "Other", "Asian")]
covariates = cbind('isFemale' = dummy.code(dat$gender, na.rm=TRUE)[, c('F')],
                   'Hispanic' = race_dummy[, 1], 
                   'Black' = race_dummy[, 2],
                   'Other' = race_dummy[, 3],
                   'Asian' = race_dummy[, 4])
educ_dummy = dummy.code(dat$pcgeduc_recode, na.rm=TRUE)[, c('College Degree', 'Masters or Professional Degree', 'Some College', 'Associates or Occupational Degree', 'High School Degree or Equivalent')]
covariates = cbind(covariates, 
                   'College' = educ_dummy[, 1],
                   'Masters' = educ_dummy[, 2],
                   'Some_College' = educ_dummy[, 3],
                   'Associates' = educ_dummy[, 4],
                   'High_School' = educ_dummy[, 5])
covariates = cbind(covariates, 
                   'Marital_No' = dummy.code(dat$HouseholdMaritalStatus, na.rm=TRUE)[, c('no')])
income_dummy = dummy.code(dat$HouseholdIncome, na.rm=TRUE)[, c('[<50K]', '[>=50K & <100K]')]
covariates = cbind(covariates, 
                   'income_lt_50k' = income_dummy[, 1],
                   'income_bw_50k_100k' = income_dummy[, 2])
dat = cbind(dat, covariates)
cols = c("cbcl_scr_syn_anxdep_r","cbcl_scr_syn_withdep_r","cbcl_scr_syn_somatic_r",
         "cbcl_scr_syn_social_r","cbcl_scr_syn_thought_r","cbcl_scr_syn_attention_r",
         "cbcl_scr_syn_rulebreak_r","cbcl_scr_syn_aggressive_r")

# drop missing values 
dat = dat[!is.na(dat$rel_family_id), ]
for (covariate in colnames(covariates)) {
    dat = dat[!is.na(dat[, covariate]), ]
}
for (covariate in cols) {
    dat = dat[!is.na(dat[, covariate]), ]
}
dat = dat[!is.na(dat[, 'acs_raked_propensity_score_scaled']), ]
print(dim(dat))

# lavaan definitions
# change the dependent in the last line of string below (current P) to generate
# results for INT/EXT factors
bifactor_formula_raw = '
P =~ cbcl_scr_syn_anxdep_r + cbcl_scr_syn_withdep_r + cbcl_scr_syn_somatic_r + cbcl_scr_syn_social_r + cbcl_scr_syn_thought_r + cbcl_scr_syn_attention_r + cbcl_scr_syn_rulebreak_r + cbcl_scr_syn_aggressive_r
int =~ cbcl_scr_syn_withdep_r + cbcl_scr_syn_somatic_r + cbcl_scr_syn_anxdep_r
ext =~ 1*cbcl_scr_syn_rulebreak_r + 1*cbcl_scr_syn_aggressive_r
P ~~ 0*int
P ~~ 0*ext
int ~~ 0*ext

P ~ isFemale + Hispanic + Black + Other + Asian + College + Masters + Some_College + Associates + High_School + Marital_No + income_lt_50k + income_bw_50k_100k + Age
'
model = cfa(bifactor_formula_raw, 
            data=dat, std.lv = TRUE) 
design = svydesign(ids=~rel_family_id,
                   strata=~site_id_l, 
                   nest=TRUE, 
                   weights=~acs_raked_propensity_score_scaled,
                   data=dat) 
survey.model = lavaan.survey(lavaan.fit = model, survey.design = design)

anova(survey.model)

write.csv(parameterEstimates(survey.model), 
          file='./table1_demographics_nonstandardized_coefficients_P.csv', quote=FALSE, row.names=FALSE)
write.csv(standardizedSolution(survey.model), 
          file='./table1_demographics_standardized_coefficients_P.csv', quote=FALSE, row.names=FALSE)
write.table(summary(survey.model),
            file='./table1_demographics_summary_P.csv', quote=FALSE, sep=',')

