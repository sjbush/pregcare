# LT 18/12/2021
# apply integrated likelihood approach to PREGCARE data

# modfied 4/01/2022
# improved numerical stability: NAs were produced for mle or for upper bound of 95% CI interval
# solved by restricting the integration interval to non underflow likelihood
# and switching to uniroot for upper CI

# modified 14/01
# new data from Anne
# add binomial CI for controls
# add upper CI when count of ALT is 0

# code for numerically integrated log likelihood
source("integratedLL.R")

# read data
pregcare = read.csv('PREGCARE-DeepSeq_04Jan22.csv')
names(pregcare)[1] = 'Index'

# unique ID is FAM.ID x gene
IDs = unique(pregcare[c("FAM.ID", "gene")])
#IDs = unique(pregcare[c("FAM.ID", "gene", "DNM_genomic_location_.GRCh38.p12.")])

# dataframe to store results
res = c()

# loop over FAM.ID x gene combination
for (iobs in 1:nrow(IDs)) {
  
  # get the data for this family x gene
  dat = subset(pregcare, FAM.ID == IDs[iobs, "FAM.ID"] & gene == IDs[iobs, "gene"])
  
  # extract counts for the controls
  controls = subset(dat, sample_type %in% paste("Control", 1:3, sep =''))
  
  # safety check, we should have 3 controls
  if (nrow(controls) != 3) {
    cat("Family ", dat$FAM.ID[1], " does not have 3 controls\n")
    next
  }
  
  ctrlREF = sum(controls$REF_Count.ALLreplicates, na.rm = TRUE)
  ctrlALT = sum(controls$ALT_Count.ALLreplicates, na.rm = TRUE)
  
  # changed dont remove the controls, just get binomial CI
  # remove controls
  #dat = subset(dat, !(sample_type %in% paste("Control", 1:3, sep ='')))
  
  # to store results for this family
  # keep only first 7 columns now 8
  famres = dat[1:8]
  famres$vaf_binom = NA
  famres$lb_binom = NA
  famres$ub_binom = NA
  famres$vaf_corrected = NA
  famres$lb_corrected = NA
  famres$ub_corrected = NA  
  
  # loop over other sample types
  for(irow in 1:nrow(dat)) {
    caseREF = dat$REF_Count.ALLreplicates[irow]
    caseALT = dat$ALT_Count.ALLreplicates[irow]
    if (!is.na(caseREF) && !is.na(caseALT) && caseREF+caseALT > 0) {
      
      # binomial test
      test_classic = binom.test(caseALT, caseREF+caseALT)
      famres$vaf_binom[irow] = test_classic$estimate
      famres$lb_binom[irow] = test_classic$conf.int[1]
      famres$ub_binom[irow] = test_classic$conf.int[2]
      
      #integrated likelihood corrected test
      # only if not a control
#      if (caseALT>0) {
      if (!(dat$sample_type[irow] %in% paste("Control", 1:3, sep =''))) {
        test_corrected = estimatep(ctrlALT, ctrlREF + ctrlALT, caseALT, caseREF + caseALT)
        famres$vaf_corrected[irow] = test_corrected[1]
        famres$lb_corrected[irow] = test_corrected[2]
        famres$ub_corrected[irow] = test_corrected[3]
      }
    }
  }
  res = rbind(res, famres)
}

write.csv(res, file='PREGCARE-DeepSeq_04Jan22_results.csv', row.names = FALSE)
