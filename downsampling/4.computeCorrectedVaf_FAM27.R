# code for numerically integrated log likelihood
source("C:/Users/sbush/Desktop/pregcare/integratedLL.R")

# read data
pregcare = read.csv('C:/Users/sbush/Desktop/pregcare/FAM27-downsampled.csv')
names(pregcare)[1] = 'Index'

# unique ID is group x avg. depth x seed
IDs = unique(pregcare[c("avg_depth_to_which_reads_downsampled", "seed")])

# dataframe to store results
res = c()

# loop over avg. depth x seed combination
for (iobs in 1:nrow(IDs)) {
  
  # get the data for this avg. depth x seed
  dat = subset(pregcare, avg_depth_to_which_reads_downsampled == IDs[iobs, "avg_depth_to_which_reads_downsampled"] & seed == IDs[iobs, "seed"])
  
  # extract counts for the controls
  controls = subset(dat, grepl("-Co",dat$group))
  #print(nrow(controls))
  
  # safety check, we should have 9 controls (3 replicates of each control sample)
  # note that for FAM27, one of the controls failed
  if (nrow(controls) != 9) {
    cat("Family ", dat$group[1], " does not have 9 controls\n")
    next
  }
  
  # the ctrlREF and ctrlALT choices are FAM-specific
  # FAM27 is a delC mutation
  ctrlREF = sum(controls$count_hq_C.ALLreplicates, na.rm = TRUE)
  ctrlALT = sum(controls$count_hq_DEL.ALLreplicates, na.rm = TRUE)
  
  # to store results for this family
  # keep only first 7 columns now 8
  famres = dat # [1:8]
  famres$vaf_binom = NA
  famres$lb_binom = NA
  famres$ub_binom = NA
  famres$vaf_corrected = NA
  famres$lb_corrected = NA
  famres$ub_corrected = NA  
  
  # loop over other sample types
  for(irow in 1:nrow(dat)) {
    
	# the caseREF and caseALT choices are FAM-specific
	# FAM27 is a delC mutation 
	caseREF = dat$count_hq_C.ALLreplicates[irow]
    caseALT = dat$count_hq_DEL.ALLreplicates[irow]
    if (!is.na(caseREF) && !is.na(caseALT) && caseREF+caseALT > 0) {
      
      # binomial test
      test_classic = binom.test(caseALT, caseREF+caseALT)
      famres$vaf_binom[irow] = test_classic$estimate
      famres$lb_binom[irow] = test_classic$conf.int[1]
      famres$ub_binom[irow] = test_classic$conf.int[2]
      
      #integrated likelihood corrected test
      # only if not a control
      if (!(dat$group[irow] %in% paste("-Co", 1:3, sep =''))) {
        test_corrected = estimatep(ctrlALT, ctrlREF + ctrlALT, caseALT, caseREF + caseALT)
        famres$vaf_corrected[irow] = test_corrected[1]
        famres$lb_corrected[irow] = test_corrected[2]
        famres$ub_corrected[irow] = test_corrected[3]
      }
    }
  }
  res = rbind(res, famres)
}

write.csv(res, file='C:/Users/sbush/Desktop/pregcare/FAM27-downsampled-corrected.csv', row.names = FALSE)