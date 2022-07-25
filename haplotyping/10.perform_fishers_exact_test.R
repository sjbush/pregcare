# calculate Fisher's exact test for (e.g.) FAM44
# data for this table is obtained from the file 2x2_tables_for_pileup_inheritance_prediction.txt, produced by script 9.pileup_phaser.pl

dat <- data.frame(
  "phasingSNP_REF" = c(10, 321),
  "phasingSNP_ALT" = c(297, 2),
  row.names = c("DNM_REF", "DNM_ALT"),
  stringsAsFactors = FALSE
)
colnames(dat) <- c("phasingSNP_REF", "phasingSNP_ALT")

dat

# this produces:
#        phasingSNP_REF phasingSNP_ALT
#DNM_REF             10            297
#DNM_ALT            321              2

fisher.test(dat)

# this produces:
#        Fisher's Exact Test for Count Data
#
#data:  dat
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
# 0.000000000 0.001026331
#sample estimates:
#  odds ratio 
#0.0002404898

# a one-line alternative:
# fam44<-rbind(c(10,321),c(297,2)); fisher.test(fam44)
