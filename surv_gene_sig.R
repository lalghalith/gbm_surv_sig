# Lia Harrington 2017 - Community Detection
# scripts/surv_gene_sig.R
#
# Usage:
# Run in command line:
#
#       Rscript scripts/surv_gene_sig.R
#
# Output:
# Produces survival curves and analyses for cell TCGA GBM data

# Set Dependencies --------------------------------------------------------

library(survival)
library(gplots)

# Prepare Expression Data -------------------------------------------------

data_file = file.path("data", "data_RNA_Seq_v2_mRNA_median_Zscores.txt")
expression <- read.table(data_file, sep = "\t", 
                         header = T, row.names = 1)[,-1]

# transpose dataset so each patient a row
expression <- data.frame(t(expression))

# just grab genes in signature 
expression2 <- data.frame(expression$COL2A1, expression$DDIT4, expression$EGLN3, 
                          expression$FAM162A, expression$KDELR3, expression$LOXL2, 
                          expression$NDRG1, expression$P4HA1, expression$P4HA2, 
                          expression$PFKL, expression$PPP1R3C, expression$STC1, 
                          expression$VEGFA, expression$VLDLR, expression$VIT,
                          expression$CLYBL, expression$APLN, expression$ANKRD37,
                          expression$ERO1A)
# normalize expression
expression2 <- data.frame(scale(expression2))

new_row_names <- c()
current_names <- row.names(expression)
for (r in 1:length(current_names)){
  name <- current_names[r]
  new_row_names = c(new_row_names, paste("TCGA-", substr(name, 6,7), "-", 
                                         substr(name, 9, 12), sep = ""))
}

expression2$PATIENT_ID <- new_row_names

# average over those with multiple measurements 
expression_ag <- aggregate(.~PATIENT_ID, FUN=mean, data=expression2)

# Prepare Clinical Data ---------------------------------------------------

data_file = file.path("data", "data_bcr_clinical_data_patient.txt")
clinical <- read.table(data_file, header = T, 
                       row.names = 1, sep = "\t")

small_clinical <- data.frame(clinical$PATIENT_ID, clinical$OS_MONTHS, 
                             clinical$OS_STATUS)

small_clinical <- small_clinical[small_clinical$clinical.OS_STATUS != "[Not Available]", ]

# create binary variable if dead or alive 
Dead = rep(NA, dim(small_clinical)[1])
for (i in 1:length(small_clinical$clinical.OS_STATUS)){
  if (small_clinical$clinical.OS_STATUS[i] == "DECEASED"){
    Dead[i] = 1
  }
  else if (small_clinical$clinical.OS_STATUS[i] == "LIVING"){
    Dead[i] = 0
  }
}

small_clinical$Dead = Dead
small_clinical$clinical.OS_MONTHS <- as.integer(small_clinical$clinical.OS_MONTHS)

# Merge Data --------------------------------------------------------------

merged_data <- merge(small_clinical, expression_ag, by.x = "clinical.PATIENT_ID", 
                     by.y = "PATIENT_ID")

head(merged_data)

# Survival Analysis -------------------------------------------------------

coxfit <- coxph(Surv(clinical.OS_MONTHS, Dead) ~., data=merged_data[, c(2, 4:23)])

cox.zph(coxfit)

coefs <- as.matrix(coxfit$coefficients)

risk_score <- as.matrix((merged_data[, 5:23] - colMeans(merged_data[, 5:23]))) %*% coefs 

merged_data$risk_score <- risk_score

group <- rep(NA, dim(merged_data)[1])

for (i in 1:dim(merged_data)[1]){
  if (risk_score[i] > 0){
    group[i] = "H"
  }
  else if (risk_score[i] < 0){
    group[i] = "L"
  }
}

merged_data$group <- group

my_surv <- survfit(Surv(clinical.OS_MONTHS, Dead) ~ group, data = merged_data)
par(mfrow=c(1,1))
plot(my_surv, conf.int=F, lty = 1, col = c("red", "black"), 
     xlab = "Months of Overall Survival", ylab = "Percent Survival")
legend("bottomleft", legend=c("H", "L"), col=c("red","black"), lty=1, 
       horiz=F, bty='n')
png('gene_signature_survival.png', res = 300)
dev.off()

risk_low <- merged_data$risk_score[merged_data$group == 'L']
risk_high <- merged_data$risk_score[merged_data$group == 'H']

mean(risk_low)
mean(risk_high)

t.test(risk_low, risk_high)

survdiff(Surv(clinical.OS_MONTHS, Dead) ~ group, data = merged_data)

get_95_ci <- function(values){
  n = length(values)
  s = sd(values)
  error <- qnorm(.975)*(s/sqrt(n))
  ci <- mean(values) + c(-1,1)*error
  return(ci)
}

get_95_ci(risk_low)
get_95_ci(risk_high)

coxfit2 <- coxph(Surv(clinical.OS_MONTHS, Dead) ~ group, data=merged_data)

# get hazard ratio
hz <- exp(1)^coxfit2$coefficients

percen_reduct <- (1-hz)*100