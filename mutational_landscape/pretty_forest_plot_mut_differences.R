
##########################################################################################
### compare the frequency of mutations before and after radiation and make forest plot ###
##########################################################################################
library(metafor)
library(dplyr)
library(ggplot2)
library(broom)  ## The broom package takes the messy output of built-in functions in R, such as lm, nls, or t.test, and turns them into tidy tibbles


cases<-read.csv(file="~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/gistic/gistic_onesample_highest_purity_May2023/output_unmarged_ZackNormalised_singlesample_highest_purity_X_plot/costum_plots/sarcoma_mutated_cases.csv", header = TRUE)
#head(cases)
#Gene Uniradiatted.cases.mutated Uniradiated.cases.unmutated Irradiated.casese.mutated. Irrdiated.cases.unmutated
#TP53                 TP53                          6                          10                         4
#MUC16               MUC16                          4                          12                         3


rownames(cases)<-cases[,1]
cases_good<-cases[,-1]

# For each row of dataset build a matrix (i.e. a 2 x 2 contingency table) and pass this matrix to the fisher.test command.
mat_fisher<-apply(cases_good, 1, function(x) {
  tbl <- matrix(as.numeric(x[1:4]), ncol=2, byrow=T)   ### make 2*2 matirx for each row
  ft <- fisher.test(tbl, alternative="two.sided")      ### do fisher test for each matrix
  tidy(ft) %>%                                         ### tidy functions belongs to "groom" package to tide up fisher results
    dplyr::select(estimate, p.value, conf.low, conf.high) }) %>%   ### estimate is the odds ratio value
    bind_rows(.id = 'grp')      

mat_fisher_good<-data.frame(mat_fisher)


### make file for the saved table ###
final_table_forsupp<-cbind(cases_good,mat_fisher_good)
final_table_forsupp<-final_table_forsupp[,c("Uniradiatted.cases.mutated","Uniradiated.cases.unmutated", "Irradiated.casese.mutated","Irrdiated.cases.unmutated","estimate","p.value","conf.low",  "conf.high")]
names(final_table_forsupp)[5]<-"odds_ratio"

### calulate adjusted p-value
final_table_forsupp$Adj_pvalue<-p.adjust(final_table_forsupp$p.value, method = "BH")
write.table(final_table_forsupp, file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/gistic/gistic_onesample_highest_purity_May2023/output_unmarged_ZackNormalised_singlesample_highest_purity_X_plot/costum_plots/diff_mutation_ratio_fisher_test_irradiated_vs_radiated.txt", row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)


### make the forest plot 
### reference from this website: http://genometoolbox.blogspot.com/2014/06/easy-forest-plots-in-r.html
mat_fisher_good$beta <- log(mat_fisher_good$estimate)
mat_fisher_good$se   <- (log(mat_fisher_good$conf.high)-log(mat_fisher_good$conf.low))/(2*1.96)

# Assign values for plotting
labs <- mat_fisher_good$grp
yi   <- mat_fisher_good$beta
sei  <- mat_fisher_good$se

# Combine data into summary estimate
res  <- rma(yi=yi, sei=sei, method="FE")
summary(res)

forest(res, transf=exp, refline=1, xlab="Odds Ratio (95%CI)", slab=labs, mlab="Summary Estimate",col="blue")
mtext(paste("Association p-value=",summary(res)$pval),side=3, line=-1)



#### alternative, quick way of plotting with ggplot ()
### reference link: https://rpubs.com/mbounthavong/forest_plots_r

ggplot(mat_fisher_good, aes(y = grp, x = estimate)) +
  geom_point(shape = 18, size = 5) +  
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.25) +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) 
  
#############
#### END ####
#############
## reference links
#https://stackoverflow.com/questions/63528274/integrate-odds-ratio-and-ci-with-p-value-in-fisher-test-r

