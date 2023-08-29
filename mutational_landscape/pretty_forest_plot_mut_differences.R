
##########################################################################################
### compare the frequency of mutations before and after radiation and make forest plot ###
##########################################################################################
library(metafor)
library(dplyr)
library(broom)

## useful links
#https://stackoverflow.com/questions/52784396/fishers-exact-test-on-rows-in-data-frame-r
#https://stackoverflow.com/questions/63528274/integrate-odds-ratio-and-ci-with-p-value-in-fisher-test-r
#http://genometoolbox.blogspot.com/2014/06/easy-forest-plots-in-r.html

cases<-read.csv(file="~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/gistic/gistic_onesample_highest_purity_May2023/output_unmarged_ZackNormalised_singlesample_highest_purity_X_plot/costum_plots/sarcoma_mutated_cases.csv", header = TRUE)

#Gene Uniradiatted.cases.mutated Uniradiated.cases.unmutated Irradiated.casese.mutated. Irrdiated.cases.unmutated
#TP53                 TP53                          6                          10                         4
#MUC16               MUC16                          4                          12                         3
#TNC                   TNC                          2                          14                         3

rownames(cases)<-cases[,1]
cases_good<-cases[,-1]

#genes<-c("TP53","MUC16","TNC","CLTC","NBEA","ATRX","COL2A1","RBM10")
#apply(cases_good, 1, 
#      function(x) {
#        tbl <- matrix(as.numeric(x[1:4]), ncol=2, byrow=T)
#        fisher.test(tbl, alternative="two.sided")$p.value
#      })

#qq_genes<-cases_good[rownames(cases_good)%in%genes,]

qq<-apply(cases_good, 1, function(x) {
  tbl <- matrix(as.numeric(x[1:4]), ncol=2, byrow=T)
  ft <- fisher.test(tbl, alternative="two.sided")
  tidy(ft) %>% 
    dplyr::select(estimate, p.value, conf.low, conf.high) }) %>% 
  bind_rows(.id = 'grp')      

qq_good<-data.frame(qq)

#qq_good[qq_good$grp=="ATRX","conf.high"]<-9.53573
#qq_good[qq_good$grp=="COL2A1","conf.high"]<-9.53573
#qq_good[qq_good$grp=="RBM10","conf.high"]<-9.53573

### make file for the saved table ###
final_table_forsupp<-cbind(cases_good,qq_good)
final_table_forsupp<-final_table_forsupp[,c("Uniradiatted.cases.mutated","Uniradiated.cases.unmutated", "Irradiated.casese.mutated","Irrdiated.cases.unmutated","estimate","p.value","conf.low",  "conf.high")]
names(final_table_forsupp)[5]<-"odds_ratio"
write.table(final_table_forsupp, file = "~/Desktop/diff_mut_info_irradiated_vs_radiated.txt", row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)

qq_good$beta <- log(qq_good$estimate)
qq_good$se   <- (log(qq_good$conf.high)-log(qq_good$conf.low))/(2*1.96)


# Assign values for plotting
labs <- qq_good$grp
yi   <- qq_good$beta
sei  <- qq_good$se

# Combine data into summary estimate
res  <- rma(yi=yi, sei=sei, method="FE")
summary(res)

forest(res, transf=exp, refline=1, xlab="Odds Ratio (95%CI)", slab=labs, mlab="Summary Estimate",col="blue")
mtext(paste("Association p-value=",summary(res)$pval),side=3, line=-1)





