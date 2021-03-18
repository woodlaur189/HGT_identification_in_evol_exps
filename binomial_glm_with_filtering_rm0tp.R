#install.packages("readxl")
#install.packages("na.tools")

library(readxl)
library(na.tools)
library(MASS)

for (rep in c("partAb-"))
{
  my_data<-NULL
  my_data <- read_excel("/Users/lwoo0005/Documents/An_H_pylori/Shami_stuff/partMCE3_allvars_selection_coefs.xlsx", sheet=paste("MCE3_",rep,"_allvars_as_counts", sep=""))
  my_data <- my_data[rowSums(is.na(my_data)) != ncol(my_data), ]
  out_file <- paste("/Users/lwoo0005/Documents/An_H_pylori/Shami_stuff/partMCE3_allvars_",rep,"_rm1tp_selection_coef_all_tps_logit_link.tsv", sep="")
  rep_data<-NULL
  all_ref_data_cols = grep("_ref", colnames(my_data))
  ref_coms=grep("_ref", colnames(my_data))
  all_alt_data_cols = grep("_alt", colnames(my_data))
  all_ref_data = my_data[,all_ref_data_cols]
  for (i in 1:(length(row.names(my_data))))
  {
    a <- NULL
    x_gens<-as.numeric(unlist(strsplit(colnames(my_data[i,][ref_coms]), split="_ref")))
    x_gens<-as.numeric(unlist(strsplit(colnames(my_data[i,][ref_coms]), split="_ref")))
    df_ref<-unlist(my_data[i,][(ref_coms)])
    df_alt<-unlist(my_data[i,][(ref_coms)+1])
    df_ref_list <- as.numeric(my_data[i,][(ref_coms)])
    df_alt_list <- as.numeric(my_data[i,][(ref_coms)+1])
    poses <- c(grep('^0$',df_ref_list),grep('^0$',df_alt_list))
    num_0s <- length(poses)
    df_all<-c(df_alt,df_ref)
    if(num_0s < 1){
      gen_names<-paste0(x_gens, sep = "", collapse=", ")
      test_glm<-glm(cbind(df_alt,df_ref)~x_gens,family=binomial(link="logit"))
      if ((coef(summary(test_glm)))[2,2]>1){
        a<-data.frame(my_data[i,][1],gen_names,my_data[i,][2],my_data[i,][3],my_data[i,][4],my_data[i,][6],my_data[i,][7], NA,NA,NA,NA,NA,NA)
        colnames(a) <- c('Type','Gens', 'Gene',"Description", "Position", "S/NS", "Mutation_annotation","Coefficient","Standard_error","p-value","Lower_CI_bound","Upper_CI_bound","Model_chisqr_test")
      }
      else {
        a<-data.frame(my_data[i,][1],gen_names,my_data[i,][2],my_data[i,][3],my_data[i,][4],my_data[i,][6],my_data[i,][7],coef(summary(test_glm))[2,1],coef(summary(test_glm))[2,2],coef(summary(test_glm))[2,4],confint(test_glm)[2,1],confint(test_glm)[2,2],anova(test_glm,test="Chisq")$"Pr(>Chi)"[2])
        colnames(a) <- c('Type','Gens', 'Gene',"Description", "Position", "S/NS", "Mutation_annotation","Coefficient","Standard_error","p-value","Lower_CI_bound","Upper_CI_bound","Model_chisqr_test")
      }
    }
    if(num_0s == 1){
      x_gens <- x_gens[-poses]
      df_ref <- df_ref[-poses]
      df_alt <- df_alt[-poses]
      gen_names<-paste0(x_gens, sep = "", collapse=", ")
      test_glm<-glm(cbind(df_alt,df_ref)~x_gens,family=binomial(link="logit"))
      if (coef(summary(test_glm)) == 0 || coef(summary(test_glm)) == 1){
        a<-data.frame(my_data[i,][1],gen_names,my_data[i,][2],my_data[i,][3],my_data[i,][4],my_data[i,][6],my_data[i,][7], NA,NA,NA,NA,NA,NA)
        colnames(a) <- c('Type','Gens', 'Gene',"Description", "Position", "S/NS", "Mutation_annotation","Coefficient","Standard_error","p-value","Lower_CI_bound","Upper_CI_bound","Model_chisqr_test")
      } 
      else {
        a<-data.frame(my_data[i,][1],gen_names,my_data[i,][2],my_data[i,][3],my_data[i,][4],my_data[i,][6],my_data[i,][7],coef(summary(test_glm))[2,1],coef(summary(test_glm))[2,2],coef(summary(test_glm))[2,4],confint(test_glm)[2,1],confint(test_glm)[2,2],anova(test_glm,test="Chisq")$"Pr(>Chi)"[2])
        colnames(a) <- c('Type','Gens', 'Gene',"Description", "Position", "S/NS", "Mutation_annotation","Coefficient","Standard_error","p-value","Lower_CI_bound","Upper_CI_bound","Model_chisqr_test")
      }
    }
    if(num_0s > 1){
      gen_names<-NA
      a<-data.frame(my_data[i,][1],gen_names,my_data[i,][2],my_data[i,][3],my_data[i,][4],my_data[i,][6],my_data[i,][7], NA,NA,NA,NA,NA,NA)
      colnames(a) <- c('Type','Gens', 'Gene',"Description", "Position", "S/NS", "Mutation_annotation","Coefficient","Standard_error","p-value","Lower_CI_bound","Upper_CI_bound","Model_chisqr_test")
    }
    rep_data <- rbind(rep_data, a)
  }
  print("done")
  write.table(rep_data, out_file, sep="\t", row.names=F)
}
