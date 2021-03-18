#install.packages("proxy")
#install.packages("spaa")
library(proxy)
bin_data_file='/Users/lwoo0005/Documents/Laura_stuff/H_py_An/Concatanated_maxHGT_p12_genomes/numberup_pangenomes/22-8-20_roary_90_evolved_50up_w_parents_split/gene_presence_absence.Rtab'
a<-read.table(bin_data_file, header=TRUE)
b<-as.data.frame.matrix(a) 
d <- b[,-1]
rownames(d) <- b[,1]
pairwise_list <- dist2list(proxy::dist(d, by_rows = FALSE, method = "Euclidean"))
write.csv(pairwise_list, "/Users/lwoo0005/Documents/Laura_stuff/H_py_An/Concatanated_maxHGT_p12_genomes/numberup_pangenomes/22-8-20_roary_90_evolved_50up_w_parents_split/gene_presence_absence_euclidean.csv")
