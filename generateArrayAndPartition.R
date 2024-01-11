#!/usr/bin/env Rscript  

library(reshape2)
library(energy)
library(dplyr)

#INPUTS <----- ENTER USER DATA

#ENTER PARTITION NO - automatically adds partition no to the output name to avoid overwriting
partition_no <- 1

#ENTER WORKING DIRECTORY
setwd("/Users/shindio/Documents/Cambridge/SEVDA Lab/Network/R-script-test")

#ENTER NAME FOR PARTITIONED INPUT DATASET FILE (csv)
input_partition <- read.csv("long.csv", stringsAsFactors = FALSE) 

#ENTER NAME FOR WHOLE DATASET FILE (csv)
input_whole <- read.csv("long.csv", stringsAsFactors = FALSE) 

#ENTER DESIRED NAME FOR OUTPUT ATTRIBUTE FILE (do not include extension)
output_att_name <- "enter_name" 

#ENTER DESIRED NAME FOR OUTPUT COEFFICIENT FILE (do not include extension)
output_coeff_name <- "enter_name2" 

#ENTER WHETHER YOU WOULD LIKE AN ATTRIBUTE FILE (only select TRUE on first partition) - select FALSE if you do not want attributes at all.
attribute_file_yes <- TRUE

#ENTER NUMBER OF YOUR EXPRESSION CATEGORIES
no_of_exp_categories <- 7 

#ENTER EXPRESSION CATEGORY NAMES (must be same number as no_of_exp_categories)
exp_category_names <- c("Cyst", "J2", "10hpi", "48hpi", "12dFemale", "12dMale", "24dFemale") 

#Function that returns the distance correlation coefficient between two nonlinear data.

dcorr <- function(data = NULL, obs,  pred, tidy = FALSE, na.rm = TRUE) {
  dcorr <- rlang::eval_tidy(
    data = data,
    rlang::quo(
      energy::dcor(x = {{obs}},
                   y = {{pred}})) )
  # Tidy = FALSE
  if (tidy == FALSE) { return(dcorr) }
  # Tidy = TRUE
  if (tidy == TRUE) { return( as.data.frame(dcorr) ) }
}

#Finds the number of genes and attributes

no_of_att <- 1
effector1 <- input_partition$Gene[1]
while (effector1 == input_partition$Gene[no_of_att]){
  no_of_att <- no_of_att + 1
}
no_of_att <- no_of_att - 1

no_of_genes <- length(input_partition$Gene)/no_of_att
no_of_genes_all <- length(input_whole$Gene)/no_of_att

#Creates a list of gene names in import order for partition input

gene_name_list <- matrix(nrow = no_of_genes, ncol = 1)
for(l in 1:length(gene_name_list[,1])){
  gene_name_list[l, 1] <- input_partition$Gene[l*no_of_att]
}

#Creates a list of gene names in import order for whole input

gene_name_list_all <- matrix(nrow = no_of_genes_all, ncol = 1)
for(l in 1:length(gene_name_list_all[,1])){
  gene_name_list_all[l, 1] <- input_whole$Gene[l*no_of_att]
}

#Generates an array of vectors from input -- each vector containing the expression values at each life cycle stage, in order.

array_of_vectors <- array(dim = c(no_of_exp_categories,no_of_genes), dimnames = list(exp_category_names,gene_name_list))
count <- 1
while (count <= no_of_genes) {
  count2 <- 1
  while (count2 <= no_of_exp_categories){
    array_of_vectors[count2,count] <- as.double(input_partition$Value[no_of_att*(count-1) + count2])
    count2 <- count2+1
  }
  count <- count+1
}

#Generates an array of vectors from whole input  -- each vector containing the expression values at each life cycle stage, in order.

array_of_vectors_all <- array(dim = c(no_of_exp_categories,no_of_genes_all), dimnames = list(exp_category_names,gene_name_list_all))
count <- 1
while (count <= no_of_genes_all) {
  count2 <- 1
  while (count2 <= no_of_exp_categories){
    array_of_vectors_all[count2,count] <- as.double(input_whole$Value[no_of_att*(count-1) + count2])
    count2 <- count2+1
  }
  count <- count+1
}


#Creates a dataframe of node attributes. 

if(attribute_file_yes == TRUE){
  if(no_of_att - no_of_exp_categories> 0){
    att_matrix <- matrix(nrow = no_of_genes_all, ncol = no_of_att - no_of_exp_categories, dimnames = list(gene_name_list_all, NULL))
    for (a in 1:no_of_genes_all){ 
      for (b in 1:(no_of_att - no_of_exp_categories)){
        att_matrix[a,b] <- input_whole$Value[no_of_att*(a-1) + no_of_exp_categories + b] 
      }
    }
  }
}


#Generates a matrix containing the regression coefficients of all possible combination of effectors

coeff <- matrix(nrow = no_of_genes, ncol = no_of_genes_all, dimnames = list(gene_name_list, gene_name_list_all))
for(i in 1:no_of_genes){
  P <- array_of_vectors[,i] 
  for(j in 1:no_of_genes_all){
    O <- array_of_vectors_all[,j]
    if(gene_name_list[i] == gene_name_list_all[j]){
      coeff[i,j] <- NA
    }
    else{
      pearsons <- cor(P,O)
      if(is.na(pearsons) == FALSE && pearsons<0){
        coeff[i,j] <- -dcorr(obs = P, pred = O)
        print(c(i,j))
      }
      else{
        coeff[i,j] <- dcorr(obs = P, pred = O)
        print(c(i,j))
      }
    }
  }
}

output_coeff <- data.frame(melt(coeff))
output_coeff <- subset(output_coeff,subset = is.na(output_coeff$value) == FALSE)
output_coeff <- output_coeff %>% arrange(Var1,Var2)

#Writes the matrixes from the previous sections in long form.

output_att <- data.frame(att_matrix)
output_coeff_name <- paste0(output_coeff_name,partition_no,".txt")
output_att_name <- paste0(output_att_name, ".txt")

if(attribute_file_yes == TRUE){
  write.table(output_att, output_att_name, row.names = TRUE, col.names = FALSE, sep = "\t",quote = FALSE)
}
write.table(output_coeff, output_coeff_name, row.names = FALSE, col.names = FALSE, sep = "\t",quote = FALSE)

