#!/usr/bin/env Rscript

#THIS SCRIPT WAS WRITTEN BY DIO S. SHIN AND EDITED BY BETH MOLLOY

library(reshape2)
library(energy)
library(dplyr)

#ENTER WORKING DIRECTORY
#setwd("PATH/WORKING_DIRECTORY")

#ENTER NAME FOR INPUT DATASET FILE (.csv) 
#NOTE INPUT FILE MUST BE A LONG FORM TABLE WITH THE EXPRESSION VALUES IN THE FIRST COLUMNS
input <- read.csv("717_putative_effectors_LF.csv", stringsAsFactors = FALSE) 

#ENTER DESIRED NAME FOR OUTPUT ATTRIBUTE FILE (do not include file extensions)
output_att_name <- "output_attribute_file" 

#ENTER DESIRED NAME FOR OUTPUT COEFFICIENT FILE (do not include file extensions)
output_coeff_name <- "output_coefficient_file" 

#ENTER NUMBER OF YOUR EXPRESSION CATEGORIES (e.g. Cyst, J2 etc.)
no_of_exp_categories <- 7 

#ENTER WHETHER YOU WOULD LIKE AN ATTRIBUTE FILE - select FALSE if you do not want attributes at all
attribute_file_yes <- TRUE

#ENTER EXPRESSION CATEGORY NAMES (must be same number as no_of_exp_categories)
exp_category_names <- c("Cyst_mean", "J2_mean", "Hours10_mean", "Hours48_mean", "Days12female_mean", "Days12male_mean", "Days24female_mean") 

#FUNCTION THAT RETURNS THE DISTANCE CORRELATION COEFFICIENT BETWEEN TWO NONLINEAR DATASETS

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

#FINDS THE NUMBER OF GENES AND ATTRIBUTES

no_of_att <- 1
effector1 <- input$Gene[1]
while (effector1 == input$Gene[no_of_att]){
  no_of_att <- no_of_att + 1
}
no_of_att <- no_of_att - 1

no_of_genes <- length(input$Gene)/no_of_att

#CREATES A LIST OF GENE NAMES IN IMPORT ORDER FOR GLAND CELL DATA

gene_name_list <- matrix(nrow = no_of_genes, ncol = 1)
for(l in 1:length(gene_name_list[,1])){
  gene_name_list[l, 1] <- input$Gene[l*no_of_att]
}

#GENERATES AN ARRAY OF VECTORS FROM INPUT - each vector containing the expression values at each life cycle stage, in order

#N.B. Ignore the following error message
##Warning messages:
##1: In as.double(input$Value[no_of_att * (count - 1) + count2]) :
##NAs introduced by coercion

array_of_vectors <- array(dim = c(no_of_exp_categories,no_of_genes), dimnames = list(exp_category_names,gene_name_list))
count <- 1
while (count <= no_of_genes) {
  count2 <- 1
  while (count2 <= no_of_exp_categories){
    array_of_vectors[count2,count] <- as.double(input$Value[no_of_att*(count-1) + count2])
    count2 <- count2+1
  }
  count <- count+1
}

#CREATES A DATAFRAME OF NODE ATTRIBUTES

if(attribute_file_yes == TRUE){
  if(no_of_att - no_of_exp_categories> 0){
    att_matrix <- matrix(nrow = no_of_genes, ncol = no_of_att - no_of_exp_categories, dimnames = list(gene_name_list, NULL))
    for (a in 1:no_of_genes){ 
      for (b in 1:(no_of_att - no_of_exp_categories)){
        att_matrix[a,b] <- input$Value[no_of_att*(a-1) + no_of_exp_categories + b] 
      }
    }
  }}

#GENERATES A MATRIX CONTAINING THE REGRESSION COEFFICIENTS OF ALL POSSIBLE COMBINATIONS OF EFFECTORS

#N.B. Ignore the following error message
##Error in .dcov(x, y, index) : NA/NaN/Inf in foreign function call (arg 2)
##In addition: Warning message:
##In .arg2dist.matrix(y) : missing values not supported

coeff <- matrix(nrow = no_of_genes, ncol = no_of_genes, dimnames = list(gene_name_list, gene_name_list))
for(i in 1:no_of_genes){
  print(i)
  P <- array_of_vectors[,i] 
  for(j in i:no_of_genes){
    O <- array_of_vectors[,j]
    if(gene_name_list[i] == gene_name_list[j]){
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

#WRITES THE MATRICES FROM THE PREVIOUS SECTION IN LONG FORM

if(attribute_file_yes == TRUE){
  output_att <- data.frame(att_matrix)
  output_att_name <- paste0(output_att_name, ".txt")
  write.table(output_att, output_att_name, row.names = TRUE, col.names = FALSE, sep = "\t",quote = FALSE)
}
output_coeff_name <- paste0(output_coeff_name,".txt")


write.table(output_coeff, output_coeff_name, row.names = FALSE, col.names = FALSE, sep = "\t",quote = FALSE)
