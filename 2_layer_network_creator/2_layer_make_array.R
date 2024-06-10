#!/usr/bin/env Rscript

#THIS SCRIPT WAS WRITTEN BY JONATHAN LONG AND BUILDS ON THE SCRIPT BY DIO S. SHIN, EDITED BY BETH MOLLOY

library(reshape2)
library(energy)
library(tidyverse)

#TAKES TWO .csv FILES CONTAINING GENES, EXPRESSION DATA AND ATTRIBUTES FOR A "MAJOR" AND "MINOR" SET OF GENES
#BOTH THE MAJOR AND MINOR INPUT FILES SHOULD HAVE THE SAME COLUMN HEADERS, IN THE SAME ORDER

#UNLIKE THE 1 LAYER NETWORK CREATOR, INPUT FILES CAN BE IN WIDE FORM INSTEAD OF LONG FORM
#THE USER CAN ALSO SPECIFY IF THEY WOULD LIKE THE FINAL NETWORK TO COLOUR DIFFERENT GENE TYPES (i.e. nematode TFs and effectors) BY DIFFERENT ATTRIBUTES

#ENTER WORKING DIRECTORY
#setwd("PATH/WORKING_DIRECTORY")

#TAKE A .csv FILE CONTAINING ALL GENES OF THE "MAJOR" TYPE I.E. INPUT 1 (These will share edges between one another)

Input_1_wide <- read.csv("putative_effectors.csv", stringsAsFactors = FALSE)

#TAKE A .csv FILE CONTAINING ALL GENES OF THE "MINOR" TYPE I.E. INPUT 2 (These will only share edges with genes of the major type) 

Input_2_wide <- read.csv("plant_genes_immune_attribute.csv", stringsAsFactors = FALSE)

#SPECIFY THE ATTRIBUTE THAT EACH GENE TYPE WILL BE COLOURED BY

colourby_1 <- "Supercluster"
colourby_2 <- "Immunity"  

#ENTER DESIRED NAME FOR OUTPUT ATTRIBUTE FILE (do not include file ext)
output_att_name <- "output_attribute_file_0.975" 

#ENTER DESIRED NAME FOR OUTPUT COEFFICIENT FILE (do not include file ext)
output_coeff_name <- "output_coefficient_file_0.975" 

#ENTER NUMBER OF YOUR EXPRESSION CATEGORIES
no_of_exp_categories <- 5 

#ENTER WHETHER YOU WOULD LIKE AN ATTRIBUTE FILE - select FALSE if you do not want attributes at all.
attribute_file_yes <- TRUE

#ENTER EXPRESSION CATEGORY NAMES (must be same number as no_of_exp_categories)
exp_category_names <- c("Hours10_mean", "Hours48_mean", "Days12female_mean", "Days12male_mean", "Days24female_mean")

#ENTER EDGE THRESHOLD
threshold <- 0.975 

#CREATES THE colour_by ATTRIBUTE, WHICH CAN BE USED TO COLOUR DIFFERENT GENE TYPES IN THE FINAL NETWORK BASED ON THEIR ATTRIBUTES
#FINDS THE COLUMN NUMBER OF THE CHOSEN ATTRIBUTE

column_no_1 <- match(colourby_1,colnames(Input_1_wide))

#ADDS THE COLUMN OF THE CHOSEN ATTRIBUTE TO THE TABLE UNDER THE HEADING colour_by

Input_1_wide <- Input_1_wide %>% mutate(Colour_by = Input_1_wide[,column_no_1])

#REPEAT FOR THE MINOR GENE TYPE

column_no_2 <- match(colourby_2,colnames(Input_2_wide))

Input_2_wide <- Input_2_wide %>% mutate(Colour_by = Input_2_wide[,column_no_2])

#CONVERT WIDE FORM DATA TO LONG FORM

Input_1 <- melt(Input_1_wide, na.rm = FALSE, id.vars="Gene_ID", variable.name = "Attribute", value.name = "Value") %>% 
  arrange(Gene_ID)
Input_2 <- melt(Input_2_wide, na.rm = FALSE, id.vars="Gene_ID", variable.name = "Attribute", value.name = "Value") %>% 
  arrange(Gene_ID)

#FUNCTION THAT RETURNS THE DISTANCE CORRELATION COEFFICIENT BETWEEN TWO NONLINEAR DATA 

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

#COMPUTE THE NUMBER OF GENES FOR BOTH INPUTS

#FIRST FOR INPUT 1

no_of_att <- 1
transcript <- Input_1$Gene[1]
while (transcript == Input_1$Gene[no_of_att]){
  no_of_att <- no_of_att + 1
}
no_of_att <- no_of_att - 1

no_of_genes1 <- length(Input_1$Gene)/no_of_att

no_of_genes1

#THEN FOR INPUT 2

no_of_genes2 <- length(Input_2$Gene)/no_of_att

no_of_genes2

#2 SEPERATE ARRAYS OF VECTORS GENERATED

#GENE NAME LIST FOR INPUT 1

gene_name_list1 <- matrix(nrow = no_of_genes1, ncol = 1)
for(l in 1:length(gene_name_list1[,1])){
  gene_name_list1[l, 1] <- Input_1$Gene[l*no_of_att]
}

#GENE NAME LIST FOR INPUT 2

gene_name_list2 <- matrix(nrow = no_of_genes2, ncol = 1)
for(l in 1:length(gene_name_list2[,1])){
  gene_name_list2[l, 1] <- Input_2$Gene[l*no_of_att]
}

#ARRAY OF VECTORS REPRESENTING INPUT 1

aov1 <- array(dim = c(no_of_exp_categories,no_of_genes1), dimnames = list(exp_category_names,gene_name_list1))
count <- 1
while (count <= no_of_genes1) {
  count2 <- 1
  while (count2 <= no_of_exp_categories){
    aov1[count2,count] <- as.double(Input_1$Value[no_of_att*(count-1) + count2])
    count2 <- count2+1
  }
  count <- count+1
}

#ARRAY OF VECTORS REPRESENTING INPUT 2

aov2 <- array(dim = c(no_of_exp_categories,no_of_genes2), dimnames = list(exp_category_names,gene_name_list2))
count <- 1
while (count <= no_of_genes2) {
  count2 <- 1
  while (count2 <= no_of_exp_categories){
    aov2[count2,count] <- as.double(Input_2$Value[no_of_att*(count-1) + count2])
    count2 <- count2+1
  }
  count <- count+1
}

#CREATE 2 DATAFRAMES OF NODE ATTRIBUTES

if(attribute_file_yes == TRUE){
  if(no_of_att - no_of_exp_categories> 0){
    att_matrix1 <- matrix(nrow = no_of_genes1, ncol = no_of_att - no_of_exp_categories, dimnames = list(gene_name_list1, NULL))
    for (a in 1:no_of_genes1){ 
      for (b in 1:(no_of_att - no_of_exp_categories)){
        att_matrix1[a,b] <- Input_1$Value[no_of_att*(a-1) + no_of_exp_categories + b] 
      }
    }
    att_matrix2 <- matrix(nrow = no_of_genes2, ncol = no_of_att - no_of_exp_categories, dimnames = list(gene_name_list2, NULL))
    for (a in 1:no_of_genes2){ 
      for (b in 1:(no_of_att - no_of_exp_categories)){
        att_matrix2[a,b] <- Input_2$Value[no_of_att*(a-1) + no_of_exp_categories + b] 
      }
    }
  }}

#ITERATE THROUGH INPUT 1, COMPUTING ALL INTERNAL EDGES

coeff1s <- matrix(nrow = no_of_genes1, ncol = no_of_genes1, dimnames = list(gene_name_list1, gene_name_list1))
print(no_of_genes1)

for(i in 1:no_of_genes1){
  print(i)
  P <- aov1[,i] 
  for(j in i:no_of_genes1){
    O <- aov1[,j]
    if(gene_name_list1[i] == gene_name_list1[j]){
      coeff1s[i,j] <- NA  #The network has no self edges
    }
    else{
      pearsons <- cor(P,O)
      if(is.na(pearsons) == FALSE && pearsons<0){
        coeff1s[i,j] <- -dcorr(obs = P, pred = O)
      }
      else{
        coeff1s[i,j] <- dcorr(obs = P, pred = O)
      }
    }
  }
}

#ITERATE THROUGH INPUT 2, COMPUTING ALL EDGES BETWEEN EACH ELEMENT OF INPUT 2 AND ALL ELEMENT OF INPUT 1
#IN THE RESULTING MATRIX, THE ROWS ARE THE GENES OF "MINOR" TYPE AND THE COLUMNS ARE THE GENES OF THE "MAJOR" TYPE

coeff12s <- matrix(nrow = no_of_genes2, ncol = no_of_genes1, dimnames = list(gene_name_list2, gene_name_list1))
print(no_of_genes2)

for(i in 1:no_of_genes2){
  print(i)
  P <- aov2[,i] 
  for(j in 1:no_of_genes1){
    O <- aov1[,j]
    pearsons <- cor(P,O)
    if(is.na(pearsons) == FALSE && pearsons<0){
      coeff12s[i,j] <- -dcorr(obs = P, pred = O)
    }
    else{
      coeff12s[i,j] <- dcorr(obs = P, pred = O)
    }
  }
}

#CONVERT INTO LONGFORM AND FILTER BY THE CUTOFF THRESHOLD

output_coeff1s <- data.frame(melt(coeff1s))
output_coeff1s <- subset(output_coeff1s,subset = is.na(output_coeff1s$value) == FALSE) %>% arrange(Var1,Var2) %>% 
  filter(-threshold > value | value > threshold)


output_coeff12s <- data.frame(melt(coeff12s))
output_coeff12s <- subset(output_coeff12s,subset = is.na(output_coeff12s$value) == FALSE) %>% arrange(Var1,Var2) %>% 
  filter(-threshold > value | value > threshold)

#WRITE THE FILE HERE IF NEEDED FOR STATISTICS

write.csv(output_coeff12s, "stats.csv")

#MAKE 2 NEW NODE LISTS, KEEPING ONLY GENES WITH AT LEAST ONE EDGE

#ALL THE COLUMNS CONTAINING "MAJOR" GENES INVOLVED IN EDGES 

nodesa <- data.frame(output_coeff1s$Var1) %>%  rename("gene" = output_coeff1s.Var1)
nodesb <- data.frame(output_coeff1s$Var2) %>%  rename("gene" = output_coeff1s.Var2)
nodesc <- data.frame(output_coeff12s$Var2) %>%  rename("gene" = output_coeff12s.Var2)

#THE COLUMN CONTAINING "MINOR" GENES INVOLVED IN EDGES

nodesd <- data.frame(output_coeff12s$Var1) %>%  rename("gene" = output_coeff12s.Var1)

#COUNT THE NUMBER OF UNIQUE "MAJOR" GENES CONNECTED TO THE LEAST ONE EDGE

new_gene_name_list1 <- unique(rbind(nodesa,nodesb,nodesc))
print(no_of_genes1)
print(length(new_gene_name_list1$gene))
new_no_of_genes1 <- length(new_gene_name_list1$gene)

#COUNT THE NUMBER OF UNIQUE "MINOR" GENE CONNECTED TO AT LEAST ONE EDGE

new_gene_name_list2 <- unique(nodesd)
print(no_of_genes2)
print(length(new_gene_name_list2$gene))
new_no_of_genes2 <- length(new_gene_name_list2$gene)

#COMBINE EDGE TABLES INTO A SINGLE TABLE WITH ALL THE EDGES

output_coeff_total <- rbind(output_coeff1s, output_coeff12s)

#REMOVE ALL REMOVED NOTED FROM THE ATTRIBUTE MATRIX

att_data1 <- data.frame(att_matrix1) %>%  rownames_to_column("gene")
att_data2 <- data.frame(att_matrix2) %>%  rownames_to_column("gene")


#TAKE THE ATTRIBUTE TABLES AND KEEP ONLY THE ENTRIES CORRESPONDING TO GENES THAT HAVE NOT BEEN REMOVED

att_data1 <- left_join(new_gene_name_list1, att_data1)
att_data2 <- left_join(new_gene_name_list2, att_data2)


#COMBINE ALL THE REMAINING NODES AND ATTRIBUTES INTO A SINGLE TABLE CONTAINING THEM ALL

output_att_total <- rbind(att_data1, att_data2) %>% column_to_rownames("gene")


#SAVE OUTPUT

if(attribute_file_yes == TRUE){
  output_att_name <- paste0(output_att_name, ".txt")
  write.table(output_att_total, output_att_name, row.names = TRUE, col.names = FALSE, sep = "\t",quote = FALSE)
}

output_coeff_name <- paste0(output_coeff_name,".txt")
write.table(output_coeff_total, output_coeff_name, row.names = FALSE, col.names = FALSE, sep = "\t",quote = FALSE)
