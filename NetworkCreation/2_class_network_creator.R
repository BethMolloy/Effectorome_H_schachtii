library(reshape2)
library(energy)
library(tidyverse)

#Takes two CSV files containing genes, expression data and attributes
#Every expression or attribute column in one input file must have a corresponding one at the same position in the other file 

#This code is updated so the data should be input in wide form instead of longform 
#The user can also specify if they would like the final network to colour different gene types by different attributes 


#Take a csv file containing all genes of the "major" type - Input 1 (These will have edges with each other)
Input_1_wide <- read.csv("../data_for_networks/Plant_data_with_attributes/Nematode_effectors.csv", stringsAsFactors = FALSE)
#Take a csv file containing all genes of the "minor" type - Input 2 (These will only have edges with the genes of major type)
Input_2_wide <- read.csv("../data_for_networks/Plant_data_with_attributes/Annotated_difference_of_means.csv", stringsAsFactors = FALSE)

#ENTER THE FINAL ATTRIBUTE FOR THE GENE TYPE TO BE COLOURED BY
colourby_1 <- "Cluster"
colourby_2 <- "Defense_immunity_gene"  

#ENTER A LIST OF GENE NAMES TO FILTER BY
#This should be a csv file consisting of a single column of gene names with no header. If this and the section on line 58 are uncommented, only genes in that list will be included in the network
#A filter can be applied to either the major genes or minor genes
#Filter_2.1 <- read.csv("../data_for_networks/CSV_files/TFs_expressed_in_gland_cells.csv", stringsAsFactors = FALSE, header = FALSE) %>% rename("Gene.ID" = "V1")
#Filter_2.2 <- read.csv("../data_for_networks/Anika_network/Anika_TFs_of_interest.csv", stringsAsFactors = FALSE, header = FALSE) %>% rename("Gene_ID" = "V1")

#ENTER DESIRED NAME FOR OUTPUT ATTRIBUTE FILE (do not include file ext)
output_att_name <- "output_attribute_file_Neg_immune_def_new_0.975" 

#ENTER DESIRED NAME FOR OUTPUT COEFFICIENT FILE (do not include file ext)
output_coeff_name <- "output_coefficient_file_Neg_immune_def_new_0.975" 

#ENTER NUMBER OF YOUR EXPRESSION CATEGORIES
no_of_exp_categories <- 5 

#ENTER WHETHER YOU WOULD LIKE AN ATTRIBUTE FILE - select FALSE if you do not want attributes at all.
attribute_file_yes <- TRUE

#ENTER EXPRESSION CATEGORY NAMES (must be same number as no_of_exp_categories)
exp_category_names <- c("Hours_10", "Hours_48", "Days_12_Female", "Days_12_Male", "Days_24_Female")

#ENTER EDGE THRESHOLD
threshold <- 0.975 

#Creates the colour_by attribute, which can be used to colour different gene types in the final network by different attributes
#Finds the column number of the chosen attribute
column_no_1 <- match(colourby_1,colnames(Input_1_wide))
  #Adds that column to the table under the heading colour_by
Input_1_wide <- Input_1_wide %>% mutate(Colour_by = Input_1_wide[,column_no_1])
  #Repeat for the minor gene type
column_no_2 <- match(colourby_2,colnames(Input_2_wide))
Input_2_wide <- Input_2_wide %>% mutate(Colour_by = Input_2_wide[,column_no_2])

#Convert wide form data to long form
Input_1 <- melt(Input_1_wide, na.rm = FALSE, id.vars="Gene_ID", variable.name = "Attribute", value.name = "Value") %>% 
  arrange(Gene_ID)
Input_2 <- melt(Input_2_wide, na.rm = FALSE, id.vars="Gene_ID", variable.name = "Attribute", value.name = "Value") %>% 
  arrange(Gene_ID)

#Applies the given filters
#Input_2 <- right_join(Input_2,Filter_2.1, by = Gene_ID) 
#Input_2 <- right_join(Input_2,Filter_2.2, by = Gene_ID)



#Preparation section taken from Dio's code
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

#Number of genes computed for both inputs

#First for input 1
no_of_att <- 1
transcript <- Input_1$Gene[1]
while (transcript == Input_1$Gene[no_of_att]){
  no_of_att <- no_of_att + 1
}
no_of_att <- no_of_att - 1

no_of_genes1 <- length(Input_1$Gene)/no_of_att

#Then for input 2

no_of_genes2 <- length(Input_2$Gene)/no_of_att

#2 seperate arrays of vectors generated

#Gene name list for input 1
gene_name_list1 <- matrix(nrow = no_of_genes1, ncol = 1)
for(l in 1:length(gene_name_list1[,1])){
  gene_name_list1[l, 1] <- Input_1$Gene[l*no_of_att]
}

#Gene name list for input 2
gene_name_list2 <- matrix(nrow = no_of_genes2, ncol = 1)
for(l in 1:length(gene_name_list2[,1])){
  gene_name_list2[l, 1] <- Input_2$Gene[l*no_of_att]
}

#Array of vectors representing input 1
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

#Array of vectors representing input 2
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

#Create 2 dataframes of node attributes

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

#Iterate through input 1, computing all internal edges
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

#Iterate through input 2, computing all edges between each element of input 2 and every element of input 1
  #In the resulting matrix, the rows are the genes of "minor" type and the columns are genes of "major" type
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

#Convert into longform and filter by the cutoff threshold
output_coeff1s <- data.frame(melt(coeff1s))
output_coeff1s <- subset(output_coeff1s,subset = is.na(output_coeff1s$value) == FALSE) %>% arrange(Var1,Var2) %>% 
  filter(-threshold > value | value > threshold)


output_coeff12s <- data.frame(melt(coeff12s))
output_coeff12s <- subset(output_coeff12s,subset = is.na(output_coeff12s$value) == FALSE) %>% arrange(Var1,Var2) %>% 
  filter(-threshold > value | value > threshold)

#Write this file here if needed for statistics
write.csv(output_coeff12s, "Output_files/RBS_plant-effector_coefficients.csv")

#Make 2 new node lists, keeping only genes with at least one edge

#All the columns containing "major" genes involved in edges
nodesa <- data.frame(output_coeff1s$Var1) %>%  rename("gene" = output_coeff1s.Var1)
nodesb <- data.frame(output_coeff1s$Var2) %>%  rename("gene" = output_coeff1s.Var2)
nodesc <- data.frame(output_coeff12s$Var2) %>%  rename("gene" = output_coeff12s.Var2)

#The column containing "minor" genes involved in edges
nodesd <- data.frame(output_coeff12s$Var1) %>%  rename("gene" = output_coeff12s.Var1)

#Count the number of unique "major" genes connected to at least one edge
new_gene_name_list1 <- unique(rbind(nodesa,nodesb,nodesc))
print(no_of_genes1)
print(length(new_gene_name_list1$gene))
new_no_of_genes1 <- length(new_gene_name_list1$gene)

#Count the number of unique "minor" genes connected to at least one edge
new_gene_name_list2 <- unique(nodesd)
print(no_of_genes2)
print(length(new_gene_name_list2$gene))
new_no_of_genes2 <- length(new_gene_name_list2$gene)

#Combine edge tables into a single table with all the edges
output_coeff_total <- rbind(output_coeff1s, output_coeff12s)

#Remove all removed nodes from the attribute matrix
att_data1 <- data.frame(att_matrix1) %>%  rownames_to_column("gene")
att_data2 <- data.frame(att_matrix2) %>%  rownames_to_column("gene")

#Takes the attribute tables, and keeps only the entries corresponding to genes that have not been removed
att_data1 <- left_join(new_gene_name_list1, att_data1)
att_data2 <- left_join(new_gene_name_list2, att_data2)

#Combine all remaining nodes and attributes into a single table containing them all
output_att_total <- rbind(att_data1, att_data2) %>% column_to_rownames("gene")


#Save outputs

if(attribute_file_yes == TRUE){
  output_att_name <- paste0("Output_files/",output_att_name, ".txt")
  write.table(output_att_total, output_att_name, row.names = TRUE, col.names = FALSE, sep = "\t",quote = FALSE)
}

output_coeff_name <- paste0("Output_files/",output_coeff_name,".txt")
write.table(output_coeff_total, output_coeff_name, row.names = FALSE, col.names = FALSE, sep = "\t",quote = FALSE)




