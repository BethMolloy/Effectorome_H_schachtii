#IMPORTANT: Use only after "GO_slim_attribute_applier_1.R"
#Code for adding the discrete presense or absense of a GO_slim term as an attribute

slim_term <- "flower development"  #Input your GO_slim term
slim_term <- gsub(" ","_",slim_term)

annotation <- read.csv(paste0("../GO_data/GO_slim/",slim_term,".csv"), header=TRUE) %>% rename(Gene_ID=x) #Automatically read the file outputed for this GO_slim by the previous program

Plant_genes <- read.csv("../data_for_networks/Plant_data_with_attributes/Annotated_difference_of_means.csv") #Take the document to be updated
Plant_genes[, slim_term] <- 0 #Create a new column full of 0s named for the slim term

col_no<- which(colnames(Plant_genes) == slim_term)[[1]]

#Loop through the list of genes with the term, adding a 1 to the new column of the existing document if that gene is within it
for (i in annotation$Gene_ID){
  Plant_genes[,col_no]<- case_match(Plant_genes$Gene_ID, i~1, .default = Plant_genes[,col_no])
}

#Output
write.csv(Plant_genes,"../data_for_networks/Plant_data_with_attributes/Annotated_difference_of_means.csv",row.names = FALSE)
