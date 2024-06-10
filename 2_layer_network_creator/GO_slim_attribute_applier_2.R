#ADDS THE PRESENCE OR ABSENCE OF A GO SLIM TERM AS AN ATTRIBUTE TO A LIST OF GENES

slim_term <- "immune"  #Input your GO_slim term in plain text 
slim_term <- gsub(" ","_",slim_term)

#READ THE FILE OUTPUT FROM GO_slim_attribute_applier_2.R

annotation <- read.csv(paste0("Go_slim_",slim_term,".csv"), header=TRUE) %>% rename(Gene_ID=x)

Plant_genes <- read.csv("plant_genes_for_plant_superclustered_only_difference_of_means_GO_slim.csv") #Take the document to be updated
Plant_genes[, slim_term] <- 0 #Create a new column full of 0s named for the slim term

col_no<- which(colnames(Plant_genes) == slim_term)[[1]]

#LOOP THROUGH THE LIST OF GENES WITH THE SPECIFIED TERM, ADDING A 1 TO THE NEW COLUMN IF THE GENE IS PRESENT

for (i in annotation$Gene_ID){
  Plant_genes[,col_no]<- case_match(Plant_genes$Gene_ID, i~1, .default = Plant_genes[,col_no])
}

#OUTPUT
write.csv(Plant_genes,"plant_genes_for_plant_superclustered_only_difference_of_means_GO_slim.csv",row.names = FALSE)
