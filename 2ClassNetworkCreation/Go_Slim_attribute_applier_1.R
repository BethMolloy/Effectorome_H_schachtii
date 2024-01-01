#This code should take a csv file containing a long list of genes with corresponding GO terms and GO slim terms, and output a csv file containing only the subset of those genes annotated with a particular term of interest

slim_term <- "example GO term" #Input your GO slim term

#Take in a csv file containing loci, their corresponding GO terms, and the GO slims assosciated with each

GOs <- read.csv("../raw_data/All_plant_GO_terms_0.975_csv.csv", header=TRUE) %>% 
  rename(GO_slim= GO.Slim.s.) %>% select(Locus,GO.term,GO.ID,category,GO_slim)  
#Bi_proc <- GOs %>% filter(category == "proc") #Include any filters desired

subset<- Bi_proc[grep(slim_term, Bi_proc$GO_slim),] #Search for your slim term of interest
unique_subset <- unique(subset$Locus) #Include each locus in the list only once

#Output
write.csv(unique_subset, paste0("../GO_data/GO_slim/",gsub(" ","_",slim_term),".csv"), row.names = FALSE)
