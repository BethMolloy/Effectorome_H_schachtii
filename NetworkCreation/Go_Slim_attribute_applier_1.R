#Use before GO_slim_attribute_applier

slim_term <- "flower development" #Input your GO slim term


#Take in a csv file containing loci, their corresponding GO terms, and the GO slims assosciated with each
#This table is also included with these files, but contains ONLY THE GENES ACTUALLY IN THE NETWORK
#So will produce incorrect or invalid results if new plant genes are added to the dataset, or if the threshold is tweaked

GOs <- read.csv("../raw_data/All_plant_GO_terms_0.975_csv.csv", header=TRUE) %>% 
  rename(GO_slim= GO.Slim.s.) %>% select(Locus,GO.term,GO.ID,category,GO_slim)  
#Bi_proc <- GOs %>% filter(category == "proc") #Include only biological processes (maybe redundant)

subset<- Bi_proc[grep(slim_term, Bi_proc$GO_slim),] #Search for your slim term of interest
unique_subset <- unique(subset$Locus) #Include each locus in the list only once

#Output
write.csv(unique_subset, paste0("../GO_data/GO_slim/",gsub(" ","_",slim_term),".csv"), row.names = FALSE)
