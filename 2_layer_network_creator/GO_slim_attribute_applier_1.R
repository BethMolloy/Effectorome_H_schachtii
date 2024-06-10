#INPUT .csv GO SLIM BULK ANALYSIS OUTPUT FILE (https://v2.arabidopsis.org/tools/bulk/go/index.jsp)
#OUTPUT ONLY THE SUBSET OF GENES WITH THE SPECIFIED TERM

library(tidyverse)

#INPUT GO SLIM TERM

slim_term <- "immune" #Input your GO slim term

#INPUT GO SLIM BULK ANALYSIS OUTPUT FILE

GOs <- read.csv("v2_GO_plant_genes.csv", header=TRUE) %>% 
  rename(GO_slim= GO.Slim.s.) %>% select(Locus,GO.term,GO.ID,category,GO_slim)  
Bi_proc <- GOs %>% filter(category == "proc") #Include any filters desired

subset<- Bi_proc[grep(slim_term, Bi_proc$GO_slim),] #Search for your GO slim term of interest
unique_subset <- unique(subset$Locus) #Include each locus in the list only once

#OUPUT
write.csv(unique_subset, paste0("Go_slim_",gsub(" ","_",slim_term),".csv"), row.names = FALSE)
