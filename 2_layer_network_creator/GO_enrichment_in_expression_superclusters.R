#THIS SCRIPT ANALYSES THE ENRICHMENT OR DEPLETION OF GENES WITH A PARTICULAR BINARY ATTRIBUTE IN DIFFERENT PARTS OF A NETWORK

#IN SUMMARY, THIS SCRIPT:
##CALCULATES THE ACTUAL NUMBER OF EDGES BETWEEN GENES IN EACH SUPERCLUSTER AND GENES WITH THE DESIRED ATTRIBUTE (i.e. defence/immunity GO slim)
##REASIGNS THE ATTRIBUTE TO PLANT GENES AT RANDOM AND RECALCULATES THE NUMBER OF EDGES WITH BETWEEN GENES IN EACH SUPERCLUSTER AND GENES WITH THE ATTRIBUTE
##REPEATS THIS 1,000 TIMES
##PLOTS A HISTOGRAM FOR EACH SUPERCLUSTER
##THE X-AXIS REPRESENTS THE NUMBER OF EDGES BETWEEN GENES IN A PARTICULAR SUPERCLUSTER AND GENES WITH THE ATTRIBUTE
##THE Y-AXIS REPRESENTS THE FREQUENCY OF OUTCOMES OUT OF 1000 REPEATS
##THE DOTTED LINES DENOTE 95 % CONFIDENCE INTERVALS
##THE RED LINES SHOW THE TRUE VALUES

library(tidyverse)
library(BiocGenerics)
library(svglite)

#Take all the edges in the network, rename the columns appropriately and drop weights which are not needed
#Takes the coefficient file outputted by "2_layer_make_array.R"
edges <- read.table("output_coefficient_file_0.975.txt", sep = "\t", header = FALSE, quote = "", dec=".")  %>% 
  rename("Plant_genes"=V1,"Nematode_genes"=V2) %>% select(-V3)

#Extract only the edges between a plant gene and an effector
edges<- edges[grep("AT", edges$Plant_genes),] 

#Read in the data containing all the nematode effectors, and keep the locus name and cluster
Input_nematode <- read_csv("putative_effectors_for_plant_immune_attribute.csv") %>% select(Gene_ID,Supercluster)
#Having NAs in the data breaks things, so these are replaced with "none"
Input_nematode$Supercluster[is.na(Input_nematode$Supercluster)] <- "none" 

#Take the  dataset containing the annotated plant genes in the network, and keeps the gene name and presence/absence of the GO slim term of interest
Input_plant <- read_csv("plant_genes_in_network_expression_and_attributes.csv") %>% select(Gene_ID,Immunity) %>% 
  rename(GO_slim = Immunity)

#Count all the unique entries in each category and sort into ascending order
category_list_plant <- unique(Input_plant$GO_slim) %>% sort() #Should be either 1 or 0

#Extract all the nematode clusters and sort into ascending alphabetical order
category_list_nematode <- unique(Input_nematode$Supercluster) %>% sort()
Number_of_categories_plant <- length(category_list_plant)
Number_of_categories_nematode <- length(category_list_nematode)

#Create a matrix that will hold all the degrees and fill it with 0s. 
#Rows represent the presence or absence of the GO slim term, first absence, then presence
#Columns represent the effector cluster in alphabetical order
rbs_stress_degrees_observed_actual <- matrix(nrow=Number_of_categories_plant, ncol=Number_of_categories_nematode, dimnames = list(category_list_plant,category_list_nematode),data=0)

#Add data to edges to indicate the cluster from which edge originates, and whether the corresponding plant gene has the GO slim term of interest
edges<- edges %>% left_join(Input_plant, by=join_by(Plant_genes==Gene_ID)) %>% left_join(Input_nematode, by=join_by(Nematode_genes==Gene_ID))

#Iterate through the whole list of edges
for(i in c(1:length(edges$Plant_genes))){
  #For each edge record the cluster
  cluster <- edges$Supercluster[i]
  #And the presence or absence of the GO attribute
  GO <- edges$GO_slim[i]
  #Finds the position of that cluster within the nematode category list
  cluster_no<- match(cluster,category_list_nematode)
  #Finds the position of that presence/absence within the plant category list
  GO_no <- match(GO, category_list_plant)
  #For each edge with a particular presence/cluster pair, the corresponding count in the degrees observed matrix is incremented by one
  rbs_stress_degrees_observed_actual[GO_no,cluster_no]<- rbs_stress_degrees_observed_actual[GO_no,cluster_no] +1
}

#This section only acts as a quick sanity check for the section after

#Take all the edges connected to a plant gene that do not show the GO slim term
zeroes <- edges %>% filter(GO_slim==0) 
#Calculate the proportion of edges connected to a plant genes that do not show the GO slim term - 
prop_0 <- length(zeroes$Plant_genes)/length(edges$Plant_genes) 

#Create the same matrix as for holding the observed values
rbs_stress_degrees_expected <- matrix(nrow=Number_of_categories_plant, ncol=Number_of_categories_nematode, dimnames = list(category_list_plant,category_list_nematode),data=0)

#Fill with the number of edges that would have been expected if they were distributed uniformly
#This is not an entirely accurate statistic, due to the discrete nature of edges, and the fact that not every plant node is connected to a uniform number of effectors and vice versa

#Iterate through all the clusters
for(i in c(1:Number_of_categories_nematode)){
  category <- category_list_nematode[i]
  #Count the number of edges originating from that cluster
  sub <- edges %>% filter(Supercluster==category)
  freq<- length(sub$Plant_genes)
  #Calculate the number of those edges that would be expected to connect to a plant gene that does not show the GO slim term, if these were distributed uniformly
  rbs_stress_degrees_expected[1,i] <- freq*(prop_0)
  rbs_stress_degrees_expected[2,i] <- freq*(1-prop_0)
}

#End of the first section of code, calculated the actual frequency of each pairing, as well as a hypothetical distribution that the null hypothesis might predict

view(rbs_stress_degrees_expected)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#This section of code bootstraps this data
#The presence or absence of the GO term on plant genes is re-assigned at random, and the frequency of each pairing recalculated a number of times equal to the repeats
#This helps to account for the bias that might otherwise be introduced by differences in connectivity
#Finally, a histogram is plotted for each cluster showing the frequency of connections to a positive plant gene in the bootstraps, as well as a line showing the true value, and lines containing 95% of the histogram data

repeats<-1000

#Create a 3 dimensional array to hold the resulting data
#Can be thought of as 1000 copies of the matrix from the first part stacked on top of one another
dimensions <- c(Number_of_categories_plant,Number_of_categories_nematode,repeats)
rbs_stress_degrees_generated <- rep(0, Number_of_categories_plant*Number_of_categories_nematode*repeats) #Create a vector containing as many 0s as there are holes in the 3 dimensional array you will generate
dim(rbs_stress_degrees_generated)<- dimensions #Convert the vector into the array

#Read the list of edges again 
edges <- read.table("output_coefficient_file_0.975_immune.txt", sep = "\t", header = FALSE, quote = "", dec=".")  %>% 
  rename("Plant_genes"=V1,"Nematode_genes"=V2) %>% select(-V3)

#Extract only the edges between a plant gene and an effector
edges<- edges[grep("AT", edges$Plant_genes),]

#Code is repeated for each repeat requested
for(i in c(1:repeats)){
  print(i) 
  #Extracts the number of genes in the network
  Plant_gene_number <- length(Input_plant$Gene_ID)
  #Extracts the number of plant genes in the network in which the GO slim term is positive
  GO_slim_number <- Input_plant %>% filter(GO_slim==1)
  GO_slim_number <- length(GO_slim_number$GO_slim)
  
  #Randomly chooses a subset of all the plant genes in the network equal in size to the number of genes in which the GO slim term is positive
  randomised <- sample(c(1:Plant_gene_number),GO_slim_number,replace = FALSE) #Without replacement, each gene can be picked only once
  bootstrap<-Input_plant
  
  #The existing GO_slim data is erased
  bootstrap$GO_slim <-0
  #And the new random data is added
  for(j in randomised){  #I'm pretty sure it would be more efficient to do this in a vectorised way instead of with a for loop, but I couldn't figure it out
    bootstrap[j,2]<-1
  }
  
  #The edges are associated with clusters and the data that has just been generated
  bootstrap_edges<- edges %>% left_join(bootstrap, by=join_by(Plant_genes==Gene_ID)) %>% left_join(Input_nematode, by=join_by(Nematode_genes==Gene_ID))
  
  #Iterate through the list of edges
  for(k in c(1:length(bootstrap_edges$GO_slim))){
    #Record the cluster of the effector of each edge
    cluster <- bootstrap_edges$Supercluster[k]
    #Record the GO slim status of the plant gene of each edge
    GO <- bootstrap_edges$GO_slim[k]
    #Record the position of that cluster within the nematode category list
    cluster_no<- match(cluster,category_list_nematode)
    #Record the position of that cluster within the plant category list
    GO_no <- match(GO, category_list_plant)
    #The corresponding count within one slice of the array is incremented by 1
    rbs_stress_degrees_generated[GO_no,cluster_no,i]<- rbs_stress_degrees_generated[GO_no,cluster_no,i] +1
  }
}

#Iterate through every cluster
for (i in c(1:27)){
  #Record the number of ones present for that cluster across all repeats
  ones<-tibble(rbs_stress_degrees_generated[2,i,])
  ones<- rename(ones, "bootstrapped_values"=colnames(ones[1])) #This is not an informative header
  
  #Plot:
  ggplot(data=ones, aes(x=bootstrapped_values)) + 
    #A histogram showing the distribution of ones across the repeats
    geom_histogram() +
    #A vertical red line showing what was actually observed
    geom_vline(xintercept=rbs_stress_degrees_observed_actual[2,i], color="red") +
    #Two vertical black dashed lines enclosing 95% of the data
    geom_vline(xintercept = quantile(ones$bootstrapped_values,probs=c(0.975)), linetype="dashed") + 
    geom_vline(xintercept = quantile(ones$bootstrapped_values,probs=c(0.025)), linetype="dashed") +
    #The title of the graph
    ggtitle(paste0("",colnames(rbs_stress_degrees_observed_actual)[i])) #This is not automatic, the title must be changed manually
  ggsave(paste0("",colnames(rbs_stress_degrees_observed_actual)[i],".svg")) #This folder must also be created and the name changed manually
}
