library(reshape2)
library(energy)
library(tidyverse)


#As Dio's gexf creator but skipping the threshold and removal - these have been inorporated into the updated network creator

#!/usr/bin/env Rscript

#this requires all nodes to have edges in the input file (though not in the output file)

#INPUTS <----- ENTER USER DATA

#ENTER WHETHER YOU WOULD LIKE NODE ATTRIBUTES ON YOUR NETWORK
attribute_file_yes <- TRUE


#ENTER ATTRIBUTES INPUT
att <- read.table("Output_files/output_attribute_file_Neg_immune_def_new_0.975.txt", sep = "\t")

#ENTER COEFFICIENTS INPUT
coeff <- read.table("Output_files/output_coefficient_file_Neg_immune_def_new_0.975.txt", sep = "\t", header = FALSE, quote = "", dec=".")


#ENTER DESIRED OUTPUT NAME (exclude file ext)
final_gexf <- "Neg_immune_network"

#ENTER NODE ATTRIBUTE NAMES
node_att_names <- c("Cluster", "Kingdom", "Combined_kingdom_cluster","TF","3d_calc","Response_to_biotic_stimulus",	"Response_to_stress",	"RBS_and_stress",	"anatomical_structure_development",	"multicellular_organism_development",	"post.embryonic_development",	"embryo_development",	"flower_development",	"any_development",	"Annotations_of_interest","AUC","Proportion_negative","Neg_class","Defense_immunity_gene","Colour_by")

#Set Header
header_opener <- '<?xml version="1.0" encoding="UTF-8"?>\n<gexf xmlns:viz="http:///www.gexf.net/1.1draft/viz" version="1.1" xmlns="http://www.gexf.net/1.1draft">\n<meta lastmodifieddate="2010-03-03,23:44">\n<creator>Gephi 0.7</creator>\n</meta>\n<graph defaultedgetype="undirected" idtype="string" type="static">\n<attributes class = "node">'

#Set Node Attributes
if(attribute_file_yes == TRUE){
  no_of_node_att <- length(node_att_names)
  set_node_atts <- ""
  if(no_of_node_att >0){
    for(i in 1:no_of_node_att){
      set_node_atts <- paste0(set_node_atts,'\n<attribute id="',i-1,'" title="',node_att_names[i],'" type="string"/>')
    }
  }
}

#Set Edge Attributes
  
edge_att_names <- c("negpos")
no_of_edge_att <- length(edge_att_names)
set_edge_att <- '\n</attributes>\n<attributes class = "edge">'
  
for(i in 1:no_of_edge_att){
  if(i == no_of_edge_att){
    set_edge_att <- paste0(set_edge_att , '\n<attribute id="',i-1,'" title="',edge_att_names[i],'" type="string"/>\n</attributes>')
  }
  else{
    set_edge_att <- paste0(set_edge_att , '\n<attribute id="',i-1,'" title="',edge_att_names[i],'" type="string"/>')
  }
}

#Set Nodes
no_of_nodes <- length(att[,1])
nodes_header <- paste0('\n<nodes count="',no_of_nodes,'">')
nodes_footer <- "\n</nodes>"

node_combined <- vector(length = no_of_nodes)

for (i in 1:no_of_nodes) {
  node_header <- paste0('\n<node id="',att[i,1],'" label="',att[i,1],'">\n<attvalues>')
  node_atts <- ""
  if(attribute_file_yes == TRUE){
    for(j in 1:no_of_node_att){
      node_atts <- paste0(node_atts , '\n<attvalue for="',j-1,'" value="',att[i,j+1],'"/>')
    }
  }
  node_combined[i] <- paste0(node_header , node_atts , '\n</attvalues>\n</node>')
}

nodes_all_combined <- paste0(grep(pattern = "<node id",x = node_combined,value=TRUE),sep = "",collapse = "")
nodes_all_combined_write <- paste0(nodes_header,nodes_all_combined,nodes_footer)

print("nodes done")

no_of_edges <- length(coeff[,3])


#Set Edges (only works for 1 edge attribute - can be modified to suit more)

edges_header <- paste0('\n<edges count="',no_of_edges,'">')
edges_vector <- vector(length = no_of_edges)


for(i in 1:no_of_edges){
  #print(i)
  edges_vector[i] <- paste0('\n<edge id="',i-1,'" source="',coeff[i,1],'" target="',coeff[i,2],'" weight="',abs(coeff[i,3]),'">\n<attvalues>\n<attvalue for="0" value="',sign(coeff[i,3]),'"/>\n</attvalues>\n</edge>') 
}
print("all edges written")

edges_all_combined <- paste0(grep(pattern = "<edge id",x = edges_vector,value=TRUE),sep = "",collapse = "")
print("edges combined")

edges_footer <- '\n</edges>\n</graph>\n</gexf>'
edges_combined_write <- paste0(edges_header, edges_all_combined, edges_footer)

#Combines all into a single file
combine_all <- paste0(header_opener,set_node_atts,set_edge_att,nodes_all_combined,edges_combined_write)

write(combine_all,paste0(final_gexf,".gexf"))
print("all done")


