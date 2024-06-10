#!/usr/bin/env Rscript

#THIS SCRIPT WAS WRITTEN BY DIO S. SHIN AND EDITED BY BETH MOLLOY

#THIS SCRIPT REQUIRES ALL INPUT NODES TO HAVE AN EDGE

#ENTER WORKING DIRECTORY
#setwd("PATH/WORKING_DIRECTORY")

#ENTER WHETHER OR NOT YOU WOULD LIKE NODE ATTRIBUTES IN THE NETWORK
attribute_file_yes <- TRUE

#ENTER COEFFICIENTS INPUT
coeff <- read.table("output_coefficient_file.txt", sep = "\t", header = FALSE, quote = "", dec=".")

#ENTER ATTRIBUTES INPUT
att <- read.table("output_attribute_file.txt", sep = "\t")

#ENTER DESIRED NETWORK THRESHOLD
threshold <- 0.950

#ENTER DESIRED OUTPUT NAME (exclude file extension)
final_gexf <- "717_putative_effectors_network_0.950"

#ENTER NODE ATTRIBUTE NAMES
node_att_names <- c("Known_or_putative_effector", "ppJ2", "pJ2", "pJ3", "Nucleotide sequence", "Amino acid sequence", "Supercluster", "Effector_family", "Expected gland updated", "HGT", "CWDE", "Annotation", "Orthogroup (OrthoMCL analysis)", "When evolved", "Secreted in Globo")

#SET HEADER

header_opener <- '<?xml version="1.0" encoding="UTF-8"?>\n<gexf xmlns:viz="http:///www.gexf.net/1.1draft/viz" version="1.1" xmlns="http://www.gexf.net/1.1draft">\n<meta lastmodifieddate="2010-03-03,23:44">\n<creator>Gephi 0.7</creator>\n</meta>\n<graph defaultedgetype="undirected" idtype="string" type="static">\n<attributes class = "node">'

#SET NODE ATTRIBUTES

if(attribute_file_yes == TRUE){
  no_of_node_att <- length(node_att_names)
  set_node_atts <- ""
  if(no_of_node_att >0){
    for(i in 1:no_of_node_att){
      set_node_atts <- paste0(set_node_atts,'\n<attribute id="',i-1,'" title="',node_att_names[i],'" type="string"/>')
    }
  }
}

#SET EDGE ATTRIBUTES

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

#SET NODES

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


#REMOVE EDGES WITH VALUES BELOW THE THRESHOLD - (Considers absolute values (magnitude))

removed_coeff <- coeff[abs(coeff[,3]) > threshold,]
no_of_edges <- length(removed_coeff[,3])

print("remove edges done")

#SET EDGES (currently only works for 1 edge attribute)

edges_header <- paste0('\n<edges count="',no_of_edges,'">')
edges_vector <- vector(length = no_of_edges)


for(i in 1:no_of_edges){
  print(i)
  edges_vector[i] <- paste0('\n<edge id="',i-1,'" source="',removed_coeff[i,1],'" target="',removed_coeff[i,2],'" weight="',abs(removed_coeff[i,3]),'">\n<attvalues>\n<attvalue for="0" value="',sign(removed_coeff[i,3]),'"/>\n</attvalues>\n</edge>') 
}
print("all edges written")

edges_all_combined <- paste0(grep(pattern = "<edge id",x = edges_vector,value=TRUE),sep = "",collapse = "")
print("edges combined")

edges_footer <- '\n</edges>\n</graph>\n</gexf>'
edges_combined_write <- paste0(edges_header, edges_all_combined, edges_footer)

#COMBINE ANALYSES INTO A SINGLE FILE

combine_all <- paste0(header_opener,set_node_atts,set_edge_att,nodes_all_combined,edges_combined_write)

write(combine_all,paste0(final_gexf,".gexf"))
print("all done")
