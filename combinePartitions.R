#!/usr/bin/env Rscript

listfile <- list.files("/home/ss2945",pattern = "txt",full.names = T, recursive = FALSE)

#extract the  files with folder name Allgenes_text
listfile_a01 <- listfile[grep("ss2945",listfile)]

#inspect the file names
head(listfile_a01)
length(listfile_a01)

#combined all the text files in listfile_a01 and store in dataframe 'Data'
for (i in 1:length(listfile_a01)){
  print(i)
  if(i==1){
    assign(paste0("data"), read.table(listfile[i],header = FALSE))
  }
  if(!i==1){
    assign(paste0("file",i), read.table(listfile[i],header = FALSE))
    data <- rbind(data,get(paste0("file",i)))
    print("complete")
  }
}

write.table(data, file = "compiled.txt", row.names = FALSE, col.names = FALSE, sep = "\t",quote = FALSE)