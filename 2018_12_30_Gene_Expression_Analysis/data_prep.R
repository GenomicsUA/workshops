dir="/Users/yuliya/Desktop/Bioinf_course/htscounts/"

myfiles <- list.files(path = dir,
                      pattern = ".counts",
                      all.files = TRUE,
                      recursive = FALSE,
                      ignore.case = FALSE,
                      include.dirs = FALSE)
DT <- list()
for (i in 1:length(myfiles) ) {
  infile = paste(dir, myfiles[i], sep = "/")
  DT[[myfiles[i]]] <- read.table(infile, header = F, stringsAsFactors = FALSE)
  cnts <- gsub("(.*)_all_counts.txt", "\\1", myfiles[i])
  colnames(DT[[myfiles[i]]]) <- c("ID", cnts)
}
data <- DT[[myfiles[1]]]
# we now add each other table with the ID column as key
for (i in 2:length(myfiles)) {
  y <- DT[[myfiles[i]]]
  z <- merge(data, y, by = c("ID"))
  data <- z
}

# ID column becomes rownames
rownames(data) <- data$ID
data <- data[,-1]

## add total counts per sample
data <- rbind(data, tot.counts=colSums(data))
data.all <- data[grep("^ENS", rownames(data), perl=TRUE, invert=FALSE), ]
write.csv(data.all, file = "Cancer.csv")