rm(list = ls())

find.bydate <- function(path, patt, ind, rank){
  ## function to get latest file
  ## path = the folder to look in, pattern = the pattern to match to
  ## ind = where the date is stored in the name if you were to strsplit by '_'
  ## rank = which file do you want? most recent = 1, second most recent = 2, etc.
  ## will print warning if rank is greater than the # of files, but will give you the oldest one
  
  files <- list.files(path, pattern = patt)
  
  if (rank > length(files)){
    print("Warning: rank greater than # of files, returning oldest file")
    rank <- length(files)
  }
  
  time_stamp <- rep (NA, length(files))
  
  for (i in 1:length(files)){
    time_stamp[i] <- unlist(strsplit(files[i], "_"))[ind]
  }
  
  time_stamp <- substr(as.character(time_stamp), 1, 8)
  diff <- as.numeric(Sys.Date() - as.Date(time_stamp, "%Y%m%d"))
  files <- files[order(diff)]
  retrieved <- read.csv(paste0(path, files[rank]), fileEncoding = 'latin1')
  print(paste0(path, files[rank]))
  return(retrieved)
}

scopus_recent <- find.bydate(path = "lit_review/", patt = "scopus", ind = 2, rank = 1)
scopus_compare <- find.bydate(path = "lit_review/", patt = "scopus", ind = 2, rank = 2)

pubmed_recent <- find.bydate(path = "lit_review/", patt = "pubmed", ind = 2, rank = 1)
pubmed_compare <- find.bydate(path = "lit_review/", patt = "pubmed", ind = 2, rank = 2)

## remove periods from pubmed Titles
pubmed_recent$Title <- gsub(".$", "", pubmed_recent$Title)
pubmed_compare$Title <- gsub(".$", "", pubmed_compare$Title)

## merge
master <- merge(scopus_compare, pubmed_compare, by = "Title",  all = TRUE)
toCheck <- merge(scopus_recent, pubmed_recent, by = "Title",  all = TRUE)
toCheck <- toCheck[!(toCheck$Title %in% master$Title), ]
write.csv(toCheck, "lit_review/toCheck.csv")
write.csv(master, "lit_review/master.csv")
