# Nitesh Turaga
# Make Peak Files Unique

# Set working Directory
###############################################################################
my.path ="~/TestRun/COPD-LungCancer-PeaksAnalysis/data/processed_unique"
sample.set = "Combined_2.1M_COPD"
setwd(my.path)
###############################################################################

# Set Path of Files to Make unique and Select Sample set
###############################################################################
make.unique.main = function(my.path,sample.set){
    # List Peaks Mappings file names
    mappings = list.files(file.path(my.path,sample.set),pattern = ".txt",full.names=TRUE,recursive=T)
    for (i in 1:length(mappings)) { 
        dat = read.table(file = mappings[i],sep = "\t",stringsAsFactors=FALSE,header = TRUE,fill=TRUE)
        dat = dat[order(dat$Peaks,decreasing=T),]
        dat.reduce = dat[!is.na(dat$ncbi_gene_id),]
        dat.reduce = dat.reduce[! is.na(dat.reduce$Name), ]

        new.name = paste0(strsplit(mappings[i],".txt")[1],".csv")
        write.csv(dat.reduce,file = new.name,row.names = FALSE)
    }
}

# Function call

make.unique.main(my.path,"Combined_2.1M_COPD")
rm(list = ls(all = TRUE)) 




# Unused code to make unique.
################################################################################
################################################################################

# # Make uniq based on all the columns
# make.uniq = function(filename) {
#     con = read.csv(filename)
#     uniq = con[!duplicated(con[,c(1:7)]),]
#     new.name = paste0(strsplit(filename,".csv")[1],"_uniq.csv")
#     write.csv(uniq,file = new.name,row.names = FALSE)
# }
# 
# # Make uniq based on all the name and region
# 
# make.uniq2 = function(filename) {
#     con = read.csv(filename)
#     uniq = con[!duplicated(con[,c(6:7)]),]
#     new.name = paste0(strsplit(filename,"_uniq.csv")[1],"_uniq2.csv")
#     write.csv(uniq,file = new.name,row.names = FALSE)
# }
# 
# # Make uniq based on all the Name
# 
# make.uniq3 = function(filename) {
#     con = read.csv(filename)
#     uniq = con[!duplicated(con[,c(6)]),]
#     new.name = paste0(strsplit(filename,"_uniq2.csv")[1],"_uniq3.csv")
#     write.csv(uniq,file = new.name,row.names = FALSE)
# }
# 
# 
# # Functionc call
# 
# 
# csvFiles = list.files("data/processed_unique/LungCancer_2.1M/",pattern = ".csv",full.names=T)
# 
# for (i in 1:length(csvFiles)) {
#     make.uniq(csvFiles[i])
#     make.uniq2(csvFiles[i])
#     make.uniq3(csvFiles[i])
# }
# 

