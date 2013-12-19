# Nitesh Turaga
# Dec 18th 2013
# Find Gain and loss of methylation in COPD samples

################################################################################
# Set working Directory

rm(list = ls(all=T))
setwd("~/TestRun/COPD-LungCancer-PeaksAnalysis/")
################################################################################


################################################################################
# Initializing step
# Preprocess for easy file reading. Make data frame
################################################################################

mypath = "data/processed_unique/Combined_2.1M_COPD"

files.normal = list.files(file.path(mypath,"Normal"),pattern = ".csv",full.names=T)
files.cancer = list.files(file.path(mypath,"Cancer"),pattern = ".csv",full.names=T)

files.to.analyze = data.frame(filenames = as.character(files.normal),status = "Normal")
files.to.analyze = rbind(data.frame(filenames = as.character(files.cancer),status = "Cancer"),files.to.analyze)

# STATUS or CASE CONTROL
T1 = "Cancer"
T2 = "Normal"

t1.files = files.to.analyze[which(files.to.analyze$status == T1),]
t2.files = files.to.analyze[which(files.to.analyze$status == T2),]
################################################################################


################################################################################
# Step 0: Mean Subtraction code: 
################################################################################
mean.substract = function(files.to.analyze) {
    # files to analyze is a data.frame
    files.mean = vector(length=nrow(files.to.analyze))
    for (i in 1:nrow(files.to.analyze)) {
        f = read.csv(as.character(files.to.analyze$filenames[i]))
        peaks.mean = mean(as.numeric(f$Peaks),na.rm=TRUE)
        #f$Peaks.mean.substracted = f$Peaks - peaks.mean
        files.mean[i] = peaks.mean
        #write.csv(f,file = as.character(files.to.analyze$filenames[i]))
        }
    files.to.analyze$mean = files.mean
    save(files.to.analyze,file = "objs/COPD_lung_files.to.analyze_dec19_00.rda")
}

# Function Call
mean.substract(files.to.analyze)
load("objs/COPD_lung_files.to.analyze_dec19_00.rda")

################################################################################






################################################################################
# READ FILES AND ORDER By Length of gene list
################################################################################
# Step 1 and 2:

file.rda.above.cuttoff = vector(length = nrow(files.to.analyze))
genelist.length = vector(length = nrow(files.to.analyze))
read.and.order = function(files.to.analyze) {    
    #Input is the data frame for files.to.analyze
    for (i in 1:nrow(files.to.analyze)) { 
        #Read Files
        f = read.csv(as.character(files.to.analyze$filenames[i]))
        #Get length of gene list
        f.genelist.length = length(!is.na(f$Name))
        genelist.length[i] = f.genelist.length        
        
#        f.order = order(f$Peaks)
#        file.ordered = f[f.order,]
        #Choose cut off
        cutoff = files.to.analyze$mean[i]
        file.above.cutoff = f[f$Peaks>cutoff,] 
        filename = paste0(strsplit(as.character(files.to.analyze$filenames[i]),".csv")[1],"_above_cutoff.rda")
        file.rda.above.cuttoff[i] = filename
        save(file.above.cutoff,file = filename )        
    }    
    files.to.analyze$file.rda.above.cutoff = file.rda.above.cuttoff
    files.to.analyze$genelist.length = genelist.length
    
    order.by.length = order(files.to.analyze$genelist.length,decreasing=TRUE)
    files.to.analyze = files.to.analyze[order.by.length,]    
    save(files.to.analyze,file = "objs/COPD_lung_files.to.analyze_dec19_01.rda")
}

#Function Call 
read.and.order(files.to.analyze)
load("objs/COPD_lung_files.to.analyze_dec19_01.rda")


################################################################################
# Step 3: Intersect genes
################################################################################
# Intersection step

intersect.gene.names = function(files.to.analyze,label.intersect = "") {
    library(BiocGenerics)
    if(nrow(files.to.analyze) <=1){
        stop("Length of given input cannot be intersected")
    }
    first = load(files.to.analyze$file.rda.above.cutoff[1])
    first.rda = file.above.cutoff$Name
    second = load(files.to.analyze$file.rda.above.cutoff[2])
    second.rda = file.above.cutoff$Name
    gene_intersect =  BiocGenerics::intersect(first.rda,second.rda)
    print(length(gene_intersect))
    filename = paste0("objs/gene_intersect_",label.intersect,"_2.rda")
    save(gene_intersect,file = filename)
    for (i in 3:nrow(files.to.analyze)) {
        next.file = load(files.to.analyze$file.rda.above.cutoff[i])
        next.rda = file.above.cutoff
        curr_intersect = BiocGenerics::intersect(gene_intersect,next.rda$Name) 
        gene_intersect = curr_intersect  
        print(length(gene_intersect))
        filename = paste0("objs/gene_intersect_",label.intersect,"_",i,".rda")
        save(gene_intersect,file = filename)
    }
    return(gene_intersect)
}

# debug(intersect.gene.names)
#Function Call
t1.files = files.to.analyze[files.to.analyze$status == "Cancer",]
genes_list_cancer = intersect.gene.names(t1.files,"Cancer")
t2.files = files.to.analyze[files.to.analyze$status == "Normal",]
genes_list_normal = intersect.gene.names(t2.files)

genes_list_normal_unique = genes_list_normal[!duplicated(genes_list_normal)]
genes_list_cancer_unique = genes_list_cancer[!duplicated(genes_list_cancer)]




################################################################################
# Step 4: Intersect genes with occurence percentage
################################################################################
# Intersection step

intersect.gene.names. = function(files.to.analyze) {
    library(BiocGenerics)
    if(nrow(files.to.analyze) <=1){
        stop("Length of given input cannot be intersected")
    }
    first = load(files.to.analyze$file.rda.above.cutoff[1])
    first.rda = file.above.cutoff$Name
    second = load(files.to.analyze$file.rda.above.cutoff[2])
    second.rda = file.above.cutoff$Name
    gene_intersect =  BiocGenerics::intersect(first.rda,second.rda)
    print(length(gene_intersect))
    for (i in 3:nrow(files.to.analyze)) {
        next.file = load(files.to.analyze$file.rda.above.cutoff[i])
        next.rda = file.above.cutoff
        curr_intersect = BiocGenerics::intersect(gene_intersect,next.rda$Name) 
        gene_intersect = curr_intersect  
        print(length(gene_intersect))
    }
    return(gene_intersect)
}

