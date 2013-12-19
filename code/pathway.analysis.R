# Nitesh Turaga
# Pathway analysis

# README: Just call functionCall after sourcing the file.

library(plyr)

##########################################
# Code chunk 2: Process the pathways files
##########################################

processPathways = function(pathwaysFile) { 
    ###################
    #test case
    #pathwaysFile = "Complete_Homeobox_(HOX)_Genes.csv"
    ###############
    a.pathway = read.csv(pathwaysFile)
    if(grep("PCR.Array.Catalog...",colnames(a.pathway))) {
        a.pathway$PCR.Array.Catalog... = NULL
    }
    gene.symb.index = which(a.pathway[1,] == "Symbol" | colnames(a.pathway) == "Symbol")
    pathway.genes = a.pathway[,gene.symb.index]
    pathway.genes = data.frame(as.character(pathway.genes))
    if (pathway.genes[1,1] == "Symbol") {
        pathway.genes = data.frame(Symbol = pathway.genes[-1,])
    }
    pathway.genes = data.frame(name = pathway.genes[which(pathway.genes$Symbol != ""),])
    #    pathway.file.name = paste0(strsplit(pathwaysFile,split=".csv")[1],"_Genes.csv")
    #    write.csv(pathway.genes,file = pathway.file.name)
    return(pathway.genes)
    
}


##########################################
# Code chunk 3: Compare with Data set with DMRS
##########################################
# Input : Dataset with a column descriptor which has the gene symbols
# It also requires the path to location of the pathways
# dataName = All the intersection files will have this name


intersectPathways = function(dataSet,dataName,pathToSave,fileType = ".csv") {
    # Default path to pathways
    pathToPathways = "~/TestRun/Gene_Lists/Gene_Lists_CSV"
    pathways = list.files(pathToPathways,pattern = ".csv")
    
    if(fileType == ".csv"){
        tab = read.csv(dataSet)
        data = tab
    }
            
    for (i in 1:length(pathways)) {
        genes = processPathways(pathways[i])        
        genes_present = intersect(genes$name,as.character(data$name))    
        toWrite = join(data.frame(name = genes_present),data,by = "name",type = "inner",match = "all")                
        fileName = paste0(pathToSave,strsplit(pathways[i],split=".csv")[1],"_",dataName,".csv")
        write.csv(toWrite,file = fileName,row.names = FALSE)
        
    }
}

# Function Calls
setwd("~/TestRun/Gene_Lists/Gene_Lists_CSV")
pathToSave = "~/TestRun/COPD-LungCancer-PeaksAnalysis/misc/"
intersectPathways("~/TestRun/COPD-LungCancer-PeaksAnalysis/misc/bumps_tcga_lusc_intersect_lungCancer.csv",dataName="tcga_LUSC_lungCancer2.1_local",pathToSave)


undebug(intersectPathways)

