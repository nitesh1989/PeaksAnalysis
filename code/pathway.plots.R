# Nitesh Turaga
# nturaga1@jhmi.edu, Nov 21
# Make plots for DMR for pathway analysis
##################################################################################
# NOTE: This code does not plot for files which are annotated by partek 
# This is pathway analysis only for files which are in the bumphunting format.

break;

#######################
# Load the required files
####################### 
load("objs/TCGA_LUSC_objects/Beta_analysis.rda")
load("objs/TCGA_LUSC_objects/pd_LUSC.rda")
load("objs/TCGA_LUSC_objects/object_LUSC.rda")
load("objs/TCGA_LUSC_objects/M_analysis.rda")
load("objs/TCGA_LUSC_objects/currentSession_allobjects.rda")
source("~/TestRun/Gene_Lists/Gene_Lists_CSV/FUNCTIONS_2012.R")

#Model Matrix
pd=pData(object)
T1="normal";T2="cancer"
keep=pd$phenotype%in%c(T1,T2)
tt=factor(pd$phenotype[keep],c(T1,T2))
X=model.matrix(~tt)
design=model.matrix(~tt)
pos = start(object)


library(charm)
library(minfi)
library(limma)

##################################################################################

modifyOverlappingFeatures = function(tab) {
    contain.percentage.index =grep("with [0-9]+.[0-9]+%",tab$Overlapping.Features)
    
    #duplicate to not modify the original    
    x= tab
    x$Overlapping.Features = as.character(x$Overlapping.Features)
    for (i in contain.percentage.index) {
        x$Overlapping.Features[i] = gsub("[0-9]+.[0-9]+% of","",x$Overlapping.Features[i])
    }
    tab = x
    return(tab)
}

##################################################################################


#Path here is path to file.
pathwayPlot = function(path,T1,T2) {
    
    path = mypath        
    files.to.plot = list.files(path,pattern = ".csv")    
    setwd(path)
    for (i in 1:length(files.to.plot)) { 
        #         i=1
        tab = read.csv(files.to.plot[i])        
        if (nrow(tab) == 1){
            next;
        }
        y=M
        plotName = paste0(strsplit(files.to.plot[i],".csv"),"_dmrPlot.pdf")
        pdf(file=file.path(path,plotName))
        for(i in 1:nrow(tab)){
            #                     i=1    
            par(1,1)
            Index =tab$indexStart[i]:tab$indexEnd[i]
            #                     Index
            x=pos[Index]
            yy=y[Index,]
            
            #                     x
            #                     yy
            matplot(jitter(x),ilogit2(yy),col=as.fumeric(tt),
                    pch=16,cex=0.75, main=paste(tab$region[i], sep=" : ", tab$name[i]),xlab=paste("location on",tab$chr[i]),ylab="Beta")
            
            tmpIndexes=split(1:ncol(yy),tt)
            for(j in seq(along=tmpIndexes)){
                yyy=rowMeans(yy[,tmpIndexes[[j]]])
                lfit=loess(yyy~x,span=0.75)
                lines(x,ilogit2(lfit$fitted),col=j)
            }
            legend("bottomleft",levels(tt),col=1:2,lty=1)
        }
        dev.off()
        
    }
    
}




# Function call with Cancer
mypath = "~/TestRun/COPD-LungCancer-PeaksAnalysis/misc/PathwayAnalysis_geneintersect2"
setwd(mypath)
pathwayPlot(mypath,"normal","cancer")


