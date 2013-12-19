
library(BiocGenerics)
library(minfi)
library(limma)
source("code/FUNCTIONS_2012.R")

tcga_lusc = read.csv("data/TCGA_LUSC/bumps_analysis.csv")
genes_LU_2.1M = load("objs/gene_intersect_Cancer_2.rda") #The object inside is called "gene_intersect"


common_genes = BiocGenerics::intersect(tcga_lusc$name,gene_intersect)

indices = BiocGenerics::match(common_genes,tcga_lusc$name)

tcga_common = tcga_lusc[indices,]
write.csv(tcga_common,file = "misc/bumps_tcga_lusc_intersect_lungCancer2.1.csv")

#####################################
#Load bumphunter and TCGA information for plotting
load("objs/TCGA_LUSC_objects/Beta_analysis.rda")
load("objs/TCGA_LUSC_objects/pd_LUSC.rda")
load("objs/TCGA_LUSC_objects/object_LUSC.rda")
load("objs/TCGA_LUSC_objects/M_analysis.rda")
load("objs/TCGA_LUSC_objects/currentSession_allobjects.rda")

#####################################

#Model Matrix
pd=pData(object)
T1="normal";T2="cancer"
keep=pd$phenotype%in%c(T1,T2)
tt=factor(pd$phenotype[keep],c(T1,T2))
X=model.matrix(~tt)
design=model.matrix(~tt)
pos = start(object)


dmrPlotLocal = function(tab,plotName,save.path,M,tt,pos) {
    
    index <- tab$L>4
    tab <-tab[index,]
    dim(tab)
    y=M
    plotName = paste0(plotName,"_dmrPlot.pdf")
    pdf(file=file.path(save.path,plotName))
    for(i in 1:nrow(tab)){
#         i=1    
        par(1,1)
        Index =tab$indexStart[i]:tab$indexEnd[i]
#         Index
        x=pos[Index]
        yy=y[Index,]
#         x
#         yy
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

#Function call
tab = tcga_common
save.path = "figs"
plotName = "TCGA_LUSC_and_lungCancer2.1Local_common_genes_Lover4"

dmrPlotLocal(tab,plotName,save.path,M,tt,pos)


