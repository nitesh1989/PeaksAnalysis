# Find Gain and loss of methylation in COPD samples
rm(list = ls(all=T))
# Load pd file
setwd("~/TestRun/COPD-LungCancer-PeaksAnalysis/")

for (i in 1:nrow(files.to.analyze)) { 
    labels[i] = (strsplit(as.character(files.to.analyze$filenames[i]),"_")[[1]])[1]    
}


mypath = "data/processed_unique/Combined_2.1M_COPD"

files.normal = list.files(file.path(mypath,"Normal"),pattern = "_uniq3.csv",full.names=T)
files.cancer = list.files(file.path(mypath,"Cancer"),pattern = "_uniq3.csv",full.names=T)

files.to.analyze = data.frame(filenames = as.character(files.normal),status = "Normal")
files.to.analyze = rbind(data.frame(filenames = as.character(files.cancer),status = "Cancer"),files.to.analyze)


# STATUS or CASE CONTROL
T1 = "Cancer"
T2 = "Normal"

t1.files = files.to.analyze[which(files.to.analyze$status == T1),]
t2.files = files.to.analyze[which(files.to.analyze$status == T2),]


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
    save(files.to.analyze,file = "objs/COPD_lung_files.to.analyze00.rda")
}

# Function Call
mean.substract(files.to.analyze)
load("objs/COPD_lung_files.to.analyze00.rda")

# Step 1 and 2:

file.rda.above.cuttoff = vector(length = nrow(files.to.analyze))

read.and.order = function(files.to.analyze) {    
    #Input is the data frame for files.to.analyze
    for (i in 1:nrow(files.to.analyze)) {
 
        f = read.csv(as.character(files.to.analyze$filenames[i]))
        f.order = order(f$Peaks)
        file.ordered = f[f.order,]
        #Choose cut off
        cutoff = files.to.analyze$mean[i]
        file.above.cutoff = file.ordered[file.ordered$Peaks>cutoff,] 
        filename = paste0(strsplit(as.character(files.to.analyze$filenames[i]),".csv")[1],"_above_cutoff.rda")
        file.rda.above.cuttoff[i] = filename
        save(file.above.cutoff,file = filename )
    }    
    files.to.analyze$file.rda.above.cutoff = file.rda.above.cuttoff
    save(files.to.analyze,file = "objs/COPD_lung_files.to.analyze01.rda")
}

#Function Call 
read.and.order(files.to.analyze)
load("objs/COPD_lung_files.to.analyze01.rda")


# Intersection step

intersect.gene.names = function(files.to.analyze) {
    library(BiocGenerics)
    if(nrow(files.to.analyze) <=1){
        stop("Length of given input cannot be intersected")
    }
    first = load(files.to.analyze$file.rda.above.cutoff[1])
    first.rda = file.above.cutoff$Name
    second = load(files.to.analyze$file.rda.above.cutoff[2])
    second.rda = file.above.cutoff$Name
    gene_intersect =  BiocGenerics::intersect(first.rda,second.rda)
    for (i in 3:nrow(files.to.analyze)) {
        next.file = load(files.to.analyze$file.rda.above.cutoff[i])
        next.rda = file.above.cutoff
        curr_intersect = BiocGenerics::intersect(gene_intersect,next.rda$Name) 
        gene_intersect = curr_intersect        
    }
    return(gene_intersect)
}

# debug(intersect.gene.names)
#Function Call
t1.files = files.to.analyze[files.to.analyze$status == "Cancer",]
genes_list_cancer = intersect.gene.names(t1.files)
t2.files = files.to.analyze[files.to.analyze$status == "Normal",]
genes_list_normal = intersect.gene.names(t2.files)

## CUTOFF = 4

#cutoff = as.numeric(4.0)
# pfna_ordered_above_cutoff = pfna_ordered[pfna_ordered$Peaks > cutoff,]
# pfoa_ordered_above_cutoff = pfoa_ordered[pfoa_ordered$Peaks > cutoff,]

genes_exposed_intersect = BiocGenerics::intersect(pfna_ordered_above_cutoff$Name,pfoa_ordered_above_cutoff$Name)

genes_exposed_intersect = as.data.frame(genes_exposed_intersect)
colnames(genes_exposed_intersect) = "Name"

genes_exposed_intersect_pfna = merge(genes_exposed_intersect,pfna_ordered_above_cutoff,by = "Name")
genes_exposed_intersect_pfna = genes_exposed_intersect_pfna[,c("CHROMOSOME","DATA_START","DATA_END","Peaks","Name","ncbi_gene_id","Region")]
write.csv(genes_exposed_intersect_pfna,file = "genes_PFNAExposed_above_cutoff_step2.csv",row.names = FALSE)

genes_exposed_intersect_pfoa = merge(genes_exposed_intersect,pfoa_ordered_above_cutoff,by = "Name")
genes_exposed_intersect_pfoa = genes_exposed_intersect_pfoa[,c("CHROMOSOME","DATA_START","DATA_END","Peaks","Name","ncbi_gene_id","Region")]
write.csv(genes_exposed_intersect_pfoa,file = "genes_PFOAExposed_above_cutoff_step2.csv",row.names = FALSE)



###############################################################################

# 3. List with intersection of genes with methylation peaks higher than 4 among 
# control pfna and control pfoa

ctr_pfna_ordered_above_cutoff = ctr_pfna_ordered[ctr_pfna_ordered$Peaks > cutoff,]
ctr_pfoa_ordered_above_cutoff = ctr_pfoa_ordered[ctr_pfoa_ordered$Peaks > cutoff,]

genes_control_intersect = BiocGenerics::intersect(ctr_pfna_ordered_above_cutoff$Name,ctr_pfoa_ordered_above_cutoff$Name)
write.csv(genes_control_intersect,file = "genes_allControl_above_cutoff_step3.csv")

genes_control_intersect = as.data.frame(genes_control_intersect)
colnames(genes_control_intersect) = "Name"

genes_control_intersect_pfna = merge(genes_control_intersect,ctr_pfna_ordered_above_cutoff,by = "Name")
genes_control_intersect_pfna = genes_control_intersect_pfna[,c("CHROMOSOME","DATA_START","DATA_END","Peaks","Name","ncbi_gene_id","Region")]
write.csv(genes_control_intersect_pfna,file = "genes_PFNAcontrol_above_cutoff_step3.csv",row.names = FALSE)

genes_control_intersect_pfoa = merge(genes_control_intersect,ctr_pfoa_ordered_above_cutoff,by = "Name")
genes_control_intersect_pfoa = genes_control_intersect_pfoa[,c("CHROMOSOME","DATA_START","DATA_END","Peaks","Name","ncbi_gene_id","Region")]
write.csv(genes_control_intersect_pfoa,file = "genes_PFOAcontrol_above_cutoff_step3.csv",row.names = FALSE)





###############################################################################
# 4. Eliminate genes common to genes_exposed_intersect and genes_control_intersect

common_genes = BiocGenerics::intersect(genes_control_intersect$Name,genes_exposed_intersect$Name)

control_index_independent = which(!is.element(genes_control_intersect$Name,common_genes))
exposed_index_independent = which(!is.element(genes_exposed_intersect$Name,common_genes))

genes_control_independant = genes_control_intersect[control_index_independent,]
genes_exposed_independant = genes_exposed_intersect[exposed_index_independent,]


write.csv(genes_control_independant,file = "genes_control_above_cutoff_step4.csv")
write.csv(genes_exposed_independant,file = "genes_exposed_above_cutoff_step4.csv")


genes_control_independant  = as.data.frame(genes_control_independant)
colnames(genes_control_independant) = "Name"

genes_control_independant_pfna = merge(genes_control_independant,ctr_pfna_ordered_above_cutoff,by = "Name")
genes_control_independant_pfna = genes_control_independant_pfna[,c("CHROMOSOME","DATA_START","DATA_END","Peaks","Name","ncbi_gene_id","Region")]
write.csv(genes_control_independant_pfna,file = "genes_PFNAcontrol_above_cutoff_step4.csv",row.names = FALSE)

genes_control_independant_pfoa = merge(genes_control_independant,ctr_pfoa_ordered_above_cutoff,by = "Name")
genes_control_independant_pfoa = genes_control_independant_pfoa[,c("CHROMOSOME","DATA_START","DATA_END","Peaks","Name","ncbi_gene_id","Region")]
write.csv(genes_control_independant_pfoa,file = "genes_PFOAcontrol_above_cutoff_step4.csv",row.names = FALSE)




genes_exposed_independant  = as.data.frame(genes_exposed_independant)
colnames(genes_exposed_independant) = "Name"

genes_exposed_independant_pfna = merge(genes_exposed_independant,pfna_ordered_above_cutoff,by = "Name")
genes_exposed_independant_pfna = genes_exposed_independant_pfna[,c("CHROMOSOME","DATA_START","DATA_END","Peaks","Name","ncbi_gene_id","Region")]
write.csv(genes_exposed_independant_pfna,file = "genes_PFNAexposed_above_cutoff_step4.csv",row.names = FALSE)

genes_exposed_independant_pfoa = merge(genes_exposed_independant,pfoa_ordered_above_cutoff,by = "Name")
genes_exposed_independant_pfoa = genes_exposed_independant_pfoa[,c("CHROMOSOME","DATA_START","DATA_END","Peaks","Name","ncbi_gene_id","Region")]
write.csv(genes_exposed_independant_pfoa,file = "genes_PFOAexposed_above_cutoff_step4.csv",row.names = FALSE)


###############################################################################
# 5. List with genes present in exposed but not in controls.
# 6. List with genes present in controls but not in exposed.
###############################################################################

genes_exposed = BiocGenerics::intersect(pfna_ordered$Name, pfoa_ordered$Name)
genes_control = BiocGenerics::intersect(ctr_pfna_ordered$Name, ctr_pfoa_ordered$Name)

commons = BiocGenerics::intersect(genes_control,genes_exposed)

control_index = which(!is.element(genes_control,commons))
exposed_index = which(!is.element(genes_exposed,commons))

genes_only_control = genes_control[control_index]
genes_only_exposed = genes_exposed[exposed_index]

write.csv(genes_only_control,file = "genes_lost_methylation_step6.csv")
write.csv(genes_only_exposed,file = "genes_gained_methylation_step6.csv")




genes_only_control  = as.data.frame(genes_only_control)
colnames(genes_only_control) = "Name"

genes_only_control_pfna = merge(genes_only_control,ctr_pfna_ordered,by = "Name")
genes_only_control_pfna = genes_only_control_pfna[,c("CHROMOSOME","DATA_START","DATA_END","Peaks","Name","ncbi_gene_id","Region")]
write.csv(genes_only_control_pfna,file = "genes_PFNAcontrol_lost_methylation_step6.csv",row.names = FALSE)

genes_only_control_pfoa = merge(genes_only_control,ctr_pfoa_ordered,by = "Name")
genes_only_control_pfoa = genes_only_control_pfoa[,c("CHROMOSOME","DATA_START","DATA_END","Peaks","Name","ncbi_gene_id","Region")]
write.csv(genes_only_control_pfoa,file = "genes_PFOAcontrol_lost_methylation_step6.csv",row.names = FALSE)




genes_only_exposed  = as.data.frame(genes_only_exposed)
colnames(genes_only_exposed) = "Name"

genes_only_exposed_pfna = merge(genes_only_exposed,pfna_ordered,by = "Name")
genes_only_exposed_pfna = genes_only_exposed_pfna[,c("CHROMOSOME","DATA_START","DATA_END","Peaks","Name","ncbi_gene_id","Region")]
write.csv(genes_only_exposed_pfna,file = "genes_PFNAexposed_gained_methylation_step6.csv",row.names = FALSE)

genes_only_exposed_pfoa = merge(genes_only_exposed,pfoa_ordered,by = "Name")
genes_only_exposed_pfoa = genes_only_exposed_pfoa[,c("CHROMOSOME","DATA_START","DATA_END","Peaks","Name","ncbi_gene_id","Region")]
write.csv(genes_only_exposed_pfoa,file = "genes_PFOAcontrol_gained_methylation_step6.csv",row.names = FALSE)
