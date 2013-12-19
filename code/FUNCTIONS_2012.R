##you can pass cutoff through the ...
regionFinder<-function(x,regionNames,chr,position,y=x,
                       summary=mean,ind=seq(along=x),order=TRUE,oneTable=TRUE,
                       ...){

  Indexes=getSegments(x[ind],regionNames[ind],...)
  
  res=vector("list",2)
  for(i in 1:2){
    res[[i]]=data.frame(chr=sapply(Indexes[[i]],function(Index) chr[ind[Index[1]]]),
         start=sapply(Indexes[[i]],function(Index) min(position[ind[Index]])),
         end=sapply(Indexes[[i]],function(Index) max(position[ind[Index]])),
         value=sapply(Indexes[[i]],function(Index) summary(y[ind[Index]])),
         area=sapply(Indexes[[i]],function(Index) abs(sum(y[ind[Index]]))),
         pns=sapply(Indexes[[i]],function(Index) regionNames[ind[Index]][1]),
         indexStart=sapply(Indexes[[i]],function(Index) min(ind[Index])),
         indexEnd=sapply(Indexes[[i]],function(Index) max(ind[Index])))
    res[[i]]$L=res[[i]]$indexEnd-res[[i]]$indexStart+1
  }
  names(res)=c("up","dn")
  if(order & !oneTable){
    if(nrow(res$up)>0) res$up=res$up[order(-res$up$area),]
    if(nrow(res$dn)>0) res$dn=res$dn[order(-res$dn$area),]
  }
  if(oneTable){
    res=rbind(res$up,res$dn)
    if(order & nrow(res)>0) res=res[order(-res$area),]
  }
  return(res)
}

logit=function(x) log(x)-log(1-x)
ilogit=function(x) 1/(1+exp(-x))
##nearestgene needs cisgenome to be installed
nearestgene <- function(object,up=5000,down=5000,species="human",
                        build="hg18",
                        datadir="/nexsan/bst2/cisreg/genomes/",
                        wd="./cisgenome",
                        binary="/home/bst/faculty/hji/projects/cisgenome_project/bin/refgene_getnearestgene",
                        inputFileName="tmp.cod",
                        outputFileName="tmp.txt",
                        r=0){##is the type of distance
  tmpscipen=.Options$scipen
  options(scipen=100)
  
  if(!file.exists(wd)) {
    cat("\nCreating directory:",wd,"\n")
    dir.create(wd)
  }
  
  inputfile=file.path(wd,inputFileName)
  outputfile=file.path(wd,outputFileName)
  write.table(data.frame(object$chr,object$start,object$end,rep("+",nrow(object))),file=inputfile,sep="\t",quote=FALSE,col.names=FALSE)

  thecall=paste(binary,"-d",file.path(datadir,species,build,"annotation/refFlat_sorted.txt"),"-dt 1 -s",species,"-i",inputfile,"-o",outputfile,"-r", r,"-up",up,"-down",down)
  system(thecall)
  res=read.delim(outputfile,as.is=TRUE,check.names=FALSE,comment.char="#",header=FALSE,col.names=c("seq_id","chr","start","end","strand","name","annotation","num","genestrand","TSS","TSE","CSS","CSE","ne","exonStarts","exonEnds"),fill=NA)
  res[res=="---"]<-NA
  res[res==""] <- NA
  options(scipen=tmpscipen)
  res
}

## This is the new regionMatch function as of 2/18/11:
regionMatch <- function(object1, object2, verbose=TRUE) {
    require(IRanges)
    m1 = object1; m2 = object2 #keep arguments named object1,object2 for backward compatibility
    stopifnot(all(c("chr","start","end")%in%colnames(m1)))
    stopifnot(all(c("chr","start","end")%in%colnames(m2)))
    if(any(is.na(m1[,c("chr","start","end")])) | any(is.na(m2[,c("chr","start","end")]))) stop("No missing values allowed in chr, start, or end columns.")
    m1$chr = as.character(m1$chr)
    m2$chr = as.character(m2$chr)
    stopifnot(is.numeric(m1$start) & is.numeric(m1$end))
    stopifnot(is.numeric(m2$start) & is.numeric(m2$end))
    stopifnot(all(m1$end>=m1$start))
    stopifnot(all(m2$end>=m2$start))
              
    ret = matrix(NA,nrow=nrow(m1),ncol=7)
    colnames(ret) = c("dist","matchIndex","type","amountOverlap","insideDist","size1","size2")
    sp1 = split(1:nrow(m1), m1$chr)
    sp2 = split(1:nrow(m2), m2$chr)          
    for(i in names(sp1)) {
        if(verbose) cat(i," ")
        inds1 = sp1[[i]]
        if(i %in% names(sp2)) {
            inds2 = sp2[[i]]
            m1IR = IRanges(start=m1$start[inds1],end=m1$end[inds1])
            m2IR = IRanges(start=m2$start[inds2],end=m2$end[inds2])
            ret[inds1,"matchIndex"] = inds2[ nearest(m1IR,m2IR) ]
        } else ret[inds1,"matchIndex"] = NA
    }
    if(verbose) cat("\n")
    m2b = m2[ret[,"matchIndex"],]
    ret = data.frame(ret, stringsAsFactors=FALSE)
    
    ret$type[(m1$start>m2b$start & m1$end<=m2b$end) |
             (m1$start>=m2b$start & m1$end<m2b$end)] = "inside"
    ret$type[m1$start<=m2b$start & m1$end>=m2b$end] = "cover"
    ret$type[m1$start>m2b$end] = "disjointR"
    ret$type[m1$end<m2b$start] = "disjointL"
    ret$type[is.na(ret$matchIndex)] = "disjoint"
    ret$type[(m1$start>m2b$start & m1$start<=m2b$end) & m1$end>m2b$end] = "overlapR"    
    ret$type[m1$start<m2b$start & (m1$end>=m2b$start & m1$end<m2b$end)] = "overlapL"

    ret$dist = 0
    ret$dist[ret$type=="disjoint"]  = NA
    ret$dist[ret$type=="disjointR"] = m2b$end[ret$type=="disjointR"] - m1$start[ret$type=="disjointR"]
    ret$dist[ret$type=="disjointL"] = m2b$start[ret$type=="disjointL"] - m1$end[ret$type=="disjointL"]
    ret$amountOverlap[ret$type=="overlapR"] = -1*(m2b$end[ret$type=="overlapR"]-m1$start[ret$type=="overlapR"]+1)
    ret$amountOverlap[ret$type=="overlapL"] = m1$end[ret$type=="overlapL"]-m2b$start[ret$type=="overlapL"]+1
    ret$type[ret$type%in%c("disjointR","disjointL")] = "disjoint"
    ret$type[ret$type%in%c("overlapR","overlapL")] = "overlap"

    ## insideDist column:
    insideIndex = ret$type=="inside" #no missing ret$type at this point
    tmp0 = cbind(m1$end[insideIndex]  -m2b$end[insideIndex],
                 m1$start[insideIndex]-m2b$start[insideIndex])
    tmp = apply(abs(tmp0),1,which.min)
    tmpinsidedist = tmp0[,1]
    tmpIndex = tmp==2
    tmpinsidedist[tmpIndex] = tmp0[tmpIndex,2]
    ret$insideDist[insideIndex] = tmpinsidedist

    ## size1 and size2 columns:
    ret$size1 = m1$end -m1$start +1
    ret$size2 = m2b$end-m2b$start+1
    
    ret
}

###COMMENTED OUT RAFA's Original. Put in Peter's
## regionMatch <- function(object1,object2,verbose=TRUE){
##   Indexes1=split(seq(along=as.character(object1$chr)),as.character(object1$chr))
##   Indexes2=split(seq(along=as.character(object2$chr)),as.character(object2$chr))
##   res=matrix(NA,nrow=nrow(object1),7)
##   for(i in seq(along=Indexes1)){
##     CHR=names(Indexes1)[i]
##     if(verbose) cat(CHR,",")
##     Index1=Indexes1[[CHR]]
##     Index2=Indexes2[[CHR]]
##     if(length(Index1)>0 & length(Index2)>0){

##       tmpmat1=object1[Index1,c("start","end")]
##       tmpmat2=object2[Index2,c("start","end")]

##       tmp11=fuzzy.match2(tmpmat1[,1],tmpmat2)
##       tmp12=fuzzy.match2(tmpmat1[,2],tmpmat2)

##       tmp21=fuzzy.match2(tmpmat2[,1],tmpmat1)
##       tmp22=fuzzy.match2(tmpmat2[,2],tmpmat1)

##       tmpIndex = which(tmp21[,1]==0 & tmp22[,1]==0)
##       if(length(tmpIndex)>0){
##         tmpIndex2 = tmpIndex[ tmp21[tmpIndex,2]==tmp22[tmpIndex,2] ]
##         if(length(tmpIndex2)>0) insideIndex2= tmp21[tmpIndex2,2] else insideIndex2=c()
##       } else insideIndex2=c()
      
##       tmpIndex = which(tmp11[,1]==0 & tmp12[,1]==0)
##       if(length(tmpIndex)>0){
##         insideIndex=tmpIndex[
##         (object1$start[Index1[tmpIndex]] >=
##          object2$start[Index2[tmp11[tmpIndex,2]]] &
##          object1$end[Index1[tmpIndex]] <=
##          object2$end[Index2[tmp11[tmpIndex,2]]] ) |
##         (object1$start[Index1[tmpIndex]] >=
##          object2$start[Index2[tmp12[tmpIndex,2]]] &
##          object1$end[Index1[tmpIndex]] <=
##          object2$end[Index2[tmp12[tmpIndex,2]]] )
##           ]          
##       } else insideIndex=c()
      
##       overlapIndex =  setdiff(which(tmp11[,1]==0 | tmp12[,1]==0),insideIndex)
        
##       coverIndex = which((sign(tmp11[,1])==1 & sign(tmp12[,1])==-1 &
##         tmpmat2[tmp11[,2],1] <= tmpmat2[tmp12[,2],1]))
      
##       coverIndex = union(coverIndex,insideIndex2)

## ###Minimum distance and Index of minimum
##       tmp3=cbind(tmp11[,1],tmp12[,1])
##       tmp=apply(abs(tmp3),1,which.min)
##       tmpdist=tmp11[,1]
##       tmpIndex=tmp==2
##       tmpdist[tmpIndex]=tmp12[tmpIndex,1]
##       tmpdist[coverIndex]=0

##       tmpMatch1=Index2[tmp11[,2]]
##       tmpMatch2=Index2[tmp12[,2]]
##       tmpMatch=tmpMatch1
##       tmpMatch[tmpIndex]=tmpMatch2[tmpIndex]
      
##       ###overlap
##       overlap1 = (tmp11[overlapIndex,1]==0)*
##         (object1$start[Index1[overlapIndex]] -
##          object2$end[tmpMatch1[overlapIndex]])
                    

##       overlap2 = (tmp12[overlapIndex,1]==0)*
##         (object1$end[Index1[overlapIndex]]-
##          object2$start[tmpMatch2[overlapIndex]])

##       overlap = pmax(abs(overlap1),abs(overlap2))
##       ##where inside..
##       tmp3=cbind(object1$end[Index1[insideIndex]] - object2$end[tmpMatch[insideIndex]],object1$start[Index1[insideIndex]] - object2$start[tmpMatch[insideIndex]])
##       tmp=apply(abs(tmp3),1,which.min)
##       tmpinsidedist=tmp3[,1]
##       tmpIndex=tmp==2
##       tmpinsidedist[tmpIndex]=tmp3[tmpIndex,2]
    
##       res[Index1,1]<-tmpdist
##       res[Index1,2]<-tmpMatch
##       res[Index1,3]<-0;if(length(insideIndex)>0)res[Index1[insideIndex],3] <- 1
##       if(any(insideIndex)){
##         res[Index1[insideIndex],4] <- tmpinsidedist
##       }
##       res[Index1,5]<-0;if(length(overlapIndex)>0)res[Index1[overlapIndex],5] <- 1
##       if(any(overlapIndex)){
##         res[Index1[overlapIndex],6]<-overlap
##       }
##       res[Index1,7]<-0;if(length(coverIndex)>0)res[Index1[coverIndex],7] <- 1
##     }
##   }
##   type=rep("disjoint",nrow(res))
##   type[res[,3]==1] <- "inside"
##   type[res[,5]==1] <- "overlap"
##   type[res[,7]==1] <- "cover"
##   type=factor(type,levels=c("disjoint","overlap","inside","cover"))
##   size1=object1$end-object1$start+1
##   size2=object2$end[res[,2]] - object2$start[res[,2]] + 1
##   res=data.frame(dist=res[,1],matchIndex=res[,2],type=type,amountOverlap=res[,6],insideDist=res[,4],size1=size1,size2=size2)
##   res
## }

pointMatch <- function(object1,object2,col1=2,col2=2,verbose=TRUE){
  Indexes1=split(seq(along=as.character(object1$chr)),as.character(object1$chr))
  Indexes2=split(seq(along=as.character(object2$chr)),as.character(object2$chr))
  res=matrix(NA,nrow=nrow(object1),2)
  for(i in seq(along=Indexes1)){
    CHR=names(Indexes1)[i]
    if(verbose) cat(CHR,",")
    Index1=Indexes1[[CHR]]
    Index2=Indexes2[[CHR]]
    if(length(Index1)>0 & length(Index2)>0){
      tmpmat1=object1[Index1,col1]
      tmpmat2=object2[Index2,col2]
      res[Index1,]=fuzzy.match(tmpmat1,tmpmat2)
    }
  }
  colnames(res)<-c("dist","matchIndex")
  return(res)
}


getSegments <- function(x,factor,cutoff=quantile(abs(x),0.99),verbose=TRUE){

  Indexes=split(seq(along=x),factor)
  regionID=vector("numeric",length(x))
  LAST = 0
  
  segmentation = vector("numeric", length(x))
  type = vector("numeric", length(x))

  for (i in seq(along = Indexes)) {
    if (verbose) if (i%%1000 == 0) cat(".")
    Index = Indexes[[i]]
    y = x[Index]
    z = sign(y) * as.numeric(abs(y) > cutoff)
    w = cumsum(c(1, diff(z) != 0)) + LAST
    segmentation[Index] = w
    type[Index] = z
    LAST = max(w)
  }
  ##add a vector of the pns
  res=list(upIndex=split(which(type>0),segmentation[type>0]),
    dnIndex=split(which(type<0),segmentation[type<0]),
    zeroIndex=split(which(type==0),segmentation[type==0]))
  names(res[[1]])<-NULL
  names(res[[2]])<-NULL
  names(res[[3]])<-NULL
  return(res)
}

matchGenes <- function(object,promoterDist=2500,verbose=TRUE,species="human",
                       build="hg18",
                       r=0,up=50000000,down=50000000){
  if(verbose) cat("Matching regions to genes.\n")
  genes=nearestgene(object,up=up,down=down,species=species,build=build,r=r)
  type=rep("",nrow(genes))
  subtype=rep("",nrow(genes))
  ctype=rep("",nrow(genes))
  dist=rep(0,nrow(genes))
  insidedistance<-rep(NA,nrow(genes))
  exonnumber<-rep(NA,nrow(genes))
  nexons<-rep(NA,nrow(genes))
  
  geneL=rep(0,nrow(genes))
  codingL=rep(0,nrow(genes))
  subdist=rep(0,nrow(genes))
  if(verbose) cat("Annotating results")
  for(i in 1:nrow(genes)){
    if(verbose & i%%1000==0)  cat(".")

    if(!is.na(genes[i,10])){
    TS = genes[i,10]
    TE = genes[i,11]
    geneL[i]=TE-TS
    CS = genes[i,12]
    CE = genes[i,13]
    codingL[i]=CE-CS
    ES = as.numeric(strsplit(genes[i,15],",")[[1]])
    EE = as.numeric(strsplit(genes[i,16],",")[[1]])
    Exons= cbind(ES,EE)
    nexons[i]=nrow(Exons)
    Strand= ifelse(genes[i,9]=="+",1,-1)
    S = genes[i,3]
    E = genes[i,4]
    
    type[i]=""
    if(genes[i,3] <= TS & genes[i,4] >= TE){
      type[i]="covers"
    } else{
      if(genes[i,3] < TS & genes[i,4] < TS){
        if(Strand==1){
          type[i]="upstream" 
          dist[i]=TS-E
        } else{
          type[i]="downstream"
          dist[i]=TE-E
        }
      }
      if(genes[i,3] > TE & genes[i,4] > TE){
        if(Strand==-1){
          type[i]="upstream"
          dist[i]=S-TE
        }  else{
          type[i]="downstream"
          dist[i]=S-TS
        }
      }
      if (genes[i,3]>= TS & genes[i,4] <= TE){
        type[i]="inside"
        if(Strand==-1) dist[i]=TE-E  else dist[i]=S-TS
      }
      
      ##overlaps
      if(type[i]==""){
        if(genes[i,3] < TS & genes[i,4] <=TE){
          ##OVERLAP FRONT
          if(Strand==1) type[i]="overlaps 5'" else{
            type[i]="overlaps 3'"
            dist[i]=TE-E
          }
          S=TS
        }
        if (genes[i,3] >= TS & genes[i,4] > TE){
          ##OVERAL BACK
          if(Strand==-1) type[i]="overlaps 5'" else{
            type[i]="overlaps 3'"
            dist[i]=S-TS
          }
          E=TE
        }
      }
    }
    m1=NA;m2=NA
    if(S >= TS & E <= TE){
      
      ##INSIDE
      tmp1=fuzzy.match2(S,Exons)
      tmp2=fuzzy.match2(E,Exons)
      
      m1=tmp1[1,1]
      m2=tmp2[1,1]
      exon1=tmp1[1,2]
      exon2=tmp2[1,2]
      m1m2Index=which.min(abs(c(m1,m2)))
      
      if(exon1==exon2 & m1==0 & m2==0){
        subtype[i]="inside exon"
        insidedistance[i]=0
        exonnumber[i]=c(exon1,exon2)[m1m2Index]
      }
      
      if( (sign(m1)==sign(m2) & (m1!=0 & m2!=0) & exon1==exon2) |
         (sign(m1)==-1 & sign(m2)==1 & exon2-exon1==1) ){
        subtype[i]="inside intron"

        insidedistance[i]=c(m1,m2)[m1m2Index]
        exonnumber[i]=c(exon1,exon2)[m1m2Index]
        
      }
      if( (exon2-exon1 > 1) |
         ( (exon2-exon1 == 1) & 
          ((sign(m1)==-1 & sign(m2)==-1) |
           (sign(m1)==1 & sign(m2)==-1) |
           (sign(m1)==1 & sign(m2)==1) |
           (sign(m1)==0 & sign(m2)==-1)|
           (sign(m1)==1 & sign(m2)==0))) |
         (exon2==exon1 & sign(m1)==1 & sign(m2)==-1)){
        subtype[i]="covers exon(s)"
        insidedistance[i]=0
        exonnumber[i]=c(exon1,exon2)[m1m2Index]
      }
      
      if( (exon2-exon1 == 1 & sign(m1)==-1 & sign(m2)==0) |
         (exon2==exon1 & sign(m1)==1 & sign(m2)==0)){
        if(Strand==1) subtype[i]="overlaps exon upstream" else subtype[i]="overlaps exon downstream"
        insidedistance[i]=0
        exonnumber[i]=c(exon1,exon2)[m1m2Index]
      }
      if( (exon2-exon1 == 1 & sign(m1)==0 & sign(m2)==1) |
         (exon2==exon1 & sign(m1)==0 & sign(m2)==-1)){
        if(Strand==-1) subtype[i]="overlaps exon upstream" else subtype[i]="overlaps exon downstream"
        insidedistance[i]=0
        exonnumber[i]=c(exon1,exon2)[m1m2Index]
      }
      if( exon2-exon1 == 1 & sign(m1)==0 & sign(m2)==0){
        subtype[i]="overlaps two exons"
        insidedistance[i]=0
        exonnumber[i]=c(exon1,exon2)[m1m2Index]
      }

      if(Strand!=1){
        insidedistance[i]= -insidedistance[i]
        exonnumber[i] = nrow(Exons) - exonnumber[i] + 1
      }
      
      ctype[i]="inside transcription region"
      
      if(S<CS & E<CS){
        if(Strand==1) ctype[i]="5' UTR" else ctype[i]="3'UTR"
      }
      
      if(S>CE & E>CE){
        if(Strand==-1) ctype[i]="5' UTR" else ctype[i]="3'UTR"
      }
      if(S<CS & E>CE){
        ctype[i]="covers transcription region"
      }
      if(S<CS & E>CS){
        if(Strand==1) ctype[i]="overlaps 5' UTR" else ctype[i]="overlaps 3'UTR"
      }
      if(S<CE & E>CE){
        if(Strand==-1) ctype[i]="overlaps 5' UTR" else ctype[i]="overlaps 3'UTR"
      }
    }
    if(FALSE){##graphical check
      plot(0,0,ylim=c(0,1),xlim=range(genes[i,c(3:4,10:11)]),xlab=paste("inside distance=",insidedistance[i],m1,m2,exonnumber[i]))
      polygon(c(TS,TE,TE,TS),c(0,0,0.5,0.5),density=0,col=2)
      polygon(c(CS,CE,CE,CS),c(0.1,0.1,0.4,0.4),density=0,col=3)
      abline(h=0.25,lwd=2)
      apply(Exons,1,function(x) polygon(c(x[1],x[2],x[2],x[1]),
                                        c(0.2,0.2,0.3,0.3),col=4))
      polygon(as.vector(genes[i,c(3,4,4,3)]),c(0.75,0.75,0.85,0.85),col=5)
      lines(c(TS,TS+1000),c(0.65,0.65),lwd=3)
      title(paste(i,Strand,type[i],subtype[i],ctype[i],dist[i],sep=":"))

      
    }
  }
  }
  if(verbose) cat("Done.\n")
  type[dist<=promoterDist & type=="upstream"] <- "promoter"
  type[dist<=promoterDist & type=="downstream"] <- "close to 3'"

  description=type
  tmpIndex=which(description=="inside")
  description[tmpIndex] <- subtype[tmpIndex]
  return(data.frame(name=I(genes[,6]),
                    annotation=I(genes[,7]),
                    description=factor(description,levels=c("upstream","promoter","overlaps 5'","inside intron","inside exon","covers exon(s)","overlaps exon upstream","overlaps exon downstream","overlaps two exons","overlaps 3'","close to 3'","downstream","covers")),
                    region=factor(type,levels=c("upstream","promoter","overlaps 5'","inside","overlaps 3'","close to 3'","downstream","covers")),
                    distance=dist,
                    subregion=factor(subtype,levels=c("inside intron","inside exon","covers exon(s)","overlaps exon upstream","overlaps exon downstream","overlaps two exons")),
                    insidedistance=insidedistance,
                    exonnumber=exonnumber,
                    nexons=nexons,
                    UTR=factor(ctype,levels=c("inside transcription region","5' UTR","overlaps 5' UTR","3'UTR","overlaps 3'UTR","covers transcription region")),
                    strand=genes[,9],
                    geneL=geneL,
                    codingL=codingL))
}

dogo <- function(names,universe,species="human"){
  if(species=="human"){
    golib="org.Hs.eg.db"
    library(golib,character.only=TRUE)
    gomap= org.Hs.egREFSEQ2EG
  } else  {
    golib="org.Mm.eg.db"
    library(golib,character.only=TRUE)
    gomap= org.Mm.egREFSEQ2EG
  }
  require(GOstats)
  x=unlist(mget(as.character(names), gomap,ifnotfound = NA))
  x=x[!is.na(x)]
  Universe=unlist(mget(as.character(universe),gomap,ifnotfound = NA))
  
  Universe=unique(c(Universe[!is.na(Universe)],unique(x)))
  
  params <- new("GOHyperGParams", geneIds = unique(x),
                universeGeneIds = Universe,
                annotation = golib,
                ontology = "BP", pvalueCutoff = 0.01, conditional = TRUE,
                testDirection="over")
  ht=hyperGTest(params)
  tab=summary(ht)
  tmp1=geneIdsByCategory(ht)
  tmp1=tmp1[tab[,1]]
  tab$IDs=sapply(tmp1,function(y) paste(names(x)[x%in%y],collapse=";"))
  return(tab)
}




affyGO <- function(names,universe=NULL,platform="hgu133plus2",species="human",pvalueCutoff=1){
  library(paste(platform,"db",sep="."),character.only=TRUE)
  gomap=get(paste(platform,"ENTREZID",sep=""))
  symmap=get(paste(platform,"SYMBOL",sep=""))
  if(species=="human"){
    golib="org.Hs.eg.db"
    library(golib,character.only=TRUE)
  } else  {
    golib="org.Mm.eg.db"
    library(golib,character.only=TRUE)
    gomap= org.Mm.egREFSEQ2EG
  }
  require(GOstats)
  x=unlist(mget(as.character(names), gomap,ifnotfound = NA))
  x=x[!is.na(x)]

  if(is.null(universe)) universe=ls(gomap)
  Universe=unlist(mget(universe,gomap,ifnotfound = NA))
  Universe=unique(c(Universe[!is.na(Universe)],unique(x)))
  
  params <- new("GOHyperGParams", geneIds = unique(x),
                universeGeneIds = Universe,
                annotation = golib,
                ontology = "BP", pvalueCutoff = pvalueCutoff,
                conditional = TRUE,
                testDirection="over")
  ht=hyperGTest(params)
  tab=summary(ht)
  tmp1=geneIdsByCategory(ht)
  tmp1=tmp1[tab[,1]]
  tab$IDs= sapply(tmp1,function(y) paste(unlist(mget(names(x)[x%in%y],symmap,ifnotfound=NA)),collapse=";"))
  return(tab)
}

myfilter2 <- function(x,filter,...){
###unlike myfilter, myfilter2 returns NAs in the edges.
  L=dim(x)[1]
  if(L>length(filter)) res=filter(x,filter,...) else{res = t(colMeans(x) %*% t(rep(1,dim(x)[1])))}

  return(res)

}


getDesc<-function(x){
 ##get gene description
  require(org.Hs.eg.db)
  genenames=sapply(mget(as.character(x),org.Hs.egREFSEQ2EG,ifnotfound = NA),function(x) x[1])
  genenames[!is.na(genenames)]= sapply(mget(genenames[!is.na(genenames)],org.Hs.egGENENAME,ifnotfound=NA),function(x) x[1])
  genenames
}


clusterMaker <- function(chr,pos,order.it=TRUE,maxGap=300){
  nonaIndex=which(!is.na(chr) & !is.na(pos))
  Indexes=split(nonaIndex,chr[nonaIndex])
  clusterIDs=rep(NA,length(chr))
  LAST=0
  for(i in seq(along=Indexes)){
    Index=Indexes[[i]]
    x=pos[Index]

    if(order.it){ Index=Index[order(x)];x=pos[Index] }

    y=as.numeric(diff(x)>maxGap)
    z=cumsum(c(1,y))
    clusterIDs[Index]=z+LAST
    LAST=max(z)+LAST
  }
  clusterIDs
}  

midpoints <- function(x) (x[1:(length(x)-1)]+x[2:length(x)])/2

pointsplot <- function(l,jitter=TRUE,factor=1,col=rep(1,length(l)),pch=21,...){
  if(!is.list(l)) stop("l must be a list\n")
  y=unlist(l)
  x=rep(seq(along=l),sapply(l,length))
  col=rep(col,sapply(l,length))
  if(jitter) x=jitter(x,factor)
  if(pch!=21) plot(x,y,col=col,pch=pch,...) else plot(x,y,pch=pch,bg=col,...)
}

tplot <- function(x,xlab="",ylab="",...){
  plot(as.numeric(names(x)),x,xlab=xlab,ylab=ylab,...)
}

scatterBox <- function(x,y,cuts=100,...) boxplot(split(y,cut(x,unique(quantile(x,seq(0,1,length=cuts+1))),include.lowest=TRUE)),...)

  
as.fumeric<- function(x,...) as.numeric(as.factor(x,...))
first<- function(x,k=1) x[k]
###CHECK THIS ~pmurakam/feinberg/CharmFiles/functions.R
boyorgirl <- function(A,xIndex,yIndex,plot=FALSE,id=1:ncol(A)){
  ##A is average of log int
  x=Biobase::rowMedians(t(A[xIndex,]),na.rm=TRUE)
  y=Biobase::rowMedians(t(A[yIndex,]),na.rm=TRUE)
  k=kmeans(y-x,c(min(y-x),max(y-x)))
  sex=factor(ifelse(k$cluster==which.min(k$centers),"F","M"),levels=c("M","F"))
  if(plot){
    plot(x,y,type="n")
    text(x,y,id,col=as.numeric(sex))
    legend("bottomleft",c("M","F"),col=c(1,2),pch=1)
  }
  return(sex)
}

splitit <- function(x) split(seq(along=x),x)
myplclust <- function( hclust, lab=hclust$labels, lab.col=rep(1,length(hclust$labels)), hang=0.1, ... )
{
  ## modifiction of plclust for plotting hclust objects *in colour*!
  ## Copyright Eva KF Chan 2009
  ## Arguments:
  ##    hclust:    hclust object
  ##    lab:        a character vector of labels of the leaves of the tree
  ##    lab.col:    colour for the labels; NA=default device foreground colour
  ##    hang:     as in hclust & plclust
  ## Side effect:
  ##    A display of hierarchical cluster with coloured leaf labels.
  y <- rep(hclust$height,2)
  x <- as.numeric(hclust$merge)
  y <- y[which(x<0)]
  x <- x[which(x<0)]
  x <- abs(x)
  y <- y[order(x)]
  x <- x[order(x)]
  plot( hclust, labels=FALSE, hang=hang, ... )
    text( x=x, y=y[hclust$order]-(max(hclust$height)*hang), labels=lab[hclust$order], col=lab.col[hclust$order], srt=90, adj=c(1,0.5), xpd=NA, ... )
  
}


rowMads <- function(x,constant = 1.4826)constant*Biobase::rowMedians(abs(x-Biobase::rowMedians(x)))

getSegments <- function(x,factor,cutoff=quantile(abs(x),0.99),verbose=TRUE){

  Indexes=split(seq(along=x),factor)
  regionID=vector("numeric",length(x))
  LAST = 0
  
  segmentation = vector("numeric", length(x))
  type = vector("numeric", length(x))

  for (i in seq(along = Indexes)) {
    if (verbose) if (i%%1000 == 0) cat(".")
    Index = Indexes[[i]]
    y = x[Index]
    z = sign(y) * as.numeric(abs(y) > cutoff)
    w = cumsum(c(1, diff(z) != 0)) + LAST
    segmentation[Index] = w
    type[Index] = z
    LAST = max(w)
  }
  if(verbose) cat("\n")
  ##add a vector of the pns
  res=list(upIndex=split(which(type>0),segmentation[type>0]),
    dnIndex=split(which(type<0),segmentation[type<0]),
    zeroIndex=split(which(type==0),segmentation[type==0]))
  names(res[[1]])<-NULL
  names(res[[2]])<-NULL
  names(res[[3]])<-NULL
  return(res)
}

##you can pass cutoff through the ...
regionFinder<-function(x,regionNames,chr,position,y=x,
                       summary=mean,ind=seq(along=x),order=TRUE,oneTable=TRUE,
                       ...){

  Indexes=getSegments(x[ind],regionNames[ind],...)
  
  res=vector("list",2)
  for(i in 1:2){
    res[[i]]=data.frame(chr=sapply(Indexes[[i]],function(Index) chr[ind[Index[1]]]),
         start=sapply(Indexes[[i]],function(Index) min(position[ind[Index]])),
         end=sapply(Indexes[[i]],function(Index) max(position[ind[Index]])),
         value=sapply(Indexes[[i]],function(Index) summary(y[ind[Index]])),
         area=sapply(Indexes[[i]],function(Index) abs(sum(y[ind[Index]]))),
         pns=sapply(Indexes[[i]],function(Index) regionNames[ind[Index]][1]),
         indexStart=sapply(Indexes[[i]],function(Index) min(ind[Index])),
         indexEnd=sapply(Indexes[[i]],function(Index) max(ind[Index])))
    res[[i]]$L=res[[i]]$indexEnd-res[[i]]$indexStart+1
  }
  names(res)=c("up","dn")
  if(order & !oneTable){
    if(nrow(res$up)>0) res$up=res$up[order(-res$up$area),]
    if(nrow(res$dn)>0) res$dn=res$dn[order(-res$dn$area),]
  }
  if(oneTable){
    res=rbind(res$up,res$dn)
    if(order & nrow(res)>0) res=res[order(-res$area),]
  }
  return(res)
}

logit=function(x) log(x)-log(1-x)
ilogit=function(x) 1/(1+exp(-x))
##nearestgene needs cisgenome to be installed
nearestgene <- function(object,up=5000,down=5000,species="human",
                        build="hg18",
                        datadir="/dexter/disk4/gcode/hji/cisreg/genomes/",
                        wd="./cisgenome",
                        binary="/home/bst/faculty/hji/projects/cisgenome_project/bin/refgene_getnearestgene",
                        inputFileName="tmp.cod",
                        outputFileName="tmp.txt",
                        r=0){##is the type of distance
  tmpscipen=.Options$scipen
  options(scipen=100)
  
  if(!file.exists(wd)) {
    cat("\nCreating directory:",wd,"\n")
    dir.create(wd)
  }
  
  inputfile=file.path(wd,inputFileName)
  outputfile=file.path(wd,outputFileName)
  write.table(data.frame(object$chr,object$start,object$end,rep("+",nrow(object))),file=inputfile,sep="\t",quote=FALSE,col.names=FALSE)

  thecall=paste(binary,"-d",file.path(datadir,species,build,"annotation/refFlat_sorted.txt"),"-dt 1 -s",species,"-i",inputfile,"-o",outputfile,"-r", r,"-up",up,"-down",down)
  system(thecall)
  res=read.delim(outputfile,as.is=TRUE,check.names=FALSE,comment.char="#",header=FALSE,col.names=c("seq_id","chr","start","end","strand","name","annotation","num","genestrand","TSS","TSE","CSS","CSE","ne","exonStarts","exonEnds"),fill=NA)
  res[res=="---"]<-NA
  res[res==""] <- NA
  options(scipen=tmpscipen)
  res
}

## This is the new regionMatch function as of 2/18/11:
regionMatch <- function(object1, object2, verbose=TRUE) {
    require(IRanges)
    m1 = object1; m2 = object2 #keep arguments named object1,object2 for backward compatibility
    stopifnot(all(c("chr","start","end")%in%colnames(m1)))
    stopifnot(all(c("chr","start","end")%in%colnames(m2)))
    if(any(is.na(m1[,c("chr","start","end")])) | any(is.na(m2[,c("chr","start","end")]))) stop("No missing values allowed in chr, start, or end columns.")
    m1$chr = as.character(m1$chr)
    m2$chr = as.character(m2$chr)
    stopifnot(is.numeric(m1$start) & is.numeric(m1$end))
    stopifnot(is.numeric(m2$start) & is.numeric(m2$end))
    stopifnot(all(m1$end>=m1$start))
    stopifnot(all(m2$end>=m2$start))
              
    ret = matrix(NA,nrow=nrow(m1),ncol=7)
    colnames(ret) = c("dist","matchIndex","type","amountOverlap","insideDist","size1","size2")
    sp1 = split(1:nrow(m1), m1$chr)
    sp2 = split(1:nrow(m2), m2$chr)          
    for(i in names(sp1)) {
        if(verbose) cat(i," ")
        inds1 = sp1[[i]]
        if(i %in% names(sp2)) {
            inds2 = sp2[[i]]
            m1IR = IRanges(start=m1$start[inds1],end=m1$end[inds1])
            m2IR = IRanges(start=m2$start[inds2],end=m2$end[inds2])
            ret[inds1,"matchIndex"] = inds2[ nearest(m1IR,m2IR) ]
        } else ret[inds1,"matchIndex"] = NA
    }
    if(verbose) cat("\n")
    m2b = m2[ret[,"matchIndex"],]
    ret = data.frame(ret, stringsAsFactors=FALSE)
    
    ret$type[(m1$start>m2b$start & m1$end<=m2b$end) |
             (m1$start>=m2b$start & m1$end<m2b$end)] = "inside"
    ret$type[m1$start<=m2b$start & m1$end>=m2b$end] = "cover"
    ret$type[m1$start>m2b$end] = "disjointR"
    ret$type[m1$end<m2b$start] = "disjointL"
    ret$type[is.na(ret$matchIndex)] = "disjoint"
    ret$type[(m1$start>m2b$start & m1$start<=m2b$end) & m1$end>m2b$end] = "overlapR"    
    ret$type[m1$start<m2b$start & (m1$end>=m2b$start & m1$end<m2b$end)] = "overlapL"

    ret$dist = 0
    ret$dist[ret$type=="disjoint"]  = NA
    ret$dist[ret$type=="disjointR"] = m2b$end[ret$type=="disjointR"] - m1$start[ret$type=="disjointR"]
    ret$dist[ret$type=="disjointL"] = m2b$start[ret$type=="disjointL"] - m1$end[ret$type=="disjointL"]
    ret$amountOverlap[ret$type=="overlapR"] = -1*(m2b$end[ret$type=="overlapR"]-m1$start[ret$type=="overlapR"]+1)
    ret$amountOverlap[ret$type=="overlapL"] = m1$end[ret$type=="overlapL"]-m2b$start[ret$type=="overlapL"]+1
    ret$type[ret$type%in%c("disjointR","disjointL")] = "disjoint"
    ret$type[ret$type%in%c("overlapR","overlapL")] = "overlap"

    ## insideDist column:
    insideIndex = ret$type=="inside" #no missing ret$type at this point
    tmp0 = cbind(m1$end[insideIndex]  -m2b$end[insideIndex],
                 m1$start[insideIndex]-m2b$start[insideIndex])
    tmp = apply(abs(tmp0),1,which.min)
    tmpinsidedist = tmp0[,1]
    tmpIndex = tmp==2
    tmpinsidedist[tmpIndex] = tmp0[tmpIndex,2]
    ret$insideDist[insideIndex] = tmpinsidedist

    ## size1 and size2 columns:
    ret$size1 = m1$end -m1$start +1
    ret$size2 = m2b$end-m2b$start+1
    
    ret
}

###COMMENTED OUT RAFA's Original. Put in Peter's
## regionMatch <- function(object1,object2,verbose=TRUE){
##   Indexes1=split(seq(along=as.character(object1$chr)),as.character(object1$chr))
##   Indexes2=split(seq(along=as.character(object2$chr)),as.character(object2$chr))
##   res=matrix(NA,nrow=nrow(object1),7)
##   for(i in seq(along=Indexes1)){
##     CHR=names(Indexes1)[i]
##     if(verbose) cat(CHR,",")
##     Index1=Indexes1[[CHR]]
##     Index2=Indexes2[[CHR]]
##     if(length(Index1)>0 & length(Index2)>0){

##       tmpmat1=object1[Index1,c("start","end")]
##       tmpmat2=object2[Index2,c("start","end")]

##       tmp11=fuzzy.match2(tmpmat1[,1],tmpmat2)
##       tmp12=fuzzy.match2(tmpmat1[,2],tmpmat2)

##       tmp21=fuzzy.match2(tmpmat2[,1],tmpmat1)
##       tmp22=fuzzy.match2(tmpmat2[,2],tmpmat1)

##       tmpIndex = which(tmp21[,1]==0 & tmp22[,1]==0)
##       if(length(tmpIndex)>0){
##         tmpIndex2 = tmpIndex[ tmp21[tmpIndex,2]==tmp22[tmpIndex,2] ]
##         if(length(tmpIndex2)>0) insideIndex2= tmp21[tmpIndex2,2] else insideIndex2=c()
##       } else insideIndex2=c()
      
##       tmpIndex = which(tmp11[,1]==0 & tmp12[,1]==0)
##       if(length(tmpIndex)>0){
##         insideIndex=tmpIndex[
##         (object1$start[Index1[tmpIndex]] >=
##          object2$start[Index2[tmp11[tmpIndex,2]]] &
##          object1$end[Index1[tmpIndex]] <=
##          object2$end[Index2[tmp11[tmpIndex,2]]] ) |
##         (object1$start[Index1[tmpIndex]] >=
##          object2$start[Index2[tmp12[tmpIndex,2]]] &
##          object1$end[Index1[tmpIndex]] <=
##          object2$end[Index2[tmp12[tmpIndex,2]]] )
##           ]          
##       } else insideIndex=c()
      
##       overlapIndex =  setdiff(which(tmp11[,1]==0 | tmp12[,1]==0),insideIndex)
        
##       coverIndex = which((sign(tmp11[,1])==1 & sign(tmp12[,1])==-1 &
##         tmpmat2[tmp11[,2],1] <= tmpmat2[tmp12[,2],1]))
      
##       coverIndex = union(coverIndex,insideIndex2)

## ###Minimum distance and Index of minimum
##       tmp3=cbind(tmp11[,1],tmp12[,1])
##       tmp=apply(abs(tmp3),1,which.min)
##       tmpdist=tmp11[,1]
##       tmpIndex=tmp==2
##       tmpdist[tmpIndex]=tmp12[tmpIndex,1]
##       tmpdist[coverIndex]=0

##       tmpMatch1=Index2[tmp11[,2]]
##       tmpMatch2=Index2[tmp12[,2]]
##       tmpMatch=tmpMatch1
##       tmpMatch[tmpIndex]=tmpMatch2[tmpIndex]
      
##       ###overlap
##       overlap1 = (tmp11[overlapIndex,1]==0)*
##         (object1$start[Index1[overlapIndex]] -
##          object2$end[tmpMatch1[overlapIndex]])
                    

##       overlap2 = (tmp12[overlapIndex,1]==0)*
##         (object1$end[Index1[overlapIndex]]-
##          object2$start[tmpMatch2[overlapIndex]])

##       overlap = pmax(abs(overlap1),abs(overlap2))
##       ##where inside..
##       tmp3=cbind(object1$end[Index1[insideIndex]] - object2$end[tmpMatch[insideIndex]],object1$start[Index1[insideIndex]] - object2$start[tmpMatch[insideIndex]])
##       tmp=apply(abs(tmp3),1,which.min)
##       tmpinsidedist=tmp3[,1]
##       tmpIndex=tmp==2
##       tmpinsidedist[tmpIndex]=tmp3[tmpIndex,2]
    
##       res[Index1,1]<-tmpdist
##       res[Index1,2]<-tmpMatch
##       res[Index1,3]<-0;if(length(insideIndex)>0)res[Index1[insideIndex],3] <- 1
##       if(any(insideIndex)){
##         res[Index1[insideIndex],4] <- tmpinsidedist
##       }
##       res[Index1,5]<-0;if(length(overlapIndex)>0)res[Index1[overlapIndex],5] <- 1
##       if(any(overlapIndex)){
##         res[Index1[overlapIndex],6]<-overlap
##       }
##       res[Index1,7]<-0;if(length(coverIndex)>0)res[Index1[coverIndex],7] <- 1
##     }
##   }
##   type=rep("disjoint",nrow(res))
##   type[res[,3]==1] <- "inside"
##   type[res[,5]==1] <- "overlap"
##   type[res[,7]==1] <- "cover"
##   type=factor(type,levels=c("disjoint","overlap","inside","cover"))
##   size1=object1$end-object1$start+1
##   size2=object2$end[res[,2]] - object2$start[res[,2]] + 1
##   res=data.frame(dist=res[,1],matchIndex=res[,2],type=type,amountOverlap=res[,6],insideDist=res[,4],size1=size1,size2=size2)
##   res
## }


pointMatch = function(object1,object2,col1=2,col2=2,verbose=TRUE){
  Indexes1=split(seq(along=as.character(object1$chr)),as.character(object1$chr))
  Indexes2=split(seq(along=as.character(object2$chr)),as.character(object2$chr))
  res=matrix(NA,nrow=nrow(object1),2)
  for(i in seq(along=Indexes1)){
    CHR=names(Indexes1)[i]
    if(verbose) cat(CHR,",")
    Index1=Indexes1[[CHR]]
    Index2=Indexes2[[CHR]]
    if(length(Index1)>0 & length(Index2)>0){
      tmpmat1=object1[Index1,col1]
      tmpmat2=object2[Index2,col2]
      tmp=fuzzy.match(tmpmat1,tmpmat2)
      tmp[,2] = Index2[tmp[,2]]
	  
      res[Index1,]=tmp    
     }
  }
  colnames(res)<-c("dist","matchIndex")
  return(res)
}


##old pointMatch. Jaffe fixed a bug
## pointMatch <- function(object1,object2,col1=2,col2=2,verbose=TRUE){
##   Indexes1=split(seq(along=as.character(object1$chr)),as.character(object1$chr))
##   Indexes2=split(seq(along=as.character(object2$chr)),as.character(object2$chr))
##   res=matrix(NA,nrow=nrow(object1),2)
##   for(i in seq(along=Indexes1)){
##     CHR=names(Indexes1)[i]
##     if(verbose) cat(CHR,",")
##     Index1=Indexes1[[CHR]]
##     Index2=Indexes2[[CHR]]
##     if(length(Index1)>0 & length(Index2)>0){
##       tmpmat1=object1[Index1,col1]
##       tmpmat2=object2[Index2,col2]
##       res[Index1,]=fuzzy.match(tmpmat1,tmpmat2)
##     }
##   }
##   colnames(res)<-c("dist","matchIndex")
##   return(res)
## }


getSegments <- function(x,factor,cutoff=quantile(abs(x),0.99),verbose=TRUE){

  Indexes=split(seq(along=x),factor)
  regionID=vector("numeric",length(x))
  LAST = 0
  
  segmentation = vector("numeric", length(x))
  type = vector("numeric", length(x))

  for (i in seq(along = Indexes)) {
    if (verbose) if (i%%1000 == 0) cat(".")
    Index = Indexes[[i]]
    y = x[Index]
    z = sign(y) * as.numeric(abs(y) > cutoff)
    w = cumsum(c(1, diff(z) != 0)) + LAST
    segmentation[Index] = w
    type[Index] = z
    LAST = max(w)
  }
  ##add a vector of the pns
  res=list(upIndex=split(which(type>0),segmentation[type>0]),
    dnIndex=split(which(type<0),segmentation[type<0]),
    zeroIndex=split(which(type==0),segmentation[type==0]))
  names(res[[1]])<-NULL
  names(res[[2]])<-NULL
  names(res[[3]])<-NULL
  return(res)
}

matchGenes <- function(object,promoterDist=2500,verbose=TRUE,species="human",
                       build="hg18",
                       r=0,up=50000000,down=50000000){
  if(verbose) cat("Matching regions to genes.\n")
  genes=nearestgene(object,up=up,down=down,species=species,build=build,r=r)
  type=rep("",nrow(genes))
  subtype=rep("",nrow(genes))
  ctype=rep("",nrow(genes))
  dist=rep(0,nrow(genes))
  insidedistance<-rep(NA,nrow(genes))
  exonnumber<-rep(NA,nrow(genes))
  nexons<-rep(NA,nrow(genes))
  
  geneL=rep(0,nrow(genes))
  codingL=rep(0,nrow(genes))
  subdist=rep(0,nrow(genes))
  if(verbose) cat("Annotating results")
  for(i in 1:nrow(genes)){
    if(verbose & i%%1000==0)  cat(".")

    if(!is.na(genes[i,10])){
    TS = genes[i,10]
    TE = genes[i,11]
    geneL[i]=TE-TS
    CS = genes[i,12]
    CE = genes[i,13]
    codingL[i]=CE-CS
    ES = as.numeric(strsplit(genes[i,15],",")[[1]])
    EE = as.numeric(strsplit(genes[i,16],",")[[1]])
    Exons= cbind(ES,EE)
    nexons[i]=nrow(Exons)
    Strand= ifelse(genes[i,9]=="+",1,-1)
    S = genes[i,3]
    E = genes[i,4]
    
    type[i]=""
    if(genes[i,3] <= TS & genes[i,4] >= TE){
      type[i]="covers"
    } else{
      if(genes[i,3] < TS & genes[i,4] < TS){
        if(Strand==1){
          type[i]="upstream" 
          dist[i]=TS-E
        } else{
          type[i]="downstream"
          dist[i]=TE-E
        }
      }
      if(genes[i,3] > TE & genes[i,4] > TE){
        if(Strand==-1){
          type[i]="upstream"
          dist[i]=S-TE
        }  else{
          type[i]="downstream"
          dist[i]=S-TS
        }
      }
      if (genes[i,3]>= TS & genes[i,4] <= TE){
        type[i]="inside"
        if(Strand==-1) dist[i]=TE-E  else dist[i]=S-TS
      }
      
      ##overlaps
      if(type[i]==""){
        if(genes[i,3] < TS & genes[i,4] <=TE){
          ##OVERLAP FRONT
          if(Strand==1) type[i]="overlaps 5'" else{
            type[i]="overlaps 3'"
            dist[i]=TE-E
          }
          S=TS
        }
        if (genes[i,3] >= TS & genes[i,4] > TE){
          ##OVERAL BACK
          if(Strand==-1) type[i]="overlaps 5'" else{
            type[i]="overlaps 3'"
            dist[i]=S-TS
          }
          E=TE
        }
      }
    }
    m1=NA;m2=NA
    if(S >= TS & E <= TE){
      
      ##INSIDE
      tmp1=fuzzy.match2(S,Exons)
      tmp2=fuzzy.match2(E,Exons)
      
      m1=tmp1[1,1]
      m2=tmp2[1,1]
      exon1=tmp1[1,2]
      exon2=tmp2[1,2]
      m1m2Index=which.min(abs(c(m1,m2)))
      
      if(exon1==exon2 & m1==0 & m2==0){
        subtype[i]="inside exon"
        insidedistance[i]=0
        exonnumber[i]=c(exon1,exon2)[m1m2Index]
      }
      
      if( (sign(m1)==sign(m2) & (m1!=0 & m2!=0) & exon1==exon2) |
         (sign(m1)==-1 & sign(m2)==1 & exon2-exon1==1) ){
        subtype[i]="inside intron"

        insidedistance[i]=c(m1,m2)[m1m2Index]
        exonnumber[i]=c(exon1,exon2)[m1m2Index]
        
      }
      if( (exon2-exon1 > 1) |
         ( (exon2-exon1 == 1) & 
          ((sign(m1)==-1 & sign(m2)==-1) |
           (sign(m1)==1 & sign(m2)==-1) |
           (sign(m1)==1 & sign(m2)==1) |
           (sign(m1)==0 & sign(m2)==-1)|
           (sign(m1)==1 & sign(m2)==0))) |
         (exon2==exon1 & sign(m1)==1 & sign(m2)==-1)){
        subtype[i]="covers exon(s)"
        insidedistance[i]=0
        exonnumber[i]=c(exon1,exon2)[m1m2Index]
      }
      
      if( (exon2-exon1 == 1 & sign(m1)==-1 & sign(m2)==0) |
         (exon2==exon1 & sign(m1)==1 & sign(m2)==0)){
        if(Strand==1) subtype[i]="overlaps exon upstream" else subtype[i]="overlaps exon downstream"
        insidedistance[i]=0
        exonnumber[i]=c(exon1,exon2)[m1m2Index]
      }
      if( (exon2-exon1 == 1 & sign(m1)==0 & sign(m2)==1) |
         (exon2==exon1 & sign(m1)==0 & sign(m2)==-1)){
        if(Strand==-1) subtype[i]="overlaps exon upstream" else subtype[i]="overlaps exon downstream"
        insidedistance[i]=0
        exonnumber[i]=c(exon1,exon2)[m1m2Index]
      }
      if( exon2-exon1 == 1 & sign(m1)==0 & sign(m2)==0){
        subtype[i]="overlaps two exons"
        insidedistance[i]=0
        exonnumber[i]=c(exon1,exon2)[m1m2Index]
      }

      if(Strand!=1){
        insidedistance[i]= -insidedistance[i]
        exonnumber[i] = nrow(Exons) - exonnumber[i] + 1
      }
      
      ctype[i]="inside transcription region"
      
      if(S<CS & E<CS){
        if(Strand==1) ctype[i]="5' UTR" else ctype[i]="3'UTR"
      }
      
      if(S>CE & E>CE){
        if(Strand==-1) ctype[i]="5' UTR" else ctype[i]="3'UTR"
      }
      if(S<CS & E>CE){
        ctype[i]="covers transcription region"
      }
      if(S<CS & E>CS){
        if(Strand==1) ctype[i]="overlaps 5' UTR" else ctype[i]="overlaps 3'UTR"
      }
      if(S<CE & E>CE){
        if(Strand==-1) ctype[i]="overlaps 5' UTR" else ctype[i]="overlaps 3'UTR"
      }
    }
    if(FALSE){##graphical check
      plot(0,0,ylim=c(0,1),xlim=range(genes[i,c(3:4,10:11)]),xlab=paste("inside distance=",insidedistance[i],m1,m2,exonnumber[i]))
      polygon(c(TS,TE,TE,TS),c(0,0,0.5,0.5),density=0,col=2)
      polygon(c(CS,CE,CE,CS),c(0.1,0.1,0.4,0.4),density=0,col=3)
      abline(h=0.25,lwd=2)
      apply(Exons,1,function(x) polygon(c(x[1],x[2],x[2],x[1]),
                                        c(0.2,0.2,0.3,0.3),col=4))
      polygon(as.vector(genes[i,c(3,4,4,3)]),c(0.75,0.75,0.85,0.85),col=5)
      lines(c(TS,TS+1000),c(0.65,0.65),lwd=3)
      title(paste(i,Strand,type[i],subtype[i],ctype[i],dist[i],sep=":"))

      
    }
  }
  }
  if(verbose) cat("Done.\n")
  type[dist<=promoterDist & type=="upstream"] <- "promoter"
  type[dist<=promoterDist & type=="downstream"] <- "close to 3'"

  description=type
  tmpIndex=which(description=="inside")
  description[tmpIndex] <- subtype[tmpIndex]
  return(data.frame(name=I(genes[,6]),
                    annotation=I(genes[,7]),
                    description=factor(description,levels=c("upstream","promoter","overlaps 5'","inside intron","inside exon","covers exon(s)","overlaps exon upstream","overlaps exon downstream","overlaps two exons","overlaps 3'","close to 3'","downstream","covers")),
                    region=factor(type,levels=c("upstream","promoter","overlaps 5'","inside","overlaps 3'","close to 3'","downstream","covers")),
                    distance=dist,
                    subregion=factor(subtype,levels=c("inside intron","inside exon","covers exon(s)","overlaps exon upstream","overlaps exon downstream","overlaps two exons")),
                    insidedistance=insidedistance,
                    exonnumber=exonnumber,
                    nexons=nexons,
                    UTR=factor(ctype,levels=c("inside transcription region","5' UTR","overlaps 5' UTR","3'UTR","overlaps 3'UTR","covers transcription region")),
                    strand=genes[,9],
                    geneL=geneL,
                    codingL=codingL))
}

dogo <- function(names,universe,species="human"){
  if(species=="human"){
    golib="org.Hs.eg.db"
    library(golib,character.only=TRUE)
    gomap= org.Hs.egREFSEQ2EG
  } else  {
    golib="org.Mm.eg.db"
    library(golib,character.only=TRUE)
    gomap= org.Mm.egREFSEQ2EG
  }
  require(GOstats)
  x=unlist(mget(as.character(names), gomap,ifnotfound = NA))
  x=x[!is.na(x)]
  Universe=unlist(mget(as.character(universe),gomap,ifnotfound = NA))
  
  Universe=unique(c(Universe[!is.na(Universe)],unique(x)))
  
  params <- new("GOHyperGParams", geneIds = unique(x),
                universeGeneIds = Universe,
                annotation = golib,
                ontology = "BP", pvalueCutoff = 0.01, conditional = TRUE,
                testDirection="over")
  ht=hyperGTest(params)
  tab=summary(ht)
  tmp1=geneIdsByCategory(ht)
  tmp1=tmp1[tab[,1]]
  tab$IDs=sapply(tmp1,function(y) paste(names(x)[x%in%y],collapse=";"))
  return(tab)
}




affyGO <- function(names,universe=NULL,platform="hgu133plus2",species="human",pvalueCutoff=1){
  library(paste(platform,"db",sep="."),character.only=TRUE)
  gomap=get(paste(platform,"ENTREZID",sep=""))
  symmap=get(paste(platform,"SYMBOL",sep=""))
  if(species=="human"){
    golib="org.Hs.eg.db"
    library(golib,character.only=TRUE)
  } else  {
    golib="org.Mm.eg.db"
    library(golib,character.only=TRUE)
    gomap= org.Mm.egREFSEQ2EG
  }
  require(GOstats)
  x=unlist(mget(as.character(names), gomap,ifnotfound = NA))
  x=x[!is.na(x)]

  if(is.null(universe)) universe=ls(gomap)
  Universe=unlist(mget(universe,gomap,ifnotfound = NA))
  Universe=unique(c(Universe[!is.na(Universe)],unique(x)))
  
  params <- new("GOHyperGParams", geneIds = unique(x),
                universeGeneIds = Universe,
                annotation = golib,
                ontology = "BP", pvalueCutoff = pvalueCutoff,
                conditional = TRUE,
                testDirection="over")
  ht=hyperGTest(params)
  tab=summary(ht)
  tmp1=geneIdsByCategory(ht)
  tmp1=tmp1[tab[,1]]
  tab$IDs= sapply(tmp1,function(y) paste(unlist(mget(names(x)[x%in%y],symmap,ifnotfound=NA)),collapse=";"))
  return(tab)
}

myfilter2 <- function(x,filter,...){
###unlike myfilter, myfilter2 returns NAs in the edges.
  L=dim(x)[1]
  if(L>length(filter)) res=filter(x,filter,...) else{res = t(colMeans(x) %*% t(rep(1,dim(x)[1])))}

  return(res)

}


getDesc<-function(x){
 ##get gene description
  require(org.Hs.eg.db)
  genenames=sapply(mget(as.character(x),org.Hs.egREFSEQ2EG,ifnotfound = NA),function(x) x[1])
  genenames[!is.na(genenames)]= sapply(mget(genenames[!is.na(genenames)],org.Hs.egGENENAME,ifnotfound=NA),function(x) x[1])
  genenames
}


clusterMaker <- function(chr,pos,order.it=TRUE,maxGap=300){
  nonaIndex=which(!is.na(chr) & !is.na(pos))
  Indexes=split(nonaIndex,chr[nonaIndex])
  clusterIDs=rep(NA,length(chr))
  LAST=0
  for(i in seq(along=Indexes)){
    Index=Indexes[[i]]
    x=pos[Index]

    if(order.it){ Index=Index[order(x)];x=pos[Index] }

    y=as.numeric(diff(x)>maxGap)
    z=cumsum(c(1,y))
    clusterIDs[Index]=z+LAST
    LAST=max(z)+LAST
  }
  clusterIDs
}  

midpoints <- function(x) (x[1:(length(x)-1)]+x[2:length(x)])/2

pointsplot <- function(l,jitter=TRUE,factor=1,col=rep(1,length(l)),pch=21,...){
  if(!is.list(l)) stop("l must be a list\n")
  y=unlist(l)
  x=rep(seq(along=l),sapply(l,length))
  col=rep(col,sapply(l,length))
  if(jitter) x=jitter(x,factor)
  if(pch!=21) plot(x,y,col=col,pch=pch,...) else plot(x,y,pch=pch,bg=col,...)
}

tplot <- function(x,xlab="",ylab="",...){
  plot(as.numeric(names(x)),x,xlab=xlab,ylab=ylab,...)
}

scatterBox <- function(x,y,cuts=100,...) boxplot(split(y,cut(x,unique(quantile(x,seq(0,1,length=cuts+1))),include.lowest=TRUE)),...)


maplot <- function(x,y,...){ A=(x+y)/2;M=y-x;plot(A,M,...)}
                             
                             

fuzzy.match = function(x, y, x.sorted = FALSE, y.sorted = FALSE) {

	l.x = length(x)
	# space for results
	r = array(0,c(l.x,2))

	if (! x.sorted) {
		# sort x
		perm = order(x)
		X = x[perm]
	}

	l.y = length(y)
	if (! y.sorted) {
		Y = array(0,1+l.y)
		# sort y
		perm.y = order(y)
		Y[1:l.y] = y[perm.y]
		Y[1+l.y] = +Inf
	}
	else
		y[1+l.y] = +Inf

	# set up
	v = -Inf
	j = 1
	if (y.sorted)
		u = y[j]
	else
		u = Y[j]
	# run over X, or x if x.sorted
	for (i in 1:l.x) {
		if (x.sorted)
			z = x[i]
		else
			z = X[i]
		e = z - v
		d = u - z
		while (d < 0) {
			v = u
			e = z - v
			j = j + 1
			if (y.sorted)
				u = y[j]
			else
				u = Y[j]
			d = u - z
		}
		# result
		if (x.sorted)
			k = i
		else
			k = perm[i]		
		if (y.sorted)
			r[k,2] = j
		else
			r[k,2] = perm.y[j]
		if (d == 0)
			r[k,1] = 0
		else if (e < d) {
			r[k,1] = -e
			if (y.sorted)
				r[k,2] = j-1
			else
				r[k,2] = perm.y[j-1]
		}
		else # e >= d
			r[k,1] = d
	}
	r
}

fuzzy.match2 = function(x, z, x.sorted = FALSE, z.sorted = FALSE) {
	# z is an Nx2 matrix whose rows are separate intervals
	# which play the role of the y-elements in fuzzy.match
	if (z.sorted)
		y = t(z)
	else {
		perm = order(z[,1])
		y = t(z[perm,])
	}
	f = fuzzy.match(x,y,x.sorted=x.sorted,y.sorted=TRUE)
	for (i in 1:length(x)) {
		d = f[i,1]	# distance from as.vector(y)
		j = f[i,2]	# nearest index of as.vector(y)
		jmod2 = j%%2
		if (z.sorted)
			f[i,2] = (j+jmod2)/2
		else
			f[i,2] = perm[(j+jmod2)/2]
		if ((d>0 && jmod2==0) || (d<0 && jmod2 == 1))
			f[i,1] = 0
	}
	f
}
#######################################################
#####
#####  Bump hunting stuff
#####
#######################################################

loessByGroup <- function(stat ,group, pos, se=rep(1,length(stat)),
                                             numProbes = 5) {
  require(limma)
  sstat = stat
  Indexes = split(seq(along=group),group)
  
  for(i in seq(along=Indexes)) {
    if(i %% 1e4 == 0) cat(".")
    Index = Indexes[[i]]
    ss = numProbes/length(Index) # make a span
    if(ss > 1) ss = 0.75
    
    if(length(Index) > 4) {
      
      sstat[Index] <- loessFit(stat[Index],pos[Index], 
                               span = ss, weights = 1/se[Index])$fitted
    } else sstat[Index] = median(stat[Index])
  }
  return(sstat)
}


