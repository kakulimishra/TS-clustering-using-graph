##This is the final copy to extract the temporal patterns. The other copies are the older copies.
###The file initially finds the temporal patterns and saves into folders for respective datasets and then computes the
##support values. The results for support values are reported in the paper in tale 3.
##The function  'plotVals' to plot the support values was also created but later unused in the paper.


library(stringr)
library(TSclust)
# library(clusteval)
library(qualV)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(plyr)
library(reshape2)



plotVals=function(supportVals)
{
  ## Although the files saved have two measures- patternLength and support. I have used only the support for the
  ## comparative analysis of the three different graph based representation.
  
  supportVals=data.frame(supportVals)
  supportVals=supportVals[!(supportVals$temp.file=="Temporal patterns weighted"),]
  
  sv=as.data.frame(supportVals[,-c(3,6)])
  #sv$patternLength=round(as.numeric(as.character(sv$patternLength)),2)
  sv$support=round(as.numeric(as.character(sv$support)),2)
  sv$window=as.numeric(as.character(sv$window))
  #copy.supportVals=melt(sv,c(1,2,5))
  copy.supportVals=melt(sv,c(1,2,4))
  
  
  myColor=c("darkgreen","blue","deeppink4","red")
  plot_labeller <- c("Support")
  
  
  myLabels=c("Deg-WG","Ei-UWG","Deg-UWG","P") ##This line should be checked each time the code is changed.
  names(plot_labeller) <- unique(as.character(copy.supportVals$variable))
  
  # p1=ggplot(data=copy.supportVals)+geom_line(data=copy.supportVals,aes(x=window,y=value,
  #              color=factor(temp.file)))+theme_bw()+ geom_point(data=copy.supportVals,aes(x=window,y=value,
  #             shape=factor(temp.file),color=factor(temp.file)),stroke=1)+
  #   scale_color_manual(values=myColor,labels=myLabels)+ scale_shape_manual(values=c(5,6,7,8),labels=myLabels)+
  #   facet_grid(rows = vars(dl), cols = vars(variable), 
  #     scales="free",labeller=labeller(variable=plot_labeller) )+
  #   theme(plot.title = element_text(hjust = 0.5,color="black",size=10),axis.ticks.y = element_blank(),
  #         legend.title = element_blank(),legend.position="bottom",panel.grid.major.y =element_blank(),
  #         legend.text=element_text(size=10),axis.title.x = element_text(size=15,color="black"),
  #         axis.title.y= element_blank(),
  #         panel.grid.major =   element_blank() , axis.text= element_text(size=15,color="black"),
  #         panel.background = element_rect(color="black",size=1.5), 
  #         strip.background =element_rect(color="black",fill="lightblue"),
  #         axis.ticks.x = element_blank(),plot.caption = element_text(hjust=0.5, size=rel(1.2)),
  #         strip.text = element_text(size = 12, colour = "black"),legend.key.width = unit(2,"cm"))+
  # facet_wrap(scales="free",.~variable ) 
  
  
  p1=ggplot(data=copy.supportVals)+geom_line(data=copy.supportVals,aes(x=window,y=value,
               color=factor(temp.file)))+theme_bw()+ geom_point(data=copy.supportVals,aes(x=window,y=value,
                                                                                                                                                  shape=factor(temp.file),color=factor(temp.file)),stroke=1)+
    scale_color_manual(values=myColor,labels=myLabels)+ scale_shape_manual(values=c(5,6,7,8),labels=myLabels)+
    theme(axis.ticks.y = element_blank(),
          legend.title = element_blank(),legend.position="bottom",panel.grid.major.y =element_blank(),
          legend.text=element_text(size=10),axis.title.x = element_text(size=15,color="black"),
          axis.title.y= element_blank(),
          panel.grid.major =   element_blank() , axis.text= element_text(size=15,color="black"),
          panel.background = element_rect(color="black",size=1.5), 
          axis.ticks.x = element_blank(),plot.caption = element_text(hjust=0.5, size=rel(1.2)),
          strip.text = element_text(size = 12, colour = "black"),legend.key.width = unit(2,"cm"))+
    scale_x_continuous(breaks=seq(0,max(copy.supportVals$window),12))
  
  
  
  return (p1)
  
}



summarization=function(dl, temp.file)
{
  ### The function is meant to summarize the temporal patterns and obtain the plots
  
  allFiles=list.files(file.path(p,dl,"Output non-overlap",temp.file))
  
  maxSupport=NULL
  for (z in 1:length(allFiles))
  {
    window=strsplit(allFiles[z],"win")[[1]][2]
    window=gsub(window,pattern=".csv",replacement="")
    print (paste0("Window-",window))
    
    
    tempDta=read.csv(file.path(p,dl,"Output non-overlap",temp.file,allFiles[z]),header=T)
    
    fileToRead=file.path(p,dl,"Output non-overlap",temp.file,allFiles[z])
    
    maxSupport=rbind(maxSupport,cbind(dl,temp.file,"patternLength"=mean(tempDta$patternLength),
                        "mean_support"=mean(tempDta$support),"min_support"=min(tempDta$support),"max_support"=max(tempDta$support),
                        "window"=window,"basename"=basename(fileToRead)))
  }
  
  print (maxSupport)
  
  return (maxSupport)
}




findPattern=function(dta,featureDist,symb,pathList,predLabel,window)
{
  temporalPatt=NULL
  temporalPatterns=NULL
  #loop for every cluster label
  for (k in 1:max(predLabel))
  {
    #extract the rows corresponding to "k"
    print (k)
    temp.symb=data.frame(symb[which(predLabel==k),],row.names=NULL)
    
    if (nrow(temp.symb)>1)
    {
      commonSegments=list()
      
      matches<-NULL
      z<-1
      apply(temp.symb[1:nrow(temp.symb),],1,function(outt) {
        
        outt=unlist(outt)
        commonSegments=list()
        l.lcs=LCS(outt,as.character(unlist(temp.symb[z,])))   ##Initialization, (trivial match)
        l.lcs=outt[l.lcs$va]
        commonSegments[[1]]=l.lcs
        for(i in 2:(nrow(temp.symb)))
        {
          l.lcs=LCS(unlist(commonSegments[[i-1]]),as.character(unlist(temp.symb[i,])))
          commonSegments[[i]]=commonSegments[[i-1]][l.lcs$va]
          
        }
        commonSegments=commonSegments[-1] ##The trivial match is removed, i.e, the initialization part
        
        lengths=(unlist(lapply(commonSegments,length)))
        
        if (any(lengths>0))
        {
          
          minLength=min(lengths[lengths>0])
          num.match=length(which(lapply(commonSegments,length)>0))
          
          temp.pattern=(unlist(lapply(commonSegments,length),use.names=F))
          temp.pattern=as.data.frame(t(commonSegments[[which.min(temp.pattern[temp.pattern>0])]]))
          colnames(temp.pattern)=seq(1,ncol(temp.pattern),1)
          
          pattLength=minLength/ncol(temp.symb)
          support=num.match/nrow(temp.symb)
          
          matches<<-rbind.fill(matches,cbind(z,"patternLength"=pattLength,"support"=support,temp.pattern))
        }

        z<<-z+1
        
      })
      
      ##Extract row which gives the highest pattern length in the dataframe matches.
      matches=matches[which.max(matches$patternLength),2:ncol(matches)]

      temporalPatterns=rbind.fill(temporalPatterns,cbind("cluster"=k, matches))

     }
    
  }
  
  
  return (temporalPatterns)
}



dataList=c("Ausgrid Data")
fileToRead=c("AusgridData_forGraph.csv")


###below to obtain the temporal patterns and save them

for (dl in 1:length(fileToRead))
{
  trainData=read.csv(file.path(getwd(),"data",fileToRead[dl]),header=T)
  dataNm=strsplit(dataList[dl]," ")[[1]][1]
  print (dataNm)
  dta=trainData[,-c(1,2,3)]
  
  
  pathList=list.files(file.path(getwd(),"Extract temporal patterns"),pattern="chosen")
  #pathList=pathList[!is.na(str_extract(string=pathList,pattern=dataFreq))]
  pathList1=(sapply(pathList,function(x) 
  {
    x=strsplit(x,"win")[[1]][2]
    gsub(x,pattern=".csv",replacement="")
  }))
  
  featureDist=list.files(file.path(getwd(),"Extract temporal patterns"),pattern="featureSpace")
  #featureDist=featureDist[!is.na(str_extract(string=featureDist,pattern=dataFreq))]
  featureDist1=(sapply(featureDist,function(x) 
  {
    x=strsplit(x,"win")[[1]][2]
    gsub(x,pattern=".csv",replacement="")
  }))
  
  
  symbolicRep=list.files(file.path(getwd(),"Extract temporal patterns"),pattern="symbolicRep")
  #symbolicRep=symbolicRep[!is.na(str_extract(string=symbolicRep,pattern=dataFreq))]
  symbolicRep1=(sapply(symbolicRep,function(x) 
  {
    x=strsplit(x,"win")[[1]][2]
    gsub(x,pattern=".csv",replacement="")
  }))
  
  
  predLabels=list.files(file.path(getwd(),"Extract temporal patterns"),pattern="ProposedLabels")
  #predLabels=predLabels[!is.na(str_extract(string=predLabels,pattern=dataFreq))]
  predLabels1=(sapply(predLabels,function(x) 
  {
    x=strsplit(x,"win")[[1]][2]
    gsub(x,pattern=".csv",replacement="")
  }))
  
  for(q in 1:length(pathList))
  {
    
    window=strsplit(pathList[q],"win")[[1]][2]
    window=gsub(window,pattern=".csv",replacement="")
    print (paste0("Window-",window))
    
    feature=read.csv(file.path(getwd(),"Extract temporal patterns",featureDist[which(featureDist1==window)]),header=T)
    listofPath=read.csv(file.path(getwd(),"Extract temporal patterns",pathList[which(pathList1==window)]),header=T)
    symb=read.csv(file.path(getwd(),"Extract temporal patterns",symbolicRep[which(symbolicRep1==window)]),header=T)
    predLabel=read.csv(file.path(getwd(),"Extract temporal patterns",predLabels[which(predLabels1==window)]),header=T)
    
    temporalPatt=findPattern(dta,feature,symb,listofPath,predLabel,window)
    
    write.csv(data.frame(temporalPatt),file.path(getwd(),"Extract temporal patterns",paste0("temporalPat_","win",window,".csv")),row.names = F)
    
  }
  
}




