###The file is used to plot Fig. 7 in the paper.  The function where the plots are obtained is "plot.patt2"

library(stringr)
library(TSclust)
library(clusteval)
library(qualV)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(plyr)
library(reshape2)


###The function "nonOverlap" is used hee to segment the dates and obtain the actual date and time values of the segments.
nonOverlap=function (timeseries, swSize) 
{
  ####  The function is used for obtaining the non-overlapping segments
  
  SW<-NULL
  timeseries=timeseries[1:(floor(length(timeseries)/swSize)*swSize)]
  for (x in seq(1,length(timeseries) ,swSize))
  {
    
    windowCol=(timeseries[x:(x+swSize-1)])
    SW <- cbind(SW, as.character(windowCol), deparse.level = 0)
  }
  SW=t(SW)
  SW=na.omit(SW)
  row.names(SW)<-NULL
  return (SW)
}





plot.patt2=function(df,llcs,winLength)
{

  
  plotList=NULL

  if (length(llcs)>1)
  {
    df=data.frame(df[llcs,])
  }

  else
  {
    df=data.frame(t(df[llcs,]))
  }
  
  data_to_plot  <-  NULL
  df$id=seq(1,nrow(df))
  apply(df,1,function(f)
  {
    ti=as.numeric(unlist(f[1:(length(f)-1)]))
    data_to_plot<<-rbind(data_to_plot,cbind(ti,as.numeric(rep(f[length(f)], length(ti)))))
  })
  
  
  # my_colors=c("black","blue","darkgreen","red","goldenrod4","gold","darkorchid1","palegreen1","sienna1",
  #             "bisque4","darkslategrey","chocolate4","deeppink2","deepskyblue","chartreuse2","darkmagenta","brown1")
  
  my_colors = c(brewer.pal(name="Dark2", n = 8), brewer.pal(name="Paired", n = 6), brewer.pal(name='Accent', n=8),
                brewer.pal(name="Set1", n = 8), brewer.pal(name="Set2", n = 8), brewer.pal(name="Set3", n = 8) )
  
  data_to_plot=data.frame(data_to_plot)
  
  data_to_plot$hr=seq(1,nrow(data_to_plot),1)
  
  colnames(data_to_plot)=c("val","id","hr")
  win=ncol(df)-1
  
  if (length(unique(data_to_plot$id))>1)
    p.plot=    ggplot(data=data_to_plot,aes(x=hr,y=val,color=factor(id)))+geom_path(group=1)+ theme_void()+ 
    geom_line(size=0.8)+geom_vline(xintercept=seq(1,nrow(data_to_plot),win),linetype='dotted')+
    theme(legend.position = 'none', legend.background = element_blank())
  
  else
    p.plot=    ggplot(data=data_to_plot,aes(x=hr,y=val,color=factor(id)))+geom_path(group=1)+ theme_void()+ 
    geom_line(size=2)+
    theme(legend.position = 'none', legend.background = element_blank())
    
    
    # scale_color_manual(values = my_colors)+ theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
    # scale_size_manual(values=my_size)+
    # theme(plot.title = element_text(hjust = 0.5),legend.title = element_blank(),
    #       legend.position = 'none',axis.title.y = element_blank(),axis.ticks.y = element_blank(),
    #       axis.text = element_text(colour='black'), 
    #       panel.grid.major.y =element_blank(),panel.grid.major.x = element_blank() ,
    #       axis.title = element_blank(),plot.caption = element_text(hjust=0.5, size=rel(1.2)))
    

  
  return (p.plot)
  
}


extract_tempPatt=function(dl,fileToRead)
{
  trainData=read.csv(file.path(p,dl,fileToRead),header=T)
  dataNm=strsplit(dl," ")[[1]][1]
  print (dataNm)
  
  if (dl=="Web traffic")
    dta=trainData[,-c(1,2)]
  
    
  dta=trainData[,-c(1,2,3)]
  
  listTemporal=list.files(file.path(p,dl,"Output non-overlap","Temporal patterns Ei weighted"))
  
  listSymbols=list.files(file.path(p,dl,"Output non-overlap","Cluster labels Ei weighted"),pattern="symbolicRep")
  
  listPredLabel=list.files(file.path(p,dl,"Output non-overlap","Cluster labels Ei weighted"),pattern="ProposedLabels")
  
  if (dl=="Ausgrid Data")
    i=7   ##48hours
  if(dl=="Stock market")
    i=4   ##25 days
  if (dl=="Web traffic")
    i=5   ##42 days
  if (dl=="London Data")
    i=1   ##12 hrs
  
    temporalPatt=read.csv(file.path(p,dl,"Output non-overlap","Temporal patterns Ei weighted",listTemporal[i]),header=T)
    window=strsplit(basename(file.path(p,dl,"Output non-overlap","Temporal patterns Ei weighted",listTemporal[i])),"win")[[1]][2]
    window=as.numeric(gsub(window,pattern=".csv",replacement=""))
    
    listSymbols1=(sapply(listSymbols,function(x) 
    {
      x=strsplit(x,"win")[[1]][2]
      gsub(x,pattern=".csv",replacement="")
    }))
    
    
    symbolicRep=read.csv(file.path(p,dl,"Output non-overlap","Cluster labels Ei weighted",listSymbols[which(listSymbols1==window)]),header=T)
    
    listPredLabel1=(sapply(listPredLabel,function(x) 
    {
      x=strsplit(x,"win")[[1]][2]
      gsub(x,pattern=".csv",replacement="")
    }))
    
    predLabel=read.csv(file.path(p,dl,"Output non-overlap","Cluster labels Ei weighted",listPredLabel[which(listPredLabel1==window)]),header=T)
    
    
    time_of_occurence=NULL
    allCluster_plots=NULL
    for (l in 6:7)
    {
      
      #extract the symbolic representation and the raw data of the required cluster
      temp.sym=data.frame(symbolicRep[which(predLabel==temporalPatt$cluster[l]),],row.names = NULL)
      
      
      
      temp.dta=data.frame(dta[which(predLabel==temporalPatt$cluster[l]),],row.names = NULL)
     
      pat1=temporalPatt[l,4:ncol(temporalPatt)]
      pat1=as.character(pat1[!is.na(pat1)])
      pat1=unique(pat1)
      
      if (dl=="Stock market" & l==1)
      {
        pat1=temporalPatt[l,4:ncol(temporalPatt)]
        pat1=as.character(pat1[!is.na(pat1)])
        pat1=pat1[1:2]
      }
        
      
      if (dl=="London Data" & l==1)
        pat1=pat1[2:11]#This is done in case of london data because the first and the last two TP is not reported in Fig 5.
      
      
      if (dl=="Ausgrid Data" & l==1)
        pat1=pat1[1:8]#This is done in case of london data because the first and the last two TP is not reported in Fig 5.
      
      pList<-NULL
      z<-0
      plotList=NULL
      
      z<-1
      plotList<-NULL
      
      # if (l==5 | l==6 & dl=="Ausgrid Data")
      # {
      #   temp.sym=temp.sym[,1:180]
      # }
      time.llcs=apply(temp.sym,1,function(f)
      {
        print (z)
        pat2=as.character(unlist(f))
        llcs=LCS(pat1,pat2)
        segmentData=nonOverlap(temp.dta[z,],window)
        plotList[[z]]<<-plot.patt2(segmentData,llcs$vb,window)
        z<<-z+1
        
      })
      
      allCluster_plots[[l]]=plotList
      # 
      # if (dl=="Ausgrid Data")
      # {
      #   grid.arrange(grobs=plotList[c(6,8,9)],ncol=1) ##cluster 1
      # }
      # 
      # 
      # if (dl=="Stock market")
      # {
      #   grid.arrange(grobs=plotList[c(3,11,17)],ncol=1) ##cluster 1
      #   grid.arrange(grobs=plotList[c(21,22,41,44)],ncol=1) ##cluster 2
      #   grid.arrange(grobs=plotList[c(31,34,52,82)],ncol=1) ##cluster 3
      #   grid.arrange(grobs=plotList[c(8,14,32)],ncol=1) ##cluster 4
      #   grid.arrange(grobs=plotList[c(5,8,11)],ncol=1) ##cluster 4
      #   
      # }
      # 
      # if (dl=="Web traffic")
      # {
      #   grid.arrange(grobs=plotList[c(7,14,16)],ncol=1) ##cluster 1
      #   grid.arrange(grobs=plotList[c(43,49,52)],ncol=1)  ##cluster 2
      #   grid.arrange(grobs=plotList[c(1,2,3)],ncol=1) ##cluster 3
      #   grid.arrange(grobs=plotList[c(38,40,42)],ncol=1) ##cluster 4
      #   grid.arrange(grobs=plotList[c(31,35,38)],ncol=1) ##cluster 6
      # }
      # 
      # if (dl=="London Data")
      # {
      #    ##cluster 1 
            # x=allCluster_plots[[1]] 
            # x[[1]]
      # }

        
      
      allCluster_plots
    }
  return (allCluster_plots)
}


p="/home/kakuli/Graph Representation of TS"
dataList=c("London Data","Ausgrid Data","Stock market","Web traffic")[4]
fileToRead=c("londonData_forGraph.csv","AusgridData_forGraph.csv","stockMarketData_forGraph.csv","webTrafficData_forGraph.csv")[4]

for (j in 1:1)#length(dataList))
{
  out=extract_tempPatt(dataList[j],fileToRead[j] )
  
  
  
}

###The date sequence for the ausgrid dataset is not common for all samples (thestart time varies). Check the dataset fro details.
date_seq=(seq(from=as.Date("2013-12-10"),to=as.Date("2017-12-10"),by=1))  ##Stock

windowLength=48
date_seq1=nonOverlap((date_seq),windowLength)



