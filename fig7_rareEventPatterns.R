##This file is actually copied from the file "graphRep5_analysis.R" only with an added function" anomalyPlots().
###Purpose of the function "anomalyPlots" is to obtain the plots for the anomalous subsequences. This was not done in the 
##orignal file just to avoid it messy.

library(TSclust)
library(igraph)
library(dplyr)
library(reshape2)
library(plyr)
library(cluster)
library(dplyr)
library(clusterSim)
library(factoextra)
library(tidyverse)
library(cluster)
library(RevEcoR)




webTraffic_twoBuildings=function(rareEvents,dl,symbolicRep)
{
  ###Web Traffic has anomaly on 3 buildings. I has plotted on 2 of them but later it was 
  ##removed from the paper because of zero values at the beginning which could be the missing data.
  
  ##This function will give the anomaly plot for two buildings. Check, Writeup 7 for the plot.
  
  
  trainData=read.csv(file.path(p,dataList[dl],fileToRead[dl]),header=T)
  
  
  
  date_seq=(seq(from=as.Date("2015-07-1"),to=as.Date("2017-09-20"),by=1))  ##web traffic
  windowLength=7
  houseNames=as.character(trainData$trueLabel)
  anomaly_pages=which(as.character(trainData$Page)%in% unique(as.character(rareEvents$houseName)))
  anomaly_symbols=c("ghmm",'ghmm',"mhmo")
  allPages=NULL
  for (p in 1:length(anomaly_pages))
  {
    
    houseData=trainData[anomaly_pages[p],]
    
    segmentData=nonOverlap(unlist(houseData[,3:ncol(houseData)]),windowLength)  ##Put the value of "windowLength"
    
    segmentDates=nonOverlap(unlist(date_seq),windowLength)
    
    #data.with.date=data.frame(cbind(segmentDates[1:nrow(segmentData),1],segmentData))
    
    #data.with.date=data.with.date[which(data.with.date$X1 %in% rareEvents$firstOccurenceDate.x[1]):which(data.with.date$X1 %in% rareEvents$firstOccurenceDate.x[nrow(rareEvents)]),]
    
    data.with.symbols=data.frame(cbind(as.character(unlist(symbolicRep[anomaly_pages[p],])[1:nrow(segmentData)]),segmentData),stringsAsFactors = F)
    
    data.with.symbols$id=rep(0,nrow(data.with.symbols))
    
    
    #data.with.symbols$id[which(data.with.symbols[,1] %in% c(as.character(rareEvents$anomaly[4]), as.character(rareEvents$vertexName.x[4])))]=1
    data.with.symbols$id[which(data.with.symbols[,1] ==anomaly_symbols[p])]=p
    
    
    #ggplot(x1)+geom_line(data=x1,aes(x=id,y=val),size=4)+theme_void()
    plotData<-NULL
    apply(data.with.symbols[,2:ncol(data.with.symbols)],1,function(f)
    {
      ti=(unlist(f[1:(length(f)-1)]))
      plotData<<-rbind(plotData,cbind(ti,as.numeric(rep(f[length(f)], length(ti)))))
    })
    
    plotData=data.frame(plotData)
    plotData$date=date_seq[1:nrow(plotData)]
    plotData$page=rep(paste0("p",p),nrow(plotData))
    allPages=rbind(allPages,plotData)
    
  }
  
  
  allPages=data.frame(allPages)
  colnames(allPages)=c("val","id","date","page")
  
  allPages$val=as.numeric(as.character(allPages$val))
  allPages$id=as.numeric(as.character(allPages$id))
  
  line_color=c("pink","red")
  my_colors=c("black","blue","sienna1")
  
  my_size=c(0.3,rep(1,10))
  
  allPages$date=as.Date(allPages$date)
  
  allPages=allPages[-(which(allPages$page=='p2')),]
  
  myLabels=c("Mamamoo","Emmanuel_Macron")
  
  data_geomPoints=(seq(min(allPages$val), max(allPages$val),20))
  data_geomPoints=allPages[c(data_geomPoints),]
  data_geomPoints=na.omit(data_geomPoints)
  
  ggplot(data=allPages,aes(x=date,y=val,color=factor(id),size=factor(id),show.legend=F))+
    geom_path(group=as.character(allPages$page),show.legend = F)+
    geom_point(data=data_geomPoints,aes(x=date,y=val,shape=factor(page)),size=2,position=position_dodge(width = 1))+
    theme_bw()+
    scale_color_manual(values = my_colors)+ theme(axis.text.x = element_text( vjust = 0.5))+
    scale_size_manual(values=my_size)+scale_shape_manual(values=c(17,19),labels=myLabels)+
    theme(plot.title = element_text(hjust = 0.5),legend.title = element_blank(),
          axis.title.y = element_blank(),axis.ticks.y = element_blank(),legend.position = "bottom",
          legend.text = element_text(size=15,color='black'),
          axis.text = element_text(colour='black',size=13), 
          panel.grid.major.y =element_blank(),panel.grid.major.x = element_blank() ,
          axis.title = element_blank(),plot.caption = element_text(hjust=0.5, size=rel(1.2)))+
    
    scale_x_date(breaks=seq(as.Date(plotData$date[1]), as.Date(plotData$date[nrow(plotData)]),by=50), date_labels= "%b\n %d")+
    scale_y_sqrt(breaks=seq(min(allPages$val),max(allPages$val),50000))
}


anomalyPlots=function(rareEvents,dl,symbolicRep)
{
  ##This function is added just to get the plots for the subsequenceswhich are anomalous, as detected in the dataframe 'rareEvents'.
  
  trainData=read.csv(file.path(p,dataList[dl],fileToRead[dl]),header=T)

  
  if (dataList[dl]=="London Data")
  {
    #trainData=trainData[-34,] ##Due to some errors, the 34th house is remoed. But this is a big loophole in the work. Need to work on it.
    date_seq=(seq(from=as.POSIXct("2012-03-01"),to=as.POSIXct("2012-08-31"),by=(60*30)))  ##London
    
    i=35
    my_colors=c("black","blue","darkgreen","red","goldenrod4","gold","darkorchid1","palegreen1","sienna1")
    
    
    plotData=t(trainData[i,4:ncol(trainData)])
 
    
    plotData=data.frame(cbind(as.character(date_seq[1:nrow(plotData)]),plotData))
    plotData$id=rep(0,nrow(plotData))
    colnames(plotData)=c("date","val","id")
    plotData$id[7000:7012]=1
    plotData$id[7013:7024]=2
    plotData$id[7025:7036]=3
    plotData$id[7037:7048]=4
    plotData$id[7049:7060]=5
    plotData$id[7059:7072]=6
    plotData$id[7073:7084]=7
    plotData$id[7085:7096]=8
    
    
    plotData$val=as.numeric(as.character(plotData$val))
    plotData$id=as.numeric(as.character(plotData$id))
    #my_colors=c("blue","red")
    
    plotData$date=as.POSIXct(plotData$date)
    sub.plotData=plotData[6300:7500,]
    
    my_size=c(0.3,rep(1.2,8))
    
    rects <- data.frame(xmin = c(sub.plotData$date[1], sub.plotData$date[800]), 
                        xmax = c(sub.plotData$date[699],sub.plotData$date[nrow(sub.plotData)]),
                        ymin = c(0,0),  
                        ymax = c(1.2,1.2),
                        fill = c("green", "yellow"))
    
    
    ggplot(data=sub.plotData,aes(x=date,y=val,color=factor(id),size=factor(id)))+geom_path(group=1)+ theme_bw()+
      scale_color_manual(values = my_colors)+ theme(axis.text.x = element_text(vjust = 0.5))+
      scale_size_manual(values=my_size)+
      theme(plot.title = element_text(hjust = 0.5),legend.title = element_blank(),
            legend.position = 'none',axis.title.y = element_blank(),axis.ticks.y = element_blank(),
            axis.text = element_text(colour='black',size=12), 
            panel.grid.major.y =element_blank(),panel.grid.major.x = element_blank() ,
            axis.title = element_blank(),plot.caption = element_text(hjust=0.5, size=rel(1.2)))+
      
      scale_x_datetime(breaks=seq(as.POSIXct(sub.plotData$date[1]), as.POSIXct(sub.plotData$date[nrow(sub.plotData)]),by=60*3000), 
                       date_labels= "%b-%d\n%H:%M")+ 
      geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill), 
                # Control the shading opacity here.
                inherit.aes = FALSE, alpha = 0.15)
    
  }

    
  
  if (dataList[dl]=="Ausgrid Data")
  {
    i=which(trainData$House==unique(rareEvents$houseName))
    
    date_seq=trainData$dates[i]
    #genDate=(seq(from=as.POSIXct(date_seq[189]),to=as.POSIXct(as.Date(date_seq[189])+182),by=(60*30)))  ##Ausgrid
    genDate=(seq(from=as.POSIXct(date_seq[1]),to=as.POSIXct(as.Date(date_seq[1])+182),by=(60*30)))  ##Ausgrid
    date_seq=genDate
    
    houseData=trainData[i[1],]
    windowLength=24
    segmentData=nonOverlap(unlist(houseData[,4:ncol(houseData)]),windowLength)  ##Put the value of "windowLength"
    
    segmentDates=nonOverlap(unlist(date_seq),windowLength)
    
    data.with.date=data.frame(cbind(segmentDates[1:nrow(segmentData),1],segmentData))
    
    #data.with.date=data.with.date[which(data.with.date$X1 %in% rareEvents$firstOccurenceDate.x[1]):which(data.with.date$X1 %in% rareEvents$firstOccurenceDate.x[nrow(rareEvents)]),]
    
    data.with.symbols=data.frame(cbind(as.character(unlist(symbolicRep[i[1],])[1:nrow(segmentData)]),segmentData),stringsAsFactors = F)
    
    data.with.symbols$id=rep(0,nrow(data.with.symbols))
    
    
   data.with.symbols$id[which(data.with.symbols[,1] %in% as.character(rareEvents$anomaly))]=1

    #ggplot(x1)+geom_line(data=x1,aes(x=id,y=val),size=4)+theme_void()
    plotData<-NULL
    apply(data.with.symbols[,2:ncol(data.with.symbols)],1,function(f)
      {
      ti=(unlist(f[1:(length(f)-1)]))
      plotData<<-rbind(plotData,cbind(ti,as.numeric(rep(f[length(f)], length(ti)))))
    })
    
    plotData=data.frame(cbind(as.character(date_seq[1:nrow(plotData)]),plotData))
    colnames(plotData)=c("date","val","id")
    
    plotData$val=as.numeric(as.character(plotData$val))
    plotData$id=as.numeric(as.character(plotData$id))
    my_colors=c("black","blue")
    my_size=c(0.3,rep(1.2,8))
    samplePlot=plotData[7500:8500,]
    samplePlot$date=as.POSIXct(samplePlot$date)
    
    rects <- data.frame(xmin = c(samplePlot$date[1], samplePlot$date[230]), 
                        xmax = c(samplePlot$date[205],samplePlot$date[nrow(samplePlot)]),
                        ymin = c(0,0),  
                        ymax = c(max(samplePlot$val),max(samplePlot$val)),
                        fill = c("green", "yellow"))
    
    ggplot(data=samplePlot,aes(x=date,y=val,color=factor(id),size=factor(id)))+geom_path(group=1)+ theme_bw()+
      scale_color_manual(values = my_colors)+ theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
      scale_size_manual(values=my_size)+
      theme(plot.title = element_text(hjust = 0.5),legend.title = element_blank(),
            legend.position = 'none',axis.title.y = element_blank(),axis.ticks.y = element_blank(),
            axis.text = element_text(colour='black',size=13), 
            panel.grid.major.y =element_blank(),panel.grid.major.x = element_blank() ,
            axis.title = element_blank(),plot.caption = element_text(hjust=0.5, size=rel(1.2)))+ 
      geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill), 
                # Control the shading opacity here.
                inherit.aes = FALSE, alpha = 0.15)
      
      
      #scale_x_continuous(breaks=seq(from=(samplePlot$date[1]), to=(samplePlot$date[nrow(samplePlot)]),by=60*120))
  }
  
  
  if (dataList[dl]=="Stock market")
  {
    date_seq=(seq(from=as.Date("2013-12-10"),to=as.Date("2017-12-10"),by=1))  ##Stock
    ##No stock in weekends, so the weekeds should be removed
    date_seq <- date_seq["Date"=!weekdays(date_seq) %in% c('Saturday','Sunday')]
    #date_seq <-  data.frame("Date"=date_seq)
    window=c(5,10,15,20,25,30,35,40,45,50)[2]
    houseNames=as.character(trainData$stock)
    i=which(as.character(trainData$stock)==as.character(unique(rareEvents$houseName)))
    
    
    houseData=trainData[i,]
    windowLength=10
    segmentData=nonOverlap(unlist(houseData[,4:ncol(houseData)]),windowLength)  ##Put the value of "windowLength"
    
    segmentDates=nonOverlap(unlist(date_seq),windowLength)
    
    data.with.date=data.frame(cbind(segmentDates[1:nrow(segmentData),1],segmentData))
    
    #data.with.date=data.with.date[which(data.with.date$X1 %in% rareEvents$firstOccurenceDate.x[1]):which(data.with.date$X1 %in% rareEvents$firstOccurenceDate.x[nrow(rareEvents)]),]
    
    data.with.symbols=data.frame(cbind(as.character(unlist(symbolicRep[i[1],])[1:nrow(segmentData)]),segmentData),stringsAsFactors = F)
    
    data.with.symbols$id=rep(0,nrow(data.with.symbols))
    
    unique.anomaly=unique(as.character(rareEvents$anomaly))
    
    for (a in 1:length(unique.anomaly))
    {
      data.with.symbols$id[which(data.with.symbols[,1] %in% unique.anomaly[a])]=a
    }
    
    #ggplot(x1)+geom_line(data=x1,aes(x=id,y=val),size=4)+theme_void()
    plotData<-NULL
    apply(data.with.symbols[,2:ncol(data.with.symbols)],1,function(f)
    {
      ti=(unlist(f[1:(length(f)-1)]))
      plotData<<-rbind(plotData,cbind(ti,as.numeric(rep(f[length(f)], length(ti)))))
    })
    
    plotData=data.frame(cbind(as.character(date_seq[1:nrow(plotData)]),plotData))
    colnames(plotData)=c("date","val","id")
    
    plotData$val=as.numeric(as.character(plotData$val))
    plotData$id=as.numeric(as.character(plotData$id))
    my_colors=c("black","blue","darkgreen","red","goldenrod4","gold","darkorchid1","palegreen1","sienna1","maroon")
    my_size=c(0.3,rep(1,10))
    plotData$date=as.Date(plotData$date)
    
    ggplot(data=plotData,aes(x=date,y=val,color=factor(id),size=factor(id)))+geom_path(group=1)+ theme_bw()+
      scale_color_manual(values = my_colors)+ theme(axis.text.x = element_text( vjust = 0.5))+
      scale_size_manual(values=my_size)+
      theme(plot.title = element_text(hjust = 0.5),legend.title = element_blank(),
            legend.position = 'none',axis.title.y = element_blank(),axis.ticks.y = element_blank(),
            axis.text = element_text(colour='black',size=12), 
            panel.grid.major.y =element_blank(),panel.grid.major.x = element_blank() ,
            axis.title = element_blank(),plot.caption = element_text(hjust=0.5, size=rel(1.2)))+
      
      scale_x_date(breaks=seq(as.Date(plotData$date[1]), as.Date(plotData$date[nrow(plotData)]),by=70), date_labels= "%b\n %d")+
      scale_y_continuous(breaks=seq(min(plotData$val),max(plotData$val),5000),trans=log2_trans())
      
  }
  
  if (dataList[dl]=="Web traffic")
  {
    
    date_seq=(seq(from=as.Date("2015-07-1"),to=as.Date("2017-09-20"),by=1))  ##web traffic
    windowLength=7
    houseNames=as.character(trainData$trueLabel)
    anomaly_pages=which(as.character(trainData$Page)%in% as.character(rareEvents$houseName[[4]]))
    anomaly_symbols=c("mhmo")
    

    houseData=trainData[anomaly_pages,]
    
    segmentData=nonOverlap(unlist(houseData[,3:ncol(houseData)]),windowLength)  ##Put the value of "windowLength"
    
    segmentDates=nonOverlap(unlist(date_seq),windowLength)
    
    data.with.symbols=data.frame(cbind(as.character(unlist(symbolicRep[anomaly_pages,])[1:nrow(segmentData)]),segmentData),stringsAsFactors = F)
    
    data.with.symbols$id=rep(0,nrow(data.with.symbols))
    
    
    data.with.symbols$id[which(data.with.symbols[,1] ==anomaly_symbols)]=1
    
    
    #ggplot(x1)+geom_line(data=x1,aes(x=id,y=val),size=4)+theme_void()
    plotData<-NULL
    apply(data.with.symbols[,2:ncol(data.with.symbols)],1,function(f)
    {
      ti=(unlist(f[1:(length(f)-1)]))
      plotData<<-rbind(plotData,cbind(ti,as.numeric(rep(f[length(f)], length(ti)))))
    })
    
    plotData=data.frame(plotData)
    plotData$date=date_seq[1:nrow(plotData)]
    
    plotData$val=as.numeric(as.character(plotData$val))
    colnames(plotData)=c("val","id","date")

    my_colors=c("black","sienna1")
    
    my_size=c(0.3,rep(1,10))
    

    
    ggplot(data=plotData,aes(x=date,y=val,color=factor(id),size=factor(id)))+
      geom_path(group=1)+
        theme_bw()+
      scale_color_manual(values = my_colors)+ 
      scale_size_manual(values=my_size)+
      theme(plot.title = element_text(hjust = 0.5),legend.title = element_blank(),
            axis.title = element_blank(),axis.ticks.y = element_blank(),legend.position ="none",
            axis.text = element_text(colour='black',size=13))+
      
      scale_x_date(breaks=seq(as.Date(plotData$date[1]), as.Date(plotData$date[nrow(plotData)]),by=50), 
                   date_labels= "%b\n %d")
    
  }

}



plotRE=function(mainGraph)
{
  ###The function is designed to obtain plots for the rare events. 
  
  subG=subgraph(mainGraph,c("femk","edkl")) ##Change required here, in the list of nodes for the subgraph.
  
  ###Design the layout, to obtain a straight line
  lo <- cbind(seq(-1,1,length.out = gorder(subG)), 0)
  
  designedPlot <- as.numeric(ego(subG, gorder(subG), "femk")[[1]])  ###Change the node inside ""
  
  E(subG)$color="black"
  
  
  plot(subG, layout=lo[order(designedPlot), ], vertex.size=10,
       vertex.color="red", edge.arrow.size=0.8,vertex.cex=0,vertex.label.dist=2.5,
       vertex.label.color="blue4",asp=0,vertex.label.font=50,edge.width=2)
  
  
  
}




edgeList=function(temp)
{
  
  #Create a dataframe that has "from" and "to" nodes, i.e, vertices onnected by edges.
  graphEdges=NULL
  for (j in 1:(length(temp)-1))
  {
    graphEdges=rbind(graphEdges,cbind("from"=as.character(temp[j]),"to"=as.character(temp[j+1])))
    
  }
  return (data.frame(graphEdges,stringsAsFactors = F))
}



nonOverlap=function (timeseries, swSize) 
{
  ####  The function is used for obtaining the non-overlapping segments
  
  SW<-NULL
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


graphConstruct=function(df,date_seq,symbolicRep,houseNames,window,dataName)
{
  pathLength=NULL
  
  allEdge.info=NULL #store info for of all the edges in the graph.
  allVertex.info=NULL
  edge.list=NULL
  
  
  for (i in 1:nrow(symbolicRep))
  {
    ###****************Segment the dates based on the window length*************************************
    
    if (dataName=="Ausgrid Data")
    {
      genDate=(seq(from=as.POSIXct(date_seq[i]),to=as.POSIXct(as.Date(date_seq[i])+182),by=(60*30)))  ##Ausgrid
      date_seq=genDate
    }
    
    
    
    date_seq1=nonOverlap(date_seq,window)
    temp=t(symbolicRep[i,])
    
    temp=data.frame(cbind.data.frame(temp,date_seq1[1:length(temp),1],stringsAsFactors = F))
    
    colnames(temp)=c("segment","dateTime")
    
    temp=temp[!duplicated(temp$segment), ]  #this will capture the unique symbols for each row.
    
    pathLength=rbind(pathLength,nrow(temp)-1)
    
    #stores the unique vertices in each row and their dates of first occurence for the respective houses
    allVertex.info=data.frame(rbind(allVertex.info,cbind.data.frame(temp,"houseName"=houseNames[i])),row.names = NULL,stringsAsFactors = F)
    
    #store all the edges in the graph
    edge.list=data.frame(rbind(edge.list,cbind.data.frame(x=edgeList(temp$segment))),row.names = NULL)
    
    temp=NULL
    
  }
  
  ###*****************************************************************************************************
  colnames(allVertex.info)=c("vertexName","firstOccurenceDate","houseName")
  colnames(edge.list)=c("from","to")
  
  edge.list=edge.list[!duplicated(edge.list),]
  
  mainGraph <- graph.data.frame(edge.list, directed = T)
  
  ###The function below should be called later, after finding the 'rareEvents'
  #plotRE(mainGraph)
  
  adjMat=data.frame(as.matrix(as_adjacency_matrix(mainGraph)),stringsAsFactors = F)
  mainGraph.comp=KosarajuSCC(mainGraph)
  
  compDetails=NULL
  
  
  for (i in 1:length(mainGraph.comp))
  {
    
    if (length(unlist(mainGraph.comp[i]))==1)
    {
      ano= names(unlist(mainGraph.comp[i]))
      
      compDetails=rbind(compDetails,cbind("component"=i,"anomaly"=ano,edge.list[(E(mainGraph) [ from(ano) | (to(ano)) ]),],
                                          "numVertices"=  length(V(mainGraph))  ,"compSize"=length(unlist(mainGraph.comp[[i]]) )))

      
    }
    else
      compDetails=rbind(compDetails,cbind("component"=i,"anomaly"=NA,"from"=NA,"to"=NA,
                                          "numVertices"=length(V(mainGraph)),"compSize"=length(unlist(mainGraph.comp[[i]]) )))
    
    
  }
  
  
  ####***********************RARE EVENT DETAILS***************************************************************
  copy.compDetails=data.frame(na.omit(as.data.frame(compDetails)),row.names=NULL) ##remove components with no anomaly
  rareEvents=NULL
  
  if (nrow(copy.compDetails)>0)
  {
    for (z in 1:nrow(copy.compDetails))
    {
      anoFrom.loc=allVertex.info[which(allVertex.info$vertexName==copy.compDetails$from[z]),]
      
      anoTo.loc=allVertex.info[which(allVertex.info$vertexName==copy.compDetails$to[z] ),]
      
      re=merge(anoFrom.loc,anoTo.loc,by="houseName")
      
      rareEvents=rbind(rareEvents,cbind("component"=copy.compDetails$component[z],"anomaly"=copy.compDetails$anomaly[z],re))
      
      
    }
  }

  return (list(compDetails, rareEvents))
  
  
  
}



mainCall=function(dl,fileToRead,dta.processed)
{
  
  g.components=NULL
  g.clique=NULL
  g.rareEvents=NULL
  
  for (dl in 1:length(dataList))
  {
    
    trainData=read.csv(file.path(p,dataList[dl],fileToRead[dl]),header=T)
    
    #trainData=trainData[-34,] ##In case of London data, use this line because,  the house #34 shows the demerits of the work
    
    dateVals=trainData$dates
    houseNames=as.character(trainData$House)
    
    if (dataList[dl]=="London Data")
    {
      date_seq=(seq(from=as.POSIXct("2012-03-01"),to=as.POSIXct("2012-08-31"),by=(60*30)))  ##London
      window=12
      
    }
    
    if (dataList[dl]=="Ausgrid Data")
    {
      date_seq=dateVals  ##Ausgrid
      window=24
    }
    
    
    if (dataList[dl]=="Stock market")
    {
      date_seq=(seq(from=as.Date("2013-12-10"),to=as.Date("2017-12-10"),by=1))  ##Stock
      ##No stock in weekends, so the weekeds should be removed
      date_seq <- date_seq["Date"=!weekdays(date_seq) %in% c('Saturday','Sunday')]
      window=10
      houseNames=as.character(trainData$stock)
    }
    
    if (dataList[dl]=="Web traffic")
    {
      date_seq=(seq(from=as.Date("2015-07-1"),to=as.Date("2017-09-10"),by=1))  
      dta=trainData[,-c(1,2)]
      window=7
      houseNames=as.character(trainData$Page)
      
    }
    
    else
      dta=trainData[,-c(1,2,3)]
      
    symbolicRep=read.csv(file.path(p,dataList[dl],"Output non-overlap","Cluster labels Deg unweighted",paste0("symbolicRep","_win",window,".csv")),header=T)
    #symbolicRep=symbolicRep[-34,]  ##In case of London data, use this line because, the house #34 shows the demerits of the work
      
    out=graphConstruct(dta.processed[dl],date_seq,symbolicRep,houseNames,window,dataList[dl])
      
    g.components=rbind(g.components,cbind.data.frame("dataName"=dataList[dl],"window"=window,out[[1]],"numSubsequences"=length(unlist(symbolicRep))))
      
    if (!is.null(out[[2]]))
      g.rareEvents=rbind(g.rareEvents,cbind("dataName"=dataList[dl],"window"=window,out[[2]]))
      
    anomalyPlots(g.rareEvents,dl,symbolicRep)
      
    
  }
}


p="/home/kakuli/Graph Representation of TS"
dataList=c("London Data","Ausgrid Data","Stock market","Web traffic")[3]
fileToRead=c("londonData_forGraph.csv","AusgridData_forGraph.csv","stockMarketData_forGraph.csv",
             "webTrafficData_forGraph.csv")[3]
dta.processed=c("londonData_logTransformed.csv","ausgridData_logTransformed.csv",
                "stockMarketData_logTransformed.csv",
                "webTraffic_logTransformed.csv")[3]


out=mainCall(dataList,fileToRead , dta.processed)
#write.csv(data.frame(out[[1]]),"/home/kakuli/Graph Representation of TS/output_graphComponents.csv",row.names=F)
#write.csv(data.frame(out[[2]]),"/home/kakuli/Graph Representation of TS/output_rareEvents.csv",row.names=F)




