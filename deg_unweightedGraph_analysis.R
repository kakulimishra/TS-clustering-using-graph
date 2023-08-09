###This file was witten to obtain the rare events for eac dataset and also to plot the graph

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


##Thsi file has the component level analysis, for discovery of rare events.


plotRE=function(mainGraph)
{
  ###The function is designed to obtain plots for the rare events. 
  
  subG=subgraph(mainGraph,c("ghmm","dhoo"))  ##Change required here, in the list of nodes for the subgraph.
  
  ###Design the layout, to obtain a straight line
  lo <- cbind(seq(-1,1,length.out = gorder(subG)), 0)
  
  designedPlot <- as.numeric(ego(subG, gorder(subG), "ghmm")[[1]])  ###Change the node inside ""
  
  E(subG)$color="black"
  
  
  plot(subG, layout=lo[order(designedPlot), ], vertex.size=8,
       vertex.color="red", edge.arrow.size=0.8,vertex.cex=0.9,vertex.label.dist=2,
       vertex.label.color="blue4",asp=0,vertex.label.font=40,edge.width=2)
  
  
  
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
      
      # compDetails=rbind(compDetails,cbind("component"=i,"anomaly"=ano
      #                           ,"compSize"=length(unlist(mainGraph.comp[[i]]) )))
      
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

  
  ####*****************************CLIQUE DETAILS*****************************************************************
  
  # getClique=largest_cliques(mainGraph)
  # 
  # 
  # cliquePath= lapply(getClique,function(f)  V(mainGraph)[as.numeric(as.character(f))]$name)
  # 
  # cliquePath=do.call(rbind,cliquePath)
  # 
  # compDetails=cbind(compDetails,"lengthClique"=length(getClique))
  
  ###**************************************************************************************************************
  
 # return (list(compDetails, rareEvents,cliquePath))
  
  
  ##***********************************************************888
  ###Find the component score from the below code andfor each compoenent.
  subG=subgraph(mainGraph, V(mainGraph)[mainGraph.comp[[1]]])
  ##See the definition of comp score in the paper.
  
  numEdges.mainGraph=length(E(mainGraph))
  numEdges.subg=length(E(subG))
  
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
    
    dateVals=trainData$dates
    houseNames=as.character(trainData$House)
    
    if (dataList[dl]=="London Data")
      date_seq=(seq(from=as.POSIXct("2012-03-01"),to=as.POSIXct("2012-08-31"),by=(60*30)))  ##London
    
    if (dataList[dl]=="Ausgrid Data")
      date_seq=dateVals  ##Ausgrid
    
    
    if (dataList[dl]=="Stock market")
    {
      date_seq=(seq(from=as.Date("2013-12-10"),to=as.Date("2017-12-10"),by=1))  ##Stock
      ##No stock in weekends, so the weekeds should be removed
      date_seq <- date_seq["Date"=!weekdays(date_seq) %in% c('Saturday','Sunday')]
      #date_seq <-  data.frame("Date"=date_seq)
      window=c(5,10,15,20,25,30,35,40,45,50)[2]
      houseNames=as.character(trainData$trueLabel)
    }
    
    if (dataList[dl]=="Web traffic")
    {
      date_seq=(seq(from=as.Date("2015-07-1"),to=as.Date("2017-09-10"),by=1))  ##Stock
      dta=trainData[,-c(1,2)]
      window=c(7,14,21,28,35,42,49)[1]
      houseNames=as.character(trainData$trueLabel)
      
    }

    else
      dta=trainData[,-c(1,2,3)]
    
    for (w in 1:length(window))
    {
      print (c(window[w],dataList[dl]))
      
      symbolicRep=read.csv(file.path(p,dataList[dl],"Output non-overlap","Cluster labels unweighted",paste0("symbolicRep","_win",window[w],".csv")),header=T)
      
      out=graphConstruct(dta.processed[dl],date_seq,symbolicRep,houseNames,window[w],dataList[dl])

      g.components=rbind(g.components,cbind.data.frame("dataName"=dataList[dl],"window"=window[w],out[[1]],"numSubsequences"=length(unlist(symbolicRep))))
      
      if (!is.null(out[[2]]))
        g.rareEvents=rbind(g.rareEvents,cbind("dataName"=dataList[dl],"window"=window[w],out[[2]]))
      

    }
  }
  
  return (list(g.components,g.rareEvents,mainGraph))
}


p="/home/kakuli/Graph Representation of TS"
dataList=c("London Data","Ausgrid Data","Stock market","Web traffic")[4]
fileToRead=c("londonData_forGraph.csv","AusgridData_forGraph.csv","stockMarketData_forGraph.csv",
             "webTrafficData_forGraph.csv")[4]
dta.processed=c("londonData_logTransformed.csv","ausgridData_logTransformed.csv",
                "stockMarketData_logTransformed.csv",
                "webTraffic_logTransformed.csv")[4]


#window=c(12,18, 24,30,36,42,48,54,60)[3]
window=c(7,14,21,28,35,42,49)

out=mainCall(dataList,fileToRead , dta.processed)
#write.csv(data.frame(out[[1]]),"/home/kakuli/Graph Representation of TS/output_graphComponents.csv",row.names=F)
#write.csv(data.frame(out[[2]]),"/home/kakuli/Graph Representation of TS/output_rareEvents.csv",row.names=F)




