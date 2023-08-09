##As the name signifies, this file is to obtain graph for Ei-UWG, as given in the paper.
### This file has-- normalization and log transformation, segmentation, symbolic representation, graph representation, path extraction using DFS
### feature selection and clustering (Kmeans). 



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
library(clusteval)

getClusters=function(distMat,k)
{
  cri=NULL
  totalWSS=9999
  
  labels=rep(1,nrow(distMat))
  newLabels=list()
  newLabels[[1]]=rep(1,nrow(distMat))
  
  #compute the change in rand index
  for (i in 2:ncol(distMat))
  {
    print (i)
    cl<-kmeans(distMat[,1:i],k,nstart=25,iter.max = 500)
    
    newLabels[[i]]=cl$cluster
    temp=1-cluster_similarity(unlist(newLabels[[i]]), unlist(newLabels[[i-1]]),similarity="rand")
    cri=rbind(cri,temp)
  }
  
  return (list(which.min(cri),newLabels[[which.min(cri)]]))
}




featureSpace=function(listAllPaths,symbolicRep)
{
  
  library(qualV)
  get.dist=data.frame(apply(symbolicRep,1,function(f)
  {
    apply(listAllPaths,1,function(f1)
    {
      
      ( LCS(as.character(f),as.character(f1[!is.na(f1)]))$QSI)
    })
  }),row.names = NULL)
  
  return (get.dist)
}

redundantRemoval=function(pathList)
{
  library(qualV)
  library(rlist)
  
  cleanlist=list()
  pathList= apply(pathList,1,function(f)
  {
    unlist(f[!is.na(f)])
  })
  for (x in 1:length(pathList))
  {
    print (x)
    new=pathList[[x]]
    partialMatch=lapply(pathList, function(f) (LCS(new,f)$QSI))
    
    if (!(any(partialMatch>0.8 & partialMatch<1)))
    {
      cleanlist[[x]]=pathList[[x]]
    }
    
  }
  cleanlist=cleanlist[sapply(cleanlist,function(f) length(f)!=0)]
  cleanlist=do.call(rbind.fill,lapply(cleanlist,function(f) as.data.frame(t(f))))
  
  return (cleanlist)
}

findDfsPath=function(graph,source,adj.graph,listAllPaths,non_centerNodes,
                     symbolicRep)
{
  library (data.table)
  graph.elist=data.frame(get.edgelist(graph),stringsAsFactors = F)
  E(graph)$name=seq(1,length(E(graph)))
  
  source=V(graph)$name[source]
  for (p in 1:length(non_centerNodes))
  {
    print (c(p,nrow(listAllPaths)))
    sink=V(graph)$name[non_centerNodes[p]]
    temp.graph=graph
    pathTo_sink=NULL
    while (degree(temp.graph,mode='in')[sink])    #while there are no incoming nodes to sink
    {
      
      d<-dfs(temp.graph,root=source,neimode="out",order = TRUE, order.out = TRUE)
      dfs.in=data.frame(as.matrix(d$order),row.names = NULL)
      dfs.in=apply(dfs.in,1,function(f) V(temp.graph)$name[f])
      
      dfs.out=data.frame(as.matrix(d$order.out),row.names = NULL)
      dfs.out=apply(dfs.out,1,function(f) V(temp.graph)$name[f])
      commonNodes=intersect(dfs.in[1:(which(dfs.in==sink)-1)],dfs.out[1:(which(dfs.out==sink)-1)])
      
      temp=dfs.in[1:which(dfs.in==sink)]
      temp=temp[!(temp %in% commonNodes)]
      
      if ( length(temp)==1)
        break
      
      
      else if ( which(dfs.out==source)>which(dfs.out==sink))
      {
        paths=data.frame(t(temp))
        colnames(paths)=seq(1,ncol(paths))
        pathTo_sink=rbind.fill(pathTo_sink,paths)
        listAllPaths=rbind.fill(listAllPaths,paths)
        
      }
      
      delPath=data.frame(cbind("X1"=temp,"X2"=shift(temp, n=1, fill=NA, type="lead")),stringsAsFactors = F)
      delPath=na.omit(delPath)
      
      delPath=(as.numeric(apply(delPath,1,function(f) (get.edge.ids(temp.graph,f)))))
      
      temp.graph=delete_edges(temp.graph,delPath)
      
      E(temp.graph)$name=seq(1,length(E(temp.graph)))
    }
    
    
  }
  
  return (listAllPaths)
  
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
  for (x in seq(1,length(timeseries) - swSize+1,((swSize/2)+1)))
  {
    windowCol=timeseries[x:(x+swSize-1)]
    SW <- cbind(SW, windowCol, deparse.level = 0)
  }
  SW=t(SW)
  row.names(SW)<-NULL
  return (SW)
}

extractPaths=function(mainGraph,symbolicRep)
{
  eigenVals=eigen_centrality(mainGraph, directed = T, scale = TRUE,
                   weights = NULL, options = arpack_defaults)
  centerNode=names(which.max(eigenVals$vector))
  
  adjMat=data.frame(as.matrix(as_adjacency_matrix(mainGraph)),stringsAsFactors = F)
  
  non_centerNodes=V(mainGraph)[! V(mainGraph) %in% V(mainGraph)[centerNode]]
  
  listAllPaths=NULL
  
  start=Sys.time()
  
  source=V(mainGraph)[centerNode]
  
  allPaths=findDfsPath(graph=mainGraph,source,adj.graph=adjMat,listAllPaths,non_centerNodes,
                       symbolicRep )
  
  colnames(allPaths)=seq(1,ncol(allPaths))
  
  listAllPaths=rbind.fill(listAllPaths,allPaths)
  
  print (Sys.time()-start)
  
  return (listAllPaths)
  
}



graphConstruct=function(df,dateVals,symbolicRep,houseNames)
{
    for (i in 1:nrow(symbolicRep))
    {
      
      temp=t(symbolicRep[i,])
      #store all the edges in the graph
      edge.list=data.frame(rbind(edge.list,cbind.data.frame(x=edgeList(temp))),row.names = NULL)
      temp=NULL
      
    }
    colnames(edge.list)=c("from","to")
    ##Remove duplicate edges
    edge.list=edge.list[!duplicated(edge.list),]
    
    mainGraph <- graph.data.frame(edge.list, directed = T)
    l <- layout_with_kk(mainGraph)
    
    # plot(mainGraph, layout=l, vertex.size=3,
    #      vertex.color="red", edge.arrow.size=0.7,vertex.cex=0.55,vertex.label.dist=2,
    #      vertex.label.color="blue4",asp=0,vertex.label.font=12,edge.width=2,vertex.size2=12)
    return (mainGraph)
}




getTS_symbols=function(rawData,swSize)
{
  library(rowr)
  breakpointsTable=sapply(alpha,function(f) qnorm(0:f/f))
  temp=NULL
  breakpointsTable=lapply(breakpointsTable,function(f) temp<<-cbind.fill(temp,(f),fill=NA))
  breakpointsTable=breakpointsTable[[length(breakpointsTable)]]
  breakpointsTable=data.frame(breakpointsTable[,which(unlist(lapply(breakpointsTable, function(x) !all(is.na(x)))))])
  
  colnames(breakpointsTable)=c(seq(3,(ncol(breakpointsTable)+2),1))

  meanVals=NULL
  sdVals=NULL
  maxVals=NULL
  minVals=NULL
  
  for (i in 1:nrow(rawData))
  {
    dta=unlist(rawData[i,])
    
    segData=nonOverlap(dta,swSize)
    #segData is a dataframe where each row corresponds to a segment. 
    
    meanVals= rbind(meanVals,t(apply(segData,1,function(f) mean(unlist(f)))))
    
    sdVals=  rbind(sdVals,t(apply(segData,1,function(f) sd(unlist(f)))))
    
    maxVals=  rbind(maxVals,t(apply(segData,1,function(f) which.max(unlist(f)))))
    
    minVals=  rbind(minVals,t(apply(segData,1,function(f) which.min(unlist(f)))))
    
  }
   meanSym=apply(meanVals,2,function(f.m)
  {
    fromMedian=apply(breakpointsTable,2,function(f)
    {
      f=f[!is.infinite(f) & !is.na(f)]
      abs(abs(f)-median (f.m))
      
    })
    
    fromMedian=(lapply(fromMedian, min))
    
    numAlpha=as.numeric (names(breakpointsTable)[which.min(unlist(fromMedian))])
    
    a=(convert.to.SAX.symbol((f.m),alpha=numAlpha))
    
    a=letters[a]
    
  })
  
  sdSym=apply(sdVals,2,function(f.m)
  {
    fromMedian=apply(breakpointsTable,2,function(f)
    {
      f=f[!is.infinite(f) & !is.na(f)]
      abs(abs(f)-median (f.m))
      
    })
    
    fromMedian=(lapply(fromMedian, min))
    
    numAlpha=as.numeric (names(breakpointsTable)[which.min(unlist(fromMedian))])
    
    a=(convert.to.SAX.symbol((f.m),alpha=numAlpha))
    
    a=letters[a]
    
  })
  
  maxSym=apply(maxVals,2,function(f.m)
  {
    fromMedian=apply(breakpointsTable,2,function(f)
    {
      f=f[!is.infinite(f) & !is.na(f)]
      abs(abs(f)-median (f.m))
      
    })
    
    fromMedian=(lapply(fromMedian, min))
    
    numAlpha=as.numeric (names(breakpointsTable)[which.min(unlist(fromMedian))])
    
    a=(convert.to.SAX.symbol((f.m),alpha=numAlpha))
    
    a=letters[a]
    
  })
  
  minSym=apply(minVals,2,function(f.m)
  {
    fromMedian=apply(breakpointsTable,2,function(f)
    {
      f=f[!is.infinite(f) & !is.na(f)]
      abs(abs(f)-median (f.m))
      
    })
    
    fromMedian=(lapply(fromMedian, min))
    
    numAlpha=as.numeric (names(breakpointsTable)[which.min(unlist(fromMedian))])
    
    a=(convert.to.SAX.symbol((f.m),alpha=numAlpha))
    
    a=letters[a]
  })
  allSym= matrix( paste(meanSym,sdSym, maxSym,minSym, sep=""),nrow=nrow(meanSym))
  
  return (as.data.frame(allSym))
}



mainCall=function(fileToRead,dta.processed)
{
  trainData=read.csv(file.path(getwd(),"data",fileToRead),header=T)
  houseNames=as.character(trainData$House)
  
  dataNm=strsplit(dl," ")[[1]][1]
  print (dataNm)
  
  dta=trainData[,-c(1,2,3)]
  logT_data = paste0(strsplit(basename(file.path(getwd(),"data", fileToRead)),"_")[[1]][1], "_","logTransformed.csv")
  dta.processed=read.csv(file.path(getwd(),"data",logT_data),header=T)
  
  ##Run the function in debug mode to find the range of alphabets.
  #numAlphabets(rawData)
  if (dl=="London Data")
    alpha=seq(10,20)
  else
    alpha=seq(12,20)
  clusteringResults = NULL
  
  for (w in 1:length(window))
  {
    print (window[w])
    
    symbolicRep=getTS_symbols(dta.processed,window[w],alpha)
    mainGraph=graphConstruct(dta.processed,symbolicRep,houseNames)
    nClust_ev=numClust[which(numClust$W_UW=="EV"   & (numClust$window==window[w]) & (numClust$Method=="Proposed")) ,]$numClust
    listAllPaths=extractPaths(mainGraph,symbolicRep)

    lengthPaths=apply(listAllPaths,1,function(f)
      {
      length(f[!is.na(f)])/nrow(symbolicRep)
    })
    percentMatch=lengthPaths/nrow(symbolicRep)

    chosenPaths=listAllPaths[which(lengthPaths>0.07),]
    chosenPaths=redundantRemoval(chosenPaths)
    featureDist=featureSpace(chosenPaths,symbolicRep)
    
    numFeatures=getClusters(featureDist,nClust_ev)
    print (numFeatures[[2]])
    clusteringResults[w] = numFeatures[[2]]
  }
  
  return (clusteringResults )
  
}

dataList=c("London Data","Ausgrid Data","Web traffic")
fileToRead=c("londonData_forGraph.csv","AusgridData_forGraph.csv","trimmedWebTraffic_labelled.csv")[1]
window=seq(6,48,6)
for (dl in 1:length(dataList))
{
  out=mainCall(dataList[dl],fileToRead[dl])
  
}
