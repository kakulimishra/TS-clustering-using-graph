##As the name signifies, this file is to obtain graph for Ei-WG, as given in the paper.
### This file has-- normalization and log transformation, segmentation, symbolic representation, graph representation, path extraction using Dijkstras

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


shortPaths=function(adj.graph,source,graph,non_centerNode,minPath.length,maxPath.length)
  #Idea is to find the shortest path between source and sink that includes maximum number of vertices in the path.
  #Note that, one path should not be the subset of other.
{
  
  #take the source and sink nodes and another intermediate node in between them. Find the shortest path between the source and sink with
  #the intermediate node in between. Delete the detected shortest path discovered and increase the intermediate nodes between source and sink.
  #Update the adjacecy matrix of the original graph in each iteration.
  
  
  graph.elist=data.frame(get.edgelist(graph),stringsAsFactors = F)
  source=V(graph)$name[source]
  E(graph)$name=seq(1,length(E(graph)),1)
  s2=Sys.time()
  allPaths=NULL
  
  # temp.graph=graph
  # adj=adj.graph
  
  for (p in 1:length(non_centerNode))
  {
    #print (p)
    sink=V(graph)$name[non_centerNode[p]]
    temp.graph=graph
    adj=adj.graph
    
    E(temp.graph)$name=seq(1,length(E(temp.graph)),1)
    
    while(length(shortest_paths(temp.graph,from=source,to=sink,output="both")$epath[[1]])!=0)
    {
      pathList=shortest_paths(temp.graph,from=source,to=sink,output="both")
      
      if(length(unlist(pathList$vpath))>=minPath.length & length(unlist(pathList$vpath))<=maxPath.length)
      {
        paths=data.frame(t((V(temp.graph)[unlist(pathList$vpath)])$name))
        
        colnames(paths)=seq(1,ncol(paths))
        
        allPaths=rbind.fill(allPaths,paths)
      }
      
      #delete the edges detected in previous discovered path. This is done to include maximum number of nodes in the shortest path.
      temp.graph=delete_edges(temp.graph,pathList$epath[[1]])
    }
    
  }
  print (Sys.time()-s2) 
  return (data.frame(allPaths,stringsAsFactors = F))   #return the unique list of edges in case here are repetations.
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
      #cleanlist[[x]]=pathList[[which.max((partialMatch>0.85 & partialMatch<1))]]
      cleanlist[[x]]=pathList[[x]]
    }
    
  }
  cleanlist=cleanlist[sapply(cleanlist,function(f) length(f)!=0)]
  cleanlist=do.call(rbind.fill,lapply(cleanlist,function(f) as.data.frame(t(f))))
  
  return (cleanlist)
}



gapFunction=function(graphPath,nClust_evw,symbolicRep)
{
  nClust_evw=3
  library(qualV)
  pair.dist=unlist(apply(symbolicRep,1,function(f)
  {
    ( LCS(as.character(f),as.character(graphPath[!is.na(graphPath)]))$QSI)
  }))
  pair.dist=sort(pair.dist,decreasing = F)
  maxGap<-0
  dt<-0
  for (j in 1:(length(pair.dist)-1))
  {
    midVal=mean(pair.dist[j],pair.dist[j+1])
    clustA=pair.dist[which(pair.dist<midVal)]
    clustB=pair.dist[which(pair.dist>midVal)]
    r=length(clustA)/length(clustB)
    if (r<(1-(1/nClust_evw)) & r>1/nClust_evw)
    {
      
      mean_clustA=mean(clustA)
      mean_clustB=mean(clustB)
      std_clustA=sd(clustA)
      std_clustB=sd(clustB)
      gap=mean_clustB-std_clustB-(mean_clustA+std_clustA)
      if (gap>maxGap)
      {
        maxGap<-gap
        dt<-midVal
        
      }
      
    }
    
    
    
  }
  return (list(maxGap,dt))
  
}



shortPaths=function(adj.graph,source,graph,non_centerNodes)
  #Idea is to find the shortest path between source and sink that includes maximum number of vertices in the path.
  #Note that, one path should not be the subset of other.
{
  
  #take the source and sink nodes and another intermediate node in between them. Find the shortest path between the source and sink with
  #the intermediate node in between. Delete the detected shortest path discovered and increase the intermediate nodes between source and sink.
  #Update the adjacecy matrix of the original graph in each iteration.
  
  graph.elist=data.frame(get.edgelist(graph),stringsAsFactors = F)
  source=V(graph)$name[source]
  E(graph)$name=seq(1,length(E(graph)),1)
  s2=Sys.time()
  allPaths=NULL

  for (p in 1:length(non_centerNodes))
  {
    #print (p)
    sink=V(graph)$name[non_centerNodes[p]]
    temp.graph=graph
    adj=adj.graph
    
    E(temp.graph)$name=seq(1,length(E(temp.graph)),1)
    
    while(length(shortest_paths(temp.graph,from=source,to=sink,output="both",mode='out')$epath[[1]])!=0)
    {
      print (c(p,paste0("Num Edges:",length(E(temp.graph)), " ",nrow(allPaths))))
      
      pathList=shortest_paths(temp.graph,from=source,to=sink,output="both",mode='out')
      
      paths=data.frame(t((V(temp.graph)[unlist(pathList$vpath)])$name))
      
      colnames(paths)=seq(1,ncol(paths))
      
      allPaths=rbind.fill(allPaths,paths)
      
      #delete the edges detected in previous discovered path. This is done to include maximum number of nodes in the shortest path.
      temp.graph=delete_edges(temp.graph,pathList$epath[[1]])
    }
    
  }
  print (Sys.time()-s2) 
  return (data.frame(allPaths,stringsAsFactors = F))   #return the unique list of edges in case here are repetations.
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


edgeWeights=function(graphEdges)
{
  #Find the distance between the nodes given in "graphEdges". The distance is given as 
  #number of position wise alphabet matches.
  
  temp=NULL
  for (i in 1:nrow(graphEdges))
  {
    temp<-rbind(temp,  mapply(function(x,y) sum(x!=y),strsplit(graphEdges$from[i],""),strsplit(graphEdges$to[i],"")))
    temp<-rbind(temp,  adist(graphEdges$from[i],graphEdges$to[i]))
    
  }
  return(temp)
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
  
  allPaths=shortPaths(adjMat,source,mainGraph,non_centerNodes)
  
  
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
      distt=edgeWeights(graphEdges = edgeList(temp))
      #store all the edges in the graph
      edge.list=data.frame(rbind(edge.list,cbind.data.frame(x=edgeList(temp),"distVal"=distt)),row.names = NULL)
      temp=NULL
      
    }
    colnames(edge.list)=c("from","to","distVal")
    ##Remove duplicate edges
    edge.list=edge.list[!duplicated(edge.list),]
    mainGraph <- graph.data.frame(edge.list, directed = T)
    l <- layout_with_kk(mainGraph)
    
    # plot(mainGraph, layout=l, vertex.size=3,
    #      vertex.color="red", edge.arrow.size=0.7,vertex.cex=0.55,vertex.label.dist=2,
    #      vertex.label.color="blue4",asp=0,vertex.label.font=12,edge.width=2,vertex.size2=12)
    
    E(mainGraph)$weight=edge.list$distVal
    return (mainGraph)
}




getTS_symbols=function(rawData,swSize)
{

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
  logT_data = paste0(strsplit(basename(file.path(getwd(),"data", fileToRead)),"_")[[1]][1], "_","logTransformed.csv")
  dta.processed=read.csv(file.path(getwd(),"data",logT_data),header=T)
  dta=trainData[,-c(1,2,3)]
  ##Run the function in debug mode to find the range of alphabets.
  #numAlphabets(rawData)
  if (dl=="London Data")
    alpha=seq(10,20)
  else
    alpha=seq(12,20)
  clusteringResults = NULL
  for (w in 1:length(window))
  {
    symbolicRep=getTS_symbols(dta.processed,window[w],alpha)
    
    numClust=read.csv(file.path(p,dl,"Output non-overlap","numClust.csv"),header=T)
    nClust_evw=numClust[which(numClust$W_UW=="EV_W"   & (numClust$window==window[w]) & (numClust$Method=="Proposed")) ,]$numClust
    
    listAllPaths=extractPaths(mainGraph,symbolicRep)
    lengthPaths=apply(listAllPaths,1,function(f)
    {
      length(f[!is.na(f)])/nrow(symbolicRep)
    })
    percentMatch=lengthPaths/nrow(symbolicRep)
    chosenPaths=redundantRemoval(listAllPaths)
    featureDist=featureSpace(chosenPaths,symbolicRep)
    
   numFeatures=getClusters(featureDist,nClust_evw)
   clusteringResults[w] = numFeatures[[2]]
    # print(numFeatures[[2]])
  }
  return (clusteringResults)
}

dataList=c("London Data","Ausgrid Data","Web traffic")
fileToRead=c("londonData_forGraph.csv","AusgridData_forGraph.csv","trimmedWebTraffic_labelled.csv")[1]
window=seq(6,48,6)

for (dl in 1:length(dataList))
{
  out=mainCall(dataList[dl],fileToRead[dl])
  
}
