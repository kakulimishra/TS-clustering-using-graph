###This file is to obtain weighted graph using highest degree as the source  for non-overlapping windows. 
### This file has-- normalization and log transformation, segmentation, symbolic representation, graph representation, path extraction using Dijkstras



library(TSclust)
library(igraph)
library(dplyr)
library(reshape2)
library(plyr)
library(cluster)
library(dplyr)


getClusters=function(distMat,k)
{
  
  library(tidyverse)
  library(factoextra)
  library(cluster)
  library(clusteval)
  
  cri=NULL
  totalWSS=9999
  
  labels=rep(1,nrow(distMat))
  newLabels=list()
  newLabels[[1]]=rep(1,nrow(distMat))
  
  #compute the change in rand index
  for (i in 2:ncol(distMat))
  {
    print (i)
    if (length(unique(unlist(distMat[,1:i])))>k)
    {
      cl<-kmeans(distMat[,1:i],k,iter.max = 500)
      newLabels[[i]]=cl$cluster
      prev=which(unlist(lapply(newLabels,length))>0)
      prev=prev[prev!=i]
      temp=1-cluster_similarity(unlist(newLabels[[i]]), unlist(newLabels[prev[length(prev)]]),similarity="rand")
      cri=rbind(cri,cbind(temp,i))
    }
    else
      next
  }
  cri=as.data.frame(cri)
  return (list(cri$i[which.min(cri$temp)],newLabels[[cri$i[which.min(cri$temp)]]]))
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
      #cleanlist[[x]]=pathList[[which.max((partialMatch>0.85 & partialMatch<1))]]
      cleanlist[[x]]=pathList[[x]]
    }
    
  }
  cleanlist=cleanlist[sapply(cleanlist,function(f) length(f)!=0)]
  cleanlist=do.call(rbind.fill,lapply(cleanlist,function(f) as.data.frame(t(f))))
  
  return (cleanlist)
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
  
  # temp.graph=graph
  # adj=adj.graph
  
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


extractPaths=function(mainGraph,symbolicRep,nClust_w,edge.list)
{

  outDeg=degree(mainGraph,mode='out')
  
  outDeg=data.frame(cbind("vertexName"=names(outDeg),"deg"=outDeg),stringsAsFactors = F)
  centerNodes=which.max(outDeg$deg)
  
  adjMat=data.frame(as.matrix(as_adjacency_matrix(mainGraph)),stringsAsFactors = F)
  
  non_centerNodes=V(mainGraph)[! V(mainGraph) %in% V(mainGraph)[centerNodes]]

  start=Sys.time()
  
  source=V(mainGraph)[centerNodes]

  allPaths=shortPaths(adjMat,source,mainGraph,non_centerNodes)
  
  colnames(allPaths)=seq(1,ncol(allPaths))
  
  print (Sys.time()-start)
  
  return (allPaths)
  
}


getTS_symbols=function(rawData,swSize,alpha)
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


edgeWeights=function(graphEdges)
{
  #Find the distance between the nodes given in "graphEdges". The distance is given as 
  #number of position wise alphabet matches.
  
  temp=NULL
  for (i in 1:nrow(graphEdges))
  {
    temp<-rbind(temp,  adist(graphEdges$from[i],graphEdges$to[i]))
    
  }
  return(temp)
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

mainCall=function(dl,fileToRead)
{
  trainData=read.csv(file.path(getwd(),"data",fileToRead),header=T)

  dataNm=strsplit(dl," ")[[1]][1]
  print (dataNm)

  dta=trainData[,-c(1,2,3)]
  logT_data = paste0(strsplit(basename(file.path(getwd(),"data", fileToRead)),"_")[[1]][1], "_","logTransformed.csv")
  
  dta.processed=read.csv(file.path(getwd(),"data",logT_data),header=T)
  ##Run the function in debug mode to find the range of alphabets.
  #numAlphabets(rawData)
  
  ##alpha=seq(10,20)  ##For London
  ##alpha=seq(12,20)  ##For Ausgrid
  
  if (dl=="London Data")
    alpha=seq(10,20)
  else
    alpha=seq(12,20)
  clusteringResults = NULL
  
  for (w in 1:length(window))
  {
    print (window[w])
    symbolicRep=getTS_symbols(dta.processed,window[w],alpha)
    ### Get the #clusters using elbow rule and keep it saved in a file to refer later.
    clusterFile = paste0(strsplit(basename(file.path(getwd(),"data", fileToRead)),"_")[[1]][1], "_","numClust.csv")
    numClust=read.csv(file.path(getwd(),"data" ,clusterFile), header=T,stringsAsFactors = F)
    nClust_w=numClust[which(numClust$W_UW=="W"  &  numClust$Method=="Proposed"  &  numClust$window==window[w]),]$numClust
    
    mainGraph=graphConstruct(dta.processed,symbolicRep)
    listAllPaths=extractPaths(mainGraph,symbolicRep,nClust_w,edge.list)
     lengthPaths=apply(listAllPaths,1,function(f)
       {
       length(f[!is.na(f)])/nrow(symbolicRep)
     })

     #chosenPaths=listAllPaths[which(lengthPaths>0.15),]
     chosenPaths=redundantRemoval(listAllPaths)
    featureDist=featureSpace(chosenPaths,symbolicRep)

    numFeatures=getClusters(featureDist,nClust_w)
    # print(numFeatures[[2]])
    clusteringResults[w] = numFeatures[[2]]
    
 }
  return (clusteringResults)
}

dataList=c("London Data","Ausgrid Data","Stock market")
fileToRead=c("londonData_forGraph.csv","AusgridData_forGraph.csv","stockMarketData_forGraph.csv")

window= seq(6,48,6)
for (dl in 1:length(dataList))
{
  out=mainCall(dataList[dl],fileToRead)
  
}