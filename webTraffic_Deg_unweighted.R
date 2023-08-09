
### This file has-- normalization and log transformation, segmentation, symbolic representation, graph representation 
##with highest degree as the source node, path extraction using DFS
### feature selection and clustering. This is applicable only for the web traffic dataset. 
##The graph for Ei-WG and Ei-UWG is common for all the datasets.


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

totalWSS=9999
getClusters=function(distMat,k)
{
  
  library(tidyverse)
  library(factoextra)
  library(cluster)
  library(clusteval)
  
  cri=NULL
  
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



mapTS=function(symbolicRep,labelledPath)
{
  library(qualV)
  labels=labelledPath[,ncol(labelledPath)]
  labelledPath=labelledPath[,-ncol(labelledPath)]
  
  get.dist=data.frame(apply(symbolicRep,1,function(f)
  {
    apply(labelledPath,1,function(f1)
    {
      
      LCS(as.character(f),as.character(f1[!is.na(f1)]))$LLCS
    })
  }),row.names = NULL)
  
  #the "get.dist dataframe has LCS distance between the paths and the symbolically represented time series.
  #each row corresponds to a TS and each column corresponds to path. 
  
  
  ts.to.path=apply(get.dist,2,function(f) which(f==max(f)))    #location of the path which gives maximum distance.
  
  #extract the labels for the respective locations 
  clustLabels=lapply(ts.to.path,function(f) labels[f])
  
  #the maximum frequency of occurence of the labels is considered as the label for TS
  clustLabels=lapply(clustLabels,function(f)
  {
    temp=data.frame(table(f))
    as.numeric(as.character(temp[which.max(temp$Freq),1]))
  })
  
  predLabel=do.call(rbind,clustLabels)
  return (list(predLabel,(get.dist)))
}



distanceComputation=function(pathList)
{
  distMat=matrix(data=NA,nrow=nrow(pathList),ncol=nrow(pathList))
  
  for (i in 1:nrow(pathList))
  {
    #print (paste0("Row ",i))
    j=i+1
    distMat[i,i]=1
    x=pathList[i,]
    x=x[!is.na(x)]
    while (j<=nrow(pathList))
    {
      y=pathList[j,]
      y=y[!is.na(y)]
      
      maxLen=max(length(x),length(y))
      
      #larger diffRatio signifies more similar paths.
      diffRatio=(length(intersect(x,y)))/maxLen
      
      distMat[i,j]= diffRatio
      distMat[j,i]=diffRatio
      j=j+1
    }
  }
  
  #based on the inbuilt clustering techniques in R, the smaller distance values are merged to form a cluster. 
  # In this case, larger distance signifies more similar paths. hence
  #I have used 1- (the distance value), to convert it to a smaller value
  
  distMat= 1-distMat
  colnames(distMat)=seq(1,nrow(pathList))
  rownames(distMat)=seq(1,nrow(pathList))
  return (distMat)
}



shortPaths=function(adj.graph,source,graph,non_centerNodes,minPath.length,maxPath.length)
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



gapFunction=function(graphPath,nClust,symbolicRep)
{
  nClust=3
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
    if (r<(1-(1/nClust)) & r>1/nClust)
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



findDfsPath=function(graph,source,adj.graph,listAllPaths,non_centerNodes,nClust,
                     symbolicRep)
{
  library (data.table)
  graph.elist=data.frame(get.edgelist(graph),stringsAsFactors = F)
  E(graph)$name=seq(1,length(E(graph)))
  
  source=V(graph)$name[source]
  
  s1=Sys.time()
  
  for (p in 1:length(non_centerNodes))
  {
    print (c(p,nrow(listAllPaths)))
    s2=Sys.time()
    sink=V(graph)$name[non_centerNodes[p]]
    temp.graph=graph
    #adj=adj.graph
    
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
        
        # if(length(temp)>=minPath.length && length(temp)<=maxPath.length && which(dfs.out==source)>which(dfs.out==sink))
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


edgeWeights=function(graphEdges)
{
  #Find the distance between the nodes given in "graphEdges". The distance is given as 
  #number of position wise alphabet matches.
  
  temp=NULL
  for (i in 1:nrow(graphEdges))
  {
    #distance between two symbols is how many alphabets mismatches are present
    temp<-rbind(temp,  mapply(function(x,y) sum(x!=y),strsplit(graphEdges$from[i],""),strsplit(graphEdges$to[i],"")))
    
  }
  return(temp)
}


nonOverlap=function (timeseries, swSize) 
{
  ####  The function is used for obtaining the non-overlapping segments
  
  SW<-NULL
  for (x in seq(1,length(timeseries) ,swSize))
  {
    windowCol=timeseries[x:(x+swSize-1)]
    SW <- cbind(SW, windowCol, deparse.level = 0)
  }
  SW=t(SW)
  SW=na.omit(SW)
  row.names(SW)<-NULL
  return (SW)
}



segmentDates=function(temp,window)
{
  #temp=data.frame(as.numeric(df))
  numRow=nrow(temp)/window
  cnt=data.frame(rep(1:numRow,each=window))
  temp=cbind.data.frame(cnt[1:nrow(temp),],(temp))
  colnames(temp)=c("cnt","date")
  
  temp1<-temp %>%
    group_by(cnt) %>%
    dplyr::mutate(
      first = dplyr::first(date))
  temp1=unique(temp1$first)
  
  return (data.frame("DateTime"=temp1,stringsAsFactors = F))
  
}


extractPaths=function(mainGraph,symbolicRep,minPath.length,maxPath.length,nClust)
{
  outDeg=degree(mainGraph,mode='out')
  
  outDeg=data.frame(cbind("vertexName"=names(outDeg),"deg"=outDeg),stringsAsFactors = F)
  centerNodes=which.max(outDeg$deg)
  
  adjMat=data.frame(as.matrix(as_adjacency_matrix(mainGraph)),stringsAsFactors = F)
  
  non_centerNodes=V(mainGraph)[! V(mainGraph) %in% V(mainGraph)[centerNodes]]
  
  listAllPaths=NULL
  
  start=Sys.time()
  
  source=V(mainGraph)[centerNodes]
  allPaths=findDfsPath(graph=mainGraph,source,adj.graph=adjMat,listAllPaths,non_centerNodes,
                       nClust,symbolicRep )
  #allPaths=shortPaths(adjMat,source,mainGraph,non_centerNodes,minPath.length,maxPath.length)
  
  colnames(allPaths)=seq(1,ncol(allPaths))
  
  listAllPaths=rbind.fill(listAllPaths,allPaths)
  
  print (Sys.time()-start)
  
  return (listAllPaths)
  
}



graphConstruct=function(df,symbolicRep,houseNames)
{
  pathLength=NULL
  
  allEdge.info=NULL #store info for of all the edges in the graph.
  allVertex.info=NULL
  edge.list=NULL
  
  
  for (i in 1:nrow(symbolicRep))
  {
    
    temp=t(symbolicRep[i,])
    
    #temp=temp[!duplicated(temp), ]  #this will capture the unique symbols foreach row.
    
    distt=edgeWeights(graphEdges = edgeList(temp))
    #print (simMeasure)
    
    #pathLength=rbind(pathLength,nrow(temp)-1)
    
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
  #      vertex.color="red", edge.arrow.size=0.2,vertex.cex=0.55,vertex.label.dist=2,
  #      vertex.label.color="blue4",asp=0,vertex.label.font=12,edge.width=2,vertex.size2=12)
  
  
  ##The subgraph plot below is used in the model diagram
  # g=subgraph(mainGraph,V(mainGraph)[1:12])
  # V(g)$label.cex=1.5  #adjust the font size of vertex labels
  # plot(g, layout=layout_with_fr, vertex.size=3,
  #      vertex.color="red", edge.arrow.size=0.2,vertex.cex=14,vertex.label.dist=3,
  #      edge.color="darkgreen",vertex.label.color="blue4",asp=0,vertex.label.font=12,edge.width=2)
  
  
  E(mainGraph)$weight=edge.list$distVal
  
  return (mainGraph)
  
}


numAlphabets=function(rawData)
{
  ###Range of alphabets are decided seeing the data distribution. The percent of values below 0.25 and above 0.95 are less than 0.1%
  ##of the raw data, hence in order to keep those values in the same window, the number of alphabets ranges from 2 to 20.
  ##This function is run in the debug mode.
  
  alpha=seq(2,20)
  breakpointsTable=sapply(alpha,function(f) qnorm(0:f/f))
  temp=NULL
  breakpointsTable=lapply(breakpointsTable,function(f) temp<<-cbind.fill(temp,(f),fill=NA))
  breakpointsTable=breakpointsTable[[length(breakpointsTable)]]
  breakpointsTable=data.frame(breakpointsTable[,which(unlist(lapply(breakpointsTable, function(x) !all(is.na(x)))))])
  
  histData=hist(unlist(rawData))
  View(cbind(histData$breaks,histData$counts/length(unlist(rawData))))
  
}



getTS_symbols=function(rawData,swSize)
{
  library(rowr)
  
  ##Run the function in debug mode to find the range of alphabets.
  #numAlphabets(rawData)
  

  alpha=seq(7,19)
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
    
    # print (i)
    
    dta=unlist(rawData[i,])
    
    segData=nonOverlap(dta,swSize)
    #segData is a dataframe where each row corresponds to a segment. 
    
    meanVals= rbind(meanVals,t(apply(segData,1,function(f) mean(unlist(f)))))
    
    sdVals=  rbind(sdVals,t(apply(segData,1,function(f) sd(unlist(f)))))
    
    maxVals=  rbind(maxVals,t(apply(segData,1,function(f) which.max(unlist(f)))))
    
    minVals=  rbind(minVals,t(apply(segData,1,function(f) which.min(unlist(f)))))
    
  }
  
  
  #Ths number of alphabets is decided for each segment- seg1, seg2,....of all time series. eg: to decide the number of alphabets for mean,
  #the median of the mean values is obtained across the segments of all time series. (seg1, seg2, and so on).
  #The alphabet which gives minimum difference between median and the beta value, is chosen as #alphabets.
  #time series in the dataset
  
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
  
  #allSym= matrix( paste(meanSym,sdSym, maxSym,minSym, sep=""),nrow=nrow(meanSym))
  allSym= matrix( paste(meanSym,sdSym, maxSym,minSym, sep=""),nrow=nrow(meanSym))
  
  
  return (as.data.frame(allSym))
}

preprocess=function(TS)
{
  ### To convert the data to a normal distribution, log10 transform is done.
  ### But as log10 transform should not have '0', all the zeros are replaced by smaller 
  #values -- min(x[x!=0])/2.
  
  # TS=data.frame(apply(TS, 2,function(x) replace(x, x == 0, min(x[x!=0])/2)))
  # 
  # 
  # ##log transform
  # TS.logTransform=log10(TS)
  # 
  # ##scaling
  # TS.scaled=scale(TS.logTransform)    #return (TS.scaled)
  
  ###~~~~~~~~~~~~~~~~~~~~~~~~~~...........OR..........~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  TS=data.frame(apply(TS, 2,function(x) replace(x, x == 0, min(x[x!=0])/2)))
  start=Sys.time()
  TS.logTransform=log10(TS)
  print (Sys.time()-start)
  ##Normalize
  
  TS.norm=apply(TS.logTransform,2,function(f)
  {
    (f-min(f))/(max(f)-min(f))
  })
  
  #TS.norm=(TS.logTransform-min(TS.logTransform))/(max(TS.logTransform)-min(TS.logTransform))
  
  return(as.data.frame(TS.norm))
}



mainCall=function(fileToRead)
{
  
  trainData=read.csv(file.path(p,"Processed data",fileToRead),header=T)
  
  dta=trainData[,-c(1,2,3)]
  
  #dta.processed=preprocess(dta)
  #write.csv(as.data.frame(dta.processed),file.path(p,"Output non-overlap","webTraffic_logTransformed.csv"),row.names=F)
  
  dta.processed=read.csv(file.path(p,"Output non-overlap","webTraffic_logTransformed.csv"),header=T)
  for (w in 1:length(window))
  {
    # print (window[w])
    # start=Sys.time()
    # symbolicRep=getTS_symbols(dta.processed,window[w])
    # print (Sys.time()-start)
    # write.csv(data.frame(symbolicRep),file.path(p,"Output non-overlap","Cluster labels Deg unweighted",paste0("symbolicRep_","win",window[w],".csv")),row.names=F)
    # symbolicRep=read.csv(file.path(p,"Output non-overlap","Cluster labels Deg unweighted",paste0("symbolicRep_","win",window[w],".csv")),header=T)
    # 
    # 
    # #mainGraph=graphConstruct(dta.processed,symbolicRep,houseNames)
    # 
    numClust=read.csv(file.path(p,"Output non-overlap","numClust.csv"),header=T,stringsAsFactors = F)
    nClust=numClust[which(numClust$W_UW=="UW"  &   numClust$Method=="Proposed"  &  numClust$window==window[w]),]$numClust
    # 
    # start=Sys.time()
    # listAllPaths=extractPaths(mainGraph,symbolicRep,nClust)
    # write.csv(data.frame(listAllPaths),file.path(p,"Output non-overlap","Cluster labels Deg unweighted",paste0("allPaths_","win",window[w],".csv")),row.names=F)
    # listAllPaths=read.csv(file.path(p,"Output non-overlap","Cluster labels Deg unweighted",paste0("allPaths_","win",window[w],".csv")),header=T)
    # 
    # lengthPaths=apply(listAllPaths,1,function(f)
    #   {
    #   length(f[!is.na(f)])/nrow(symbolicRep)
    # })
    # percentMatch=lengthPaths/nrow(symbolicRep)
    # 
    # chosenPaths=listAllPaths[which(lengthPaths>=0.004),]
    # chosenPaths=redundantRemoval(chosenPaths)
    # write.csv(data.frame(chosenPaths),file.path(p,"Output non-overlap","Cluster labels Deg unweighted",paste0("chosenPaths_","win",window[w],".csv")),row.names=F)
    # 
    # chosenPaths=read.csv(file.path(p,"Output non-overlap","Cluster labels Deg unweighted",paste0("chosenPaths_","win",window[w],".csv")),header=T)
    # 
    # 
    # 
    # timeStart=Sys.time()
    # featureDist=featureSpace(chosenPaths,symbolicRep)
    # timeEnd=Sys.time()
    # 
    # print (paste0("Time taken for finding paths::", difftime(timeEnd, timeStart, units='mins')))
    # 
    # 
    # write.csv(data.frame(t(featureDist)),file.path(p,"Output non-overlap","Cluster labels Deg unweighted",paste0("featureSpace_","win",window[w],".csv")),row.names=F)
    featureDist=read.csv(file.path(p,"Output non-overlap","Cluster labels Deg unweighted",paste0("featureSpace_","win",window[w],".csv")),header=T)
    
    
    # ### Get the #clusters using elbow rule and keep it saved in a file to refer later. Comment out the lines when not in use
    # fviz_nbclust((featureDist), kmeans, method = "wss",k.max=floor(nrow(dta)/8))
    # a= fviz_nbclust((featureDist), kmeans, method = "wss",k.max=floor(nrow(dta)/8))
    # View(a$data)
    # 
    # # startTime=Sys.time()
    #numFeatures=getClusters(featureDist,nClust)
    # #print(numFeatures[[2]])
    # # endTime=Sys.time()
    # # print (paste0("Time taken for finding paths::", difftime(startTime, endTime, units='mins')))
    # 
    # 
    pc<-pam(featureDist,diss=F,k=nClust,metric='euclidean')
    label=pc$clustering

    write.csv(data.frame(label),file.path(p,"Output non-overlap","Cluster labels Deg unweighted",paste0("PAMLabels_","win",window[w],".csv")),row.names = F)
    
    
    # #  cl<-kmeans((featureDist[,1:numFeatures]),nClust,iter.max = 500,nstart=25)
    # write.csv(data.frame(numFeatures[[2]]),file.path(p,dl,"Output non-overlap","Cluster labels Deg unweighted",paste0("ProposedLabels_","win",window[w],".csv")),row.names = F)
    # #write.csv(data.frame(label),file.path(p,dl,"Output non-overlap","Cluster labels Deg unweighted",paste0("ProposedLabels_",dataFreq,"_wl",window[w],".csv")),row.names = F)
    
  }
  
  return (numFeatures)
}

p="/home/kakuli/Graph Representation of TS/Web traffic"
fileToRead="webTrafficData_forGraph.csv"
window=c(7,14,21,28,35,42,49)

out=mainCall(fileToRead)
  