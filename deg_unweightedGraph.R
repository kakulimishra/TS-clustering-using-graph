###This is the last and the final version of unweighted graph representation for non-overlapping windows. The other saved copies of graphRep are lder copies.
### This file has-- normalization and log transformation, segmentation, symbolic representation, graph representation, ath extraction using DFS
### feature selection and clustering

# library(TSclust)
# library(igraph)
# library(dplyr)
# library(reshape2)
# library(plyr)
# library(cluster)
# library(dplyr)
# library(clusterSim)
# library(factoextra)
# library(tidyverse)
# library(cluster)
library(clusteval)
totalWSS=9999



getClusters=function(distMat,k)
{
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
      cleanlist[[x]]=pathList[[x]]
    }
    
  }
  cleanlist=cleanlist[sapply(cleanlist,function(f) length(f)!=0)]
  cleanlist=do.call(rbind.fill,lapply(cleanlist,function(f) as.data.frame(t(f))))
  
  return (cleanlist)
}



findDfsPath=function(graph,source,adj.graph,listAllPaths,non_centerNodes,nClust_uw,
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
      
      ####The "commonNodes are removed . This helps in finding out the longer length paths between the nodes.
      commonNodes=intersect(dfs.in[1:(which(dfs.in==sink)-1)],dfs.out[1:(which(dfs.out==sink)-1)])
      
      temp=dfs.in[1:which(dfs.in==sink)]
      temp=temp[!(temp %in% commonNodes)]
      
      if (length(commonNodes)>0)
        print (commonNodes)
      
      if ( length(temp)==1)
        break
      
      
      else if ( which(dfs.out==source)>which(dfs.out==sink))  ##this condition ensures that the path exists b/w source and the sink
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


extractPaths=function(mainGraph,symbolicRep,minPath.length,maxPath.length,nClust_uw)
{
  outDeg=degree(mainGraph,mode='out')
  
  outDeg=data.frame(cbind("vertexName"=names(outDeg),"deg"=outDeg),stringsAsFactors = F)
  centerNodes=which.max(outDeg$deg)
  
  adjMat=data.frame(as.matrix(as_adjacency_matrix(mainGraph)),stringsAsFactors = F)
  
  non_centerNodes=V(mainGraph)[! V(mainGraph) %in% V(mainGraph)[centerNodes]]
  
  listAllPaths=NULL
  
  source=V(mainGraph)[centerNodes]
  allPaths=findDfsPath(graph=mainGraph,source,adj.graph=adjMat,listAllPaths,non_centerNodes,
                       nClust_uw,symbolicRep )
  
  ##Another method, to explore
  #allPaths=shortPaths(adjMat,source,mainGraph,non_centerNodes,minPath.length,maxPath.length)

  colnames(allPaths)=seq(1,ncol(allPaths))
  
  listAllPaths=rbind.fill(listAllPaths,allPaths)

  ##**********************************graph plot*************************************************
  
  E(mainGraph)$name=seq(1,length(E(mainGraph)))
  mySize=rep(2,length(V(mainGraph)))
  mySize[centerNodes]=12
  
  vColors=rep("#281EE8",length(V(mainGraph)))
  vColors[centerNodes]="#BF2C2C"
  
  eWidth=rep(0.2,length(E(mainGraph)))
  eWidth[E(mainGraph)[from(centerNodes)]$name]=1
  
  eColor=rep("grey66",length(E(mainGraph)))
  eColor[E(mainGraph)[from(centerNodes)]$name]= "#EC2055" #"#3C804F"

  p<-ggraph(mainGraph,layout = "focus",focus=centerNodes) + 
    geom_edge_link0(edge_colour = eColor, edge_width = eWidth, edge_alpha = 0.5, 
                    arrow = arrow(angle = 30, length = unit(0.10, "inches"), ends = "last", type = "closed")) + 
    geom_node_point(fill = vColors, colour = "#000000", 
                    size =mySize , stroke = 0.8, shape = 21) + 
    theme_graph() + geom_node_label(aes(filter = mySize>=10, label = V(mainGraph)$name[centerNodes]),fontface="bold",label.size = 0.5)+
    theme(legend.position = "none")
  
##***************************************************************  
  
  return (listAllPaths)
  
}



graphConstruct=function(df,symbolicRep)
{

  for (i in 1:nrow(symbolicRep))
  {
    temp=t(symbolicRep[i,])
    #store all the edges in the graph
    edge.list=data.frame(rbind(edge.list,x=edgeList(temp)),row.names = NULL)
    temp=NULL
    
  }
  colnames(edge.list)=c("from","to")
  ##Remove duplicate edges
  edge.list=edge.list[!duplicated(edge.list),]
  mainGraph <- graph.data.frame(edge.list, directed = T)
  l <- layout_with_kk(mainGraph)
  
  # plot(mainGraph, layout=l, vertex.size=3,
  #      vertex.color="red", edge.arrow.size=0.2,vertex.cex=0.55,vertex.label.dist=2,
  #      vertex.label.color="blue4",asp=0,vertex.label.font=12,edge.width=2,vertex.size2=12)
  

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



getTS_symbols=function(rawData,swSize,alpha)
{
  # library(rowr)
  
  ##Run the function in debug mode to find the range of alphabets.
  #numAlphabets(rawData)
  
  ##alpha=seq(10,20)  ##For London
  ##alpha=seq(12,20)  ##For Ausgrid
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

mainCall=function(dl,fileToRead)
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

    symbolicRep=getTS_symbols(dta.processed,window[w],alpha)
    
    ### Get the #clusters using elbow rule and keep it saved in a file to refer later.

    clusterFile = paste0(strsplit(basename(file.path(getwd(),"data", fileToRead)),"_")[[1]][1], "_","numClust.csv")
    
    numClust=read.csv(file.path(getwd(),"data" ,clusterFile), header=T,stringsAsFactors = F)
    nClust_uw=numClust[which(numClust$W_UW=="UW"  &  numClust$Method=="Proposed"  &  numClust$window==window[w]),]$numClust

    mainGraph=graphConstruct(dta.processed,symbolicRep)

    listAllPaths=extractPaths(mainGraph,symbolicRep,nClust_uw)
  
    lengthPaths=apply(listAllPaths,1,function(f)
      {
      length(f[!is.na(f)])/nrow(symbolicRep)
    })
    percentMatch=lengthPaths/nrow(symbolicRep)

    chosenPaths=listAllPaths[which(lengthPaths > 0.15),]
    chosenPaths=redundantRemoval(chosenPaths)

    featureDist=featureSpace(chosenPaths,symbolicRep)

    numFeatures=getClusters(featureDist,nClust_uw)
    ##Save numFeatures[[2]] as the cluster labels
    # print (numFeatures[[2]])
    clusteringResults[w] = numFeatures[[2]]
  }
  return (clusteringResults)
}

dataList=c("London Data","Ausgrid Data","Stock market")
fileToRead=c("londonData_forGraph.csv","AusgridData_forGraph.csv","stockMarketData_forGraph.csv")

window=seq(6,48,6)
for (dl in 1:length(dataList))
{
  out=mainCall(dataList[dl],fileToRead[dl])
  
}
