### This file has-- normalization and log transformation, segmentation, symbolic representation, weighted graph representation 
##with highest degree as the source node, path extraction using Dijkstras
### feature selection and clustering. This is applicable only for the web traffic dataset. 
##The graph for Ei-WG and Ei-UWG is common for all the datasets.



library(TSclust)
library(igraph)
library(dplyr)
library(reshape2)
library(plyr)
library(cluster)
library(dplyr)
library(factoextra)


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





mapTS=function(symbolicRep,labelled.path)
{
  library(qualV)
  labels=labelled.path[,ncol(labelled.path)]
  labelled.path=labelled.path[,-ncol(labelled.path)]
  
  get.dist=data.frame(apply(symbolicRep,2,function(f)
  {
    apply(labelled.path,1,function(f1)
    {
      
      LCS(as.character(f),as.character(f1[!is.na(f1)]))$LLCS
    })
  }),row.names = NULL)
  
  print (get.dist)
  
  ts.to.path=apply(get.dist,2,function(f) which(f==max(f)))    #location of the path which gives maximum distance.
  clustLabels=lapply(ts.to.path,function(f) labels[f])
  clustLabels=plyr::ldply(clustLabels, rbind)
  write.csv(data.frame(clustLabels),file.path(p,"Results","clustLabels.csv"),row.names = F)
  
}



clusterFormation=function(pathList,distMat)
{
  library(clusteval)
  
  clustVal=NULL
  for (i in 2:(nrow(pathList)/2))
  {
    hclust.complete<-hclust(distMat,method="complete")
    plot(hclust.complete)
    cut_com <- cutree(hclust.complete, k = i)
    
    str(si <- silhouette(cut_com,distMat))
    (ssi <- summary(si))
    clustVal=data.frame(rbind(clustVal,cbind("num"=i,"val"=ssi$avg.width)),stringsAsFactors = F)
    
  }
  
  #clustVal=data.frame(melt(clustVal[,2:3],1),row.names = NULL)
  colnames(sil_score)=c("id","variable","value")
  sil_score$id=sil_score$id+1
  plot(ggplot(clustVal)+geom_line(data=clustVal,aes(x=num,y=val))+ xlab("#Clusters")+ylab("Score")+
         theme_bw())+ scale_x_continuous(breaks = seq(0,max(clustVal$num),5))+theme(panel.background = element_rect(colour="black"),
                                                                                    axis.text = element_text(color="black"),
                                                                                    axis.ticks = element_blank())
  
  
  #elbow rule
  #num=clustVal[which(clustVal$num==min(clustVal$num)),]$num
  
  num=5
  hclust.complete<-hclust(distMat,method="complete")
  plot(hclust.complete)
  cut_com <- cutree(hclust.complete, k = num)
  
  labelled.path <- mutate(pathList, cluster = cut_com)
  
  return (labelled.path)
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



editDist=function(graphEdges)
{
  temp=NULL
  for (i in 1:nrow(graphEdges))
  {
    temp<-rbind(temp,  mapply(function(x,y) sum(x==y),strsplit(graphEdges$from[i],""),strsplit(graphEdges$to[i],"")))
    
  }
  return(temp)
}



toSegments=function(temp,window)
{
  numRow=length(temp)/window
  cnt=data.frame(rep(1:numRow,each=window))
  temp=data.frame(cbind(cnt[1:length(temp),],unlist(temp)))
  colnames(temp)=c("cnt","vals")
  
  return (temp)
}




arrangeData=function(dta,window)
{
  ids=rep(1:window,nrow(dta))
  dta=(cbind("id"=ids[1:nrow(dta)],dta))
  dta=reshape(dta, idvar = "winVal", timevar = "id", direction = "wide")
  return(dta)
}




extractPaths=function(mainGraph,symbolicRep,nClust,edge.list)
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




edgeWeights=function(graphEdges)
{
  #Find the distance between the nodes given in "graphEdges". The distance is given as 
  #number of position wise alphabet matches.
  
  temp=NULL
  for (i in 1:nrow(graphEdges))
  {
    #distance between two symbols is how many alphabets mismatches are present
    #temp<-rbind(temp,  mapply(function(x,y) sum(x!=y),strsplit(graphEdges$from[i],""),strsplit(graphEdges$to[i],"")))
    temp<-rbind(temp,  adist(graphEdges$from[i],graphEdges$to[i]))
    
  }
  return(temp)
}





graphConstruct=function(df,dateVals,symbolicRep,houseNames)
{
  pathLength=NULL
  
  allEdge.info=NULL #store info for of all the edges in the graph.
  allVertex.info=NULL
  edge.list=NULL
  
  
  for (i in 1:nrow(symbolicRep))
  {
    
    pathLength=NULL
    
    allEdge.info=NULL #store info for of all the edges in the graph.
    allVertex.info=NULL
    edge.list=NULL
    
    
    for (i in 1:nrow(symbolicRep))
    {
      
      temp=t(symbolicRep[i,])
      
      distt=edgeWeights(graphEdges = edgeList(temp))
      #print (simMeasure)
      
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
    
    
    ##The subgraph plot below is used in the model diagram
    # g=subgraph(mainGraph,V(mainGraph)[1:12])
    # V(g)$label.cex=1.5  #adjust the font size of vertex labels
    # plot(g, layout=layout_with_fr, vertex.size=3,
    #    vertex.color="red", edge.arrow.size=0.7,vertex.cex=14,vertex.label.dist=3,
    #      edge.color="darkgreen",vertex.label.color="blue4",asp=0,vertex.label.font=12,edge.width=2)
    
    
    E(mainGraph)$weight=edge.list$distVal
    return (mainGraph)
  }
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
  
  TS=(apply(TS, 2,function(x) { as.data.frame(((replace(x, x == 0, min(x[x!=0])/2)))) }))
  
  TS=do.call(cbind,TS)
  TS.logTransform=log10(TS)
  
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
  trainData=read.csv(file.path(p,fileToRead),header=T)
  
  dateVals=(trainData[,2])
  houseNames=as.character(trainData$houseId)
  
  
  dataFreq=strsplit(basename(file.path(p,fileToRead)),"_")[[1]][3]
  dataFreq=gsub(dataFreq,pattern=".csv",replacement="")
  
  dta=trainData[,-c(1,2,3)]
  
  #dta.processed=preprocess(dta)
  dta.processed=read.csv(file.path(p,"Output non-overlap","webTraffic_logTransformed.csv"),header=T)
  
  
  for (w in 1:length(window))
  {
    print (window[w])
    
    ## The symbolic rpresentation is same as that of the unweighted part.hence the symbolic representation as not been 
    ### performed here again. It is directly read.
    
    symbolicRep=read.csv(file.path(p,"Output non-overlap","Cluster labels Deg weighted",paste0("symbolicRep_","win",window[w],".csv")),header=T)
    numClust=read.csv(file.path(p,"Output non-overlap","numClust.csv"),header=T)
    nClust=numClust[which(numClust$W_UW=="W"  &  numClust$Method=="Proposed"  &  numClust$window==window[w]),]$numClust
    
    #  mainGraph=graphConstruct(dta.processed,dateVals,symbolicRep,houseNames)
    # 
    # timeStart=Sys.time()
    # listAllPaths=extractPaths(mainGraph,symbolicRep,nClust,edge.list)
    # 
    #  timeEnd=Sys.time()
    # 
    #  print (paste0("Time taken for finding paths::", difftime(timeEnd, timeStart, units='mins')))
    # 
    #  write.csv(data.frame(listAllPaths),file.path(p,"Output non-overlap","Cluster labels Deg weighted",paste0("EDweighted_allPaths_","win",window[w],".csv")),row.names=F)
    # listAllPaths=read.csv(file.path(p,"Output non-overlap","Cluster labels Deg weighted",paste0("EDweighted_allPaths_","win",window[w],".csv")),header=T)
    # 
    # lengthPaths=apply(listAllPaths,1,function(f)
    #   {
    #   length(f[!is.na(f)])/nrow(symbolicRep)
    # })
    # 
    # chosenPaths=listAllPaths[which(lengthPaths>0.05),]
    # chosenPaths=redundantRemoval(listAllPaths)
    # write.csv(data.frame(chosenPaths),file.path(p,"Output non-overlap","Cluster labels Deg weighted",paste0("EDweighted_chosenPaths_","win",window[w],".csv")),row.names=F)
    # 
    #  chosenPaths=read.csv(file.path(p,"Output non-overlap","Cluster labels Deg weighted",paste0("EDweighted_chosenPaths_","win",window[w],".csv")),header=T)
    # 
    # featureDist=featureSpace(chosenPaths,symbolicRep)
    # write.csv(data.frame(t(featureDist)),file.path(p,"Output non-overlap","Cluster labels Deg weighted",paste0("EDweighted_featureSpace_","win",window[w],".csv")),row.names=F)
    featureDist=read.csv(file.path(p,"Output non-overlap","Cluster labels Deg weighted",paste0("EDweighted_featureSpace_","win",window[w],".csv")),header=T)
    
    
    ### Get the #clusters using elbow rule and keep it saved in a file to refer later. Comment out the lines when not in use
    #fviz_nbclust((featureDist), kmeans, method = "wss",k.max=floor(nrow(dta)/8))
    #a= fviz_nbclust((featureDist), kmeans, method = "wss",k.max=floor(nrow(dta)/8))
    #View(a$data)
    
    # startTime=Sys.time()
    
    numFeatures=getClusters(featureDist,nClust)
    # print(numFeatures[[2]])
    # endTime=Sys.time()
    # print (paste0("Time taken for finding paths::", difftime(startTime, endTime, units='mins')))
    
    
    # pc<-pam(featureDist,diss=F,k=nClust,metric='euclidean')
    # label=pc$clustering
    
    
    
    write.csv(data.frame(numFeatures[[2]]),file.path(p,"Output non-overlap","Cluster labels Deg weighted",paste0("EDweighted_ProposedLabels_","win",window[w],".csv")),row.names = F)
    #write.csv(data.frame(label),file.path(p,"Output non-overlap","Cluster labels Deg weighted",paste0("ProposedLabels_",dataFreq,"_wl",window[w],".csv")),row.names = F)
    
    #return (numFeatures)
    
  }
  
}
p="/home/kakuli/Graph Representation of TS/Web traffic"
fileToRead="webTrafficData_forGraph.csv"
window=c(7,14,21,28,35,42,49)
out=mainCall(fileToRead)
