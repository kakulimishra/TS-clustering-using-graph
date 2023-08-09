###Fig 5 in the paper has graph for each dataset and patterns. This file is used to obtain the graph. 
###The file is executed in the JU system because it needs a package 'ggraph' which could not be installed in my laptop.


library (igraph)
library (stringr)
library(ggraph)
plot_the_graph=function(mainGraph,temporalPatt,edge.list)
{
  my_colorList=c("blue","darkgreen","red","goldenrod4","Pacific blue","darkorchid1","Wild watermelon")
  temporalPatt=temporalPatt[,-c(2,3)]
  

  # tempList<-NULL
  # tempList<-apply(temporalPatt,1,function(f){
  #   clustNo=f[1]
  #   f=f[2:length(f)]
  #   f=unique(as.character(f[!is.na(f)])
  #   tempList<<-rbind(tempList,cbind(rep(clustNo,length(f)),(f)))
  # })
  

    colorNodes=temporalPatt[1,-1]
    #colorNodes=unique(as.character(colorNodes[!is.na(colorNodes)]))[-c(1,13)]      #London
    #colorNodes=unique(as.character(colorNodes[!is.na(colorNodes)]))     #Ausgrid
    #colorNodes=unique(as.character(colorNodes[!is.na(colorNodes)]))     #Stock
    colorNodes=unique(as.character(colorNodes[!is.na(colorNodes)]))     #Web traffic
    
    list1=which(V(mainGraph)$name %in% colorNodes)
    list2=which(!(V(mainGraph)$name %in% colorNodes))
    
    
    
    subVertices=V(mainGraph)[V(mainGraph)[c(list1,list2[1:40])]]
    subG=subgraph(mainGraph,subVertices)
    
    colorEdges=data.frame(cbind("from"=colorNodes[1:(length(colorNodes)-1)],"to"=colorNodes[2:length(colorNodes)]))
    #colorEdges=na.omit(colorEdges)  #stock
    
    getEdge=apply(colorEdges,1,function(colr)
      {
      E(subG) [from(colr[1]) & (to(colr[2])) ]
    })
    eColor=rep("grey66",length(E(subG)))
    eColor[unlist(getEdge)]= "#EC2055" #"#3C804F"
    

    
    

    mySize=rep(2,length(V(subG)))
    mySize[which(V(subG)$name %in% as.character(unlist(colorEdges)))]=8
    startVetex.size=10
    mySize[which(V(subG)$name ==colorEdges$from[1])]= startVetex.size
    
    vColor=rep("black",length(V(subG)))
    vColor[which(V(subG)$name %in% as.character(unlist(colorEdges)))]=my_colorList[6]
    

     eWidth=rep(0.2,length(E(subG)))
     eWidth[unlist(getEdge)]=1
     
     arrowLength=rep(0.15,length(E(subG)))
     arrowLength[unlist(getEdge)]=0.3
     

    ##5C5151
    # p<-ggraph(subG,layout = "dh") +
    #   geom_edge_link0(edge_colour = eColor, edge_width = eWidth, edge_alpha = 0.5,
    #                   arrow = arrow(angle = 30, length = unit(0.10, "inches"), ends = "last", type = "closed")) +
    #   geom_node_point(fill = vColor, colour = "#000000",
    #                   size =mySize , stroke = 0.8, shape = 21) +
    #   theme_graph() + geom_node_label(aes(filter = mySize>=10, label = V(subG)$name[colorNodes]),fontface="bold",label.size = 0.5)+
    #   theme(legend.position = "none")+scale_fill_brewer(palette = "Set1")
     
     
     p<-ggraph(subG,layout = "lgl") +
       geom_edge_link0(edge_colour = eColor, edge_width = eWidth, edge_alpha = 0.5,
                       arrow = arrow(angle = 30, length = unit(arrowLength, "inches"), ends = "last", type = "closed")) +
          geom_node_point(fill = vColor, colour = "#000000",
                          size =mySize , stroke = 0.8, shape = 21) + 
       geom_node_label(aes(filter = (mySize== startVetex.size), label = V(subG)$name[which(V(subG)$name==colorEdges$from[1])]),fontface="bold",label.size = 0.5)+
       theme_graph() 

    
  

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




graphConstruct=function(df,symbolicRep)
{
  pathLength=NULL
  
  allEdge.info=NULL #store info for of all the edges in the graph.
  allVertex.info=NULL
  edge.list=NULL
  
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
  return (list(mainGraph,edge.list))
  
}



mainCall=function(dl,fileToRead,temporalPatt,window)
{
  # dta.processed=read.csv(file.path(p,dl,"Output non-overlap","webTraffic_logTransformed.csv"),header=T)
  # symbolicRep=read.csv(file.path(p,dl,"Output non-overlap","Cluster labels Ei weighted",paste0("symbolicRep_","win",window,".csv")),header=T)
  #   
  #  out=graphConstruct(dta.processed,symbolicRep)
  # mainGraph=out[[1]]
  # edge.list=out[[2]]
  plot_the_graph(mainGraph,temporalPatt,edge.list)

  
  
  return (numFeatures)
}

p="/home/ujjwal/Kakuli/Graph Representation of TS"
dataList=c("London Data","Ausgrid Data","Stock market","Web traffic")[3]
fileToRead="ausgridData_forGraph.csv"
temporalFiles=list.files(file.path(p,dataList[dl],"Output non-overlap","Temporal patterns Ei weighted"))

window=12  #London
window=48 #ausgrid
window=25#stock
window=42#web traffic

for (dl in 1:length(dataList))
{
  vars=which (str_detect(temporalFiles,as.character(paste0("win",window,".csv")))==T)
  temporalPatt=read.csv(file.path(p,dataList[dl],"Output non-overlap","Temporal patterns Ei weighted",temporalFiles[vars]),header=T)
  out=mainCall(dataList[dl],fileToRead,temporalPatt,window)
  
}