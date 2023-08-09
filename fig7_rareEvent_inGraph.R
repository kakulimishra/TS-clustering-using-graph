library (igraph)
library (stringr)
library(ggraph)

# 
# plot_for_london=function(mainGraph,anEvents,edge.list)
# {
#   my_colorList=c("blue","darkgreen","red","goldenrod4","gold","darkorchid1","palegreen1","sienna1")
#   colorNodes=unique(as.character((anoEvents$anomaly)))
#   
#   
#   list1=which(V(mainGraph)$name %in% colorNodes)
#   list2=which(!(V(mainGraph)$name %in% colorNodes))
#   
#   
#   
#   #subVertices=V(mainGraph)[V(mainGraph)[c(list1,list2[1:40])]]
#   #subG=induced_subgraph(mainGraph,subVertices)
#   
#   mySize=rep(2,length(V(mainGraph)))
#   mySize[which(V(mainGraph)$name %in% colorNodes)]=13
#   
#   vColor=rep("cyan4",length(V(mainGraph)))
#   vColor[which(V(mainGraph)$name %in% colorNodes)]=my_colorList
#   
#   
#   ##append edges to mainGraph
#   E(mainGraph)$name=seq(1, length(E(mainGraph)))
#   appendEdges=c(colorNodes[1],rep(colorNodes[2:7],each=2), colorNodes[length(colorNodes)])
#   mainGraph=add_edges(mainGraph,appendEdges,directed=T)
#   eColor=rep("grey66",length(E(mainGraph)))
#   eColor[get.edge.ids(mainGraph,appendEdges)]=my_colorList[1:7]
#   
#   eWidth=rep(0.2,length(E(mainGraph)))
#   eWidth[get.edge.ids(mainGraph,appendEdges)]=3
#   
#   my_arrowSize=rep(0.15,length(E(mainGraph)))
#   my_arrowSize[get.edge.ids(mainGraph,appendEdges)]=0.3
#   
#   
#   p<-ggraph(mainGraph,layout = "lgl") +
#     geom_edge_link0(edge_colour = eColor, edge_width = eWidth, edge_alpha = 0.5,
#                     arrow = arrow(angle = 30, length = unit(my_arrowSize, "inches"), ends = "last", type = "closed")) +
#     geom_node_point(fill = vColor, colour = "#000000",
#                     size =mySize , stroke = 0.8, shape = 21) +theme_graph()+
#     geom_node_text(aes(filter = mySize == 13, label = V(mainGraph)$name),family = "serif")
#   
# }
# 

plot_webTraffic=function(mainGraph,anoEvents,edge.list)
{
  
  my_colorList=c("blue","sienna1")
  
  colorNodes=as.character(anoEvents$anomaly[nrow(anoEvents)])
  #colorNodes=unique(as.character(anoEvents$anomaly))
  
  list1=which(V(mainGraph)$name %in% colorNodes)
  list2=which(!(V(mainGraph)$name %in% colorNodes))
  
  
  
  #subVertices=V(mainGraph)[V(mainGraph)[c(list1,list2[1:40])]]
  #subG=induced_subgraph(mainGraph,subVertices)
  
  mySize=rep(2,length(V(mainGraph)))
  mySize[which(V(mainGraph)$name %in% colorNodes)]=13
  
  vColor=rep("cyan4",length(V(mainGraph)))
  vColor[which(V(mainGraph)$name %in% colorNodes)]=my_colorList[2]
  
  eColor=rep("grey66",length(E(mainGraph)))
  eWidth=rep(0.2,length(E(mainGraph)))
  #eWidth[get.edge.ids(mainGraph,colorNodes)]=3
  
    my_arrowSize=rep(0.15,length(E(mainGraph)))
  

  p<-ggraph(mainGraph,layout = "lgl") +
    geom_edge_link0(edge_colour = eColor, edge_width = eWidth, edge_alpha = 0.5,
                    arrow = arrow(angle = 30, length = unit(my_arrowSize, "inches"), ends = "last", type = "closed")) +
    geom_node_point(fill = vColor, colour = "#000000",
                    size =mySize , stroke = 0.8, shape = 21) +theme_graph()+
    geom_node_text(aes(filter = mySize == 13, label = V(mainGraph)$name),family = "serif")

  
}


plot_stock=function(mainGraph,anoEvents,edge.list)
{
  my_colorList=c("blue","darkgreen","red","goldenrod4","gold","darkorchid1","palegreen1","sienna1")
  colorNodes=unique(as.character((anoEvents$anomaly)))   
  

  list1=which(V(mainGraph)$name %in% colorNodes)
  list2=which(!(V(mainGraph)$name %in% colorNodes))
  
  
  
  #subVertices=V(mainGraph)[V(mainGraph)[c(list1,list2[1:40])]]
  #subG=induced_subgraph(mainGraph,subVertices)

  mySize=rep(2,length(V(mainGraph)))
  mySize[which(V(mainGraph)$name %in% colorNodes)]=13
  
  vColor=rep("cyan4",length(V(mainGraph)))
  vColor[which(V(mainGraph)$name %in% colorNodes)]=my_colorList
  
  
  getEdges=c(colorNodes[1],rep(colorNodes[2:4],each=2), colorNodes[length(colorNodes)]) 
  
  eColor=rep("grey66",length(E(mainGraph)))
  eColor[get.edge.ids(mainGraph,getEdges)]=my_colorList[1:4]
  
  eWidth=rep(0.2,length(E(mainGraph)))
  eWidth[get.edge.ids(mainGraph,getEdges)]=3
  
  
    my_arrowSize=rep(0.15,length(E(mainGraph)))
    my_arrowSize[get.edge.ids(mainGraph,getEdges)]=0.3

  p<-ggraph(mainGraph,layout = "lgl") +
    geom_edge_link0(edge_colour = eColor, edge_width = eWidth, edge_alpha = 0.5,
                    arrow = arrow(angle = 30, length = unit(my_arrowSize, "inches"), ends = "last", type = "closed")) +
    geom_node_point(fill = vColor, colour = "#000000",
                    size =mySize , stroke = 0.8, shape = 21) +theme_graph()+
    geom_node_text(aes(filter = mySize == 13, label = V(mainGraph)$name),family = "serif")

  
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



mainCall=function(dl,dta.processed,anoEvents,window)
{
  
  symbolicRep=read.csv(file.path(p,dl,"Output non-overlap","Cluster labels EV weighted",paste0("symbolicRep_","win",window,".csv")),header=T)

   out=graphConstruct(dta.processed,symbolicRep)
  mainGraph=out[[1]]
  edge.list=out[[2]]

    #plot_for_london(mainGraph,anoEvents,edge.list)
    plot_webTraffic(mainGraph,anoEvents,edge.list)
    #plot_stock(mainGraph,anoEvents,edge.list)
    
  
}

p="/home/ujjwal/Kakuli/Graph Representation of TS"
dataList=c("London Data","Ausgrid Data","Stock market","Web traffic")[4]
processedData=c("londonData_logTransformed.csv","ausgridData_logTransformed.csv","stockMarketData_logTransformed.csv",
                "webTraffic_logTransformed.csv")[4]
anomalyEvents=read.csv("/home/ujjwal/Kakuli/Graph Representation of TS/output_rareEvents.csv",header=T)

#window=c(12,18, 24,30,36,42,48,54,60)[3]

#window=c(5,10,15,20,25,30,35,40,45,50)[2]
window=c(7,14,21,28,35,42,49)[1]


for (dl in 1:length(dataList))
{
  dta.processed=read.csv(file.path(p,dataList[dl],"Output non-overlap",processedData[dl]),header=T)
  for (i in 1:length(window))
  {
  
    anoEvents=anomalyEvents[which(anomalyEvents$dataName==dataList[dl] &  anomalyEvents$window==window[i]),]
    out=mainCall(dataList[dl],dta.processed,anoEvents,window[i])
    
  }

  
}