###This file is to obtain the final RI values, reported in the paper in Table 2. Although other measures are computed but due to poor results,
####those are not reported in the paper. The results are saved in the local directory as: "tableResults_nonOverlap.csv"


library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(scales)
library(clusteval)
library(DescTools)


plotRI=function(resultTable)
{

  
  #The Entropy, NMI and FM shows poor result in the proposed method, hence only the RI is plotted
  
  resultTable=resultTable[!(resultTable$clustering=="weighted_ProposedLabels"),]
  
  #resultTable=resultTable[!(resultTable$clustering=="PAM" | resultTable$clustering=="FastGreedy" |resultTable$clustering=="KM"|resultTable$clustering=="Walktrap"  ),]
  resultTable=resultTable[(resultTable$clustering == "ProposedLabels") |
                            (resultTable$clustering == "EDweighted_ProposedLabels") |
                            (resultTable$clustering == "EV_ProposedLabels") |
                            resultTable$clustering=='EV_weighted_ProposedLabels',]
  
  
  
  plot_list=NULL
  
  count=0
  dataSets=unique(resultTable$data)
  dataFreq=unique(resultTable$dataFreq)
  for (i in 1:length(dataSets))
  {

      count=count+1
      
      subsetResult=resultTable[resultTable$data==dataSets[i],]
      subsetResult=subsetResult[,c(2,7,8)]
     
      
      
      subsetResult$rIndex=as.numeric(subsetResult$rIndex)
      subsetResult$window=as.numeric(subsetResult$window)
      subsetResult$clustering=as.character(subsetResult$clustering)
      mylevels=c("EDweighted_ProposedLabels",  "ProposedLabels" ,"EV_weighted_ProposedLabels", "EV_ProposedLabels"  )
      # meltData=melt(subsetResult,c(2,4))
      # 
      myLabels=c("Deg-WG","Deg-UWG","Ei-WG", "Ei-UWG") ###Check this everytime the code is changed. This is important to label correctly
      myColor=c("blue","red","deeppink4","darkgreen")
      # 
      # metricNames <- c("RI")
      # names(metricNames) <- unique(as.character(meltData$variable))
      # 
      p1=ggplot(data=subsetResult)+geom_line(data=subsetResult,aes(x=window,y=rIndex,
                                                           color=factor(clustering,levels=mylevels)),size=0.8)+xlab("Window")+theme_bw()+
        geom_point(data=subsetResult,aes(x=window,y=rIndex,
                          shape=factor(clustering, levels=mylevels),color=factor(clustering, levels=mylevels))
                   ,stroke=1,size=2)+
        scale_color_manual(values=myColor, labels=myLabels)+scale_shape_manual(values=c(5,6,7,8),labels=myLabels)+
        theme(plot.title = element_text(hjust = 0.5,color="black",size=15),axis.ticks.y = element_blank(),
              legend.title = element_blank(),legend.position="bottom",panel.grid.major.y =element_blank(),
              legend.text=element_text(size=10),axis.title.x = element_text(size=15,color="black"),
              axis.title.y= element_blank(),
              panel.grid.major =   element_blank() , axis.text= element_text(size=15,color="black"),
              panel.background = element_rect(color="black",size=1.5), 
              strip.background =element_rect(color="black",fill="lightgreen"),
              axis.ticks.x = element_blank(),plot.caption = element_text(hjust=0.5, size=rel(1.2)),
              strip.text = element_text(size = 12, colour = "black"),legend.key.width = unit(2,"cm"))+
        ggtitle(dataSets[i])+ scale_x_continuous(breaks=seq(7,max(subsetResult$window),7))
      
      
      plot_list[[count]]=p1
      
    
    
  }
  
  return (plot_list)
  
}




plotAllMeasures=function(resultTable)
{
  
  library(data.table)
#resultTable=resultTable[!(resultTable$clustering=="PAM" | resultTable$clustering=="FastGreedy" |resultTable$clustering=="KM"|resultTable$clustering=="Walktrap"  ),]
 resultTable=resultTable[(resultTable$clustering %like% "Proposed") |(resultTable$clustering %like% "weighted_Proposed") |
                             (resultTable$clustering %like% "SAX") | (resultTable$clustering %like% "EV"),]
  
  plot_list=NULL
  
  myColor=c("darkgreen","blue","deeppink4","red")
  
  count=0
  dataSets=unique(resultTable$data)
  dataFreq=unique(resultTable$dataFreq)
  for (i in 1:length(dataSets))
  {

      count=count+1
      
      subsetResult=resultTable[resultTable$data==dataSets[i] ,]
      subsetResult=subsetResult[,c(1,2,6,7,8,9)]
      meltData=melt(subsetResult,c(3,4,5,6))
      
     
      meltData$window=as.numeric(meltData$window)
      meltData$value=as.numeric(meltData$value)
      
      #myLabels=levels(factor(meltData$clustering))
      myLabels=c("EV","UWG","SAX","WG")
      metricNames <- c("En","RI")
      names(metricNames) <- unique(as.character(meltData$variable))
      
      p1=ggplot(data=meltData)+geom_line(data=meltData,aes(x=window,y=value,
                                  color=factor(clustering)))+xlab("Window")+ylab("Values")+theme_bw()+
        geom_point(data=meltData,aes(x=window,y=value,
                                     shape=factor(clustering),color=factor(clustering)),stroke=1)+
        scale_color_manual(values=myColor,labels=myLabels)+scale_shape_manual(values=c(5,6,7,8),labels=myLabels)+
        theme(plot.title = element_text(hjust = 0.5,color="black",size=15),axis.ticks.y = element_blank(),
              legend.title = element_blank(),legend.position="bottom",panel.grid.major.y =element_blank(),
              legend.text=element_text(size=10),axis.title.x = element_text(size=15,color="black"),
              axis.title.y= element_text(size=15,color="black"),
              panel.grid.major =   element_blank() , axis.text= element_text(size=15,color="black"),
              panel.background = element_rect(color="black",size=1.5), 
              strip.background =element_rect(color="black",fill="lightgreen"),
              axis.ticks.x = element_blank(),plot.caption = element_text(hjust=0.5, size=rel(1.2)),
              strip.text = element_text(size = 12, colour = "black"),legend.key.width = unit(2,"cm"))+
        
        facet_wrap(scales="free_y",.~variable  ,labeller=labeller(variable=metricNames) ) +
        ggtitle(dataSets[i])
        #facet_grid(rows = vars(dataFreq), cols = vars(variable),space="free_x")
      
      
      plot_list[[count]]=p1
      
    
    
  }
  
  return (plot_list)
  
}


clustValidity=function(tp,tn,fp,fn,numSamples)
{
  numSamples=(numSamples*(numSamples-1))/2
  acc= round((tp+tn)/numSamples,2)
  pprecision=round(tp/(tp+fp),2)
  jc=round(tp/(tp+fp+fn),2)
  
  rIndex=round((tp+tn)/(tp+tn+fp+fn),2)
  return (list(acc,jc,pprecision,rIndex))
}



fp_tn=function(trueLabel,predLabel)
{
  pairedSamples=NULL
  start=Sys.time()
  trueLabel=lapply(unique(trueLabel), function(x) which(trueLabel%in% x))
  
  temp=NULL
  i=1
  j=i+1
  #trueLabel=trueLabel[1:2]
  temp<<-lapply(trueLabel,function(x) {
    while(i!=length(trueLabel) && j<=length(trueLabel))
    {
      temp<<-data.frame(rbind(temp, expand.grid(x,trueLabel[[j]])),stringsAsFactors = F)
      j<<-j+1
      temp
    }
    i<<-i+1
    j<<-i+1
  })
  
  fpCount<-0
  fpCount<<-apply(temp,1,function(f)
  {
    if(predLabel[f[1]]==predLabel[f[2]])
      fpCount<<-fpCount+1
  })
  
  tnCount<-nrow(temp)-fpCount
  
  
  # fnCount<<-apply(temp,1,function(f)
  # {
  #   if(predLabel[f[1]]!=predLabel[f[2]])
  #     fnCount<<-fnCount+1
  # })
  # 
  
  return (list(fpCount,tnCount))
  
}



tp_fn=function(trueLabel,predLabel)
{
  pairedSamples=NULL
  start=Sys.time()
  trueLabel=lapply(unique(trueLabel), function(x) which(trueLabel%in% x))
  
  temp=NULL
  temp<<-lapply(trueLabel,function(x) {
    for(i in 1:(length(x)-1))
    {
      temp<<-data.frame(rbind(temp, expand.grid(x[i],x[(i+1):length(x)])),stringsAsFactors = F)
      temp
    }
  })
  
  tpCount<-0
  tpCount<<-apply(temp,1,function(f)
  {
    if(predLabel[f[1]]==predLabel[f[2]])
      tpCount<<-tpCount+1
  })
  
  fnCount=nrow(temp)-tpCount
  
  return (list(tpCount,fnCount))
}


entrpy=function(trueLabel,predLabel)
{
  length(unique(trueLabel))
  
  entrpyList=unlist(lapply(unique(predLabel), function(x) 
  {
    pos= which(predLabel%in% x)
    pi=table(trueLabel[pos])/length(pos)
    pi=pi[pi>0]
    logPi=log(pi,2)
    #sum(a*b)*(-1/log(length(unique(trueLabel))))
    -(sum(pi*logPi))
  }))
  
  x=table(predLabel)/length(predLabel)
  avgEntropy=sum(x*entrpyList)
  
  return (round(avgEntropy,2))
  
  
}


ffmeasure=function(trueLabel,predLabel)
{
  #mergeLabels=cbind(trueLabel,predLabel)
  pprecision=NULL
  pprecision<-lapply(unique(predLabel), function(x) 
  {
    pos= which(predLabel%in% x)
    pprecision<<-rbind(pprecision,data.frame(cbind(x, unique(trueLabel[pos]), (table(trueLabel[pos])/sum(predLabel[pos])))))
  })
  pprecision=pprecision[[length(pprecision)]]
  colnames(pprecision)=c("clusters","classes","vals")
  
  
  temp=NULL
  
  temp=lapply(unique(predLabel),function(x)
  {
    pos= which(predLabel%in% x)
    val=as.data.frame(table(trueLabel[pos]))
    temp<<-rbind(temp,val)
  })
  temp=temp[[length(temp)]]
  colnames(temp)=c("trueLabel","Freq")
  temp=aggregate(Freq~trueLabel,temp,sum)
  
  rrecall=NULL
  rrecall=lapply(unique(predLabel), function(x) 
  {
    pos= which(predLabel%in% x)
    temp1=as.data.frame(table(trueLabel[pos]))
    rrecall<<-rbind(rrecall,data.frame(cbind(x, unique(trueLabel[pos]), (table(trueLabel[pos])/temp$Freq[temp1$Var1]))))
  })
  
  rrecall=rrecall[[length(rrecall)]]
  
  colnames(rrecall)=c("clusters","classes","vals")
  pr=data.frame(merge(pprecision,rrecall,by=c("clusters","classes")),stringsAsFactors = F)
  x= 2* (pr$vals.x * pr$vals.y)
  y=(pr$vals.x+pr$vals.y)
  pr$fMeasure=(x/y)
  
  pr=aggregate(fMeasure~classes,pr,max)
  colnames(pr)=c("class","val")
  
  #fMeasure.weights=  (table(predLabel))/table(trueLabel)
  #fMeasure=sum(fMeasure.weights * pr$val)
  
  fMeasure=sum(pr$val*temp$Freq)/sum(temp$Freq)
  
  return (round(fMeasure,2))
}



qualityMeasures=function(trueLabel,predLabel)
{
  
  ffmeasure=ffmeasure(trueLabel,predLabel)
  entrpy=entrpy(trueLabel,predLabel)
  
  out=tp_fn(trueLabel,predLabel)
  tp=out[[1]]
  fn=out[[2]]
  
  out=fp_tn(trueLabel,predLabel)
  fp=out[[1]]
  tn=out[[2]]
  
  #predLabel=pamx(dist_mat,n_clust)
  a=clustValidity(tp,tn,fp,fn,length(trueLabel))
  acc=a[[1]]
  
  pprecision=a[2]
  pprecision=round(pprecision[[1]],2)
  #pprecision=MLmetrics::Precision(trueLabel,predLabel)
  jcc=a[[3]]
  
  #rIndex=a[[4]]
  rIndex=round(cluster_similarity(predLabel, trueLabel,similarity="rand"),2)
  
  p=funtimes::purity(trueLabel,predLabel)$pur
  nmi=aricode::NMI(trueLabel,predLabel)
  nmi=round(nmi,2)
  #clust_accuracy=rbind(clust_accuracy,data.frame(cbind(jcc,ffmeasure,entrpy,rIndex,"fileName"=basename(fileList[files]),"distName"=basename(dirList[dirs])),stringsAsFactors = F))
  #resultTable=rbind(resultTable,data.frame(cbind(ffmeasure, entrpy, rIndex),stringsAsFactors = F))
  
  return (list(entrpy,rIndex,p,nmi,ffmeasure)) 
}




getTableResults=function()
{
  tableResults=NULL

  
  for (t in 1:length(trueLabel.files))
  {

      trueLabel=read.csv(file.path(getwd(),dataName[t],trueLabel.files[t]),header=T)$trueLabel
      
      p=c("Cluster labels Ei weighted","Cluster labels Ei unweighted","Cluster labels Deg unweighted","Cluster labels Deg weighted")
     
      w_uw=c("W","UW","EiV-UW","Ei")
      
      
      for (j in 1:length(p))
      {
        listFiles=  list.files(file.path(getwd(),dataName[t],"Output non-overlap",p[j]),pattern="Labels_")
        
        featureMat= list.files(file.path(getwd(),dataName[t],"Output non-overlap",p[j]),pattern="featureSpace_")

        for (z in 1:length(listFiles))
        {
          predLabel=read.csv(file.path(getwd(),dataName[t],"Output non-overlap",p[j],listFiles[z]),header=T)
          dataToCluster=read.csv(file.path(getwd(),dataName[t],"Output non-overlap",p[j],featureMat[z]),header=T)
          #print (basename(file.path(getwd(),dataName[t],"Output non-overlap",p[j],listFiles[z])))
          print (listFiles[z])
          #clustType=strsplit(listFiles[z],"Labels")[[1]][1]
          bsname=  (basename(file.path(getwd(),dataName[t],"Output non-overlap",p[j],listFiles[z])))
          bsname=sub("_win.*", "", bsname)
        
          
          
          window=strsplit(listFiles[z],"win")[[1]][2]
          window=gsub(window,pattern=".csv",replacement="")
          library(stringr)
          
          out=qualityMeasures(unlist(trueLabel,use.names = F),unlist(predLabel,use.names = F))

          tableResults=rbind(tableResults,data.frame(cbind("entropy"=out[[1]], "rIndex"=out[[2]],"purity"=out[[3]],"nmi"=out[[4]],"Fmeasure"=out[[5]],"data"=dataName[t],
                      "clustering"=bsname,"window"=window, "W_UW"=w_uw[j]),stringsAsFactors = F))
          
        }
        
      }
    
  }
  
  
  return (tableResults)
}


dataName=c("London Data","Ausgrid Data","Stock market","Web traffic")
trueLabel.files=c("londonData_forGraph.csv", "AusgridData_forGraph.csv","stockMarketData_forGraph.csv","webTrafficData_forGraph.csv")
setwd("/home/kakuli/Graph Representation of TS")
#tableResults=getTableResults()
#plot_list=plotAllMeasures(tableResults)

tableResults=read.csv("/home/kakuli/Graph Representation of TS/tableResults_nonOverlap.csv",header=T)
plot_list=plotRI(tableResults)

