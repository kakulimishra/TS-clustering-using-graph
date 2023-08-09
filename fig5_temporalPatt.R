###The file is used to plot for the patterns shown in Fig. 5 in the paper.  Every dataset has 2 functions ,for example, "webTrafficData" and "getPlots_webTraffic"
###The first function prepares the data for plots and the second function is to get the plots.

###Run the file in debug mode


library(stringr)
library(TSclust)
library(clusteval)
library(qualV)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(plyr)
library(reshape2)


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

getPlots_webTraffic=function(df,locs,date_seq)
{
  
  myColors=rep(c("tomato","blue","red","deepskyblue3","goldenrod4","deeppink4","darkslategray"),100)
  
  ts.all=NULL
  df=df[101:105,1:798] ###make it divisible by thw window length,i.e, 42 for the web traffic data
  
  for (i in 1:nrow(df))
  {
    temp=nonOverlap(df[i,],42)
    temp=data.frame(temp)
    temp$id=0
    temp$id[locs]=1
    #temp$id[locs]=1
    tsData<-NULL
    
    apply(temp,1,function(f)
    {
      ti=(unlist(f[1:(length(f)-1)]))
      tsData<<-rbind(tsData,cbind(ti,as.numeric(rep(f[length(f)], length(ti)))))
    })
    tsData=data.frame(tsData)
    tsData$date=date_seq[1:nrow(tsData)]
    tsData=tsData[580:800,]
    tsData$ts=rep(paste("T",i),nrow(tsData))
    
    
    ts.all=rbind(ts.all,tsData)
  }
  
  my_colors=c("black","darkorchid1")
  my_size=c(0.3,1)
  
  ts.all$ti=as.numeric(as.character(ts.all$ti))
  
  ggplot(ts.all)+ geom_line(aes(x=date,y=ti,color=factor(V2),group=ts,size=factor(V2)))+ theme_bw()+
    scale_color_manual(values=my_colors)+scale_size_manual(values=my_size)+
    theme(plot.title = element_text(hjust = 0.5),legend.title = element_blank(),
          legend.position = 'none',axis.title.y = element_blank(),axis.ticks.y = element_blank(),
          axis.text = element_text(colour='black'), 
          panel.grid.major.y =element_blank(),panel.grid.major.x = element_blank() ,
          axis.title = element_blank(), 
          axis.text.x = element_text(angle = 45, vjust = 0.5, color='black'), axis.text.y = element_text(color='black'))+
    scale_x_date(breaks=seq(as.Date(date_seq[1]), as.Date(date_seq[length(date_seq)]),by=20), date_labels= "%b-%d\n%H:%M")
  
  
}



webTrafficData=function(dl,fileToRead)
{
  trainData=read.csv(file.path(p,dl,fileToRead),header=T)
  dataNm=strsplit(dl," ")[[1]][1]
  print (dataNm)
  
  dta=trainData[,-c(1,2)]
  
  listTemporal=list.files(file.path(p,dl,"Output non-overlap","Temporal patterns EV weighted"))
  
  listSymbols=list.files(file.path(p,dl,"Output non-overlap","Cluster labels Ei weighted"),pattern="symbolicRep")
  
  listPredLabel=list.files(file.path(p,dl,"Output non-overlap","Cluster labels Ei weighted"),pattern="ProposedLabels")
  
  date_seq=(seq(from=as.Date("2015-07-1"),to=as.Date("2017-09-20"),by=1))  ##web Traffic
  
  #for the saved copies of Temporal patterns detected in different window lengths.
  for (i in 5:5 )  ##window 42 for the web traffic dataset
    ###Get plots for the temporal patterns windows which shows the highest average support. 
  {
    
    temporalPatt=read.csv(file.path(p,dl,"Output non-overlap","Temporal patterns EV weighted",listTemporal[i]),header=T)
    window=strsplit(basename(file.path(p,dl,"Output non-overlap","Temporal patterns EV weighted",listTemporal[i])),"win")[[1]][2]
    window=as.numeric(gsub(window,pattern=".csv",replacement=""))
    
    listSymbols1=(sapply(listSymbols,function(x) 
    {
      x=strsplit(x,"win")[[1]][2]
      gsub(x,pattern=".csv",replacement="")
    }))
    
    
    symbolicRep=read.csv(file.path(p,dl,"Output non-overlap","Cluster labels Ei weighted",listSymbols[which(listSymbols1==window)]),header=T)
    
    plot_list=NULL
    time_of_occurence=NULL
    
    listPredLabel1=(sapply(listPredLabel,function(x) 
    {
      x=strsplit(x,"win")[[1]][2]
      gsub(x,pattern=".csv",replacement="")
    }))
    
    predLabel=read.csv(file.path(p,dl,"Output non-overlap","Cluster labels Ei weighted",listPredLabel[which(listPredLabel1==window)]),header=T)
    
    
    for (l in 1:1)  ##plot only the first cluster of ausgrid dataset
    {
      
      #extract the symbolic representation and the raw data of the required cluster
      temp.sym=data.frame(symbolicRep[which(predLabel==temporalPatt$cluster[l]),],row.names = NULL)
      
      temp.dta=data.frame(dta[which(predLabel==temporalPatt$cluster[l]),],row.names = NULL)
      
      pat1=temporalPatt[l,4:ncol(temporalPatt)]
      pat1=as.character(pat1[!is.na(pat1)])
      #pat1=unique(pat1)
      
      
      ####Time of occurrence of the temporal patterns.
      time.llcs=apply(temp.sym,1,function(f)
      {
        pat2=as.character(unlist(f))
        llcs=LCS(pat1,pat2)
        llcs$vb
      })
      
      time_of_occurence[[l]]=time.llcs
      
      x=data.frame(time_of_occurence[[1]])
      
      locs=c(15, 18)  ##for the web traffic dataset
      
      ts_for_plots=which(apply(x,2,function(f) all(f %in% locs ))==TRUE)
      
      getPlots_webTraffic(temp.dta[ts_for_plots,],locs,date_seq)
      
    }
  }
  
  
}


getPlots_stock=function(df,locs,date_seq)
{   
  
  ts.all=NULL
  df=df[-c(3,4,9,13),1:1015]  ##Make the number of columns a multiple of window length, 35

  for (i in 4:nrow(df))
  {
    temp=nonOverlap(df[i,],35)
    temp=data.frame(temp)
    temp$id=0
    temp$id[locs[1:3]]=1
    #temp$id[locs]=1
    tsData<-NULL
    
    apply(temp,1,function(f)
    {
      ti=(unlist(f[1:(length(f)-1)]))
      tsData<<-rbind(tsData,cbind(ti,as.numeric(rep(f[length(f)], length(ti)))))
    })
    tsData=data.frame(tsData)
    tsData$date=date_seq[1:nrow(tsData)]
    tsData=tsData[600:nrow(tsData),]
    tsData$ts=rep(paste("T",i),nrow(tsData))
    
    
    ts.all=rbind(ts.all,tsData)
  }
  
  my_colors=c("black","darkorchid1")
  my_size=c(0.5,1.5)
  
  ts.all$ti=as.numeric(as.character(ts.all$ti))
  
  ggplot(ts.all)+ geom_line(aes(x=date,y=ti,color=factor(V2),group=ts, size=factor(V2)))+ theme_bw()+
    scale_color_manual(values=my_colors)+scale_size_manual(values=my_size)+
    theme(plot.title = element_text(hjust = 0.5),legend.title = element_blank(),
          legend.position = 'none',axis.title.y = element_blank(),axis.ticks.y = element_blank(),
          axis.text = element_text(colour='black'), 
          panel.grid.major.y =element_blank(),panel.grid.major.x = element_blank() ,
          axis.title = element_blank(), 
          axis.text.x = element_text( vjust = 0.5,color = "black"),
          axis.text.y = element_text(color='black'))+
    scale_x_date(breaks=seq(as.Date(date_seq[1]), as.Date(date_seq[length(date_seq)]),by=20), date_labels= "%b\n %d")+
    scale_y_continuous(breaks=seq(min(ts.all$ti),32900,1000),trans=log2_trans())
  
}


stockMarketData=function(dl,fileToRead)
{
  
  trainData=read.csv(file.path(p,dl,fileToRead),header=T)
  dataNm=strsplit(dl," ")[[1]][1]
  print (dataNm)
  
  dta=trainData[,-c(1,2,3)]
  
  listTemporal=list.files(file.path(p,dl,"Output non-overlap","Temporal patterns EV weighted"))
  
  listSymbols=list.files(file.path(p,dl,"Output non-overlap","Cluster labels Ei weighted"),pattern="symbolicRep")
  
  listPredLabel=list.files(file.path(p,dl,"Output non-overlap","Cluster labels Ei weighted"),pattern="ProposedLabels")
  
  date_seq=(seq(from=as.Date("2013-12-10"),to=as.Date("2017-12-10"),by=1))  ##Stock
  date_seq <- date_seq["Date"=!weekdays(date_seq) %in% c('Saturday','Sunday')]
  
  #for the saved copies of Temporal patterns detected in different window lengths.
  for (i in 6:6 )  ##window 35stock market for the ausgrid dataset
    ###Get plots for the temporal patterns windows which shows the highest average support. 
  {
    
    temporalPatt=read.csv(file.path(p,dl,"Output non-overlap","Temporal patterns EV weighted",listTemporal[i]),header=T)
    window=strsplit(basename(file.path(p,dl,"Output non-overlap","Temporal patterns EV weighted",listTemporal[i])),"win")[[1]][2]
    window=as.numeric(gsub(window,pattern=".csv",replacement=""))
    
    listSymbols1=(sapply(listSymbols,function(x) 
    {
      x=strsplit(x,"win")[[1]][2]
      gsub(x,pattern=".csv",replacement="")
    }))
    
    
    symbolicRep=read.csv(file.path(p,dl,"Output non-overlap","Cluster labels Ei weighted",listSymbols[which(listSymbols1==window)]),header=T)

    plot_list=NULL
    time_of_occurence=NULL
    
    listPredLabel1=(sapply(listPredLabel,function(x) 
    {
      x=strsplit(x,"win")[[1]][2]
      gsub(x,pattern=".csv",replacement="")
    }))
    
    predLabel=read.csv(file.path(p,dl,"Output non-overlap","Cluster labels Ei weighted",listPredLabel[which(listPredLabel1==window)]),header=T)
    
    
    for (l in 1:1)  ##plot only the first cluster of ausgrid dataset
    {
      
      #extract the symbolic representation and the raw data of the required cluster
      temp.sym=data.frame(symbolicRep[which(predLabel==temporalPatt$cluster[l]),],row.names = NULL)
      
      temp.dta=data.frame(dta[which(predLabel==temporalPatt$cluster[l]),],row.names = NULL)
      
      pat1=temporalPatt[l,4:ncol(temporalPatt)]
      pat1=as.character(pat1[!is.na(pat1)])
      #pat1=unique(pat1)
      
      
      ####Time of occurrence of the temporal patterns.
      time.llcs=apply(temp.sym,1,function(f)
      {
        pat2=as.character(unlist(f))
        llcs=LCS(pat1,pat2)
        llcs$vb
      })
      
      time_of_occurence[[l]]=time.llcs
      
      x=data.frame(time_of_occurence[[1]])
      
      locs=c(27, 29)  ##for the ausgrid dataset
      
      ts_for_plots=which(apply(x,2,function(f) all(f %in% locs ))==TRUE)
      
      getPlots_stock(temp.dta[ts_for_plots,],locs,date_seq)
      
    }
  }
  
  
}





getPlots_ausgrid=function(df,locs,date_seq)
{
  myColors=rep(c("tomato","blue","red","deepskyblue3","goldenrod4","deeppink4","darkslategray"),100)
  
  ts.all=NULL
  df=df[1:5,]
  for (i in 1:nrow(df))
  {
    temp=nonOverlap(df[i,],48)
    temp=data.frame(temp)
    temp$id=0
    temp$id[locs[1:3]]=1
    tsData<-NULL
    
    apply(temp,1,function(f)
    {
      ti=(unlist(f[1:(length(f)-1)]))
      tsData<<-rbind(tsData,cbind(ti,as.numeric(rep(f[length(f)], length(ti)))))
    })
    tsData=data.frame(tsData)
    tsData=tsData[2600:3050,]
    tsData$date=as.POSIXct(date_seq[2600:3050])
    tsData$ts=rep(paste("T",i),nrow(tsData))
    
    
    ts.all=rbind(ts.all,tsData)
  }
  
  my_colors=c("black","darkorchid1")
  my_size=c(0.3,1)
  
  ts.all$ti=as.numeric(as.character(ts.all$ti))

  ggplot(ts.all)+ geom_line(aes(x=date,y=ti,color=factor(V2),group=ts,size=factor(V2)))+ theme_bw()+
    scale_color_manual(values=my_colors)+scale_size_manual(values=my_size)+
    theme(plot.title = element_text(hjust = 0.5),legend.title = element_blank(),
          legend.position = 'none',axis.title.y = element_blank(),axis.ticks.y = element_blank(),
          axis.text = element_text(colour='black'), 
          panel.grid.major.y =element_blank(),panel.grid.major.x = element_blank() ,
          axis.title = element_blank(), 
          axis.text.x = element_text( vjust = 0.5))+
    scale_x_datetime(breaks=seq(as.POSIXct(date_seq[1]), as.POSIXct(date_seq[length(date_seq)]),by=60*800), date_labels= "%b-%d\n%H:%M")
  
    
    
    #scale_x_continuous(breaks=seq(0,length(unique(ts.all$hr)),48))+
    
    

}


ausgridData=function(dl,fileToRead)
{
  trainData=read.csv(file.path(p,dl,fileToRead),header=T)
  dataNm=strsplit(dl," ")[[1]][1]
  print (dataNm)
  
  if (dl=="Web traffic")
    dta=trainData[,-c(1,2)]
  
  
  dta=trainData[,-c(1,2,3)]
  
  listTemporal=list.files(file.path(p,dl,"Output non-overlap","Temporal patterns EV weighted"))
  
  listSymbols=list.files(file.path(p,dl,"Output non-overlap","Cluster labels Ei weighted"),pattern="symbolicRep")
  
  listPredLabel=list.files(file.path(p,dl,"Output non-overlap","Cluster labels Ei weighted"),pattern="ProposedLabels")
  
  
  #for the saved copies of Temporal patterns detected in different window lengths.
  for (i in 7:  length(listTemporal))  ##window 48 for the ausgrid dataset
    ###Get plots for the temporal patterns windows which shows the highest average support. 
  {
    
    temporalPatt=read.csv(file.path(p,dl,"Output non-overlap","Temporal patterns EV weighted",listTemporal[i]),header=T)
    window=strsplit(basename(file.path(p,dl,"Output non-overlap","Temporal patterns EV weighted",listTemporal[i])),"win")[[1]][2]
    window=as.numeric(gsub(window,pattern=".csv",replacement=""))
    
    listSymbols1=(sapply(listSymbols,function(x) 
    {
      x=strsplit(x,"win")[[1]][2]
      gsub(x,pattern=".csv",replacement="")
    }))
    
    
    symbolicRep=read.csv(file.path(p,dl,"Output non-overlap","Cluster labels Ei weighted",listSymbols[which(listSymbols1==window)]),header=T)
    
    plot_list=NULL
    time_of_occurence=NULL
    
    listPredLabel1=(sapply(listPredLabel,function(x) 
    {
      x=strsplit(x,"win")[[1]][2]
      gsub(x,pattern=".csv",replacement="")
    }))
    
    predLabel=read.csv(file.path(p,dl,"Output non-overlap","Cluster labels Ei weighted",listPredLabel[which(listPredLabel1==window)]),header=T)
    
    date_seq=trainData$dates[which(predLabel==1)[1]]
    
    genDate=(seq(from=as.POSIXct(date_seq),to=as.POSIXct(as.Date(date_seq)+182),by=(60*30)))  ##Ausgrid
    
    
    for (l in 1:1)  ##plot only the first cluster of ausgrid dataset)
    {
      
      #extract the symbolic representation and the raw data of the required cluster
      temp.sym=data.frame(symbolicRep[which(predLabel==temporalPatt$cluster[l]),],row.names = NULL)
      
      temp.dta=data.frame(dta[which(predLabel==temporalPatt$cluster[l]),],row.names = NULL)
      
      pat1=temporalPatt[l,4:ncol(temporalPatt)]
      pat1=as.character(pat1[!is.na(pat1)])
      pat1=unique(pat1)

      
      ####Time of occurrence of the temporal patterns.
      time.llcs=apply(temp.sym,1,function(f)
      {
        pat2=as.character(unlist(f))
        llcs=LCS(pat1,pat2)
        llcs$vb
      })
      
      time_of_occurence[[l]]=time.llcs
      
      x=data.frame(time_of_occurence[[1]])
      
      locs=c(56,59,63,67,87,101,120,176,180)  ##for the ausgrid dataset
      
      ts_for_plots=which(apply(x,2,function(f) all(f %in% locs ))==TRUE)
      
      getPlots_ausgrid(temp.dta[ts_for_plots,],locs,genDate)
      
    }
  }
}


getPlots_london=function(df,locs,date_seq)
{
  
  myColors=rep(c("tomato","blue","red","deepskyblue3","goldenrod4","deeppink4","darkslategray"),100)
  
  ts.all=NULL
  df=df[,1:8784]  ##make it divisible by 12(window length)
  for (i in 1:nrow(df))
  {
    temp=nonOverlap(df[i,],12)  #window length 12 for the london dataset
    temp=data.frame(temp)
    temp$id=0
    temp$id[locs]=1
    #temp$id[locs]=1
    tsData<-NULL
    
    apply(temp,1,function(f)
    {
      ti=(unlist(f[1:(length(f)-1)]))
      tsData<<-rbind(tsData,cbind(ti,as.numeric(rep(f[length(f)], length(ti)))))
    })
    tsData=data.frame(tsData)
    tsData=tsData[4500:4900,]
    tsData$date=as.POSIXct(date_seq[4500:4900])
    tsData$ts=rep(paste("T",i),nrow(tsData))
    
    
    ts.all=rbind(ts.all,tsData)
  }
  
  my_colors=c("black","darkorchid1")
  my_size=c(0.3,1)
  
  ts.all$ti=as.numeric(as.character(ts.all$ti))
  
  ggplot(ts.all)+ geom_line(aes(x=date,y=ti,color=factor(V2),group=ts,size=factor(V2)))+ theme_bw()+
    scale_color_manual(values=my_colors)+scale_size_manual(values=my_size)+
    theme(plot.title = element_text(hjust = 0.5),legend.title = element_blank(),
          legend.position = 'none',axis.title.y = element_blank(),axis.ticks.y = element_blank(),
          axis.text = element_text(colour='black'), 
          panel.grid.major.y =element_blank(),panel.grid.major.x = element_blank() ,
          axis.title = element_blank(), 
          axis.text.x = element_text( vjust = 0.5))+
    scale_x_datetime(breaks=seq(as.POSIXct(date_seq[1]), as.POSIXct(date_seq[length(date_seq)]),by=60*800), date_labels= "%b-%d\n%H:%M")
  
}


londonData=function(dl,fileToRead)
{
  
  trainData=read.csv(file.path(p,dl,fileToRead),header=T)
  dataNm=strsplit(dl," ")[[1]][1]
  print (dataNm)
  dta=trainData[,-c(1,2,3)]
  
  listTemporal=list.files(file.path(p,dl,"Output non-overlap","Temporal patterns EV weighted"))
  
  listSymbols=list.files(file.path(p,dl,"Output non-overlap","Cluster labels Ei weighted"),pattern="symbolicRep")
  
  listPredLabel=list.files(file.path(p,dl,"Output non-overlap","Cluster labels Ei weighted"),pattern="ProposedLabels")
  
  
  #for the saved copies of Temporal patterns detected in different window lengths.
  for (i in 1:  1)  ##window 12 for the london dataset
    ###Get plots for the temporal patterns windows which shows the highest average support. 
  {
    
    temporalPatt=read.csv(file.path(p,dl,"Output non-overlap","Temporal patterns EV weighted",listTemporal[i]),header=T)
    window=strsplit(basename(file.path(p,dl,"Output non-overlap","Temporal patterns EV weighted",listTemporal[i])),"win")[[1]][2]
    window=as.numeric(gsub(window,pattern=".csv",replacement=""))
    
    listSymbols1=(sapply(listSymbols,function(x) 
    {
      x=strsplit(x,"win")[[1]][2]
      gsub(x,pattern=".csv",replacement="")
    }))
    
    
    symbolicRep=read.csv(file.path(p,dl,"Output non-overlap","Cluster labels Ei weighted",listSymbols[which(listSymbols1==window)]),header=T)
    
    plot_list=NULL
    time_of_occurence=NULL
    
    listPredLabel1=(sapply(listPredLabel,function(x) 
    {
      x=strsplit(x,"win")[[1]][2]
      gsub(x,pattern=".csv",replacement="")
    }))
    
    predLabel=read.csv(file.path(p,dl,"Output non-overlap","Cluster labels Ei weighted",listPredLabel[which(listPredLabel1==window)]),header=T)
    
    date_seq=(seq(from=as.POSIXct("2012-03-01"),to=as.POSIXct("2012-08-31"),by=(60*30)))  ##London
    
    
    for (l in 1:1)  ##plot only the first cluster of ausgrid dataset)
    {
      
      #extract the symbolic representation and the raw data of the required cluster
      temp.sym=data.frame(symbolicRep[which(predLabel==temporalPatt$cluster[l]),],row.names = NULL)
      
      temp.dta=data.frame(dta[which(predLabel==temporalPatt$cluster[l]),],row.names = NULL)
      
      pat1=temporalPatt[l,4:ncol(temporalPatt)]
      pat1=as.character(pat1[!is.na(pat1)])
      pat1=unique(pat1)[2:10]
      
      
      ####Time of occurrence of the temporal patterns.
      time.llcs=apply(temp.sym,1,function(f)
      {
        pat2=as.character(unlist(f))
        llcs=LCS(pat1,pat2)
        llcs$vb
      })
      
      time_of_occurence[[l]]=time.llcs
      
      x=data.frame(time_of_occurence[[1]])
      
      #locs=c(384, 402, 404, 438, 458, 462, 466, 471, 543)  ##for the london dataset
      locs=c(384,398,404)
      
      ts_for_plots=which(apply(x[c(1,2,3),],2,function(f) all(f %in% locs ))==TRUE)
      
      getPlots_london(temp.dta[ts_for_plots,],locs,date_seq)
      
    }
  }
  
  
}


p="/home/kakuli/Graph Representation of TS"
dataList=c("London Data","Ausgrid Data","Stock market","Web traffic")[2]
fileToRead=c("londonData_forGraph.csv","AusgridData_forGraph.csv","stockMarketData_forGraph.csv","webTrafficData_forGraph.csv")[2]

for (j in 1:1)#length(dataList))
{
  
  #out=londonData(dataList[j],fileToRead[j] )
  out=ausgridData(dataList[j],fileToRead[j] )
  #out=stockMarketData(dataList[j],fileToRead[j] )
  #out=webTrafficData(dataList[j],fileToRead[j] )
}

