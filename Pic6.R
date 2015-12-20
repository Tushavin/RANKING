library(ggplot2)
library(scales)
library(readxl)
library(extrafont)

# Подготовка данных
xls.data<-read_excel("datafile.xlsx")
mydata<-xls.data[,c(1,3,9,15)]
mydata[is.na(mydata[,2]),2]<-0
mydata[is.na(mydata[,3]),3]<-0
mydata[,3]<-as.factor(mydata[,3])
max.idx<-nrow(mydata)
mydata<-mydata[order(mydata$created),]
lag<-as.numeric(diff(mydata$created))
graph<-data.frame(x=lag[1:(max.idx-2)],y=lag[2:(max.idx-1)], hour=as.numeric(format(mydata$created[2:(max.idx-1)],"%H")),
                  wd=as.numeric(format(mydata$created[2:(max.idx-1)],"%u")))
test<-subset(graph,wd>1 & wd<6 & hour>11 & hour<15)
UCL<-function(x,alpha=0.05) log(alpha/2)/log(1-1/x)-1
LCL<-function(x,alpha=0.05) log(1-alpha/2)/log(1-1/x)
# Построение графика
#font_import()
g<-ggplot(test,aes(x=x,y=y))+
  geom_point(col="grey")+
  scale_x_log10(breaks=c(1,60,3600,3600*24),
                labels=c("1 cек.","1 мин.",  "1 час", "1 сутки"))+
  scale_y_log10(breaks=c(1,60,3600,3600*24),
                labels=c("1 cек.","1 мин.", "1 час", "1 сутки"))+
  stat_density2d(col="black")+
  geom_hline(aes(yintercept=c(LCL(mean(test$x)),UCL(mean(test$x)))))+
  geom_vline(aes(xintercept=c(LCL(mean(test$x)),UCL(mean(test$x)))))+
  labs(x="Интервал до обращения",y="Интервал после обращения") 
g+theme_bw()+theme(text=element_text(size=18,family="Times New Roman"))

ggsave("Pic6.png",width=4.5,height=4.5,dpi=300)

