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

# Построение графика
#font_import()
#loadfonts("win")

ggplot(graph,aes(x=1,y=x))+geom_violin()+
  geom_boxplot(width=0.1,fill="black")+
  theme(axis.title.x=element_blank())+
  scale_y_log10(breaks=c(1,60,3600,3600*24),labels=c("1 cек.","1 мин.", "1 час", "1 сутки"))+
  stat_summary(fun.y=median,geom="point",fill="white",shape=21,size=2.5)+
  scale_x_continuous(breaks=NULL)+
  labs(x="",y="Интервал до обращения")+theme_bw()+
  theme(text=element_text(size=18,family="Times New Roman"))
  
ggsave("Pic1.png",width=4.5,height=4.5,dpi=300)
