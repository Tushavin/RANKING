doe<-readRDS("doe.RDs")
WD<-getwd()
setwd("Pictures/Artcl03")
library(VennDiagram)
png("Pic03_01.png",width = 6, height = 6,units ="in" ,res=300)
plot.venn<-function(doe) {
  grid.newpage()
  draw.quad.venn(area1 = nrow(doe), 
                 area2 = nrow(doe),
                 area3 = nrow(doe),
                 area4=nrow(doe),
                 n12 = sum(doe$SF),
                 n13 = sum(doe$TrueS),
                 n14 = sum(doe$SLP),
                 n23 = sum(doe$TrueF),
                 n24 = sum(doe$FLP),
                 n34 = sum(doe$TrueLP),
                 n123 = sum(doe$Sch == doe$Ken & doe$Ken==doe$TS),
                 n124 = sum(doe$Sch == doe$Ken & doe$Ken==doe$LP),
                 n134 = sum(doe$Sch == doe$TS & doe$TS==doe$LP),
                 n234 = sum(doe$Ken == doe$TS & doe$TS==doe$LP),
                 n1234 = sum(doe$Ken == doe$TS & doe$TS==doe$LP & doe$Ken==doe$Sch), 
                 category = c("Шульце", "Кемени", "Исходные данные", "ЛП"),
                 cex = 2,cat.cex = 1.5,
                 fill = c("skyblue", "pink1", "mediumorchid","white"))
}
plot.venn(doe)
dev.off()

png("Pic03_02.png",width = 6, height = 6,units ="in" ,res=300)
plot.venn(subset(doe,val>3.5 & experts>12))
dev.off()

png("Pic03_03.png",width = 6, height = 4,units ="in" ,res=300)
library(rpart)
library(rpart.plot)
set.seed(2015)
idx<-sample(1:nrow(doe),3800)
test<-doe[idx,]
verify<-doe[-idx,]
mytree <- rpart(as.factor(TrueS) ~ val + experts + W, data = test)
prp(mytree,type = 2,extra = 104,fallen.leaves=TRUE)
dev.off()

setwd(WD)
