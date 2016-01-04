doe<-readRDS("doe.RDs")
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
plot.venn(subset(doe,val>3.5 & experts>11))
dev.off()

