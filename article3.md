# Сравнение подходов к ранжированию
Тушавин В. А.  
4 января 2016 г.  

Целью настоящего исследования является сравнение алгоритмов ранжирования объектов с точки зрения их эффективности.
Сравниваются три алгоритма. Два из них были ранее описаны в статье Тушавин В. А. Ранжирование показателей качества с использованием методов Кемени-Янга и Шульце //Экономика и менеджмент систем управления. 2005. № 4.4.
Третий - быстрый алгоритм нахождения медианы Кемени из пакета `ConsRank`.

### Определение библиотек и функций


```r
library(gtools)
library(ConsRank)
```

```
## Loading required package: MASS
## Loading required package: proxy
## 
## Attaching package: 'proxy'
## 
## Следующие объекты скрыты от 'package:stats':
## 
##     as.dist, dist
## 
## Следующий объект скрыт от 'package:base':
## 
##     as.matrix
## 
## Loading required package: rgl
## 
## Attaching package: 'ConsRank'
## 
## Следующий объект скрыт от 'package:base':
## 
##     labels
```

```r
library(lpSolve)
library(irr)

# Функции для нахождения медианы Кемени
# Нахождение расстояния между оценками
kendall_tau<-function(rank.a,rank.b) {
  tau<-0
  n<-length(rank.a)
  for(k in 1:ncol(z<-combn(n,2))) {
    i=z[1,k]
    j=z[2,k]
    tau<-tau+(sign(rank.a[i]-rank.a[j]) == -sign(rank.b[i]-rank.b[j]))
  } 
  return(tau)
}


# Построение графа
build_graph<-function(ranks) {
  n_voters<-nrow(ranks)
  n_candidates<-ncol(ranks)
  edge_weights<-matrix(0,nrow=n_candidates,ncol=n_candidates)
  for(k in 1:ncol(z<-combn(n_candidates,2))) {
    i=z[1,k]
    j=z[2,k]
    preference<-ranks[, i] - ranks[, j]
    h.ij <- sum(preference < 0) 
    h.ji <- sum(preference > 0)
    if(h.ij > h.ji) edge_weights[i, j] <- h.ij - h.ji else if(h.ij < h.ji) edge_weights[j, i] <- h.ji - h.ij
  }
  return(edge_weights)
}

# Нахождение медианы Кемени посредством решения задачи ЛП
rank_solve<-function(ranks) {
  tic = proc.time()[3]
  n_voters<-nrow(ranks)
  n_candidates<-ncol(ranks)
  # Строим граф
  edge_weights<-build_graph(ranks)
  # Задаем параметры. 
  # Коээфициенты при целевой функции
  objective.in<- as.vector(t(edge_weights))
  # Коэффициенты для каждой пары
  pairwise_constraints <- matrix(0,
                                 n_candidates * (n_candidates - 1) / 2, n_candidates ^ 2)
  for(k in 1:nrow(z<-combinations(n_candidates,2))) {
    i=z[k,1]
    j=z[k,2] 
    pairwise_constraints[k,c((i-1)*n_candidates+j,(j-1)*n_candidates+i)]<-1
  }
  # Коэффициенты для каждой тройки
  triangle_constraints <-matrix(0,n_candidates *
                                  (n_candidates - 1) *
                                  (n_candidates - 2), n_candidates ^ 2)
  
  for(m in 1:nrow(z<-permutations(n_candidates,3))) {
    i=z[m,1]
    j=z[m,2]
    k=z[m,3]
    triangle_constraints[m,c((i-1)*n_candidates+j,(j-1)*n_candidates+k,(k-1)*n_candidates+i)]<-1
  }
  constraints<-rbind(pairwise_constraints,triangle_constraints)
  constraint_rhs<-rep(1,nrow(pairwise_constraints)+nrow(triangle_constraints))
  constraint_signs<-c(rep("==",nrow(pairwise_constraints)),rep(">=",nrow(triangle_constraints)))
  z<-lp("min",objective.in, constraints, constraint_signs, constraint_rhs,all.int=T) 
  x<-matrix(z$solution,nrow=n_candidates,ncol=n_candidates,byrow=T)
  best_rank<-apply(x,1,sum)
  tau<-sum(apply(ranks,1,function(x){kendall_tau(x,best_rank)}))
  toc = proc.time()[3]
  eltime = toc - tic
  Consensus=best_rank+1
  names(Consensus)<-LETTERS[1:n_candidates]
  return(list(min_dist=tau,best_rank=best_rank,Consensus=Consensus,Eltime=eltime))
}
```

Поскольку предыдущая версия функции нахождения ранжирования методом Шульце требовала полное ранжирование заданное в виде последовательности номеров элементов, то необходимо её немного преобразовать для работы с таблицей рангов 


```r
# Модифицированния функция нахождения итогового
# ранжирования методом Шульце
Schulze.m<-function(ranks,Wk=NULL) {
  tic = proc.time()[3]
  if (class(ranks) == "data.frame") {
    ranks = as.matrix(ranks)
  }  
  n_voters<-nrow(ranks)
  n_candidates<-ncol(ranks)
  ranks[is.na(ranks)]<-Inf
  if (n_voters == 1) {
    consensus = ranks
  } else {
    mtx<-matrix(data=0,nrow=n_candidates,ncol=n_candidates)
    rownames(mtx)<-colnames(ranks)
    colnames(mtx)<-colnames(ranks)
    index<-combinations(n_candidates,2)
    for(i in 1:n_voters){
      temp<-matrix(data=0,nrow=n_candidates,ncol=n_candidates)
      for(idx in 1:nrow(index)) {
        x1=index[idx,1]
        x2=index[idx,2]
        if(ranks[i,x1]<ranks[i,x2]) temp[x1,x2]<-temp[x1,x2]+1
        if(ranks[i,x1]>ranks[i,x2]) temp[x2,x1]<-temp[x2,x1]+1
      }
      if(!is.null(Wk)) temp<-temp*Wk[i]
      mtx<-mtx+temp
    }
  result<-matrix(data=0,nrow=n_candidates,ncol=n_candidates)
  for(i in 1:n_candidates)
    for(j in 1:n_candidates)
      if(i!=j) result[i,j]<-ifelse(mtx[i,j] > mtx[j,i],mtx[i,j],0)
  for(i in 1:n_candidates)
    for(j in 1:n_candidates)
      if(i!=j) for(k in 1:n_candidates)
        if(i!=k & j !=k) result[j,k]<-max(result[j,k],
                                          min(result[j,i],result[i,k]))
  vec<-rep(0,n_candidates)
  for(k in 1:nrow(z<-combinations(n_candidates,2))) {
    i=z[k,1]
    j=z[k,2]
    if(result[i,j]>result[j,i]) 
      vec[j]<-vec[j]+1
    else if(result[i,j]<result[j,i])
      vec[i]<-vec[i]+1
  }
  } 
  consensus<-matrix(vec,nrow=1,ncol=n_candidates)
  colnames(consensus)<-colnames(ranks)
toc = proc.time()[3]
eltime = toc - tic
return(list(Consensus=consensus+1,Eltime=eltime))
}
```

Проведем тестирование на примере из Википедии


```r
test<-data.frame(order=c("ACBED","ADECB",
                 "BEDAC",
                 "CABED",
                 "CAEBD",
                 "CBADE",
                 "DCEBA",
                 "EBADC"),
                  wk=c(5,5,8,3,7,2,7,8))  
ranks<-matrix(0,nrow=nrow(test),ncol=5)
colnames(ranks)<-LETTERS[1:5]
for(i in 1:nrow(test)) 
  ranks[i,t(asc(as.character(test$order[i]))-64)]<-1:5
Wk<-test$wk
Schulze.m(ranks,Wk)
```

```
## $Consensus
##      A B C D E
## [1,] 2 4 3 5 1
## 
## $Eltime
## elapsed 
##   0.003
```

Результаты свопали полностью. Функция работает.
Проверим результат с помощью пакета `ConsRank`


```r
FASTcons(ranks,Wk)
```

```
## $Consensus
##      A B C D E
## [1,] 3 2 5 4 1
## 
## $Tau
##           [,1]
## [1,] 0.1555556
## 
## $Eltime
## elapsed 
##   1.197
```

Имеется расхождение, поскольку пример является несбалансированным по рангам.

Сравним результаты тестового примера из пакета `ConsRank`


```r
data(sports)
Schulze.m(sports)
```

```
## $Consensus
##      Baseball Football Basketball Tennis Cycling Swimming Jogging
## [1,]        4        6          3      5       1        2       7
## 
## $Eltime
## elapsed 
##   0.047
```

```r
FASTcons(sports,maxiter=10)
```

```
## $Consensus
##      Baseball Football Basketball Tennis Cycling Swimming Jogging
## [1,]        4        6          3      5       1        2       7
## 
## $Tau
##           [,1]
## [1,] 0.1443223
## 
## $Eltime
## elapsed 
##   0.404
```

Результаты удовлетворительны. 

### Задание функций тестирования

Пусть имеется матрица рангов n x m, где n - эксперты,  а m - ранжируемые показатели и имеется случайный вектор рангов длиной m. Пусть имеется полный конценсус и все эксперты поставили одинаковую оценку равную этому случайному вектору. Внесем искажения в матрицу, поменяв в каждой строке, кроме первой, два случайно выбранных рядом стоящих ранга местами. Попробуем отыскать исходное ранжирование с помощью метода Шульце, медианы Кемени методом ЛП и методом ветвей и границ.



```r
change<-function(vec,cnt=1) {
  idx<-1:(length(vec)-1)
  for(i in 1:cnt) {
    k<-sample(idx,1)
    x1<-which(vec==k)
    x2<-which(vec==(k+1))
    x<-vec[x1]
    vec[x1]<-vec[x2]
    vec[x2]<-x
  }  
  return(vec)
}

dotest<-function(val=5,experts=10,tests=1,correct=0) {
  result<-list()
  # Start
  for(i in 1:tests) {
    x0<-sample(1:val,val)
    names(x0)<-LETTERS[1:val]
    z<-matrix(rep(x0,experts),byrow=T,ncol=val)
    colnames(z)<-LETTERS[1:val]
    z1<-t(apply(z,1,change,cnt=1))
    if(correct>0) for(j in 1:correct) z1[j,]<-x0
    result[[i]]<-list(kendall=kendall(t(z1)),
                      TrueValue=x0,
                      Shulze=Schulze.m(z1),
                      FASTcons=QuickCons(z1),
                      LP=rank_solve(z1-1))
  }
return(result)
}
```

### Тестирование


```r
if(!file.exists("doe.RDs")) { 
  stop("Бля!")
  set.seed(2015)
  doe<-expand.grid(testn=1:100,val=3:8,experts=c(4,8,16,32,64,128,256,512))
  doe$W<-NA
  doe$Shulze.time<-NA
  doe$FASTcons.time<-NA
  doe$LP.time<-NA
  doe$TS<-""
  doe$Sch<-""
  doe$Ken<-""
  doe$LP<-""
  doe$TrueS<-FALSE
  doe$TrueF<-FALSE
  doe$TrueLP<-FALSE
  doe$SF<-FALSE
  doe$SLP<-FALSE
  doe$FLP<-FALSE
  
  pb <- txtProgressBar(min = 0, max = nrow(doe), style = 3,file = stderr())
  for(i in 1:nrow(doe)) { 
    r1<-dotest(val=doe$val[i],experts=doe$experts[i],tests=1,correct=1)
    doe$W[i]<-r1[[1]]$kendall$value
    doe$Shulze.time[i]<-r1[[1]]$Shulze$Eltime
    doe$FASTcons.time[i]<-r1[[1]]$FASTcons$Eltime
    doe$LP.time[i]<-r1[[1]]$LP$Eltime
    
    doe$TS[i]<-paste0(names(r1[[1]]$TrueValue)[order(r1[[1]]$TrueValue)],collapse="")
    doe$Sch[i]<-paste0(names(r1[[1]]$Shulze$Consensus[1,])[order(r1[[1]]$Shulze$Consensus[1,])],collapse="")
    doe$Ken[i]<-paste0(names(r1[[1]]$FASTcons$Consensus[1,])[order(r1[[1]]$FASTcons$Consensus[1,])],collapse="")
    doe$LP[i]<-paste0(names(r1[[1]]$LP$Consensus)[order(r1[[1]]$LP$Consensus)],collapse="")
    
    doe$TrueS[i]<-(doe$TS[i] == doe$Sch[i])
    doe$TrueF[i]<-(doe$TS[i] == doe$Ken[i])
    doe$TrueLP[i]<-(doe$TS[i] == doe$LP[i])
    doe$SF[i]<-(doe$Sch[i] == doe$Ken[i])
    doe$SLP[i]<-(doe$Sch[i] == doe$LP[i])
    doe$FLP[i]<-(doe$Ken[i] == doe$LP[i])
    setTxtProgressBar(pb, i)
  }
  close(pb)
  ### Сохранение данных
  saveRDS(doe,"doe.RDs")
  library(xlsx)
  write.xlsx(doe, "mydata.xlsx", sheetName="Results",row.names=FALSE)
} else doe<-readRDS("doe.RDs")
  
knitr::kable(agg.doe<-aggregate(cbind(TrueS,TrueF,TrueLP,SF,SLP,FLP,Shulze.time,FASTcons.time,LP.time)~val+experts,data=doe,sum))
```



 val   experts   TrueS   TrueF   TrueLP    SF   SLP   FLP   Shulze.time   FASTcons.time   LP.time
----  --------  ------  ------  -------  ----  ----  ----  ------------  --------------  --------
   3         4      37      36        0    67    63    64         0.082           0.819     0.190
   4         4      52      60       37    80    53    45         0.113           1.113     0.223
   5         4      67      64       52    87    57    64         0.161           1.614     0.292
   6         4      72      73       63    93    67    68         0.229           2.322     0.432
   7         4      76      77       68    95    62    65         0.309           3.347     0.696
   8         4      89      86       79    95    72    75         0.443           4.893     1.092
   3         8      26      25        0    69    74    75         0.095           0.965     0.232
   4         8      58      61       40    85    80    77         0.153           1.413     0.307
   5         8      85      85       79    92    86    86         0.237           1.965     0.391
   6         8      93      91       87    98    94    96         0.334           2.691     0.538
   7         8      93      93       92    98    97    97         0.447           3.927     0.909
   8         8      98      97       95    99    97    98         0.734           6.533     1.421
   3        16      15      12        0    83    85    88         0.204           1.519     0.376
   4        16      81      80       74    95    93    94         0.250           1.848     0.473
   5        16      99      99       98   100    99    99         0.428           3.293     0.772
   6        16      99     100       99    99   100    99         0.674           4.841     1.027
   7        16     100     100       99   100    99    99         0.909           5.447     1.191
   8        16      99     100       99    99   100    99         0.855           5.863     1.313
   3        32      10      15        0    79    90    85         0.211           1.490     0.432
   4        32      92      94       91    98    99    97         0.462           2.331     0.648
   5        32     100     100      100   100   100   100         0.698           3.666     1.004
   6        32     100     100      100   100   100   100         0.961           4.633     1.296
   7        32     100     100      100   100   100   100         1.240           5.937     1.622
   8        32     100     100      100   100   100   100         1.564           7.989     2.160
   3        64       8      15        0    89    92    85         0.430           2.746     0.814
   4        64      99      99       99   100   100   100         0.790           4.173     1.227
   5        64     100     100      100   100   100   100         0.933           3.908     1.298
   6        64     100     100      100   100   100   100         1.300           5.148     1.631
   7        64     100     100      100   100   100   100         1.963           7.143     2.351
   8        64     100     100      100   100   100   100         2.577           9.512     3.202
   3       128       7      11        0    92    93    89         0.706           4.342     1.446
   4       128     100     100      100   100   100   100         1.197           5.150     1.775
   5       128     100     100      100   100   100   100         1.755           6.262     2.332
   6       128     100     100      100   100   100   100         4.284          13.640     5.005
   7       128     100     100      100   100   100   100         5.174          15.745     5.701
   8       128     100     100      100   100   100   100         4.999          13.862     5.398
   3       256       4       3        0    95    96    97         1.493           8.693     3.024
   4       256     100     100      100   100   100   100         2.862          12.151     4.355
   5       256     100     100      100   100   100   100         4.187          13.792     5.307
   6       256     100     100      100   100   100   100         4.995          14.446     5.793
   7       256     100     100      100   100   100   100         7.514          18.975     7.993
   8       256     100     100      100   100   100   100         9.351          21.477     9.596
   3       512       1       2        0    99    99    98         3.117          17.981     6.029
   4       512     100     100      100   100   100   100         4.340          16.396     6.351
   5       512     100     100      100   100   100   100         6.641          19.912     8.273
   6       512     100     100      100   100   100   100         9.394          24.662    10.737
   7       512     100     100      100   100   100   100        15.015          35.548    15.773
   8       512     100     100      100   100   100   100        18.084          39.237    18.170

### Время выполнения скрипта в зависимости от алгоритма


```r
library(ggplot2)
library(reshape2)
agg.doe<-agg.doe[,c(1,2,9:11)]
names(agg.doe)[3:5]<-c("Шульце","Кемени","ЛП")
md<-melt(agg.doe,id=c("val","experts"))
ggplot(md,aes(x=val,y=value,fill=variable))+
  geom_bar(stat="identity",colour="black")+
  scale_fill_discrete(name = "Алгоритм")+
  labs(x="Число параметров",y="Время, сек.")+
  theme_bw()+theme(text=element_text(size=14))
```

![](article3_files/figure-html/sc.time-1.png) 

```r
ggsave("Pic03_00.png",width=6,height=4,dpi=300)
```


### Построение диаграммы Венна


```r
library(VennDiagram)
```

```
## Loading required package: grid
## Loading required package: futile.logger
## 
## Attaching package: 'futile.logger'
## 
## Следующий объект скрыт от 'package:gtools':
## 
##     scat
```

```r
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
                 fill = c("skyblue", "pink1", "mediumorchid","white"))
}
plot.venn(doe)
```

![](article3_files/figure-html/venn.first-1.png) 

```
## (polygon[GRID.polygon.133], polygon[GRID.polygon.134], polygon[GRID.polygon.135], polygon[GRID.polygon.136], polygon[GRID.polygon.137], polygon[GRID.polygon.138], polygon[GRID.polygon.139], polygon[GRID.polygon.140], text[GRID.text.141], text[GRID.text.142], text[GRID.text.143], text[GRID.text.144], text[GRID.text.145], text[GRID.text.146], text[GRID.text.147], text[GRID.text.148], text[GRID.text.149], text[GRID.text.150], text[GRID.text.151], text[GRID.text.152], text[GRID.text.153], text[GRID.text.154], text[GRID.text.155], text[GRID.text.156], text[GRID.text.157], text[GRID.text.158], text[GRID.text.159])
```


# Анализ зависимости


```r
library(rpart)
library(rpart.plot)
library(caret)
```

```
## Loading required package: lattice
```

```r
library(ROCR)
```

```
## Loading required package: gplots
## 
## Attaching package: 'gplots'
## 
## Следующий объект скрыт '.GlobalEnv':
## 
##     plot.venn
## 
## Следующий объект скрыт от 'package:stats':
## 
##     lowess
```

```r
set.seed(2015)
idx<-sample(1:nrow(doe),3800)
test<-doe[idx,]
verify<-doe[-idx,]

(mylogit <- glm(TrueS ~ val + experts + W, data = test, family=binomial(link='logit')))
```

```
## 
## Call:  glm(formula = TrueS ~ val + experts + W, family = binomial(link = "logit"), 
##     data = test)
## 
## Coefficients:
## (Intercept)          val      experts            W  
##   -5.077245     0.707930     0.003041     4.254768  
## 
## Degrees of Freedom: 3799 Total (i.e. Null);  3796 Residual
## Null Deviance:	    3710 
## Residual Deviance: 2097 	AIC: 2105
```

```r
summary(mylogit)
```

```
## 
## Call:
## glm(formula = TrueS ~ val + experts + W, family = binomial(link = "logit"), 
##     data = test)
## 
## Deviance Residuals: 
##      Min        1Q    Median        3Q       Max  
## -3.10587   0.09522   0.19664   0.46476   1.94230  
## 
## Coefficients:
##               Estimate Std. Error z value Pr(>|z|)    
## (Intercept) -5.0772454  0.2269911 -22.368  < 2e-16 ***
## val          0.7079303  0.0910490   7.775 7.53e-15 ***
## experts      0.0030411  0.0004047   7.515 5.71e-14 ***
## W            4.2547684  0.4587072   9.276  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for binomial family taken to be 1)
## 
##     Null deviance: 3709.8  on 3799  degrees of freedom
## Residual deviance: 2097.3  on 3796  degrees of freedom
## AIC: 2105.3
## 
## Number of Fisher Scoring iterations: 6
```

```r
confint(mylogit)
```

```
## Waiting for profiling to be done...
```

```
##                    2.5 %       97.5 %
## (Intercept) -5.532272880 -4.641835510
## val          0.533344617  0.890616275
## experts      0.002261516  0.003848818
## W            3.360742352  5.160190314
```

```r
fitted.results<-predict(mylogit,newdata=verify,type='response')
fitted.results <- (fitted.results > 0.5)
confusionMatrix(fitted.results,verify$TrueS)
```

```
## Confusion Matrix and Statistics
## 
##           Reference
## Prediction FALSE TRUE
##      FALSE   158   21
##      TRUE     55  766
##                                           
##                Accuracy : 0.924           
##                  95% CI : (0.9058, 0.9397)
##     No Information Rate : 0.787           
##     P-Value [Acc > NIR] : < 2.2e-16       
##                                           
##                   Kappa : 0.7593          
##  Mcnemar's Test P-Value : 0.0001535       
##                                           
##             Sensitivity : 0.7418          
##             Specificity : 0.9733          
##          Pos Pred Value : 0.8827          
##          Neg Pred Value : 0.9330          
##              Prevalence : 0.2130          
##          Detection Rate : 0.1580          
##    Detection Prevalence : 0.1790          
##       Balanced Accuracy : 0.8576          
##                                           
##        'Positive' Class : FALSE           
## 
```

```r
fitted.results<-predict(mylogit,newdata=verify,type='response')
pr <- prediction(fitted.results, verify$TrueS)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)
```

![](article3_files/figure-html/hh-1.png) 

```r
auc <- performance(pr, measure = "auc")
auc <- auc@y.values[[1]]
auc
```

```
## [1] 0.9147771
```

```r
(mytree <- rpart(as.factor(TrueS) ~ val + experts + W, data = test))
```

```
## n= 3800 
## 
## node), split, n, loss, yval, (yprob)
##       * denotes terminal node
## 
##  1) root 3800 727 TRUE (0.19131579 0.80868421)  
##    2) val< 3.5 615  87 FALSE (0.85853659 0.14146341) *
##    3) val>=3.5 3185 199 TRUE (0.06248038 0.93751962)  
##      6) experts< 12 813 173 TRUE (0.21279213 0.78720787)  
##       12) val< 4.5 162  71 TRUE (0.43827160 0.56172840)  
##         24) W>=0.7375 84  27 FALSE (0.67857143 0.32142857) *
##         25) W< 0.7375 78  14 TRUE (0.17948718 0.82051282) *
##       13) val>=4.5 651 102 TRUE (0.15668203 0.84331797) *
##      7) experts>=12 2372  26 TRUE (0.01096121 0.98903879) *
```

```r
prp(mytree,type = 2,extra = 104,fallen.leaves=TRUE)
```

![](article3_files/figure-html/hh-2.png) 

```r
fitted.results<-predict(mytree,newdata=verify,type = "class")
misClasificError <- mean(fitted.results != verify$TrueS)
print(paste('Accuracy',1-misClasificError))
```

```
## [1] "Accuracy 0.936"
```

```r
confusionMatrix(fitted.results,verify$TrueS)
```

```
## Confusion Matrix and Statistics
## 
##           Reference
## Prediction FALSE TRUE
##      FALSE   178   29
##      TRUE     35  758
##                                          
##                Accuracy : 0.936          
##                  95% CI : (0.919, 0.9504)
##     No Information Rate : 0.787          
##     P-Value [Acc > NIR] : <2e-16         
##                                          
##                   Kappa : 0.8071         
##  Mcnemar's Test P-Value : 0.532          
##                                          
##             Sensitivity : 0.8357         
##             Specificity : 0.9632         
##          Pos Pred Value : 0.8599         
##          Neg Pred Value : 0.9559         
##              Prevalence : 0.2130         
##          Detection Rate : 0.1780         
##    Detection Prevalence : 0.2070         
##       Balanced Accuracy : 0.8994         
##                                          
##        'Positive' Class : FALSE          
## 
```

```r
pr <- prediction(predict(mytree,newdata=verify,type = "vector"), verify$TrueS)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)
```

![](article3_files/figure-html/hh-3.png) 

```r
auc <- performance(pr, measure = "auc")
auc <- auc@y.values[[1]]
auc
```

```
## [1] 0.899416
```

Построим такую же диаграмму для числа оцениваемых параметров больше 3.5 и число эеспертов больше 11


```r
plot.venn(subset(doe,val>3.5 & experts>11))
```

![](article3_files/figure-html/venn.fin-1.png) 

```
## (polygon[GRID.polygon.170], polygon[GRID.polygon.171], polygon[GRID.polygon.172], polygon[GRID.polygon.173], polygon[GRID.polygon.174], polygon[GRID.polygon.175], polygon[GRID.polygon.176], polygon[GRID.polygon.177], text[GRID.text.178], text[GRID.text.179], text[GRID.text.180], text[GRID.text.181], text[GRID.text.182], text[GRID.text.183], text[GRID.text.184], text[GRID.text.185], text[GRID.text.186], text[GRID.text.187], text[GRID.text.188], text[GRID.text.189], text[GRID.text.190], text[GRID.text.191], text[GRID.text.192], text[GRID.text.193], text[GRID.text.194], text[GRID.text.195], text[GRID.text.196])
```


### Информация о параметрах R


```r
sessionInfo()
```

```
## R version 3.2.3 (2015-12-10)
## Platform: x86_64-apple-darwin13.4.0 (64-bit)
## Running under: OS X 10.11.2 (El Capitan)
## 
## locale:
## [1] ru_RU.UTF-8/ru_RU.UTF-8/ru_RU.UTF-8/C/ru_RU.UTF-8/ru_RU.UTF-8
## 
## attached base packages:
## [1] grid      stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] ROCR_1.0-7          gplots_2.17.0       caret_6.0-62       
##  [4] lattice_0.20-33     rpart.plot_1.5.3    rpart_4.1-10       
##  [7] VennDiagram_1.6.16  futile.logger_1.4.1 reshape2_1.4.1     
## [10] ggplot2_2.0.0       irr_0.84            lpSolve_5.6.13     
## [13] ConsRank_1.0.2      rgl_0.95.1435       proxy_0.4-15       
## [16] MASS_7.3-45         gtools_3.5.0       
## 
## loaded via a namespace (and not attached):
##  [1] splines_3.2.3        colorspace_1.2-6     htmltools_0.3       
##  [4] stats4_3.2.3         yaml_2.1.13          mgcv_1.8-9          
##  [7] e1071_1.6-7          nloptr_1.0.4         lambda.r_1.1.7      
## [10] foreach_1.4.3        plyr_1.8.3           stringr_1.0.0       
## [13] MatrixModels_0.4-1   munsell_0.4.2        gtable_0.1.2        
## [16] caTools_1.17.1       codetools_0.2-14     evaluate_0.8        
## [19] labeling_0.3         knitr_1.11           SparseM_1.7         
## [22] class_7.3-14         quantreg_5.19        pbkrtest_0.4-4      
## [25] parallel_3.2.3       highr_0.5.1          Rcpp_0.12.2         
## [28] KernSmooth_2.23-15   scales_0.3.0         formatR_1.2.1       
## [31] gdata_2.17.0         lme4_1.1-10          digest_0.6.8        
## [34] stringi_1.0-1        tools_3.2.3          bitops_1.0-6        
## [37] magrittr_1.5         futile.options_1.0.0 car_2.1-1           
## [40] Matrix_1.2-3         minqa_1.2.4          rmarkdown_0.9       
## [43] iterators_1.0.8      nnet_7.3-11          nlme_3.1-122
```
