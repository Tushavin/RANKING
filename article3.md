# Сравнение подходов к ранжированию
Тушавин В. А.  
4 января 2016 г.  

Целью настоящего исследования является сравнение алгоритмов ранжирования объектов с точки зрения их эффективности.

Сравниваются три алгоритма. Два из них были ранее описаны в статье Тушавин В. А. Ранжирование показателей качества с использованием методов Кемени-Янга и Шульце //Экономика и менеджмент систем управления. 2005. № 4.4.

Третий - быстрый алгоритм нахождения медианы Кемени из пакета `ConsRank`.

### Определение библиотек и функций


```r
library(ggplot2)
library(scales)
library(gtools)
library(ConsRank)
```

```
## Warning: package 'ConsRank' was built under R version 3.2.3
```

```
## Loading required package: MASS
## Loading required package: proxy
```

```
## Warning: package 'proxy' was built under R version 3.2.3
```

```
## 
## Attaching package: 'proxy'
## 
## The following objects are masked from 'package:stats':
## 
##     as.dist, dist
## 
## The following object is masked from 'package:base':
## 
##     as.matrix
## 
## Loading required package: rgl
## 
## Attaching package: 'ConsRank'
## 
## The following object is masked from 'package:base':
## 
##     labels
```

```r
library(lpSolve)
library(irr)
library(reshape2)
library(VennDiagram)
```

```
## Loading required package: grid
```

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
## Warning: package 'ROCR' was built under R version 3.2.3
```

```
## Loading required package: gplots
## 
## Attaching package: 'gplots'
## 
## The following object is masked from 'package:stats':
## 
##     lowess
```

```r
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
rank_solve<-function(ranks,Wk=NULL) {
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
  consensus<-matrix(best_rank+1,nrow=1,ncol=n_candidates)
  colnames(consensus)<-colnames(ranks)
  return(list(min_dist=tau,best_rank=best_rank,Consensus=consensus,Eltime=eltime))
}
```

Поскольку предыдущая версия функции ранжирования методом Шульце требовала полное ранжирование заданное в виде последовательности номеров элементов, то необходимо её немного преобразовать для работы с таблицей рангов.


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

Проведем тестирование на [примере из Википедии](https://en.wikipedia.org/wiki/Schulze_method)

number of voters | order of preference
-----------------|---------------------
5 | ACBED
5 | ADECB
8 | BEDAC
3 | CABED
7 | CAEBD
2 | CBADE
7 | DCEBA
8 | EBADC

Schulze ranking is E > A > C > B > D, and E wins. 


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
##    0.01
```

Результаты совпали полностью. Функция работает.
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
##    0.86
```

```r
QuickCons(ranks,Wk)
```

```
## $Consensus
##      A B C D E
## [1,] 3 2 5 4 1
## 
## $Tau
## [1] 0.1555556
## 
## $Eltime
## elapsed 
##    0.05
```

```r
EMCons(ranks,Wk)
```

```
## [1] "round 1"
## [1] "evaluating 1 branches"
## [1] "evaluating 3 branches"
## [1] "evaluating 13 branches"
## [1] "evaluating 33 branches"
## [1] "round 2"
## [1] "evaluating 1 branches"
## [1] "evaluating 3 branches"
## [1] "evaluating 13 branches"
## [1] "evaluating 18 branches"
```

```
## $Consensus
##      A B C D E
## [1,] 3 2 5 4 1
## 
## $Tau
## [1] 0.1555556
## 
## $Eltime
## elapsed 
##    0.11
```

Имеется расхождение, поскольку пример является несбалансированным по рангам.

Сравним результаты тестового примера из пакета `ConsRank`. Данные представляют собой ранжирование 130 студентами семи видов спорта в соответствии с их предпочтениями.
(Marden, J. I. (1996). Analyzing and modeling rank data. CRC Press.)


```r
data(sports)
colnames(sports)
```

```
## [1] "Baseball"   "Football"   "Basketball" "Tennis"     "Cycling"   
## [6] "Swimming"   "Jogging"
```

```r
colnames(sports)<-c("Бейсбол", "Футбол", "Баскетбол", "Теннис", "Велоспорт", "Плавание", "Бег трусцой")
dim(sports)
```

```
## [1] 130   7
```

```r
FASTcons(sports,maxiter=10)
```

```
## $Consensus
##      Бейсбол Футбол Баскетбол Теннис Велоспорт Плавание Бег трусцой
## [1,]       4      6         3      5         1        2           7
## 
## $Tau
##           [,1]
## [1,] 0.1443223
## 
## $Eltime
## elapsed 
##    0.33
```

```r
QuickCons(sports)
```

```
## $Consensus
##      Бейсбол Футбол Баскетбол Теннис Велоспорт Плавание Бег трусцой
## [1,]       4      6         3      5         1        2           7
## 
## $Tau
##           [,1]
## [1,] 0.1443223
## 
## $Eltime
## elapsed 
##    0.08
```

```r
EMCons(sports)
```

```
## [1] "round 1"
## [1] "evaluating 1 branches"
## [1] "evaluating 1 branches"
## [1] "evaluating 1 branches"
## [1] "evaluating 1 branches"
## [1] "evaluating 1 branches"
## [1] "evaluating 1 branches"
```

```
## $Consensus
##      Бейсбол Футбол Баскетбол Теннис Велоспорт Плавание Бег трусцой
## [1,]       4      6         3      5         1        2           7
## 
## $Tau
## [1] 0.1443223
## 
## $Eltime
## elapsed 
##    0.06
```

```r
Schulze.m(sports)
```

```
## $Consensus
##      Бейсбол Футбол Баскетбол Теннис Велоспорт Плавание Бег трусцой
## [1,]       4      6         3      5         1        2           7
## 
## $Eltime
## elapsed 
##    0.03
```

```r
rank_solve(sports)
```

```
## $min_dist
## [1] 1168
## 
## $best_rank
## [1] 3 5 2 4 0 1 6
## 
## $Consensus
##      Бейсбол Футбол Баскетбол Теннис Велоспорт Плавание Бег трусцой
## [1,]       4      6         3      5         1        2           7
## 
## $Eltime
## elapsed 
##    0.04
```

Результаты совпадают. 

### Сравнение времени выполнения алгоритмов


```r
#  Тест по времени
set.seed(1968)
d.len<-c()
d.n<-c()
d.mth<-c()
d.time<-c()

pb <- txtProgressBar(min = 0, max = 24, style = 3,file = stderr())
zzz<-0
for(rank_len in 3:8)
  for(n_ranks in c(5,10,15,20)) {
    ranks<-c()
    for(i in 1:n_ranks) ranks<-c(ranks,sample(1:rank_len,rank_len))
    ranks<-matrix(ranks,ncol=rank_len,byrow=T)
    
    start.time <- Sys.time()
    z<-FASTcons(ranks)
    end.time <- Sys.time()
    time.taken <-end.time - start.time
    d.len<-c(d.len,rank_len)
    d.n<-c(d.n,n_ranks)
    d.mth<-c(d.mth,"FASTcons")
    d.time<-c(d.time,as.numeric(time.taken))
    
     start.time <- Sys.time()
    z<-QuickCons(ranks)
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    d.len<-c(d.len,rank_len)
    d.n<-c(d.n,n_ranks)
    d.mth<-c(d.mth,"QuickCons")
    d.time<-c(d.time,as.numeric(time.taken))
    
    start.time <- Sys.time()
    z<-EMCons(ranks,PS=F)
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    d.len<-c(d.len,rank_len)
    d.n<-c(d.n,n_ranks)
    d.mth<-c(d.mth,"EMCons")
    d.time<-c(d.time,as.numeric(time.taken))
    
    start.time <- Sys.time()
    z<-Schulze.m(ranks)
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    d.len<-c(d.len,rank_len)
    d.n<-c(d.n,n_ranks)
    d.mth<-c(d.mth,"Schulze")
    d.time<-c(d.time,as.numeric(time.taken))
    
    start.time <- Sys.time()
    z<-rank_solve(ranks)
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    d.len<-c(d.len,rank_len)
    d.n<-c(d.n,n_ranks)
    d.mth<-c(d.mth,"LP")
    d.time<-c(d.time,as.numeric(time.taken))
    zzz<-zzz+1
    setTxtProgressBar(pb, zzz)
  }
close(pb)

mydata<-data.frame(Показателей=d.len,Экспертов=d.n,Время=d.time,Метод=d.mth)
mydata$Экспертов<-as.factor(mydata$Экспертов)
mydata$Метод<-as.factor(mydata$Метод)
g<-ggplot(aggregate(Время~Показателей+Метод,data=mydata,mean),aes(x=Показателей,y=Время,linetype=Метод,col=Метод))+
  geom_point()+
  geom_line(size=1)+
  scale_y_log10(breaks=trans_breaks("log10",function(x) 10^x),
                labels=trans_format("log10",math_format(10^.x)),
                minor_breaks=log10(5)+-4:1)+
  xlab("Количество оцениваемых показателей")+ylab("Время расчета, сек.")
g<-g+theme_bw(base_size = 16)
g  #+theme(legend.position=c(1,1),legend.justification=c(1,1))
```

![](article3_files/figure-html/alg.time-1.png) 



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
    doe$LP[i]<-paste0(names(r1[[1]]$LP$Consensus[1,])[order(r1[[1]]$LP$Consensus[1,])],collapse="")
    
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
   3         4      37      36        0    67    63    64          0.12            0.61      0.23
   4         4      52      60       37    80    53    45          0.13            0.79      0.30
   5         4      67      64       52    87    57    64          0.11            1.35      0.36
   6         4      72      73       63    93    67    68          0.21            1.80      0.55
   7         4      76      77       68    95    62    65          0.18            2.77      0.66
   8         4      89      86       79    95    72    75          0.31            3.90      0.81
   3         8      26      25        0    69    74    75          0.10            0.63      0.26
   4         8      58      61       40    85    80    77          0.21            0.95      0.27
   5         8      85      85       79    92    86    86          0.21            1.48      0.32
   6         8      93      91       87    98    94    96          0.29            2.01      0.58
   7         8      93      93       92    98    97    97          0.35            2.68      0.83
   8         8      98      97       95    99    97    98          0.49            4.05      0.96
   3        16      15      12        0    83    85    88          0.11            0.71      0.18
   4        16      81      80       74    95    93    94          0.24            1.12      0.34
   5        16      99      99       98   100    99    99          0.21            1.58      0.58
   6        16      99     100       99    99   100    99          0.34            2.42      0.68
   7        16     100     100       99   100    99    99          0.58            3.21      0.80
   8        16      99     100       99    99   100    99          0.66            4.56      1.25
   3        32      10      15        0    79    90    85          0.20            1.06      0.45
   4        32      92      94       91    98    99    97          0.29            1.55      0.51
   5        32     100     100      100   100   100   100          0.41            2.23      0.56
   6        32     100     100      100   100   100   100          0.57            2.87      0.80
   7        32     100     100      100   100   100   100          0.98            3.75      1.00
   8        32     100     100      100   100   100   100          1.15            5.42      1.41
   3        64       8      15        0    89    92    85          0.20            1.66      0.81
   4        64      99      99       99   100   100   100          0.39            2.22      0.88
   5        64     100     100      100   100   100   100          0.72            3.01      1.14
   6        64     100     100      100   100   100   100          0.92            3.91      1.42
   7        64     100     100      100   100   100   100          1.34            5.59      1.69
   8        64     100     100      100   100   100   100          2.02            6.90      2.15
   3       128       7      11        0    92    93    89          0.57            2.85      0.97
   4       128     100     100      100   100   100   100          0.81            3.72      1.32
   5       128     100     100      100   100   100   100          1.34            4.79      1.75
   6       128     100     100      100   100   100   100          1.99            6.14      2.23
   7       128     100     100      100   100   100   100          2.79            7.92      2.78
   8       128     100     100      100   100   100   100          3.64            9.94      3.61
   3       256       4       3        0    95    96    97          0.96            5.27      1.85
   4       256     100     100      100   100   100   100          1.71            6.57      2.39
   5       256     100     100      100   100   100   100          2.84            8.16      3.20
   6       256     100     100      100   100   100   100          3.91           10.53      4.09
   7       256     100     100      100   100   100   100          5.32           13.15      5.25
   8       256     100     100      100   100   100   100          6.86           16.32      6.63
   3       512       1       2        0    99    99    98          1.95            9.91      3.48
   4       512     100     100      100   100   100   100          3.38           12.48      4.52
   5       512     100     100      100   100   100   100          5.55           16.16      6.29
   6       512     100     100      100   100   100   100          7.78           20.34      8.17
   7       512     100     100      100   100   100   100         10.71           24.74     10.19
   8       512     100     100      100   100   100   100         13.65           29.16     12.42

### Время выполнения скрипта в зависимости от алгоритма


```r
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
## (polygon[GRID.polygon.149], polygon[GRID.polygon.150], polygon[GRID.polygon.151], polygon[GRID.polygon.152], polygon[GRID.polygon.153], polygon[GRID.polygon.154], polygon[GRID.polygon.155], polygon[GRID.polygon.156], text[GRID.text.157], text[GRID.text.158], text[GRID.text.159], text[GRID.text.160], text[GRID.text.161], text[GRID.text.162], text[GRID.text.163], text[GRID.text.164], text[GRID.text.165], text[GRID.text.166], text[GRID.text.167], text[GRID.text.168], text[GRID.text.169], text[GRID.text.170], text[GRID.text.171], text[GRID.text.172], text[GRID.text.173], text[GRID.text.174], text[GRID.text.175])
```


# Анализ зависимости


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
plot(prf,colorize=T)
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
plot(prf,colorize=T)
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

Построим такую же диаграмму для числа оцениваемых параметров больше 3.5 и число экспертов больше 11


```r
plot.venn(subset(doe,val>3.5 & experts>11))
```

![](article3_files/figure-html/venn.fin-1.png) 

```
## (polygon[GRID.polygon.176], polygon[GRID.polygon.177], polygon[GRID.polygon.178], polygon[GRID.polygon.179], polygon[GRID.polygon.180], polygon[GRID.polygon.181], polygon[GRID.polygon.182], polygon[GRID.polygon.183], text[GRID.text.184], text[GRID.text.185], text[GRID.text.186], text[GRID.text.187], text[GRID.text.188], text[GRID.text.189], text[GRID.text.190], text[GRID.text.191], text[GRID.text.192], text[GRID.text.193], text[GRID.text.194], text[GRID.text.195], text[GRID.text.196], text[GRID.text.197], text[GRID.text.198], text[GRID.text.199], text[GRID.text.200], text[GRID.text.201], text[GRID.text.202])
```


### Информация о параметрах R


```r
sessionInfo()
```

```
## R version 3.2.2 (2015-08-14)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 8 x64 (build 9200)
## 
## locale:
## [1] LC_COLLATE=Russian_Russia.1251  LC_CTYPE=Russian_Russia.1251   
## [3] LC_MONETARY=Russian_Russia.1251 LC_NUMERIC=C                   
## [5] LC_TIME=Russian_Russia.1251    
## 
## attached base packages:
## [1] grid      stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] ROCR_1.0-7        gplots_2.17.0     caret_6.0-62     
##  [4] lattice_0.20-33   rpart.plot_1.5.3  rpart_4.1-10     
##  [7] VennDiagram_1.6.9 reshape2_1.4.1    irr_0.84         
## [10] lpSolve_5.6.13    ConsRank_1.0.2    rgl_0.95.1367    
## [13] proxy_0.4-15      MASS_7.3-45       gtools_3.5.0     
## [16] scales_0.3.0      ggplot2_1.0.1    
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_0.12.2        highr_0.5.1        nloptr_1.0.4      
##  [4] formatR_1.2.1      plyr_1.8.3         class_7.3-14      
##  [7] bitops_1.0-6       iterators_1.0.8    tools_3.2.2       
## [10] digest_0.6.8       lme4_1.1-10        evaluate_0.8      
## [13] gtable_0.1.2       nlme_3.1-122       mgcv_1.8-9        
## [16] Matrix_1.2-3       foreach_1.4.3      parallel_3.2.2    
## [19] yaml_2.1.13        SparseM_1.7        proto_0.3-10      
## [22] e1071_1.6-7        stringr_1.0.0      knitr_1.11        
## [25] caTools_1.17.1     MatrixModels_0.4-1 stats4_3.2.2      
## [28] nnet_7.3-11        rmarkdown_0.8.1    gdata_2.17.0      
## [31] minqa_1.2.4        car_2.1-0          magrittr_1.5      
## [34] codetools_0.2-14   htmltools_0.2.6    splines_3.2.2     
## [37] pbkrtest_0.4-2     colorspace_1.2-6   labeling_0.3      
## [40] quantreg_5.19      KernSmooth_2.23-15 stringi_1.0-1     
## [43] munsell_0.4.2
```
