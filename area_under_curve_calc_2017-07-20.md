# area-under_curve
Stacey Harmer  
7/20/2017  



My goal is to automate a previously created function so that I can easily determine qPCR area-under-curve for individual genes and across unique time frames (time frames unique to each genotype and represent one circadian day).


```r
library(tidyverse)
```

```
## Loading tidyverse: ggplot2
## Loading tidyverse: tibble
## Loading tidyverse: tidyr
## Loading tidyverse: readr
## Loading tidyverse: purrr
## Loading tidyverse: dplyr
```

```
## Conflicts with tidy packages ----------------------------------------------
```

```
## filter(): dplyr, stats
## lag():    dplyr, stats
```

```r
library(MESS)
```

```
## Loading required package: geepack
```

```
## Loading required package: geeM
```

```
## Loading required package: Matrix
```

```
## 
## Attaching package: 'Matrix'
```

```
## The following object is masked from 'package:tidyr':
## 
##     expand
```

Below is the original code:


```r
Akiva.12deg<- read_csv("qRT-PCR_data_12c_6genes_n.csv") #read it in as a tibble
```

```
## Parsed with column specification:
## cols(
##   time = col_integer(),
##   expn = col_double(),
##   gene = col_character(),
##   gt = col_character(),
##   light = col_character(),
##   rep = col_character()
## )
```

```r
head(Akiva.12deg)
```

```
## # A tibble: 6 × 6
##    time      expn  gene    gt light   rep
##   <int>     <dbl> <chr> <chr> <chr> <chr>
## 1    24 0.7444322   LUX   Col    LL  rep1
## 2    27 0.3420660   LUX   Col    LL  rep1
## 3    30 0.7249853   LUX   Col    LL  rep1
## 4    33 1.5997045   LUX   Col    LL  rep1
## 5    36 1.2732694   LUX   Col    LL  rep1
## 6    39 0.8051121   LUX   Col    LL  rep1
```

```r
tail(Akiva.12deg)
```

```
## # A tibble: 6 × 6
##    time      expn  gene      gt light   rep
##   <int>     <dbl> <chr>   <chr> <chr> <chr>
## 1    42 0.9244249  PRR9 cca1lhy    LL  rep2
## 2    45 1.1748752  PRR9 cca1lhy    LL  rep2
## 3    48 0.9544120  PRR9 cca1lhy    LL  rep2
## 4    51 0.6510480  PRR9 cca1lhy    LL  rep2
## 5    54 0.5624138  PRR9 cca1lhy    LL  rep2
## 6    57 0.4736999  PRR9 cca1lhy    LL  rep2
```

Make the data wide

```r
PCR.12.wide <- spread(Akiva.12deg,key = gene, value = expn ) %>%
  mutate(gt.rep = paste(gt, rep, sep = "."))
head(PCR.12.wide)
```

```
## # A tibble: 6 × 11
##    time            gt light   rep      ELF4       LUX      PRR5       PRR7
##   <int>         <chr> <chr> <chr>     <dbl>     <dbl>     <dbl>      <dbl>
## 1    24       cca1lhy    LL  rep1 0.7614978 1.0159563 0.8754983 1.08389818
## 2    24       cca1lhy    LL  rep2 0.5325001 0.8054884 0.5434463 0.62551410
## 3    24 cca1lhyrve468    LL  rep1 0.6674427 0.8242006 0.4024810 0.06825593
## 4    24 cca1lhyrve468    LL  rep2 0.7459150 0.6513596 0.5311108 0.09719690
## 5    24           Col    LL  rep1 0.4225847 0.7444322 0.1598923 0.06971710
## 6    24           Col    LL  rep2 0.3057271 0.4399131 0.1450424 0.10571808
## # ... with 3 more variables: PRR9 <dbl>, TOC1 <dbl>, gt.rep <chr>
```

Original code:


```r
PRR5_12deg_24_48.6 <- sapply(unique(PCR.12.wide$gt.rep), function(gt.rep) {
  tmp <- PCR.12.wide[PCR.12.wide$gt.rep==gt.rep,]
  auc(tmp$time,tmp$PRR5,type="spline", from =24, to = 48.6)
})
PRR5_12deg_24_48.6
```

```
##       cca1lhy.rep1       cca1lhy.rep2 cca1lhyrve468.rep1 
##           41.80379           32.55078           19.65421 
## cca1lhyrve468.rep2           Col.rep1           Col.rep2 
##           23.24228           20.81449           26.72926 
##        rve468.rep1        rve468.rep2 
##           16.23943           18.51161
```

I want to be able to step through each of the genes.  would it be better if I used long insetad of wide dataframe?


```r
Akiva.12deg <- Akiva.12deg %>%
  mutate(gt.rep = paste(gt, rep, sep = "."))

head(Akiva.12deg)
```

```
## # A tibble: 6 × 7
##    time      expn  gene    gt light   rep   gt.rep
##   <int>     <dbl> <chr> <chr> <chr> <chr>    <chr>
## 1    24 0.7444322   LUX   Col    LL  rep1 Col.rep1
## 2    27 0.3420660   LUX   Col    LL  rep1 Col.rep1
## 3    30 0.7249853   LUX   Col    LL  rep1 Col.rep1
## 4    33 1.5997045   LUX   Col    LL  rep1 Col.rep1
## 5    36 1.2732694   LUX   Col    LL  rep1 Col.rep1
## 6    39 0.8051121   LUX   Col    LL  rep1 Col.rep1
```

```r
#all_12deg_24_48.6 <- sapply(unique(Akiva.12deg$gt.rep), function(gt.rep) {
#  gene <- Akiva.12deg[Akiva.12deg$gene==gene,],
#  tmp <- Akiva.12deg[Akiva.12deg$gt.rep==gt.rep,]
#  auc(tmp$time,tmp$gene,type="spline", from =24, to = 48.6)
#})
```
that failed. 


```r
head(PCR.12.wide)
```

```
## # A tibble: 6 × 11
##    time            gt light   rep      ELF4       LUX      PRR5       PRR7
##   <int>         <chr> <chr> <chr>     <dbl>     <dbl>     <dbl>      <dbl>
## 1    24       cca1lhy    LL  rep1 0.7614978 1.0159563 0.8754983 1.08389818
## 2    24       cca1lhy    LL  rep2 0.5325001 0.8054884 0.5434463 0.62551410
## 3    24 cca1lhyrve468    LL  rep1 0.6674427 0.8242006 0.4024810 0.06825593
## 4    24 cca1lhyrve468    LL  rep2 0.7459150 0.6513596 0.5311108 0.09719690
## 5    24           Col    LL  rep1 0.4225847 0.7444322 0.1598923 0.06971710
## 6    24           Col    LL  rep2 0.3057271 0.4399131 0.1450424 0.10571808
## # ... with 3 more variables: PRR9 <dbl>, TOC1 <dbl>, gt.rep <chr>
```

I suspect I need to incoroprate something like the below

 for(i in 5:10)
    { 
      PCR.12.wide[, i];
    }
but I can't figure it out

Yes, that would be one way to do it, but I wonder if we could use summarize instead


```r
results <- Akiva.12deg %>%
  group_by(gene,gt,rep) %>%
  summarize(auc=auc(time,expn,type="spline", from =24, to = 48.6))
```

let's compare the prr5 results from this method and the Stacey's single call

```r
results %>% filter(gene=="PRR5")
```

```
## Source: local data frame [8 x 4]
## Groups: gene, gt [4]
## 
##    gene            gt   rep      auc
##   <chr>         <chr> <chr>    <dbl>
## 1  PRR5       cca1lhy  rep1 41.80379
## 2  PRR5       cca1lhy  rep2 32.55078
## 3  PRR5 cca1lhyrve468  rep1 19.65421
## 4  PRR5 cca1lhyrve468  rep2 23.24228
## 5  PRR5           Col  rep1 20.81449
## 6  PRR5           Col  rep2 26.72926
## 7  PRR5        rve468  rep1 16.23943
## 8  PRR5        rve468  rep2 18.51161
```


```r
PRR5_12deg_24_48.6
```

```
##       cca1lhy.rep1       cca1lhy.rep2 cca1lhyrve468.rep1 
##           41.80379           32.55078           19.65421 
## cca1lhyrve468.rep2           Col.rep1           Col.rep2 
##           23.24228           20.81449           26.72926 
##        rve468.rep1        rve468.rep2 
##           16.23943           18.51161
```

Success!
