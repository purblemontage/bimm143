# class05_data_visualization
James Woolley, A16440072

## using GGPLOT

The ggplot2 package needs to be installed here as it does not come with
R “out of the box”. Used the `install.packages()` function to do it.

To use ggplot, you need to load it before using any of the functions it
comes with. Do this with the `library()` function.

``` r
head(cars)
```

      speed dist
    1     4    2
    2     4   10
    3     7    4
    4     7   22
    5     8   16
    6     9   10

``` r
library(ggplot2)
ggplot(cars) +
  aes(x=speed,y=dist) +
  geom_point() +
  labs(title="Distance to Stop at Differenrt Speeds of Cars",
       x="Speed (MPH)",
       y="Distance (Ft)",
       subtitle="This is a plot with a line of best fit",
       caption="Dataset: Cars")+
  geom_smooth(method = "lm",se=FALSE)
```

    `geom_smooth()` using formula = 'y ~ x'

![](class05_render_files/figure-commonmark/unnamed-chunk-2-1.png)

All ggplot figures have at least 3 things. 1- data 2- aesthetic mapping
3- geoms

There are lots of graphing systems in R, including one that comes with
it.

``` r
plot(cars)
```

![](class05_render_files/figure-commonmark/unnamed-chunk-3-1.png)

``` r
url <- "https://bioboot.github.io/bimm143_S20/class-material/up_down_expression.txt"
genes <- read.delim(url)
head(genes)
```

            Gene Condition1 Condition2      State
    1      A4GNT -3.6808610 -3.4401355 unchanging
    2       AAAS  4.5479580  4.3864126 unchanging
    3      AASDH  3.7190695  3.4787276 unchanging
    4       AATF  5.0784720  5.0151916 unchanging
    5       AATK  0.4711421  0.5598642 unchanging
    6 AB015752.4 -3.6808610 -3.5921390 unchanging

``` r
nrow(genes)
```

    [1] 5196

``` r
colnames(genes)
```

    [1] "Gene"       "Condition1" "Condition2" "State"     

``` r
ncol(genes)
```

    [1] 4

``` r
table(genes$State)
```


          down unchanging         up 
            72       4997        127 

``` r
round( table(genes$State)/nrow(genes) * 100, 2 )
```


          down unchanging         up 
          1.39      96.17       2.44 

``` r
library(ggplot2)
p <- ggplot(genes)+
  aes(x=Condition1,y=Condition2, col=State) +
  geom_point()+
  scale_colour_manual( values=c("orange","green","white") )+
   labs(title="Gene Expression Changes Upon Drug Treatment",
       x="Control (no drug)",
       y="Drug Treatment")
p
```

![](class05_render_files/figure-commonmark/unnamed-chunk-6-1.png)

``` r
allstates <- genes$State
nstates <- table(allstates)
nstates/nrow(genes)
```

    allstates
          down unchanging         up 
    0.01385681 0.96170131 0.02444188 
