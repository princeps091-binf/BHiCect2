# BHiCect2

<!-- badges: start -->

<!-- badges: end -->

The goal of BHiCect2 is to cluster HiC data as described in our [manuscript](). <!-- markdownlint-disable-line MD042 -->
Briefly, we decompose intra-chromosomal HiC data into nested clusters of chromosome regions across multiple resolutions
starting from the complete chromosome all the way to DNA-loops at the maximum resolution provided.

## Installation

You can install the current version of BHiCect2 from [GitHub](https://github.com/) with:

```r
# install.packages("devtools")
devtools::install_github("princeps091-binf/BHiCect2")
```

## Example

BHiCect2 offers a core function to cluster the input HiC data.
The expected input data to BHiCect2 is a list of dataframes containing the HiC data in a three columns format for the various resolution provided by the user.

```r
library(BHiCect2)
## basic example using chromosome 22 data from human GM12878 data (Rao et al. 2014)
data(chr_dat_l)

# manual entry for resolutions
res_set<-c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num<-c(1e6L,5e5L,1e5L,5e4L,1e4L,5e3L)
names(res_num)<-res_set

BHiCect_results<-BHiCect(res_set,res_num,chr_dat_l,'smpl.cl',4)
col_chr_dat_l <- make_mres_color_map(res_set,res_set,chr_dat_l)
mres_col_map <- custom_color_map_fn(res_set,'Set1')
cl_heat_dat_tbl <- produce_mres_heat_tbl("100kb_58_1649_37000000_42700000",BHiCect_results,'100kb',c('100kb','50kb','10kb','5kb'),res_num,col_chr_dat_l,4)
base_heatmap(cl_heat_dat_tbl,mres_col_map,res_num,'chr22')
```
