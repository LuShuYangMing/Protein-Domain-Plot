# Protein-Domain-Plot
Plot protein domain via querying Uniprot Database based on R

#### Package Versions and R sessionInfo
```
R version 3.5.1 (2018-07-02)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows >= 8 x64 (build 9200)

Matrix products: default

locale:
[1] LC_COLLATE=Chinese (Simplified)_China.936  LC_CTYPE=Chinese (Simplified)_China.936   
[3] LC_MONETARY=Chinese (Simplified)_China.936 LC_NUMERIC=C                              
[5] LC_TIME=Chinese (Simplified)_China.936    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] readxl_1.2.0         ggthemes_4.0.1       tidyr_0.8.2          ggplot2_3.1.0       
[5] rlang_0.3.1          dplyr_0.7.8          httr_1.4.0           RevoUtils_11.0.1    
[9] RevoUtilsMath_11.0.0

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.0       pillar_1.3.1     compiler_3.5.1   cellranger_1.1.0 plyr_1.8.4      
 [6] bindr_0.1.1      tools_3.5.1      digest_0.6.18    jsonlite_1.6     tibble_2.0.1    
[11] gtable_0.2.0     pkgconfig_2.0.2  rstudioapi_0.9.0 curl_3.3         yaml_2.2.0      
[16] bindrcpp_0.2.2   withr_2.1.2      stringr_1.3.1    grid_3.5.1       tidyselect_0.2.5
[21] glue_1.3.0       R6_2.3.0         purrr_0.2.5      magrittr_1.5     scales_1.0.0    
[26] assertthat_0.2.0 colorspace_1.3-2 labeling_0.3     stringi_1.2.4    lazyeval_0.2.1  
[31] munsell_0.5.0    crayon_1.3.4    
```

#### Instructions for Use 

Please modify the set.wd(path_to_gmtv-Rv1468c.xlsx) in snp_plot.R, then run the code within R

#### Reproducibility
The expected output is as shown in output.pdf, and the labels and colors could be further modified by Adobe Illustrator
