The purpose of this package is perform higher level statistical analyses on RNA-Seq and other datasets.


## Assembling this package
In R:
``` r
devtools::install_github("DanteBortone/housekeeping") # if needed

housekeeping::assemble_package(package_name = "binfotron", my_version = "0.3-19", 
  my_dir = "/datastore/nextgenout5/share/labs/Vincent_Lab/members/dbortone/rstudio-common/packages/binfotron", 
  should_build = FALSE)
```

## Push changes
In bash:
``` bash
cd /datastore/nextgenout5/share/labs/Vincent_Lab/members/dbortone/rstudio-common/packages/binfotron
my_comment="Changed defaulats on gene signatures. Broke up coxph N reporting into two columns and changed the spacing on differential expression gene signatures so I can separate __FDR_pValue__ from __pvalue__ ."
git commit -am "$my_comment"; git push origin binfotron_dev_dsb
git tag -a 0.3-18 -m "$my_comment"; git push -u origin --tags
```

## Install
Restart R
In R (local library, packrat library):
``` r
devtools::install_github("Benjamin-Vincent-Lab/binfotron")
```

Or for a specific version:
``` r
devtools::install_github("Benjamin-Vincent-Lab/binfotron", ref = "0.3-18")
```
