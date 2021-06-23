The purpose of this package is perform higher level statistical analyses on RNA-Seq and other datasets.


## Assembling this package
In R:
``` r
devtools::install_github("DanteBortone/housekeeping") # if needed

housekeeping::assemble_package(package_name = "binfotron", my_version = "0.5-01", 
  my_dir = "/datastore/nextgenout5/share/labs/Vincent_Lab/members/dbortone/rstudio-common/packages/binfotron", 
  should_build = FALSE)
```

## Push changes
In bash:
``` bash
cd /datastore/nextgenout5/share/labs/Vincent_Lab/members/dbortone/rstudio-common/packages/binfotron
my_comment="Updated readmes on gmts to include the version of the repo that made them."
git add .
git commit -am "$my_comment"; git push
git tag -a 0.5-01 -m "$my_comment"; git push -u origin --tags
```

## Install
Restart R
In R (local library, packrat library):
``` r
devtools::install_github("Benjamin-Vincent-Lab/binfotron")
```

Or for a specific version:
``` r
devtools::install_github("Benjamin-Vincent-Lab/binfotron", ref = "0.5-00")
```
