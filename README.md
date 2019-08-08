The purpose of this package is perform higher level statistical analyses on RNA-Seq and other datasets.


## Assembling this package
In R:
``` r
devtools::install_github("DanteBortone/housekeeping") # if needed

housekeeping::assemble_package(package_name = "binfotron", my_version = "0.2-11",
  my_dir = "/datastore/alldata/shiny-server/rstudio-common/dbortone/packages/binfotron",
  should_build = FALSE)
```

## Push changes
In bash:
``` bash
cd /datastore/alldata/shiny-server/rstudio-common/dbortone/packages/binfotron
my_comment="Bug fixing for differential express and volcano plots."
git commit -am "$my_comment"; git push origin master
git tag -a 0.2-11 -m "$my_comment"; git push -u origin --tags
```

## Install
Restart R
In R (local library, packrat library):
``` r
devtools::install_bitbucket("unc_lineberger/binfotron")
```

Or for a specific version:
``` r
devtools::install_bitbucket("unc_lineberger/binfotron", ref = "0.2-10")
```

