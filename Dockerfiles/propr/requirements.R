# data.table to manage large datasets
pkgs = c("R.utils",
         "data.table"
         )
install.packages(pkgs)

# partial correlation packages
pkgs = c("ppcor",
         "corpcor"
         )
install.packages(pkgs)

# zero imputation package
install.packages("zCompositions")

# normalization packages
BiocManager::install("edgeR")  # for TMM normalization

# propr package
devtools::install_github("tpq/propr")