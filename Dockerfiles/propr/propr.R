# partial correlation packages
pkgs = c(
    "ppcor",
    "corpcor"
)
install.packages(pkgs)

# zero imputation package
# install.packages("zCompositions")
devtools::install_github('cran/Zcompositions@1a329810b3f5c773ccca63628fc6fbbf31f38233')

# easyCODA for alr reference computation (Greenacre et al 2021)
pkgs = c(
    "vegan",
    "easyCODA"
)
install.packages(pkgs)

# propr package
devtools::install_github("tpq/propr")