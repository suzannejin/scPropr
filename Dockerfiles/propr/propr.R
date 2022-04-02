# partial correlation packages
pkgs = c(
    "ppcor",
    "corpcor"
)
install.packages(pkgs)

# zero imputation package
install.packages("zCompositions")

# easyCODA for alr reference computation (Greenacre et al 2021)
pkgs = c(
    "vegan",
    "easyCODA"
)
install.packages(pkgs)

# propr package
devtools::install_github("tpq/propr")