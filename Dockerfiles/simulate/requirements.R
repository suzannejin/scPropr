# install dependencies
pkgs = c('R.utils',
         'pscl',
         'copula',
         'Rtsne', 
         'plyr', 
         'reshape2', 
         'gridExtra',
         'ggplot2',
         'ggpubr',
         'cowplot',
         'argparse',
         'optparse'
         )
install.packages(pkgs)

# install scDesign2
options('download.file.method' = 'curl')   # somehow the recommended method 'libculr' gives me error while using devtools::install_github() on my mac
devtools::install_github("JSB-UCLA/scDesign2@554f2c4b1a7ee6cc04969a287df9b3b77d7bb2fe")