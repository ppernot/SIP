# Install packages if necessary
## CRAN packages
libs <- c(
  "xtable","mixtools","inlmisc","rlist","boot","corrplot",
  "repmis","sessioninfo","ggplot2","dplyr","grid","WRS2","WRS",
  "distillery"
)
for (lib in libs) {
  if (!require(lib, character.only = TRUE, quietly = TRUE)) {
    install.packages(
      lib,
      dependencies = TRUE,
      repos = "https://cran.univ-paris1.fr"
    )
  }
}
## Load packages and generate biblio
repmis::LoadandCite(libs,file='../article/packages.bib')

## Github package
lib = "ErrViewLib"
if(!require(lib,character.only = TRUE))
  devtools::install_github(paste0("ppernot/",lib))
library(lib,character.only = TRUE)

# Parallel options for bootstrap
options(boot.parallel = "multicore")
options(boot.ncpus = 4)

# Set random seed for reprod
set.seed(123)

# Set graphical params
gPars = list(
  cols     = rev(inlmisc::GetColors(8))[1:7],
  cols_tr  = rev(inlmisc::GetColors(8, alpha = 0.2))[1:7],
  cols_tr2 = rev(inlmisc::GetColors(8, alpha = 0.5))[1:7],
  pty      = 's',
  mar      = c(3,3,1.6,.2),
  mgp      = c(2,.75,0),
  tcl      = -0.5,
  lwd      = 4.0,
  cex      = 4.0,
  cex.leg  = 0.7,
  reso     = 1200  # (px) base resolution for png figs
)
# Expose gPars list
for (n in names(gPars))
  assign(n, rlist::list.extract(gPars, n))

# Define Data and Results repositories
dataRepo = '../data/'
figRepo  = '../results/figs/'
tabRepo  = '../results/tables/'

sink(file ='./sessionInfo.txt')
print(sessioninfo::session_info())
sink()
