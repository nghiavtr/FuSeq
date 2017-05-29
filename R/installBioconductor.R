#####
#Install some core packages of bioconductor for my docker
#These codes are not mine (Nghia), I got them from https://github.com/Bioconductor/bioc_docker/tree/master/src/core
#Note: must install libcurl4-openssl-dev libxml2-dev in ubuntu (if using) before insall GenomicFeatures because two dependencies R packages curl and xml require them
############################################################
source("http://bioconductor.org/biocLite.R")# This link works - Nghia
biocLite()
pkgs <- c(
  
  "GenomicFeatures"
)

ap.db <- available.packages(contrib.url(biocinstallRepos()))
ap <- rownames(ap.db)

pkgs_to_install <- pkgs[pkgs %in% ap]

biocLite(pkgs_to_install)

warnings()

if (!is.null(warnings()))
{
  w <- capture.output(warnings())
  if (length(grep("is not available|had non-zero exit status", w)))
    quit("no", 1L)
}