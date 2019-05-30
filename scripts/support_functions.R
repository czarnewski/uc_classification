
#Automatically install packages, if needed
inst_packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) {
    install.packages(new.pkg, dependencies = TRUE, quiet = T)
  if (!new.pkg %in% installed.packages()[, "Package"]){
    Bio.pkg <- new.pkg[!(new.pkg %in% installed.packages()[, "Package"])]
    source("http://bioconductor.org/biocLite.R",prompt.echo = F,echo=F); biocLite(Bio.pkg,suppressUpdates = T,suppressAutoUpdate = T)}}
    data.frame("Pkg_Status"=suppressMessages(sapply(pkg, require, character.only = TRUE)))}



