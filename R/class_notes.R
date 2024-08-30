#################################################################
##### HOW TO CHECK R VERSIONS AND GIT
##################################################################
print(R.version.string)
print(RStudio.Version()$version)
##in terminal, run: "git --version"

#################################################
##### CLASS NOTES
##################################################
## Create an Rstudio project
  #usethis::create_project("~/Desktop/biostat776classnotes")

## Start version controlling it
  #usethis::use_git()

## Share it via GitHub with the world
  #usethis::use_github()


## Below is an example for installing more than one package at a time:
install.packages(
  c("postcards", "usethis", "gitcreds", "here", "sessioninfo", "remotes")
  )

## installing remotes and using it to download packages from github
library(remotes)
remotes::install_github()


##using reprex to share code that has errors for help
stop("This R error is weird")
  #copy code with error
  #type reprex::reprex() in R console and run
  #clipboard has error code
  # add this in comment on Github

##using here for relative directories instead of firm ones
here::here() # check current 'here' directory

## I can now easily share code to access files from this project
## such as access to the flight.csv file saved under the data
## directory.
here::here("data", "flights.csv")


## Reproducibility information
library(sessioninfo)
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()


##markdown -- allows us to have a source doc with text + code (essentially a publication)
  #weave (human readable docs) vs tangle (machine readable text)
  #to create: File --> new file --> R Markdown








