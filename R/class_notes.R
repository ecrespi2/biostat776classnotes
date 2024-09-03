#################################################################
##### HOW TO CHECK R VERSIONS AND GIT
##################################################################
print(R.version.string)
print(RStudio.Version()$version)
##in terminal, run: "git --version"

#################################################
##### CLASS NOTES: CLASSES 1-2
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

#################################################
##### CLASS NOTES: LECTURE 6
##################################################
## options for citing rmarkdown
citation("rmarkdown")
print(<citation>, bibtex=TRUE)
toBibtex(.)

## Generate my-refs.bib file
knitr::write_bib("rmarkdown", file = "my-refs.bib")
list.files()

## Add bibliography to Rmarkdown
 # add: "bibliography: my-refs.bib" to Rmarkdown

## Add each citation
  # [@key] for single citation
  # [@key1; @key2] multiple citation can be separated by semi-colon
  # [-@key] in order to suppress author name, and just display the year
  # [see @key1 p 12; also this ref @key2] is also a valid syntax

## Specify citation style
  #Add "csl: biomed-central.csl" to Rmarkdown but pick citation from: https://www.zotero.org/styles

## Add citation to bibliography but not text
  # To Rmarkdown:
      #nocite: |
        #@item1, @item2

## Add everything from bibliography
  #Add "nocite: '@*'" to Rmarkdown file

## SciWheel
  # https://sciwheel.com/?lg
  # online reference management software

#################################################
##### CLASS NOTES: LECTURE 7
##################################################
## Get workiing directory
getwd()

## Set working directory
setwd("/Users/lizziecrespi/Desktop/biostat776classnotes")

## Use relative paths
  #setwd("..\data")

## Using here
library(here)
here::here() ##tells current here directory
here("data","team_standings") ## provides directory for team_standings in data folder

### READING AND WRITING DATA
  #Reading in data from files
    #load(): for reading in single or multiple R objects (opposite of save()) with a .Rda or .RData file format (objects must be same name)
    #readRDS(): for reading in a single object with a .Rds file format (can rename objects)
    #unserialize(): for reading single R objects in binary form
  #Writing data to files
    #save(): for saving an arbitrary number of R objects in binary format (possibly compressed) to a file.
    #saveRDS(): for saving a single object
    #serialize(): for converting an R object into a binary format for outputting to a connection (or file).
    #save.image(): short for ‘save my current workspace’; while this sounds nice, it’s not terribly useful for reproducibility (hence not suggested); it’s also what happens when you try to quit R and it asks if you want to save your work space.

## Saving data
x <- 1:5
save(x, file = here("data", "x.Rda"))
saveRDS(x, file = here("data", "x.Rds"))
list.files(path = here("data"))

## Reading RDS file with either readRDS or load
new_x1 <- readRDS(here("data", "x.Rds"))
new_x2 <- load(here("data", "x.Rda"))

## Reading CSV
standings <- read.csv(here("data","team_standings.csv"))
View(standings)

## Clean up space (remove some files)
file.remove(here("data", "x.Rda"))
file.remove(here("data", "x.Rds"))
rm(x)

## readr package
  # read_csv(): Reads comma-separated file
  # read_csv2(): Reads semicolon-separated file
  # read_tsv(): Reads tab-separated file
  # read_delim(): General function for reading delimited files
  # read_fwf(): Reads fixed width files
  # read_log(): Reads log files

library(readr)
teams <- read_csv(here("data", "team_standings.csv"))
read_csv("The first line of metadata
  The second line of metadata
  x,y,z
  1,2,3",
         skip = 2) #tell it to skip the first two lines
teams <- read_csv(here("data", "team_standings.csv"),
                  col_types = "cc") #specify column types

logs <- read_csv(here("data", "2016-07-19.csv.bz2"),
                 n_max = 10) # reads a gzip-compressed CSV file from the RStudio CRAN mirror.

logs <- read_csv(here("data", "2016-07-19.csv.bz2"),
                 col_types = "ccicccccci",
                 n_max = 10) # reads same file but fix column types

logdates <- read_csv(here("data", "2016-07-19.csv.bz2"),
                     col_types = cols_only(date = col_date()),
                     n_max = 10) #read only first column





