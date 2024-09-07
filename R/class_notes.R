#################################################################
##### HOW TO CHECK R VERSIONS AND GIT
##################################################################
print(R.version.string)
print(RStudio.Version()$version)
##in terminal, run: "git --version"

#################################################
##### CLASS NOTES: Lectures 1-5 - introduction, git/github, reprex, Rmarkdown, reproducibility
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
##### CLASS NOTES: LECTURE 6 - reference management
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
##### CLASS NOTES: LECTURE 7 - Reading and writing data
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
                     n_max = 10) #read only first


#################################################
##### CLASS NOTES: LECTURE 8 -  Managing data frames with the Tidyverse
##################################################

########     TIBBLES      ########
## Data frames vs. tibbles
  # data frames remove spaces from names and convert them to periods or "x"; tibbles do not
  # tibbles don't have row names
  # tibbles print first ten rows and columns that fit on one screen; data frames can print too much

## Load tidyverse (tibble is part of tidyverse package)
library("tidyverse")
library ("here")

## Read in RDS data as dataframe
chicago <- readRDS(here("data", "chicago.rds"))

## Get dimensions of data frame
dim(chicago)

## Get some info on what hte data frame looks like
str(chicago)

## Get some info on what hte data frame looks like as a tibble
str(as_tibble(chicago))

## Create a tibble from scratch
df <- tibble(
  a = 1:5,
  b = 6:10,
  c = 1,
  z = (a + b)^2 + c
)

## Tibbles can have column names that data frames cannot
tibble(
  `two words` = 1:5,
  `12` = "numeric",
  `:)` = "smile",
)

## Subsetting tibbles
df$z  #lists everything in column z
df[["z"]]  #lists everything in column z
df[[4]]  #lists everything in the fourth column (z)



########     DPLYR    ########
##  Some notes
  # The first argument for DPLYR commnas is always a data fram type object
  # Subsequent arguments describe what to do with the data frame - don't need $
  # Return result is a new data frame
  # Key 'verbs": select(), filter(), arrange(), rename(), mutate(), summarize(), slice_*
  # %>%: the “pipe” operator is used to connect multiple verb actions together into a pipeline

library(dplyr)
chicago <- as_tibble(chicago)
str(chicago)

## SELECT: return a subset of the columns of a data frame, using a flexible notation
subset <- select(chicago, city:dptp)   # select only columns from city-dptp
select(chicago, -(city:dptp))   # select all columns except those from city:dptp
subset <- select(chicago, ends_with("2"))   # select only columns ending in 2
subset <- select(chicago, starts_with("d"))   # select only columns starting with d

## FILTER: extract a subset of rows from a data frame based on logical conditions
chic.f <- filter(chicago, pm25tmean2 > 30)   #keep only records where pm25tmean >30
chic.f <- filter(chicago, pm25tmean2 > 30 & tmpd > 80)   #keep only records where pm25tmean >30 and tmpd > 80

## ARRANGE: reorder rows of a data frame
chicago <- arrange(chicago, date) # sort rows by date (ascending)
chicago <- arrange(chicago, desc(date)) # sort rows by date (descending)

## RENAME: rename variables in a data frame
chicago <- rename(chicago, dewpoint = dptp, pm25 = pm25tmean2) #rename dewpoint column to dptp and pm25 to pm25tmean2


## MUTATE: add new variables/columns or transform existing variables
chicago <- mutate(chicago, pm25detrend = pm25 - mean(pm25, na.rm = TRUE)) #create var pm25detrend which is pm25-mean(pm25)
chicago <- mutate(chicago, year = as.POSIXlt(date)$year + 1900) #extract year from date


## SUMMARIZE: generate summary statistics of different variables in the data frame, possibly within strata
years <- group_by(chicago, year) #use group_by to create dataframe that splits chicago df by year
summarize(years,
          pm25 = mean(pm25, na.rm = TRUE),
          o3 = max(o3tmean2, na.rm = TRUE),
          no2 = median(no2tmean2, na.rm = TRUE)) # summarize mean, max, and median of some variables in years df

## USING %>%
chicago %>%
  mutate(year = as.POSIXlt(date)$year + 1900) %>%
  group_by(year) %>%
  summarize(
    pm25 = mean(pm25, na.rm = TRUE),
    o3 = max(o3tmean2, na.rm = TRUE),
    no2 = median(no2tmean2, na.rm = TRUE)
  )

mutate(chicago, month = as.POSIXlt(date)$mon + 1) %>%
  group_by(month) %>%
  summarize(
    pm25 = mean(pm25, na.rm = TRUE),
    o3 = max(o3tmean2, na.rm = TRUE),
    no2 = median(no2tmean2, na.rm = TRUE)
  ) # compute the average pollutant level by month

## SLICE_*: shows rows of data
slice_sample(chicago, n = 10) # show 10 randomly selected rows
slice_head(chicago, n = 5) # show first 5 rows
slice_tail(chicago, n = 5) # show last 5 rows

## side note about random seelction
set.seed(20240905) #set seed so you will get same numbers every time
rnorm(n=5) # randomly generate 5 numbers



#################################################
##### CLASS NOTES: LECTURE 9 -  Tidy data and the Tidyverse
##################################################

## Load tidyverse, which includes dplyr, tidyr, readr, ggplot2
library(tidyr)
relig_income



#PIVOT_LONGER: turn wide dataset to long
relig_income %>%
  pivot_longer(-religion, names_to = "income", values_to = "respondents") %>%
  mutate(religion = factor(religion), income = factor(income))

relig_income %>%
  pivot_longer(-religion, names_to = "income", values_to = "respondents") # gather everything EXCEPT religion to tidy data

#PIVOT_WIDER: turn long dataset to wide
relig_income %>%
  pivot_longer(-religion, names_to = "income", values_to = "respondents") %>%
  mutate(religion = factor(religion), income = factor(income)) %>%
  group_by(income) %>%
  summarize(total_respondents = sum(respondents)) %>%
  pivot_wider(
    names_from = "income",
    values_from = "total_respondents"
  ) %>%
  knitr::kable() # summarize the total number of respondents per income category




