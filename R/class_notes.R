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


#################################################
##### CLASS NOTES: LECTURE 10 -  Joining data in R
##################################################

## Notes on join functions
library(knitr)
join_funcs <- data.frame(
  func = c(
    "`left_join()`",
    "`right_join()`",
    "`inner_join()`",
    "`full_join()`"
  ),
  does = c(
    "Includes all observations in the left data frame, whether or not there is a match in the right data frame",
    "Includes all observations in the right data frame, whether or not there is a match in the left data frame",
    "Includes only observations that are in both data frames",
    "Includes all observations from both data frames"
  )
)
knitr::kable(join_funcs, col.names = c("Function", "What it includes in merged data frame"))

  # left_join(): Includes all observations in the left data frame, whether or not there is a match in the right data frame
  # right_join(): Includes all observations in the right data frame, whether or not there is a match in the left data frame
  # inner_join(): Includes only observations that are in both data frames
  # full_join(): Includes all observations from both data frames

## Load tidyverse
library(tidyverse)
help(package = "tidyverse")

## Create two tibbles 'outcomes' and 'subjects'
outcomes <- tibble(
  id = rep(c("a", "b", "c"), each = 3),
  visit = rep(0:2, 3),
  outcome = rnorm(3 * 3, 3)
)

subjects <- tibble(
  id = c("a", "b", "c"),
  visit = c(0, 1, 0),
  house = c("detached", "rowhouse", "rowhouse"),
)

print(subjects)
print(outcomes)


## LEFT JOIN
left_join(outcomes, subjects, by = c("id", "visit"))


## Change table
subjects <- tibble(
  id = c("b", "c"),
  visit = c(1, 0),
  house = c("rowhouse", "rowhouse"),
)

subjects


## LEFT JOIN
left_join(x = outcomes, y = subjects, by = c("id", "visit"))


## INNER JOIN
inner_join(x = outcomes, y = subjects, by = c("id", "visit"))


## RIGHT JOIN
right_join(x = outcomes, y = subjects, by = c("id", "visit"))


#################################################
##### CLASS NOTES: LECTURE 13 - ggplot
##################################################
## Recommended books:
  # ggplot2: Elegant Graphics for Data Analysis (3e) https://ggplot2-book.org/
  # R Graphics Cookbook, 2nd edition https://r-graphics.org/

## Load libraries and data
library("tidyverse")
library("here")
maacs <- read_csv(here("data", "bmi_pm25_no2_sim.csv"),
                  col_types = "nnci"
)
maacs


## Create empty plot for logpm25 vs nocturnal symptoms
g <- ggplot(maacs, aes(
  x = logpm25,
  y = NocturnalSympt
))
summary(g)
class(g)
print(g)

## Create empty plot for logpm25 vs nocturnal symptoms in a different way
g <- maacs %>%
  ggplot(aes(logpm25, NocturnalSympt))
print(g)


## Create scatterplot of PM2.5 and days with nocturnal symptoms
g <- maacs %>%
  ggplot(aes(logpm25, NocturnalSympt))
g + geom_point()


## Create scatterplot with smoother of PM2.5 and days with nocturnal symptoms
g +
  geom_point() +
  geom_smooth()


## Create scatterplot with regression line of PM2.5 and days with nocturnal symptoms
g +
  geom_point() +
  geom_smooth(method = "lm")


## Create scatterplot with regression line of PM2.5 and days with nocturnal symptoms
library("palmerpenguins")
glimpse(penguins)

penguins %>%
  ggplot(aes(
    x = flipper_length_mm,
    y = bill_length_mm,
    color = species)) +
  geom_point() +
  geom_smooth()

## Create scatterplot with regression line and facet grid by bmicat of PM2.5 and days with nocturnal symptoms
g +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(. ~ bmicat)


## Change color, size, and transparency of point
g + geom_point(color = "steelblue", size = 4, alpha = 1 / 2)


## Make point color based on bmicat
g + geom_point(aes(color = bmicat), size = 4, alpha = 1 / 2)

## Add smooth line with specified linewidth and type
g +
  geom_point(aes(color = bmicat),
             size = 2,
             alpha = 1 / 2
  ) +
  geom_smooth(
    linewidth = 4,
    linetype = 3,
    method = "lm",
    se = FALSE #specifies if SE is shown
  )

## Specify plot theme
apropos("theme")

g +
  geom_point(aes(color = bmicat)) +
  theme_bw(base_family = "Times")


## Specify plot labels
g +
  geom_point(aes(color = bmicat)) +
  labs(title = "MAACS Cohort") +
  labs(
    x = expression("log " * PM[2.5]),
    y = "Nocturnal Symptoms"
  )


## Specify data frame to work with
testdat <- data.frame(
  x = 1:100,
  y = rnorm(100)
)
testdat[50, 2] <- 100 ## Add an Outlier!
plot(testdat$x,
     testdat$y,
     type = "l",
     ylim = c(-3, 3)
)


## Time series plot with default settings
g <- ggplot(testdat, aes(x = x, y = y))
g + geom_line()


## Time series plot with modified scale on y axis
g +
  geom_line() +
  ylim(-3, 3)


## Time series plot with restricted scale on y axis
g +
  geom_line() +
  coord_cartesian(ylim = c(-3, 3))


## Install bbplot
  ## remotes::install_github("bbc/bbplot")
library("bbplot")

## Basic ggplot2 object with our data
g <- maacs %>%
  ggplot(aes(logpm25, NocturnalSympt))

## A plot we made before, but this time without the SE lines
g +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_grid(. ~ bmicat)

## Now let's add bbplot::bbc_style()
g +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_grid(. ~ bmicat) +
  bbplot::bbc_style()


## Add labels
g +
  geom_point() +
  geom_smooth(colour = "#1380A1", method = "lm", se = FALSE) +
  facet_grid(. ~ bmicat) +
  bbplot::bbc_style() +
  labs(
    title = "Child asthma's link to air quality worsens in overweight children",
    subtitle = "Number of days with symptoms vs PM2.5 by weight group"
  )


## Get rid of points and change colors
g +
  geom_smooth(aes(colour = bmicat), method = "lm", se = FALSE, linewidth = 2) +
  scale_colour_manual(values = c("#FAAB18", "#1380A1")) + ##
  bbplot::bbc_style() +
  labs(
    title = "Child asthma's link to air quality worsens in overweight children",
    subtitle = "Number of days with symptoms vs PM2.5 by weight group"
  )


## Install ThemePark from GitHub
  ## remotes::install_github("MatthewBJane/theme_park")


## Barbie-inspired theme
g +
  geom_smooth(aes(colour = bmicat), method = "lm", se = FALSE, linewidth = 2) +
  scale_colour_manual(values = c("#FAAB18", "#1380A1")) +
  ThemePark::theme_barbie() +
  labs(
    title = "Child asthma's link to air quality worsens in overweight children",
    subtitle = "Number of days with symptoms vs PM2.5 by weight group"
  )

## Oppenheimer-inspired theme
g +
  geom_smooth(aes(colour = bmicat), method = "lm", se = FALSE, linewidth = 2) +
  scale_colour_manual(values = c("#FAAB18", "#1380A1")) +
  ThemePark::theme_oppenheimer() +
  labs(
    title = "Child asthma's link to air quality worsens in overweight children",
    subtitle = "Number of days with symptoms vs PM2.5 by weight group"
  )


## Install ggthemes from CRAN
  ## install.packages("ggthemes")

## Your favorite statistics class theme ;)
## I bet that you could fool a few people into thinking
## that you are not using R ^_^'-- MAKES STATA LOOKING CHART
g +
  geom_smooth(aes(colour = bmicat), method = "lm", se = FALSE, linewidth = 2) +
  scale_colour_manual(values = c("#FAAB18", "#1380A1")) +
  ggthemes::theme_stata() +
  labs(
    title = "Child asthma's link to air quality worsens in overweight children",
    subtitle = "Number of days with symptoms vs PM2.5 by weight group"
  )


## Save our plot into an object
g_complete <- g +
  geom_point(aes(colour = bmicat)) +
  geom_smooth(aes(colour = bmicat), method = "lm", se = FALSE, linewidth = 2) +
  scale_colour_manual(values = c("#FAAB18", "#1380A1"))

## Make it interactive with plotly::ggplotly()
library("plotly")
plotly::ggplotly((g_complete))

## Install colorblindr from GitHub
  ## remotes::install_github("clauswilke/colorblindr")


## Check how grid will look to people with colorblindness
colorblindr::cvd_grid(g_complete)


## Create cutpoints
cutpoints <- quantile(maacs$logno2_new, seq(0, 1, length = 4), na.rm = TRUE)

## create cutpoints and check levels
maacs$no2tert <- cut(maacs$logno2_new, cutpoints)
levels(maacs$no2tert)

## Setup ggplot with data frame
g <- maacs %>%
  ggplot(aes(logpm25, NocturnalSympt))

## Add layers
g + geom_point(alpha = 1 / 3) +
  facet_grid(bmicat ~ no2tert) +
  geom_smooth(method = "lm", se = FALSE, col = "steelblue") +
  theme_bw(base_family = "Avenir", base_size = 10) +
  labs(x = expression("log " * PM[2.5])) +
  labs(y = "Nocturnal Symptoms") +
  labs(title = "MAACS Cohort")




#################################################
##### CLASS NOTES: LECTURE 14 - R nuts and bolts
##################################################


## Assign 5 to vector x and print
x <- 5        ## nothing printed
x             ## auto-printing occurs
print(x)      ## explicit printing
  # The [1] shown in the output indicates that x is a vector and 5 is its first element.
  # Typically with interactive work, we do not explicitly print objects with the print() function;
  # it is much easier to just auto-print them by typing the name of the object and hitting return/enter.
  # However, when writing scripts, functions, or longer programs, there is sometimes
  #a need to explicitly print objects because auto-printing does not work in those settings.


## Assign "hello" to vector msg
msg <- "hello"


## Assign numbers 11 through 30 to vector msg and autoprint
x <- 11:30
x

## Print numbers 5-0 and -15-15
5:0
-15:15

## Use vector function to create vectors (less common)
vector(mode = "numeric", length = 4)     ## makes numeric vector with all 0's
vector(mode = "logical", length = 4)     ## makes logical vector with all FALSE
vector(mode = "character", length = 4)   ## makes character vector with all ""


## Use c() function to create vectors (more common)
x <- c(0.5, 0.6)       ## numeric
x <- c(TRUE, FALSE)    ## logical
x <- c(T, F)           ## logical
x <- c("a", "b", "c")  ## character
x <- 9:29              ## integer
x <- c(1+0i, 2+4i)     ## complex


## Check type of 4 and 4L
typeof(4)    ## double
typeof(4L)   ## integer

## Calculations
x <- sqrt(2) ^ 2
x

## Get letters and check type
letters            ## a-z
typeof(letters)    ## character

## Get 1:10 and check type
1:10            ## 1-10
typeof(1:10)    ## integer


## Store a list with a, b, and 1:10
x <- list("a", "b", 1:10)
x
length(x)           ## check length of x -- 3  (a, b, 1:10)
typeof(x)           ## check type of x -- list
attributes(x)       ## check attributes of x -- NULL
  # Lists (unlike vectors) don't have to be homogenous, which can be useful
x[1]
x[2]
x[3]

## Create and print several vectors; vectors have to be all the same type
y <- c(1.7, "a")       ## R 'coerces' data so both are characters
y
typeof(y)
y <- c(TRUE, 2)        ## R 'coerces' data so both are double
y
typeof(y)
y <- c("a", TRUE)      ## R 'coerces' data so both are characters
y
typeof(y)

## Explicit coercion
x <- 0:6
class(x)           ## vector is integer
as.numeric(x)      ## make vector numeric
as.logical(x)      ## make vector logical
as.character(x)    ## make vector character

x <- c("a", "b", "c")
as.numeric(x)     ## make vector numeric - generates NAs
as.logical(x)     ## make vector numeric - still NAs


## Create empty matrix with 2 rows and 3 columns
m <- matrix(nrow = 2, ncol = 3)
m
dim(m)            ## check matrix dimension (row, column)
attributes(m)     ## check matrix attributes (row, column)

## Create matrix with 2 rows and 3 columns, 1-6
m <- matrix(1:6, nrow = 2, ncol = 3)
m
dim(m)            ## check matrix dimension (row, column)
attributes(m)     ## check matrix attributes (row, column)

## Create matrix with numbers 1-10, auto makes
m <- 1:10     ## makes it a 1x10
m
dim(m) <- c(2, 5)    ## makes it a 2x5
m


## Make two vectors and combine into a matrix (have to be same type)
x <- 1:3
y <- 10:12
cbind(x, y)   ## makes two columns, 1-3 and 10-12
rbind(x, y)   ## makes two rows, 1-3 and 10-12

## Make another list
x <- list(1, "a", TRUE, 1 + 4i)
x

## Make a vector with five empty lists
x <- vector("list", length = 5)
x


## Make a factor which allows underlying numeric coding with labels
x <- factor(c("yes", "yes", "no", "yes", "no"))
x
table(x)
unclass(x)    ## See the underlying representation of factor

## ------------------------------------------------------------------------------------------------------------------------
x <- factor(c("yes", "yes", "no", "yes", "no"))
x  ## Levels are put in alphabetical order
x <- factor(c("yes", "yes", "no", "yes", "no"),
            levels = c("yes", "no"))
x

## Create a vector with NAs in it
x <- c(1, 2, NA, 10, 3)
is.na(x)         ## Return a logical vector indicating which elements are NA
is.nan(x)        ## Return a logical vector indicating which elements are NaN

## Now create a vector with both NA and NaN values
x <- c(1, 2, NaN, NA, 4)
is.na(x)
is.nan(x)


## Make data frame with varnames foo and bar
x <- data.frame(foo = 1:4, bar = c(T, T, F, F))
x
nrow(x)
ncol(x)
attributes(x)

## Convert data frame to matrix
data.matrix(x)
attributes(data.matrix(x))


## Load packages
library(tidyverse)
library(palmerpenguins)
penguins

## make vector 1:3 and then add names to vector
x <- 1:3
names(x)
names(x) <- c("New York", "Seattle", "Los Angeles")
x
names(x)
attributes(x)


## Make list of numbers 1:3 with names
x <- list("Los Angeles" = 1, Boston = 2, London = 3)
x
names(x)


## Make 2x2 matrix with numbers 1-4, add row and column names
m <- matrix(1:4, nrow = 2, ncol = 2)
m
dimnames(m) <- list(c("a", "b"), c("c", "d"))
m


## Change row and column names for matrix
colnames(m) <- c("h", "f")
rownames(m) <- c("x", "z")
m

##################################################
##### CLASS NOTES: LECTURE 15 - Control Structures
##################################################

## Draw random value between min and max
x <- runif(n = 1, min = 0, max = 10)
x


## check if x > 3
x > 3


## If x> 3, first condition occurs. Else, second condition occurs.

if (x > 3) {
  y <- 10
} else {
  y <- 0
}

y

## Different but equivalent way of doing if else
y <- if (x > 3) {
  10
} else {
  0
}

y

## Load packages
library(tidyverse)
library(palmerpenguins)
penguins

## Print numbers 1-10
for (i in 1:10) {
  print(i)
}

for (i in 1:10) {
  i
}


## Define a loop to iterate over and print each element
x <- c("a", "b", "c", "d")

for (i in 1:4) {
  ## Print out each element of 'x'
  print(x[i])
}

## seq_along generates integer sequence based on length of object
x
seq_along(x)

## Generate a sequence based on length of 'x'
for (i in seq_along(x)) {
  print(x[i])
}


## Can name i anything
for (babyshark in x) {
  print(babyshark)
}

## Can't use numbers though
  # for (1999 in x) {
  #  print(1999)
  # }


## Technically don't need curly braces, but use them anyway for simplicity!
for (i in 1:4) print(x[i])

## Create 2x3 matrix with numbers 1-6
x <- matrix(1:6, nrow = 2, ncol = 3)
x

## Create nested for loop: does something for every row and then every column in row
for (i in seq_len(nrow(x))) {
  for (j in seq_len(ncol(x))) {
    print(x[i, j])
  }
}

  # lots of loops can be confusing
  # can break them up with functions if desired

## While is rarely used;
  # specify value outside loop
  # have condition that uses that value
  # update something based on that

count <- 0
while (count < 10) {
  print(count)
  count <- count + 1
}


## More complicated while loop
z <- 5
set.seed(1)

while (z >= 3 && z <= 10) {
  coin <- rbinom(1, 1, 0.5)

  if (coin == 1) { ## random walk
    z <- z + 1
  } else {
    z <- z - 1
  }
}
print(z)


## Difference between one & and two &&
-2:2
((-2:2) >= 0) & ((-2:2) <= 0)

(2 >= 0) && (-2 <= 0)
(-2 >= 0) && (-2 <= 0)


## Repeat loops
  # Be careful as can create never ending loop
  #| eval: false
  ## x0 <- 1
  ## tol <- 1e-8
  ##
  ## repeat {
  ##     x1 <- computeEstimate()
  ##
  ##     if (abs(x1 - x0) < tol) { ## Close enough?
  ##         break
  ##     } else {
  ##         x0 <- x1
  ##     }
  ## }


## Using next in loops
  ## for (i in 1:100) {
  ##     if (i <= 20) {
  ##         ## Skip the first 20 iterations
  ##         next
  ##     }
  ##     ## Do something here
  ## }


## Using break in loops
  ## for (i in 1:100) {
  ##     print(i)
  ##
  ##     if (i > 20) {
  ##         ## Stop loop after 20 iterations
  ##         break
  ##     }
  ## }



##################################################
##### CLASS NOTES: LECTURE 16 - Functions
##################################################
## Create empty function and check class
f <- function() {
  ## This is an empty function
}

class(f)
f()


## Create function that prints "Hello world!"
f <- function() {
  # this is the function body
  hello <- "Hello, world!\n"
  cat(hello)
}
f()


## cat() is sometimes prefereable to print()
  # It doesn't output new lines
hello <- "Hello, world!\n"
print(hello)
cat(hello)


## Create function that prints Hello World the number of times specified
f <- function(num) {
  for (i in seq_len(num)) {
    hello <- "Hello, world!\n"
    cat(hello)
  }
}
f(3)


## Create function that prints Hellow world specified number of times and prints number of letters in all "hello"
f <- function(num) {
  hello <- "Hello, world!\n"
  for (i in seq_len(num)) {
    cat(hello)
  }
  chars <- nchar(hello) * num
  chars
}
meaningoflife <- f(3)
print(meaningoflife)

## Make same function but set default to 1 so if they don't enter anything you dont get an error message
f <- function(num = 1) {
  hello <- "Hello, world!\n"
  for (i in seq_len(num)) {
    cat(hello)
  }
  chars <- nchar(hello) * num
  chars
}

f() ## Use default value for 'num'
f(2) ## Use user-specified value


## Returns list of all formal arguments of a function
formals(f)


## Can specify function argument name if desired
f(num = 2)


## Figure out what arguments are in rnorm function & the order
str(rnorm)
mydata <- rnorm(100, 2, 1) ## Generate some data


## Play with SD function
sd(mydata)                    ## Positional match first argument, default for 'na.rm'
sd(x = mydata)                ## Specify 'x' argument by name, default for 'na.rm'
sd(x = mydata, na.rm = FALSE) ## Specify both arguments by name
sd(na.rm = FALSE, x = mydata) ## Can change position if you specify argument name
sd(na.rm = FALSE, mydata)     ## Can mix specifying arguments and not


## Check arguments in function f
args(f)


## Check arguments in function lm
args(lm)
lm(data = mydata, y ~ x, model = FALSE, 1:100)
lm(y ~ x, mydata, 1:100, model = FALSE)

## Create function which returns square of first argument
f <- function(a, b) {
  a^2
}
f(2)


## Error because no function
f <- function(a, b) {
  print(a)
  print(b)
}
f(45)


## Try to call mean function without argument
mean

## Paste function combines strings
args(paste)

paste("one", "two", "three")
paste("one", "two", "three", "four", "five", sep = "_")
paste("a", "b", sep = ":")
paste("a", "b", se = ":")

## Examples of function names
    # f()
        # Too short
    # my_awesome_function()
       # Not a verb, or descriptive
    ## impute_missing()
      # Long, but clear
    ## collapse_years()
      # Long, but clear

## snake_case vs camelCase -- PICK ONE AND STICK WITH IT
col_mins <- function(x, y) {}
rowMaxes <- function(x, y) {}   # DO NOT MIX CASES

## Good function names - help with autocompletion
  # input_select()
  # input_checkbox()
  # input_text()

## Not so good function names
  # select_input()
  # checkbox_input()
  # text_input()


# Don't do this! -- Avoid overriding existing functions
  # T <- FALSE
  # c <- 10
  # mean <- function(x) sum(x)


## Uses y set from outside function, keep this in mind!
f <- function(x) {
  x + y
}

y <- 100
f(10)

y <- 1000
f(10)


## IDK what this is doing
`+` <- function(x, y) {
  if (runif(1) < 0.1) {
    sum(x, y)
  } else {
    sum(x, y) * 1.1
  }
}
table(replicate(1000, 1 + 2))


## Remove `+` function
rm(`+`)


## ------------------------------------------------------------------------------------------------------------------------
#| eval: false
## mean(is.na(x))
##
## x / sum(x, na.rm = TRUE)




##################################################
##### CLASS NOTES: LECTURE 17 - Vectorization and loop functionals
##################################################

## Create vector inches with specific numbers
inches <- c(69, 62, 66, 70, 70, 73, 67, 73, 67, 70)


## Print vector inches * 2.54
inches * 2.54


## Print vector inches - 69
inches - 69


## Create vectors x and y 1-10 and add (x + y = 121)
x <- 1:10
y <- 1:10
x + y

## Print sqrt of vector x
sqrt(x)

## Print vector x*y
x * y


## In a lot of other languages, we would have to do:
  ## Check that x and y have the same length
  stopifnot(length(x) == length(y))

  ## Create our result object
  result <- vector(mode = "integer", length = length(x))

  ## Loop through each element of x and y, calculate the sum,

  ## then store it on 'result'
  for (i in seq_along(x)) {
    result[i] <- x[i] + y[i]
  }
  ## Check that we got the same answer
  identical(result, x + y)


## Create sum function
my_sum <- function(a, b) {
  a + b
}
my_sum (4, 100)


## Same but with an extra check to make sure that 'a' and 'b'
## have the same lengths.
my_sum <- function(a, b) {
  ## Check that a and b are of the same length
  stopifnot(length(a) == length(b))
  a + b
}

my_sum(c(1, 3, 4), c(1, 4, 5))
my_sum(c(1, 3, 4), c(1, 4))

## Documenting function with roxygen2 syntax
#' Title
#'
#' @param a
#' @param b
#'
#' @return
#' @export
#'
#' @examples
my_sum <- function(a, b) {
  ## Check that a and b are of the same length
  stopifnot(length(a) == length(b))
  a + b
}


## Documenting function with roxygen2 syntax
#' Title
#'
#' Description
#'
#' Details
#'
#' @param a What is `a`?
#' @param b What is `b`?
#'
#' @return What does the function return?
#' @export ## Do we want to share this function? yes!
#'
#' @examples
#' ## How do you use this function?
my_sum <- function(a, b) {
  ## Check that a and b are of the same length
  stopifnot(length(a) == length(b))
  a + b
}


## Documenting function with roxygen2 syntax
#' Sum two vectors
#'
#' This function does the element wise sum of two vectors.
#'
#' It really is just an example function that is powered by the `+` operator
#' from [base::Arithmetic].
#'
#' @param a An `integer()` or `numeric()` vector of length `L`.
#' @param b An `integer()` or `numeric()` vector of length `L`.
#'
#' @return An `integer()` or `numeric()` vector of length `L` with
#' the element-wise sum of `a` and `b`.
#' @export
#'
#' @examples
#' ## Generate some input data
#' x <- 1:10
#' y <- 1:10
#'
#' ## Perform the element wise sum
#' my_sum(x, y)
my_sum <- function(a, b) {
  ## Check that a and b are of the same length
  stopifnot(length(a) == length(b))
  a + b
}


## CREATING PACKAGES
library("testthat")
test_that("my_sum works", {
  x <- seq_len(10)
  expect_equal(my_sum(x, x), x + x)

  expect_error(my_sum(x, seq_len(5)))
})


## CREATING PACKAGES
#| eval: false
## ## Install biocthis if you don't have it
## if (!require("BiocManager", quietly = TRUE)) {
##     install.packages("BiocManager")
## }
##
## BiocManager::install("biocthis")
##
## ## Create an empty R package that is also an
## ## RStudio project
## usethis::create_package("~/Desktop/sum776")
##
## ## On the new RStudio window, create the
## ## scripts that will guide you into making a package
## biocthis::use_bioc_pkg_templates()


## CREATING PACKAGES
remotes::install_github("lcolladotor/sum776")


## MAPPLY()
## Check the arguments to mapply()
args(mapply)

## Apply mapply() to our function my_sum() with the inputs 'x' and 'y'
mapply(sum776::my_sum, x, y)

## Or write an anynymous function that is:
## * not documented
## * not tested
## * not shared
##
## :(
mapply(function(a, b) {
  a + b
}, x, y)


## purr (alternative to mapply)
library("purrr") ## part of tidyverse

## Check the arguments of map2_int()
args(purrr::map2_int)

## Apply our function my_sum() to our inputs
purrr::map2_int(x, y, sum776::my_sum)

## You can also use anonymous functions
purrr::map2_int(x, y, function(a, b) {
  a + b
})

## purrr even has a super short formula-like syntax
## where .x is the first input and .y is the second one
purrr::map2_int(x, y, ~ .x + .y)

## This formula syntax has nothing to do with the objects 'x' and 'y'
purrr::map2_int(1:2, 3:4, ~ .x + .y)


## lapply (another alternative but only uses lists)
lapply
x <- list(a = 1:5, b = rnorm(10))
x
lapply(x, mean) ##return mean of each list in list x


## purr returns dbl
purrr::map_dbl(x, mean)

## Return mean of each list in list x
x <- list(a = 1:4, b = rnorm(10), c = rnorm(20, 1), d = rnorm(100, 5))
lapply(x, mean)


## For each in x, return list???
x <- 1:4
lapply(x, runif)


## ??
purrr::map(x, runif)


## Can specify addition arguments in function after function
x <- 1:4
lapply(x, runif, min = 0, max = 10)
purrr::map(x, runif, min = 0, max = 10)


## Create list of matrices
x <- list(a = matrix(1:4, 2, 2), b = matrix(1:6, 3, 2))
x
apply(x, function(elt) {
  elt[, 1]
})


## ------------------------------------------------------------------------------------------------------------------------
f <- function(elt) {
  elt[, 1]
}
lapply(x, f)


## ------------------------------------------------------------------------------------------------------------------------
x <- list(a = 1:4, b = rnorm(10), c = rnorm(20, 1), d = rnorm(100, 5))
lapply(x, mean)


## sapply is simplified version of lapply
sapply(x, mean)
purrr::map(x, mean)
purrr::map_dbl(x, mean)


## split() takes vector or other object and splits it into groups
str(split)

x <- c(rnorm(10), runif(10), rnorm(10, 1))
f <- gl(3, 10) # generate factor levels
f
split(x, f) ## get list with name of each element = factor levels, and values x


## Get list with name of each element = factor levels, and values mean(x)
lapply(split(x, f), mean)


## Load dataset
library("datasets")
head(airquality)


## Split data frame by month; creates list where each element is named as month and contains a dataframe
s <- split(airquality, airquality$Month)
str(s)


## Get column means for Ozone, Solar.R and Wind in list s
lapply(s, function(x) {
  colMeans(x[, c("Ozone", "Solar.R", "Wind")])
})


## Get column means for Ozone, Solar.R and Wind in list s
sapply(s, function(x) {
  colMeans(x[, c("Ozone", "Solar.R", "Wind")])
})


## Get column means for Ozone, Solar.R and Wind in list s, exclude NA
sapply(s, function(x) {
  colMeans(x[, c("Ozone", "Solar.R", "Wind")],
           na.rm = TRUE
  )
})


## Get column means for Ozone, Solar.R and Wind in list s
purrr::map(s, function(x) {
  colMeans(x[, c("Ozone", "Solar.R", "Wind")], na.rm = TRUE)
})


## Get column means for Ozone, Solar.R and Wind in list s
purrr::map_dfc(s, function(x) {
  colMeans(x[, c("Ozone", "Solar.R", "Wind")], na.rm = TRUE)
})


## ------------------------------------------------------------------------------------------------------------------------
## Make sure we get data.frame / tibble outputs for each element
## of the list
purrr:::map(s, function(x) {
  tibble::as_tibble(colMeans(x[, c("Ozone", "Solar.R", "Wind")], na.rm = TRUE))
})

## Now we can combine them with list_cbind()
purrr:::map(s, function(x) {
  tibble::as_tibble(colMeans(x[, c("Ozone", "Solar.R", "Wind")], na.rm = TRUE))
}) %>% purrr::list_cbind()

## And we can then add the actual variable it came from with mutate()
purrr:::map(s, function(x) {
  tibble::as_tibble(colMeans(x[, c("Ozone", "Solar.R", "Wind")], na.rm = TRUE))
}) %>%
  purrr::list_cbind() %>%
  dplyr::mutate(Variable = c("Ozone", "Solar.R", "Wind"))


## ------------------------------------------------------------------------------------------------------------------------
## Sadly map_dfr() is now superseded (aka not recommended)
purrr:::map_dfr(s, function(x) {
  colMeans(x[, c("Ozone", "Solar.R", "Wind")], na.rm = TRUE)
})

## This is how we would have added back the Month variable
purrr:::map_dfr(s, function(x) {
  colMeans(x[, c("Ozone", "Solar.R", "Wind")], na.rm = TRUE)
}) %>%
  dplyr:::mutate(Month = as.integer(names(s)))


## ------------------------------------------------------------------------------------------------------------------------
## Get data.frame / tibble outputs, but with each variable as a separate
## column. Here we used the t() or transpose() function.
purrr:::map(s, function(x) {
  tibble::as_tibble(t(colMeans(x[, c("Ozone", "Solar.R", "Wind")], na.rm = TRUE)))
})

## Now we can row bind each of these data.frames / tibbles into a
## single one
purrr:::map(s, function(x) {
  tibble::as_tibble(t(colMeans(x[, c("Ozone", "Solar.R", "Wind")], na.rm = TRUE)))
}) %>% purrr::list_rbind()

## Then with mutate, we can add the Month back
purrr:::map(s, function(x) {
  tibble::as_tibble(t(colMeans(x[, c("Ozone", "Solar.R", "Wind")], na.rm = TRUE)))
}) %>%
  purrr::list_rbind() %>%
  dplyr:::mutate(Month = as.integer(names(s)))


## ------------------------------------------------------------------------------------------------------------------------
## group_by() is in a way splitting our input data.frame / tibble by
## our variable of interest. Then summarize() helps us specify how we
## want to use that data, before it's all put back together into a
## tidy tibble.
airquality %>%
  dplyr::group_by(Month) %>%
  dplyr::summarize(
    Ozone = mean(Ozone, na.rm = TRUE),
    Solar.R = mean(Solar.R, na.rm = TRUE),
    Wind = mean(Wind, na.rm = TRUE)
  )


## str(tapply): apply a function over subsets of a vector


## ------------------------------------------------------------------------------------------------------------------------
## Simulate some data
x <- c(rnorm(10), runif(10), rnorm(10, 1))
## Define some groups with a factor variable
f <- gl(3, 10)
f
tapply(x, f, mean)


## ------------------------------------------------------------------------------------------------------------------------
tapply(x, f, range)


## ------------------------------------------------------------------------------------------------------------------------
split(x, f) %>% purrr::map_dbl(mean)
split(x, f) %>% purrr::map(range)


## str(tapply): apply a function over the margins of an array
str(apply)


## ------------------------------------------------------------------------------------------------------------------------
x <- matrix(rnorm(200), 20, 10)
head(x)
apply(x, 2, mean) ## Take the mean of each column


## ------------------------------------------------------------------------------------------------------------------------
apply(x, 1, sum) ## Take the mean of each row


## ------------------------------------------------------------------------------------------------------------------------
#| eval: false
## apply(x, 2, mean)


## ------------------------------------------------------------------------------------------------------------------------
#| eval: false
## apply(x, 1, sum)


## ------------------------------------------------------------------------------------------------------------------------
array_branch(x, 2) %>% map_dbl(mean)
array_branch(x, 1) %>% map_dbl(sum)


## ------------------------------------------------------------------------------------------------------------------------
x <- matrix(rnorm(200), 20, 10)
head(x)

## Get row quantiles
apply(x, 1, quantile, probs = c(0.25, 0.75))


## ------------------------------------------------------------------------------------------------------------------------
array_branch(x, 1) %>%
  map(quantile, probs = c(0.25, 0.75)) %>%
  map(~ as.data.frame(t(.x))) %>%
  list_rbind()


## ------------------------------------------------------------------------------------------------------------------------
sumsq <- function(mu, sigma, x) {
  sum(((x - mu) / sigma)^2)
}


## ------------------------------------------------------------------------------------------------------------------------
x <- rnorm(100) ## Generate some data
sumsq(mu = 1, sigma = 1, x) ## This works (returns one value)


## ------------------------------------------------------------------------------------------------------------------------
sumsq(1:10, 1:10, x) ## This is not what we want


## ------------------------------------------------------------------------------------------------------------------------
vsumsq <- Vectorize(sumsq, c("mu", "sigma"))
vsumsq(1:10, 1:10, x)

## The details are a bit complicated though
## as we can see below
vsumsq



#################################################
##### CLASS NOTES: LECTURE 18 - Debugging R Code
##################################################

## ------------------------------------------------------------------------------------------------------------------------
#| eval: false
## remotes::install_github("jalvesaq/colorout")


## ------------------------------------------------------------------------------------------------------------------------
#| eval: false

## ## Open your .Rprofile file
## usethis::edit_r_profile()
##
## ## Copy paste the following code taken from
## ## https://lcolladotor.github.io/bioc_team_ds/config-files.html#rprofile
##
## ## Change colors
## # Source https://github.com/jalvesaq/colorout
## if (Sys.getenv("TERM") %in% c("term", "xterm-256color", "cygwin", "screen")) {
##     if (!requireNamespace("colorout", quietly = TRUE) & .Platform$OS.type != "windows") {
##         cat('To install colorout use: remotes::install_github("jalvesaq/colorout")\n')
##     }
## }


## ------------------------------------------------------------------------------------------------------------------------
#| warning: true
#| error: true
require("colorout")

## From colorout's README documentation
x <- data.frame(
  logic = c(TRUE, TRUE, FALSE),
  factor = factor(c("abc", "def", "ghi")),
  string = c("ABC", "DEF", "GHI"),
  real = c(1.23, -4.56, 7.89),
  cien.not = c(1.234e-23, -4.56 + 45, 7.89e78),
  date = as.Date(c("2012-02-21", "2013-02-12", "2014-03-04"))
)
rownames(x) <- seq_len(3)
x

summary(x[, c(1, 2, 4, 6)])

warning("This is an example of a warning.")

example.of.error

library("KernSmooth")

colorout::setOutputColors()


## ------------------------------------------------------------------------------------------------------------------------
library("reprex")


## ------------------------------------------------------------------------------------------------------------------------
#| eval: false
## (y <- 1:4)
## mean(y)


## ------------------------------------------------------------------------------------------------------------------------
#| warning: true
log(-1)


## ------------------------------------------------------------------------------------------------------------------------
print_message <- function(x) {
  if (x > 0) {
    print("x is greater than zero")
  } else {
    print("x is less than or equal to zero")
  }
  invisible(x)
}


## ------------------------------------------------------------------------------------------------------------------------
#| error: true
print_message(1)


## ------------------------------------------------------------------------------------------------------------------------
#| error: true
print_message(NA)


## ------------------------------------------------------------------------------------------------------------------------
print_message2 <- function(x) {
  if (is.na(x)) {
    print("x is a missing value!")
  } else if (x > 0) {
    print("x is greater than zero")
  } else {
    print("x is less than or equal to zero")
  }
  invisible(x)
}


## ------------------------------------------------------------------------------------------------------------------------
print_message2(NA)


## ------------------------------------------------------------------------------------------------------------------------
#| error: true
x <- log(c(-1, 2))
print_message2(x)


## ------------------------------------------------------------------------------------------------------------------------
print_message3 <- function(x) {
  if (length(x) > 1L) {
    stop("'x' has length > 1")
  }
  if (is.na(x)) {
    print("x is a missing value!")
  } else if (x > 0) {
    print("x is greater than zero")
  } else {
    print("x is less than or equal to zero")
  }
  invisible(x)
}


## ------------------------------------------------------------------------------------------------------------------------
#| error: true
print_message3(1:2)


## ------------------------------------------------------------------------------------------------------------------------
print_message3_no_call <- function(x) {
  if (length(x) > 1L) {
    stop("'x' has length > 1", call. = FALSE)
  }
  if (is.na(x)) {
    print("x is a missing value!")
  } else if (x > 0) {
    print("x is greater than zero")
  } else {
    print("x is less than or equal to zero")
  }
  invisible(x)
}


## ------------------------------------------------------------------------------------------------------------------------
#| error: true
print_message3_no_call(99:100)
print_message3(99:100)


## ------------------------------------------------------------------------------------------------------------------------
print_message3_tidyverse <- function(x) {
  if (length(x) > 1L) {
    rlang::abort("'x' has length > 1")
  }
  if (is.na(x)) {
    rlang::warn("x is a missing value!")
  } else if (x > 0) {
    rlang::inform("x is greater than zero")
  } else {
    rlang::inform("x is less than or equal to zero")
  }
  invisible(x)
}


## ------------------------------------------------------------------------------------------------------------------------
#| error: true
print_message3_tidyverse(99:100)
print_message3_tidyverse(NA)
print_message3_tidyverse(1)
print_message3_tidyverse(-1)


## ------------------------------------------------------------------------------------------------------------------------
print_message3_cli <- function(x) {
  if (length(x) > 1L) {
    len <- length(x)

    ## Avoid the print() calls from
    ## https://github.com/ComunidadBioInfo/praiseMX/blob/master/R/praise_crear_emi.R
    praise_mx_log <- capture.output({
      praise_mx <- praiseMX:::praise_bien()
    })
    cli::cli_abort(
      c(
        "This function is not vectorized:",
        "i" = "{.var x} has length {len}.",
        "x" = "{.var x} must have length 1.",
        ">" = "Try using {.code purrr::map(x, print_message3_cli)} to loop your input {.var x} on this function.",
        "v" = praise::praise(),
        "v" = praise_mx
      )
    )
  }
  if (is.na(x)) {
    rlang::warn("x is a missing value!")
  } else if (x > 0) {
    rlang::inform("x is greater than zero")
  } else {
    rlang::inform("x is less than or equal to zero")
  }
  invisible(x)
}


## ------------------------------------------------------------------------------------------------------------------------
#| error: true
set.seed(20230928)
print_message3_cli(-1:1)
purrr::map(-1:1, print_message3_cli)


## ------------------------------------------------------------------------------------------------------------------------
print_message4 <- Vectorize(print_message2)
out <- print_message4(c(-1, 2))


## ------------------------------------------------------------------------------------------------------------------------
#| error: true
lm(y ~ x)
rlang::last_error()


## ----error=TRUE----------------------------------------------------------------------------------------------------------
f <- function(a) g(a)
g <- function(b) h(b)
h <- function(c) i(c)
i <- function(d) {
  if (!is.numeric(d)) {
    stop("`d` must be numeric", call. = FALSE)
  }
  d + 10
}
f("a")


## ------------------------------------------------------------------------------------------------------------------------
options(width = 120)
sessioninfo::session_info()


## ------------------------------------------------------------------------------------------------------------------------
options(width = 120)
sessioninfo::session_info()
