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


## Create function to compute sum of squares
sumsq <- function(mu, sigma, x) {
  sum(((x - mu) / sigma)^2)
}

x <- rnorm(100) ## Generate some data
sumsq(mu = 1, sigma = 1, x) ## This works (returns one value)


## Cannot use vectors for mu or sigma with our current function
sumsq(1:10, 1:10, x) ## This is not what we want


## Create function for vectorized sum of squares
vsumsq <- Vectorize(sumsq, c("mu", "sigma")) ## Vectorize automatically creates vectorized version of your function
vsumsq(1:10, 1:10, x)

## The details are a bit complicated though, as we can see below
vsumsq

## Parallelize your functions using furrr

#################################################
##### CLASS NOTES: LECTURE 18 - Debugging R Code
##################################################

## Finding the root cause is challenging but one of the best packages is colorout -- changes colors of output to make clearer
remotes::install_github("jalvesaq/colorout")

## Make colorout load automatically when you open R
  # Open your .Rprofile file
  # usethis::edit_r_profile()
  # Copy paste the following code taken from
  # https://lcolladotor.github.io/bioc_team_ds/config-files.html#rprofile
  # Change colors
  # Source https://github.com/jalvesaq/colorout
  # if (Sys.getenv("TERM") %in% c("term", "xterm-256color", "cygwin", "screen")) {
  #     if (!requireNamespace("colorout", quietly = TRUE) & .Platform$OS.type != "windows") {
  #         cat('To install colorout use: remotes::install_github("jalvesaq/colorout")\n')
  #     }
  # }

  # Note that you can also change the default colors for colorout

## From colorout's README documentation
require("colorout")
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

## You can also Google erros when they happen

## Or you can use reprex
library("reprex")

(y <- 1:4)
mean(y)
  # And then enter reprex() in the R console and the issue has been copied to your clipboard; can post it on Github

## Produces warning message
log(-1)

## Write function that prints message depending on input
print_message <- function(x) {
  if (x > 0) {
    print("x is greater than zero")
  } else {
    print("x is less than or equal to zero")
  }
  invisible(x)  # Note:invisible return means value does not get auto-printed
}

print_message(1)     ## "x is greater than zero"
print_message(NA)    ## "Error in if (x > 0) { : missing value where TRUE/FALSE needed"


## Edit function to handle NAs
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

print_message2(NA)


## Rework function for cases with longer inputs than expected
x <- log(c(-1, 2))
print_message2(x)

print_message3 <- function(x) {
  if (length(x) > 1L) {
    stop("'x' has length > 1")    # stop creates an error message
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

print_message3(1:2)


## ------------------------------------------------------------------------------------------------------------------------
print_message3_no_call <- function(x) {
  if (length(x) > 1L) {
    stop("'x' has length > 1", call. = FALSE)  # Adding call. = FALSE creates shorter error message which has less information; might be better for googling when using custom functions
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

print_message3_no_call(99:100)
print_message3(99:100)


## Tidyverse has two functions that create error messages
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

print_message3_tidyverse(99:100)  # Abort changes how error displayed (i.e., stop)
print_message3_tidyverse(NA)      # Warning prints in diff color (i.e., warning)
print_message3_tidyverse(1)       # Inform prints message (i.e., message)
print_message3_tidyverse(-1)


## Cli can provide more detailed error messages
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
        "v" = praise::praise(),  #package gives people praise
        "v" = praise_mx          #package gives people praise in spanish
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

set.seed(20230928)
print_message3_cli(-1:1)
purrr::map(-1:1, print_message3_cli)

## Vectorize function
print_message4 <- Vectorize(print_message2)
out <- print_message4(c(-1, 2))


## Understanding where you are getting an error, on the user side rather than function builder
lm(y ~ x)
traceback()            # Shows how many layers deep you were when error occured

rlang::global_entrace()  # Need global entrace first because lm doesnt use rlang internally
lm(y~x)
rlang::last_error()    # Gives additional information about error


#? debug()
#? recover()
#? browser()
#? trace()

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



#################################################
##### CLASS NOTES: LECTURE 19 - Error Handling and Generation
##################################################

## Can't add two strings
"hello" + "world"


## Can't turn "seven" to numeric
as.numeric(c("5", "6", "seven"))


## Messages can be supressed but prints cannot -- can be better to use messages
f <- function() {
  message("This is a message.")
}

f()
suppressMessages(f())


## Error messages with stop
f2 <- function() {
  message(Sys.time(), " This is a message.")
}
f2()

stop("Something erroneous has occurred!")  ## error message


name_of_function <- function() {
  stop("Something bad happened.")
}

name_of_function()


##Stop if not (reverse of stop)
error_if_n_is_greater_than_zero <- function(n) {
  stopifnot(n <= 0)
  n
}

error_if_n_is_greater_than_zero(5)


## ------------------------------------------------------------------------------------------------------------------------
#| error: true
error_if_n_squared <- function(n) {
  ## Ok use
  stopifnot(n <= 0)

  ## Create an internal object
  n_squared <- n^2

  ## Not ok, since we are using the internal object n_squared
  stopifnot(n_squared <= 10)
  n
}

error_if_n_squared(-2)

## This generates a confusing error message to our users
error_if_n_squared(-4)


## rlang arg_match helps with close matches
fn <- function(x = c("foo", "bar")) {
  x <- rlang::arg_match(x)

  ## Known scenario 1
  if (x == "foo") {
    print("I know what to do here with 'x = foo'")
  }

  ## Known scenario 2
  if (x == "bar") {
    print("I know what to do here with 'x = bar'")
  }
}
fn("foo")
fn("zoo")


## Stop vs. Warning
  # Use stop when its really bad and it needs to stop
  # Use warning when its maybe OK to continue

## ------------------------------------------------------------------------------------------------------------------------
#| warning: true
warning("Consider yourself warned!")


## ------------------------------------------------------------------------------------------------------------------------
#| warning: true
make_NA <- function(x) {
  warning("Generating an NA.")
  NA
}

make_NA("Sodium")


## ------------------------------------------------------------------------------------------------------------------------
message("In a bottle.")


## ------------------------------------------------------------------------------------------------------------------------
as.numeric(c("5", "6", "seven"))


## ------------------------------------------------------------------------------------------------------------------------
beera <- function(expr) {
  tryCatch(expr,
           error = function(e) {
             message("An error occurred:\n", e)
           },
           warning = function(w) {
             message("A warning occured:\n", w)
           },
           finally = {
             message("Finally done!")
           }
  )
}

beera({
  2 + 2
})

beera({
  "two" + 2
})

beera({
  as.numeric(c(1, "two", 3))
})


## ----error=TRUE----------------------------------------------------------------------------------------------------------
is_even <- function(n) {
  n %% 2 == 0
}

is_even(768)

is_even("two")


## tryCatch will rerun code a bit later to see if it will work
is_even_error <- function(n) {
  tryCatch(n %% 2 == 0,
           error = function(e) {
             FALSE
           }
  )
}

is_even_error(714)

is_even_error("eight")


## ------------------------------------------------------------------------------------------------------------------------
is_even_check <- function(n) {
  is.numeric(n) && n %% 2 == 0
}

is_even_check(1876)

is_even_check("twelve")


## Microbenchmark package gives how long it takes for each function to tbe applied to the same data
library(microbenchmark)
microbenchmark(sapply(letters, is_even_check))



###############################################################
##### CLASS NOTES: LECTURE 20 - Working with dates and times
###############################################################

# date = <date>
# time = <time>
# date-time = <dttm> OR "POSIXct" OR "POSIXt"

## lubridate package to deal with dates
library("tidyverse")
library("lubridate")
help(package = "lubridate")   ## get list of all functions in lubridate


## lubridate versions
today()    ## today's date in "yyyy-mm-dd" format; class = "date"
now()      ## today's date in "yyyy-mm-dd mm:ss:ss EDT" format; class = "POSIXCT"

## base R versions
base::Sys.Date()   ## equivalent of lubridate's today()
base::Sys.time()   ## equivalent of lubridate's now()

## Store today's date in x and check class (="Date")
x <- today()
class(x)


## Indicate that a field is a date and the order of the mdy
ymd("1970-01-01")
ymd("2017-01-31")
mdy("January 31st, 2017")
dmy("31-Jan-2017")
ymd(20170131)      # still works with numbers

## Base R equivalent to lubridate's mdy/ymd/etc, but quickly becomes complicated
as.Date("1970-01-01")
as.Date("January 31st, 2017", "%B %dst, %Y")
as.Date(gsub("st,", "", "January 31st, 2017"), "%B %d %Y")




## Alternative formulations
ymd("2016-09-13") ## International standard
ymd("2016/09/13") ## Just figure it out
mdy("09-13-2016") ## Mostly U.S.
dmy("13-09-2016") ## Europe

x <- c(
  "2016-04-05",
  "2016/05/06",
  "2016,10,4"
)
ymd(x)


## Look at flights data and create YMD variable with make_date
library("nycflights13")

flights %>%
  select(year, month, day) %>%
  mutate(departure = make_date(year, month, day))


## make date-time column using make_datetime
args(make_datetime)

flights %>%
  select(year, month, day, hour, minute) %>%
  mutate(departure = make_datetime(year, month, day, hour, min = minute))


## Converting date to date-time and vice versa
today()
as_datetime(today())

now()
as_date(now())


## Create date-time objects from strings
ymd_hms("2017-01-31 20:11:59")
mdy_hm("01/31/2017 08:01")

ymd_hms("2016-09-13 14:00:00")
ymd_hms("2016-09-13 14:00:00", tz = "America/New_York") ## can also provide time zone
ymd_hms("2016-09-13 14:00:00", tz = "") ##defaults to your computers time zone


## POSIXct; # seconds since january 1 1970
x <- ymd_hm("1970-01-01 01:00")
class(x)
unclass(x)   # number of seconds
typeof(x)    # a very large number technically
attributes(x)


## POSIXt; similar as POSIXct but more info, # seconds since january 1 1970
y <- as.POSIXlt(x)
y
typeof(y)
attributes(y)  ## POSIXt is rare, POSIXct is more common


## Time Zone suck!
x <- ymd_hm("1970-01-01 01:00", tz = "")
x
attributes(x)   # stores empty time zone


## Specify time zone
attr(x, "tzone") <- "US/Pacific"
x

attr(x, "tzone") <- "US/Eastern"
x


## Artithmetic with date-times
x <- ymd("2012-01-01", tz = "") ## Midnight
y <- dmy_hms("9 Jan 2011 11:34:21", tz = "")
x  #counts as midnight
y
x - y ## this works

x < y ## this works
x > y ## this works
x == y ## this works
x + y ## what??? why does this not work? Cannot add dates

x + 3 * 60 * 60   # add 3 hours (3 hours * 60 minutes * 60 seconds)

## Can also do addition with dates, but unit is days
y <- date(y)
y
y + 1     ##adds a day


## It recognizes leap years
x <- ymd("2012-03-01")
y <- ymd("2012-02-28")
x - y

x <- ymd("2013-03-01")
y <- ymd("2013-02-28")
x - y


## Difference between time zones make it actually the same time
x <- ymd_hms("2012-10-25 01:00:00", tz = "")
y <- ymd_hms("2012-10-25 05:00:00", tz = "GMT")
y - x


## Shows years with leap seconds
.leap.seconds


## can get info about a date, like month, day of week, etc.
x <- ymd_hms(c(
  "2012-10-25 01:13:46",
  "2015-04-23 15:11:23"
), tz = "")
year(x)
month(x)
day(x)
weekdays(x)

x <- ymd_hms(c(
  "2012-10-25 01:13:46",
  "2015-04-23 15:11:23"
), tz = "")
minute(x)
second(x)
hour(x)
week(x)


## Load in storm data
library(here)
library(readr)
library(dplyr)
storm <- read_csv(here("data", "storms_2004.csv.gz"), progress = FALSE)
storm
names(storm)

## Filter out and mutate some things
storm_sub <-
  storm %>%
  select(BEGIN_DATE_TIME, EVENT_TYPE, DEATHS_DIRECT) %>%
  mutate(begin = dmy_hms(BEGIN_DATE_TIME)) %>%
  rename(type = EVENT_TYPE, deaths = DEATHS_DIRECT) %>%
  select(begin, type, deaths)
storm_sub


## Make histogram of when storm start; ggplot already knows what is useful to print
library("ggplot2")
storm_sub %>%
  ggplot(aes(x = begin)) +
  geom_histogram(bins = 20) +
  theme_bw()


## Another histogram
#| fig-width: 12
#| fig-height: 12
library(ggplot2)
storm_sub %>%
  ggplot(aes(x = begin)) +
  facet_wrap(~type) +
  geom_histogram(bins = 20) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(angle = 90))


## Scatter plot
storm_sub %>%
  ggplot(aes(x = begin, y = deaths)) +
  geom_point()


## Another scatter plot
storm_sub %>%
  filter(month(begin) == 6) %>%
  ggplot(aes(begin, deaths)) +
  geom_point()


## Another scatter plot
storm_sub %>%
  filter(month(begin) == 6, day(begin) == 16) %>%
  ggplot(aes(begin, deaths)) +
  geom_point()


## ----eval=FALSE----------------------------------------------------------------------------------------------------------
## ymd(c("2010-10-10", "bananas"))
##
## ## Compare against base R's behavior:
## as.Date(c("2010-10-10", "bananas"))


## ------------------------------------------------------------------------------------------------------------------------
unclass(today())


## ------------------------------------------------------------------------------------------------------------------------
d1 <- "January 1, 2010"
d2 <- "2015-Mar-07"
d3 <- "06-Jun-2017"
d4 <- c("August 19 (2015)", "July 1 (2015)")
d5 <- "12/30/14" # Dec 30, 20




###############################################################
##### CLASS NOTES: LECTURE 21 - Regular Expressions
###############################################################
## Regular expression (i.e., regec or regexp) is a concise language for describing patterns in character strings
  # Can be used for:
    # searching for a pattern or string
    # replacing one part of a string

## Functions
  # grepl(): takes two arguments and returns a logical (true/false)

## String basics
string1 <- "This is a string"
string2 <- 'If I want to include a "quote" inside a string, I use single quotes'
string1
string2

## vector with strings
c("one", "two", "three")


#### GREPL(): takes two arguments and returns a logical
## searching for a regular expression ("a") in a string ("Maryland")
regular_expression <- "a"
string_to_search <- "Maryland"
grepl(pattern = regular_expression, x = string_to_search) # returns "TRUE"

## searching for a regular expression ("u") in a string ("Maryland")
regular_expression <- "u"
string_to_search <- "Maryland"
grepl(pattern = regular_expression, x = string_to_search) # returns "FALSE"

## More grepl searches
grepl("land", "Maryland")
grepl("ryla", "Maryland")
grepl("Marly", "Maryland")
grepl("dany", "Maryland")

# METACHARACTERS

  # . represents any character other than a new line
      grepl(".", "Maryland")            # TRUE
      grepl(".", "*&2[0+,%<@#~|}")      # TRUE
      grepl(".", "")                    # FALSE
      grepl("a.b", c("aaa", "aab", "abb", "acadb"))  # F,T,T,T

  # + indicates one or more of preceding expression is present
      grepl("a+", "Maryland")     # Does "Maryland" contain one or more of "a"
      grepl("x+", "Maryland")     # Does "Maryland" contain one or more of "x" ?

  # * indicates zero or more of preceding expression is present
      grepl("x*", "Maryland")     # Does "Maryland" contain zero or more of "x" ?
      grepl("(xx)*", "Maryland")  # Does "Maryland" contain zero or more of "x" ?

  # ? indicates zero or 1 of preceding expression is NOT present or present AT MOST 1 time

  # {} indicates number of character or metacharacter it contains
      grepl("s{2}", "Mississippi")        # Does "Mississippi" contain exactly 2 adjacent "s" ?
      grepl("ss", "Mississippi")          # Does "Mississippi" contain exactly 2 adjacent "s" ?
      grepl("s{1,3}", "Mississippi")      # Does "Mississippi" contain between 1 and 3 adjacent "s" ?
      grepl("i{2,3}", "Mississippi")      # Does "Mississippi" contain between 2 and 3 adjacent "i" ?

      grepl("(iss){2}", "Mississippi")   # Does "Mississippi" contain between 2 adjacent "iss" ?
      grepl("(ss){2}", "Mississippi")     # Does "Mississippi" contain between 2 adjacent "ss" ?
      grepl("(i.{2}){3}", "Mississippi")  # Does "Mississippi" contain the pattern of an "i" followed by 2 of any character, with that pattern repeated three times adjacently?

  # () creates a capture group
  # Character sets
    # words: "\\w"
      grepl("\\w", "abcdefghijklmnopqrstuvwxyz0123456789")

    # digits: "\\d"
      grepl("\\d", "0123456789")

    # whitespace: "\\s"
    # NOT words: "\\W"
    # NOT digits: "\\D"
    # NOT whitespace: "\\S"
  # "\\" to indicate you want the character, not metacharacter
      grepl("\\+", "tragedy + time = humor")
      grepl("\\.", "https://publichealth.jhu.edu")
      x <- c("\'", "\"", "\\")
      x
      writeLines(x)
  # "\n" = new line

  # "\t" = tab
      x <- c("\\t", "\\n", "\u00b5")
      x
      writeLines(x)


# more grepl examples
      grepl("\\s", "\n\t   ")

      grepl("\\d", "abcdefghijklmnopqrstuvwxyz")

      grepl("\\D", "abcdefghijklmnopqrstuvwxyz")

      grepl("\\w", "\n\t   ")




## Brackets specify specfic character sets
grepl("[aeiou]", "rhythms")

## ^ inside brackets matches all characters EXCEPT the lowercase vowels
grepl("[^aeiou]", "rhythms")


## - means anything between the left and right
grepl("[a-m]", "xyz")
grepl("[a-m]", "ABC")
grepl("[a-mA-M]", "ABC")


## ^ is looking for things at the beginning
grepl("^a", c("bab", "aab"))

## $ is looking for things at the end
grepl("b$", c("bab", "aab"))

## combining some
grepl("^[ab]*$", c("bab", "aab", "abc"))


## | is OR metacharacter -- matches either the regex on left OR right
grepl("a|b", c("abc", "bcd", "cde"))      # looking for a or b
grepl("North|South", c("South Dakota", "North Carolina", "West Virginia"))

## state.name dataset with a vector of strings for each US state
head(state.name)
length(state.name)
grepl("land", state.name)

## state names that begin and end with a vowel
start_end_vowel <- "^[AEIOU]{1}.+[aeiou]{1}$"
vowel_state_lgl <- grepl(start_end_vowel, state.name)
head(vowel_state_lgl)
state.name[vowel_state_lgl]

## Create tibble that shows meaning of each metacharacter
library(knitr)

mc_tibl <- data.frame(
  Metacharacter =
    c(
      ".", "\\\\w", "\\\\W", "\\\\d", "\\\\D",
      "\\\\s", "\\\\S", "[xyz]", "[^xyz]", "[a-z]",
      "^", "$", "\\\\n", "+", "*", "?", "|", "{5}", "{2, 5}",
      "{2, }"
    ),
  Meaning =
    c(
      "Any Character", "A Word", "Not a Word", "A Digit", "Not a Digit",
      "Whitespace", "Not Whitespace", "A Set of Characters",
      "Negation of Set", "A Range of Characters",
      "Beginning of String", "End of String", "Newline",
      "One or More of Previous", "Zero or More of Previous",
      "Zero or One of Previous", "Either the Previous or the Following",
      "Exactly 5 of Previous", "Between 2 and 5 or Previous",
      "More than 2 of Previous"
    ),
  stringsAsFactors = FALSE
)
kable(mc_tibl, align = "c")


## Antoher example
grepl("[Ii]", c("Hawaii", "Illinois", "Kentucky"))


## grep() -- old fashioned version that returns indices of the vector
grep(pattern = "[Ii]", x = c("Hawaii", "Illinois", "Kentucky"))


## sub() replaces the first instance of regex found in each string
sub(pattern = "[Ii]", replacement = "1", x = c("Hawaii", "Illinois", "Kentucky"))


## gsub replaces every instance of the regex that is matched
gsub("[Ii]", "1", c("Hawaii", "Illinois", "Kentucky"))


## split up strings on patterns
two_s <- state.name[grep("ss", state.name)]
two_s
strsplit(x = two_s, split = "ss")


## stringr from tidyverse package
  # all take data and then string as arguments
library(stringr)
state_tbl <- paste(state.name, state.area, state.abb)
head(state_tbl)

## str_extract --> returns substring that matches expression
str_extract(state_tbl, "[0-9]+")

## str_detect --> equivalent to grepl() and returns true/false
str_detect(state_tbl, "[0-9]+")
grepl("[0-9]+", state_tbl)

## str_order --> numeric vector corresponding to alphabetical order of strings in vector
head(state.name)
str_order(state.name)

head(state.abb)
str_order(state.abb)

## str_replace --> equivalent to sub(); replaces regular expressions
str_replace(string = state.name, pattern = "[Aa]", replace = "B")
sub(pattern = "[Aa]", replacement = "B", x = state.name)


## str_pad --> pads strings with other characters
str_pad("Thai", width = 8, side = "left", pad = "-")
str_pad("Thai", width = 8, side = "right", pad = "-")
str_pad("Thai", width = 8, side = "both", pad = "-")


## str_to_title --> puts strings in title case
cases <- c("CAPS", "low", "Title")
str_to_title(cases)


## str_trim --> deletes white space from both sides of a string
to_trim <- c("   space", "the    ", "    final frontier  ")
str_trim(to_trim)


## str_wrap --> inserts new lines in strings to make them only so wide in characters
pasted_states <- paste(state.name[1:20], collapse = " ")
cat(str_wrap(pasted_states, width = 80))
cat(str_wrap(pasted_states, width = 30))


## gets words from a string
a_tale <- "It was the best of times it was the worst of times it was the age of wisdom it was the age of foolishness"
word(a_tale, 2)
word(a_tale, end = 3) # end = last word to extract
word(a_tale, start = 11, end = 15) # start = first word to extract


## ------------------------------------------------------------------------------------------------------------------------
head(stringr::words)
length(stringr::words)


###############################################################
##### CLASS NOTES: LECTURE 22 - Working with factors
###############################################################

## Factors are usefule when:
  # you want to include categorical variables in regression models
  # you want to plot categorical data
  # you want to display character vectors in a non-alphabetical order

## Issues with vectors

# Sorting alphabetically might not be useful
x <- c("Dec", "Apr", "Jan", "Mar")
sort(x)

# no protection from typos
x_typo <- c("Dec", "Apr", "Jam", "Mar")

## Create vector with months and create a factor with it
month_levels <- c(
  "Jan", "Feb", "Mar", "Apr", "May", "Jun",
  "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
)

y <- factor(x, levels = month_levels)
y
# Check attributes and levels of y
attributes(y)
levels(y)
class(y)

# Sort now sorts based on the order of the level
sort(y)

## Any values not in the level will be silently converted to NA
y_typo <- factor(x_typo, levels = month_levels)
y_typo


## What if new categories are added, categories change, or mess up relation between numeric and label?
## Get numeric value of factor
library(tidyverse)

x1_original <- c(10, 10, 10, 50, 60, 20, 20, 40)
x1_factor <- factor(x1_original)
attributes(x1_factor)

tibble(x1_original, x1_factor) %>%
  mutate(x1_numeric = as.numeric(x1_factor))


## as.numeric with vectors
as.numeric(c("hello"))
as.numeric(factor(c("hello")))
as.numeric(factor(c("hello", "goodbye")))  # since we didn't specify levels, R stores it alphabetically as factor

## factor silently makes a missing value if the values in the data and levels do not match
factor("a", levels = "c")

## BECAUSE OF DISPUTES OVER HOW FACTORS WORK, DEFAULT FOR DATA IMPORTS IS TO
## READ STRINGS AS STRINGS NOT FACTORS

## Make income level factor and random Y and perform linear regression
  # factors must be used if you want to be able to pick reference category
income_level <- c(
  rep("low", 10),
  rep("medium", 10),
  rep("high", 10)
)
income_level

x <- factor(income_level)
x

y <- rnorm(30) # generate some random obs from a normal dist
lm(y ~ x)      # high as reference category

## factors used to be more efficient in storage but this is less of an issue now
income_level <- c(
  rep("low", 10000),
  rep("medium", 10000),
  rep("high", 10000)
)

format(object.size(income_level), units = "Kb") # size of the character string
format(object.size(factor(income_level)), units = "Kb") # size of the factor

#### FORCATS PACKAGE
## load forcats package for working with categorical variables (part of tidvyerse)
library("forcats")
gss_cat
gss_cat %>%
  count(race)


## A bar chart showing the distribution of race.
  # There are ~2000 records with race "Other", 3000 with race "Black" and other,
  # 15,000 with race "White".
gss_cat %>%
  ggplot(aes(x = race)) +
  geom_bar()


## count number in each religion factor
gss_cat %>%
  count(relig)

attributes(gss_cat$relig)


## Scatterplot with poorly ordered attributes
#| fig-alt: >
#|   A scatterplot of with tvhours on the x-axis and religion on the y-axis.
#|   The y-axis is ordered seemingly aribtrarily making it hard to get
#|   any sense of overall pattern.
relig_summary <- gss_cat %>%
  group_by(relig) %>%
  summarise(
    tvhours = mean(tvhours, na.rm = TRUE),
    n = n()
  )

relig_summary %>%
  ggplot(aes(x = tvhours, y = relig)) +
  geom_point()


## Scatterplot ordered by increasing order ot tvhrs
#| fig-alt: >
#|   The same scatterplot as above, but now the religion is displayed in
#|   increasing order of tvhours. "Other eastern" has the fewest tvhours
#|   under 2, and "Don't know" has the highest (over 5).
relig_summary %>%
  ggplot(aes(
    x = tvhours,
    y = fct_reorder(.f = relig, .x = tvhours)
  )) +
  geom_point()


## ------------------------------------------------------------------------------------------------------------------------
relig_summary %>%
  mutate(relig = fct_reorder(relig, tvhours)) %>%
  ggplot(aes(x = tvhours, y = relig)) +
  geom_point()


## Scatterplot w/ and w/o ordering by age
#| fig-alt: >
#|   A scatterplot with age on the x-axis and income on the y-axis. Income
#|   has been reordered in order of average age which doesn't make much
#|   sense. One section of the y-axis goes from $6000-6999, then <$1000,
#|   then $8000-9999.
rincome_summary <-
  gss_cat %>%
  group_by(rincome) %>%
  summarise(
    age = mean(age, na.rm = TRUE),
    n = n()
  )

## Original rincome order
rincome_summary %>%
  ggplot(aes(x = age, y = rincome)) +
  geom_point()

## rincome re-ordered by age's values
rincome_summary %>%
  ggplot(aes(x = age, y = fct_reorder(.f = rincome, .x = age))) +
  geom_point()


## load penguins data
library(palmerpenguins)
penguins

## Reorder scatterplot again
#| fig-alt: >
#|   The same scatterplot but now "Not Applicable" is displayed at the
#|   bottom of the y-axis. Generally there is a positive association
#|   between income and age, and the income band with the highest average
#|   age is "Not applicable".
rincome_summary %>%
  ggplot(aes(age, fct_relevel(rincome, "Not applicable"))) +
  geom_point()


## line plot with rearranging of legend order via forcat
  #fct_reorder2(f,x,y) reorders factor f by y values associated with largest x values
#| layout-ncol: 2
#| fig-width: 4
#| fig-height: 2
#| fig-alt: >
#|   - A line plot with age on the x-axis and proportion on the y-axis.
#|     There is one line for each category of marital status: no answer,
#|     never married, separated, divorced, widowed, and married. It is
#|     a little hard to read the plot because the order of the legend is
#|     unrelated to the lines on the plot.
#|   - Rearranging the legend makes the plot easier to read because the
#|     legend colours now match the order of the lines on the far right
#|     of the plot. You can see some unsuprising patterns: the proportion
#|     never marred decreases with age, married forms an upside down U
#|     shape, and widowed starts off low but increases steeply after age
#|     60.
by_age <-
  gss_cat %>%
  filter(!is.na(age)) %>%
  count(age, marital) %>%
  group_by(age) %>%
  mutate(prop = n / sum(n))

by_age %>%
  ggplot(aes(age, prop, colour = marital)) +
  geom_line(na.rm = TRUE)

by_age %>%
  ggplot(aes(age, prop, colour = fct_reorder2(marital, age, prop))) +
  geom_line() +
  labs(colour = "marital")


## Bar chart ordered by frequency of marital status
#| fig-alt: >
#|   A bar char of marital status ordered in from least to most common:
#|   no answer (~0), separated (~1,000), widowed (~2,000), divorced
#|   (~3,000), never married (~5,000), married (~10,000).
gss_cat %>%
  mutate(marital = marital %>% fct_infreq() %>% fct_rev()) %>%
  ggplot(aes(marital)) +
  geom_bar()


## fct_recode for partyid factor to make clearer labels
gss_cat %>%
  count(partyid)
gss_cat %>%
  mutate(partyid = fct_recode(partyid,
                              "Republican, strong"    = "Strong republican",
                              "Republican, weak"      = "Not str republican",
                              "Independent, near rep" = "Ind,near rep",
                              "Independent, near dem" = "Ind,near dem",
                              "Democrat, weak"        = "Not str democrat",
                              "Democrat, strong"      = "Strong democrat"
  )) %>%
  count(partyid)

## combin groups using fct_recode
gss_cat %>%
  mutate(partyid = fct_recode(partyid,
                              "Republican, strong"    = "Strong republican",
                              "Republican, weak"      = "Not str republican",
                              "Independent, near rep" = "Ind,near rep",
                              "Independent, near dem" = "Ind,near dem",
                              "Democrat, weak"        = "Not str democrat",
                              "Democrat, strong"      = "Strong democrat",
                              "Other"                 = "No answer",
                              "Other"                 = "Don't know",
                              "Other"                 = "Other party"
  )) %>%
  count(partyid)


## fct_collapse to combine factors
gss_cat %>%
  mutate(partyid = fct_collapse(partyid,
                                "other" = c("No answer", "Don't know", "Other party"),
                                "rep" = c("Strong republican", "Not str republican"),
                                "ind" = c("Ind,near rep", "Independent", "Ind,near dem"),
                                "dem" = c("Not str democrat", "Strong democrat")
  )) %>%
  count(partyid)


## fct_lump_lowfreq to lump together smaller groups
gss_cat %>%
  mutate(relig = fct_lump_lowfreq(relig)) %>%
  count(relig)   ## grouped too much though


## fct_lump_n allows you to specify number of groups
gss_cat %>%
  mutate(relig = fct_lump_n(relig, n = 10)) %>%
  count(relig, sort = TRUE) %>%
  print(n = Inf)


## also fct_lump_min()
## also fct_lump_prop()


## Ordered factors imply a strict ordering and equal distance between levels
ordered(c("a", "b", "c"))


###############################################################
##### CLASS NOTES: LECTURE 23 - Tidytext and sentiment analysis
###############################################################

# tidy text format -- one-token-per-row where a token is a meaninful unit of text

# text can be stored as:
  # string
  # corpus -- contain raw strings but have some metadata
  # document-term matrix (one row for each document, on column for each term)


## ----out.width = "95%", echo = FALSE-------------------------------------------------------------------------------------
knitr::include_graphics("http://r4ds.had.co.nz/images/tidy-1.png")


## ----echo=FALSE, out.width = '90%', fig.cap="A flowchart of a typical text analysis using tidy data principles."---------
knitr::include_graphics("https://www.tidytextmining.com/images/tmwr_0101.png")


## ----eval=FALSE----------------------------------------------------------------------------------------------------------
## ?unnest_tokens

## load packages for working with tidy text data
library(tidyverse)
library(stringr)
library(tidytext) ## needs to be installed
library(janeaustenr) ## needs to be installed

## ------------------------------------------------------------------------------------------------------------------------
#| label: pengpreface
#| echo: false
#| fig-cap: 'Preface from R Programming for Data Science'
#| out-width: '90%'
knitr::include_graphics("../../images/peng_preface.png")


## each sentence in one string of vector
peng_preface <-
  c(
    "I started using R in 1998 when I was a college undergraduate working on my senior thesis.",
    "The version was 0.63.",
    "I was an applied mathematics major with a statistics concentration and I was working with Dr. Nicolas Hengartner on an analysis of word frequencies in classic texts (Shakespeare, Milton, etc.).",
    "The idea was to see if we could identify the authorship of each of the texts based on how frequently they used certain words.",
    "We downloaded the data from Project Gutenberg and used some basic linear discriminant analysis for the modeling.",
    "The work was eventually published and was my first ever peer-reviewed publication.",
    "I guess you could argue it was my first real 'data science' experience."
  )

peng_preface


## take character vector and turn it into a tibble
peng_preface_df <- tibble(
  line = 1:7,
  text = peng_preface
)
peng_preface_df


## unnest_tokens --> create tibble with a row for each word and columns for the word and the line it was from
peng_token <-
  peng_preface_df %>%
  unnest_tokens(
    output = word,
    input = text,
    token = "words"
  )

peng_token %>%
  head()

peng_token %>%
  tail()


## unnest_tokens --> create tibble with a row for each character and columns for the character and the line it was from
peng_preface_df %>%
  unnest_tokens(word,
                text,
                token = "characters"
  ) %>%
  head()


## unnest_tokens --> create tibble with a row for every 3 words and columns for the 3 words and the line it was from
peng_preface_df %>%
  unnest_tokens(word,
                text,
                token = "ngrams",
                n = 3
  ) %>%
  head()


## unnest_tokens --> create tibble with a row for every 4 characters and columns for the 4 characters and the line it was from
peng_preface_df %>%
  unnest_tokens(word,
                text,
                token = "character_shingles",
                n = 4
  ) %>%
  head()


## unnest_tokens --> create tibble with a row for every set between spaces
peng_preface_df %>%
  unnest_tokens(word,
                text,
                token = stringr::str_split,
                pattern = " "
  ) %>%
  head()


## New example text to tibble
gorman_hill_we_climb <-
  c(
    "When day comes we ask ourselves, where can we find light in this neverending shade?",
    "The loss we carry, a sea we must wade.",
    "We’ve braved the belly of the beast, we’ve learned that quiet isn’t always peace and the norms and notions of what just is, isn’t always justice.",
    "And yet the dawn is ours before we knew it, somehow we do it, somehow we’ve weathered and witnessed a nation that isn’t broken but simply unfinished."
  )

hill_df <- tibble(
  line = seq_along(gorman_hill_we_climb),
  text = gorman_hill_we_climb
)
hill_df

## unnest_tokens: token = words
hill_df %>%
  unnest_tokens(
    output = wordsforfun,   ## column name in new data where words will go
    input = text,           ## column name in original data
    token = "words"         ## desired token
  )

## unnest_tokens: token = characters
hill_df %>%
  unnest_tokens(
    output = wordsforfun,   ## column name in new data where words will go
    input = text,           ## column name in original data
    token = "characters",         ## desired token
  )

## unnest_tokens: token = sentences
hill_df %>%
  unnest_tokens(
    output = wordsforfun,   ## column name in new data where words will go
    input = text,           ## column name in original data
    token = "sentences",         ## desired token
  )

## unnest_tokens: token = ngrams
hill_df %>%
  unnest_tokens(
    output = wordsforfun,   ## column name in new data where words will go
    input = text,           ## column name in original data
    token = "ngrams",         ## desired token
    n = 3
  )

## unnest_tokens: token = character_shingles
hill_df %>%
  unnest_tokens(
    output = wordsforfun,            ## column name in new data where words will go
    input = text,                   ## column name in original data
    token = "character_shingles",   ## desired token
    n = 4
  )

## load jane austen package and put in tibble
library(janeaustenr)
head(prideprejudice, 20)

pp_book_df <- tibble(text = prideprejudice)

## unnest_tokens, token = words
pp_book_df %>%
  unnest_tokens(
    output = word,
    input = text,
    token = "words"
  )


## unnest_tokens, token = paragraphs
tmp <- pp_book_df %>%
  unnest_tokens(
    output = paragraph,
    input = text,
    token = "paragraphs"
  )
tmp
tmp[3, 1]


## unnest_tokens, token = sentences; gets tricked by Mr. and Mrs.
pp_book_df %>%
  unnest_tokens(
    output = sentence,
    input = text,
    token = "sentences"
  )


## unnest_tokens, token = paragraphs and then token = words
paragraphs <-
  pp_book_df %>%
  unnest_tokens(
    output = paragraph,
    input = text,
    token = "paragraphs"
  ) %>%
  mutate(paragraph_number = row_number())

paragraphs

paragraphs %>%
  unnest_tokens(
    output = word,
    input = paragraph
  )


## ANALYZING TEXT ONCE EXTRACTED

## stop words (i.e., common words like the, of, to, etc.)
data(stop_words)     # set of stop words
table(stop_words$lexicon)     # lexicons of stop words
stop_words %>%
  head(n = 10)

## anti_join to remove stop words
words_by_paragraph <-
  paragraphs %>%
  unnest_tokens(
    output = word,
    input = paragraph
  ) %>%
  anti_join(stop_words)   # return all rows of x without a match in y

words_by_paragraph


## Look at top six words in the book
words_by_paragraph %>%
  count(word, sort = TRUE) %>%
  head()


## make plot of word count by word for words with at least 150
words_by_paragraph %>%
  count(word, sort = TRUE) %>%
  filter(n > 150) %>%
  mutate(word = fct_reorder(word, n)) %>%
  ggplot(aes(word, n)) +
  geom_col() +
  xlab(NULL) +
  coord_flip()


## get line number and chapter to find where all the chapters are
austen_books() %>%
  head()

original_books <-
  austen_books() %>%
  group_by(book) %>%
  mutate(
    linenumber = row_number(),
    chapter = cumsum(                                  ## cumulative sum of logical vectors (where true = 1, false = 0)
      str_detect(text,
                 pattern = regex(
                   pattern = "^chapter [\\divxlc]",     ## find anything that says "Chapter" and hten a number in digits or roman numerals
                   ignore_case = TRUE
                 )
      )
    )
  ) %>%
  ungroup()

original_books


## unnest_tokens = words and remove stop words
tidy_books <- original_books %>%
  unnest_tokens(word, text) %>%
  anti_join(stop_words)

tidy_books


## most common words across all books by jane austen
tidy_books %>%
  count(word, sort = TRUE) %>%
  filter(n > 600) %>%
  mutate(word = fct_reorder(word, n)) %>%
  ggplot(aes(word, n)) +
  geom_col() +
  xlab(NULL) +
  coord_flip()


## ----echo=FALSE, out.width = '90%', fig.cap="A flowchart of a typical text analysis that uses tidytext for sentiment analysis."----
knitr::include_graphics("https://www.tidytextmining.com/images/tmwr_0201.png")


## SENTIMENT ANALYSES
  # NRC lexicon -- positive, negative, anger, anticipation, discust, fear, joy, sadness, surprise, trust
    get_sentiments("nrc")


    # bing lexicon -- positive, negative
    get_sentiments("bing")


    # afinn lexicon -- -5 to 5 based on negativity vs positivity
    get_sentiments("afinn")

## Get joy sentiments and count sentiment frequency
nrc_joy <- get_sentiments("nrc") %>%
  filter(sentiment == "joy")

tidy_books %>%
  filter(book == "Emma") %>%
  inner_join(nrc_joy) %>%
  count(word, sort = TRUE)


## get sentiment score for each word
tidy_books %>%
  inner_join(get_sentiments("bing"))


## count positive and negative words in each section of the book
tidy_books %>%
  inner_join(get_sentiments("bing")) %>%
  count(book,
        index = linenumber %/% 80,
        sentiment
  )


## pivot_wider because we multiple rows for same sentiment
jane_austen_sentiment <-
  tidy_books %>%
  inner_join(get_sentiments("bing")) %>%
  count(book,
        index = linenumber %/% 80,
        sentiment
  ) %>%
  pivot_wider(
    names_from = sentiment,
    values_from = n,
    values_fill = 0
  ) %>%
  mutate(sentiment = positive - negative)

jane_austen_sentiment


## plot sentiment across all the books
jane_austen_sentiment %>%
  ggplot(aes(x = index, y = sentiment, fill = book)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(. ~ book, ncol = 2, scales = "free_x")


## create word cloud
library(wordcloud)

tidy_books %>%
  anti_join(stop_words) %>%
  count(word) %>%
  with(wordcloud(word, n, max.words = 100))


## ----echo=FALSE, out.width = '90%', fig.cap=" A flowchart of a typical text analysis that combines tidytext with other tools and data formats, particularly the `tm` or `quanteda` packages. Here, we show how to convert back and forth between document-term matrices and tidy data frames, as well as converting from a Corpus object to a text data frame."----
knitr::include_graphics("https://www.tidytextmining.com/images/tmwr_0501.png")


## ------------------------------------------------------------------------------------------------------------------------
tidy_austen <-
  austen_books() %>%
  mutate(line = row_number()) %>%
  unnest_tokens(word, text) %>%
  anti_join(stop_words)

tidy_austen


## ------------------------------------------------------------------------------------------------------------------------
austen_sparse <- tidy_austen %>%
  count(line, word) %>%
  cast_sparse(row = line, column = word, value = n)

austen_sparse[1:10, 1:10]


## ------------------------------------------------------------------------------------------------------------------------
austen_dtm <- tidy_austen %>%
  count(line, word) %>%
  cast_dtm(document = line, term = word, value = n)

austen_dtm


## ------------------------------------------------------------------------------------------------------------------------
class(austen_dtm)
dim(austen_dtm)
as.matrix(austen_dtm[1:20, 1:10])


## ------------------------------------------------------------------------------------------------------------------------
options(width = 120)
sessioninfo::session_info()


###############################################################
##### CLASS NOTES: LECTURE 24 - Best practices for data analyses
###############################################################
## ETHICS
  # Applied ethics
  # Ethical theory
  # Metaethics

## Many data ethics case studies are around predictive algorithms
  # example: black-box predictive algorithms and facial-recognition in policing
  # How was the algorithm developed?
    # What data used to train?
    # What criteria used to tune?
    # Accuracy rates by demographic?
    # Who has access to the algorithm and data?
    # Should data be accessible to researchers? Bad actors?
  # How is the algorithm used?
    # Misuse -- inaccurately tying people to crimes
    # do people have a right to know they are in databases?
    # companies profit w/o financial responsibility for misuse?

## Best practices for sharing data
  # FAIR Principles
    # Findable -- good metadata and data so easy to find
    # Accessible -- once found, clear how to access (authentication, authorization)
    # Interoperable -- integrated with other data
    # Reusable -- optimize reuse of data

## Addressing concerns re: sharing data
  # decrease impact of novel work -- share only after publication
  # time spent on sharing data -- don't have to store internally tho
  # human subjects data -- challenging but can be done (controlled-access repositories)

## What data to share?
  # data itself
  # metadata
  # data dictionary
  # source code
  # licensing -- choosealicense.com

## When adapting code
  # Copy code
  # Version control before making any changes and
      # put permalink to original source on commit message
      # Use github coauthored message
  # Make edits
  # auto-style code using styler

## Best practices for data visualization
  # Questions to ask yourself
    # What is the question?
    # Why are we building the visualization?
    # For whom are we producing the visualization?
  # A good visualization tells a complete story in a single frame
  # What type of plot to make? How to optimize effectiveness?

  # Developing plots
    # Show comparisons
    # Show causality, mechanism, explanation
    # Show multivariate data
    # Integrate multiple modes of evidence
    # Describe and document the evidence
    # Content is king - good plots start with good questions
  # Optimizing plots
    # Maximize the data/ink ratio – if “ink” can be removed without reducing the information being communicated, then it should be removed.
    # Maximize the range of perceptual conditions – your audience’s perceptual abilities may not be fully known, so it’s best to allow for a wide range, to the extent possible (or knowable).
    # Show variation in the data, not variation in the design.
  # Bad Plots
    # Display as little information as possible.
    # Obscure what you do show (with chart junk).
    # Use pseudo-3D and color gratuitously.
    # Make a pie chart (preferably in color and 3D).
    # Use a poorly chosen scale.
    # Ignore significant figures.

  # Principles
    # Create expository graphs to tell a story (figure and caption should be self-sufficient; it’s the first thing people look at)
    # Be accurate and clear
    # Let the data speak
    # Make axes, labels and titles big
    # Make labels full names (ideally with units when appropriate)
    # Add informative legends; use space effectively
    # Show as much information as possible, taking care not to obscure the message
    # Science not sales: avoid unnecessary frills (especially gratuitous 3D)
    # In tables, every digit should be meaningful

  # Pie charts and donut plots suck! use bar plots instead
  # Avoid 3D barplots
  # Instead of bar plots, can also use box plots
  # Instead of paired bar plots, try scatter plots or line plots




library(tidyverse)


## ------------------------------------------------------------------------------------------------------------------------
#| fig-align: center
#| echo: false
#| fig-cap-location: "top"
#| fig-width: 4
knitr::include_graphics("http://upload.wikimedia.org/wikipedia/en/e/e9/John_Tukey.jpg")


## ------------------------------------------------------------------------------------------------------------------------
#| warning: false
d <- airquality %>%
  mutate(Summer = ifelse(Month %in% c(7, 8, 9), 2, 3))
with(d, {
  plot(Temp, Ozone, col = unclass(Summer), pch = 19, frame.plot = FALSE)
  legend("topleft",
         col = 2:3, pch = 19, bty = "n",
         legend = c("Summer", "Non-Summer")
  )
})


## ------------------------------------------------------------------------------------------------------------------------
#| warning: false
airquality %>%
  mutate(Summer = ifelse(Month %in% c(7, 8, 9),
                         "Summer", "Non-Summer"
  )) %>%
  ggplot(aes(Temp, Ozone)) +
  geom_point(aes(color = Summer), size = 2) +
  theme_minimal()


## ------------------------------------------------------------------------------------------------------------------------
browsers <- c(
  Chrome = 60, Safari = 14, UCBrowser = 7,
  Firefox = 5, Opera = 3, IE = 3, Noinfo = 8
)
browsers.df <- gather(
  data.frame(t(browsers)),
  "browser", "proportion"
)


## ------------------------------------------------------------------------------------------------------------------------
pie(browsers, main = "Browser Usage (July 2022)")


## ------------------------------------------------------------------------------------------------------------------------
#| eval: false
## ?pie


## ------------------------------------------------------------------------------------------------------------------------
p <- browsers.df %>%
  ggplot(aes(
    x = reorder(browser, -proportion),
    y = proportion
  )) +
  geom_bar(stat = "identity")
p


## ------------------------------------------------------------------------------------------------------------------------
#| eval: false
## ?ggplot2::theme


## ------------------------------------------------------------------------------------------------------------------------
p <- p + xlab("Browser") +
  ylab("Proportion of Users")
p


## ------------------------------------------------------------------------------------------------------------------------
p + ggtitle("Browser Usage (July 2022)")


## ------------------------------------------------------------------------------------------------------------------------
p + ggtitle("Browser Usage (July 2022)") +
  theme(plot.title = element_text(hjust = 0.5))


## ------------------------------------------------------------------------------------------------------------------------
p <- p + ggtitle("Browser Usage (July 2022)") +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 15)
  )
p


## ------------------------------------------------------------------------------------------------------------------------
p + theme_bw()


## ------------------------------------------------------------------------------------------------------------------------
p + theme_dark()


## ------------------------------------------------------------------------------------------------------------------------
p + theme_classic() # axis lines!


## ------------------------------------------------------------------------------------------------------------------------
p + ggthemes::theme_base()


## ------------------------------------------------------------------------------------------------------------------------
set.seed(1000)
dat <- data.frame(
  "Treatment" = rnorm(10, 30, sd = 4),
  "Control" = rnorm(10, 36, sd = 4)
)
gather(dat, "type", "response") %>%
  ggplot(aes(type, response)) +
  geom_boxplot() +
  geom_point(position = "jitter") +
  ggtitle("Response to drug treatment")


## ------------------------------------------------------------------------------------------------------------------------
set.seed(1000)
dat <- data.frame(
  "Treatment" = rgamma(10, 10, 1),
  "Control" = rgamma(10, 1, .01)
)
gather(dat, "type", "response") %>%
  ggplot(aes(type, response)) +
  geom_boxplot() +
  geom_point(position = "jitter")


## ------------------------------------------------------------------------------------------------------------------------
gather(dat, "type", "response") %>%
  ggplot(aes(type, response)) +
  geom_boxplot() +
  geom_point(position = "jitter") +
  scale_y_log10()


## ------------------------------------------------------------------------------------------------------------------------
set.seed(1000)
before <- runif(6, 5, 8)
after <- rnorm(6, before * 1.15, 2)
li <- range(c(before, after))
ymx <- max(abs(after - before))

par(mfrow = c(1, 2))
plot(before, after,
     xlab = "Before", ylab = "After",
     ylim = li, xlim = li
)
abline(0, 1, lty = 2, col = 1)

plot(before, after - before,
     xlab = "Before", ylim = c(-ymx, ymx),
     ylab = "Change (After - Before)", lwd = 2
)
abline(h = 0, lty = 2, col = 1)


## ------------------------------------------------------------------------------------------------------------------------
z <- rep(c(0, 1), rep(6, 2))
par(mfrow = c(1, 2))
plot(z, c(before, after),
     xaxt = "n", ylab = "Response",
     xlab = "", xlim = c(-0.5, 1.5)
)
axis(side = 1, at = c(0, 1), c("Before", "After"))
segments(rep(0, 6), before, rep(1, 6), after, col = 1)

boxplot(before, after,
        names = c("Before", "After"),
        ylab = "Response"
)


## ----message=FALSE-------------------------------------------------------------------------------------------------------
x <- read_csv("https://github.com/kbroman/Talk_Graphs/raw/master/R/fig8dat.csv") %>%
  as_tibble(.name_repair = make.names)

p <- x %>%
  gather("drug", "proportion", -log.dose) %>%
  ggplot(aes(
    x = log.dose, y = proportion,
    color = drug
  )) +
  geom_line()
p


## ------------------------------------------------------------------------------------------------------------------------
p + ggtitle("Survival proportion") +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 15)
  )


## ------------------------------------------------------------------------------------------------------------------------
p + ggtitle("Survival proportion") +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 15),
    legend.position = c(0.2, 0.3)
  )


## ------------------------------------------------------------------------------------------------------------------------
transparent_legend <- theme(
  legend.background = element_rect(fill = "transparent"),
  legend.key = element_rect(
    fill = "transparent",
    color = "transparent"
  )
)

p + ggtitle("Survival proportion") +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 15),
    legend.position = c(0.2, 0.3)
  ) +
  transparent_legend


## ------------------------------------------------------------------------------------------------------------------------
heights <- cbind(
  rnorm(8, 73, 3), rnorm(8, 73, 3), rnorm(8, 80, 3),
  rnorm(8, 78, 3), rnorm(8, 78, 3)
)
colnames(heights) <- c("SG", "PG", "C", "PF", "SF")
rownames(heights) <- paste("team", 1:8)
heights


## ------------------------------------------------------------------------------------------------------------------------
round(heights, 1)


## ----fig.cap=""----------------------------------------------------------------------------------------------------------
transparent_legend <- theme(
  legend.background = element_rect(fill = "transparent"),
  legend.key = element_rect(
    fill = "transparent",
    color = "transparent"
  )
)

p + ggtitle("Survival proportion") +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 15),
    legend.position = c(0.2, 0.3)
  ) +
  xlab("dose (mg)") +
  transparent_legend


## ------------------------------------------------------------------------------------------------------------------------
options(width = 120)
sessioninfo::session_info()

###############################################################
##### CLASS NOTES: LECTURE 25 - Python for R Users
###############################################################

## SEE PYTHON EXAMPLE RMARKDOWN

## ----repl-python, echo=FALSE, fig.cap='Using the repl_python() function', fig.align='center'-----------------------------
knitr::include_graphics("https://rstudio.github.io/reticulate/images/python_repl.png")



## ------------------------------------------------------------------------------------------------------------------------
options(width = 120)
sessioninfo::session_info()
