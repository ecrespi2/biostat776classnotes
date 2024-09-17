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
##### CLASS NOTES: LECTURE 11 - ggplot
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
##### CLASS NOTES: LECTURE 12 - R nuts and bolts
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


## ------------------------------------------------------------------------------------------------------------------------
x <- data.frame(foo = 1:4, bar = c(T, T, F, F))
x
nrow(x)
ncol(x)
attributes(x)


## ------------------------------------------------------------------------------------------------------------------------
data.matrix(x)
attributes(data.matrix(x))


## ------------------------------------------------------------------------------------------------------------------------
#| message: false
# try it yourself

library(tidyverse)
library(palmerpenguins)
penguins



## ------------------------------------------------------------------------------------------------------------------------
x <- 1:3
names(x)
names(x) <- c("New York", "Seattle", "Los Angeles")
x
names(x)
attributes(x)


## ------------------------------------------------------------------------------------------------------------------------
x <- list("Los Angeles" = 1, Boston = 2, London = 3)
x
names(x)


## ------------------------------------------------------------------------------------------------------------------------
m <- matrix(1:4, nrow = 2, ncol = 2)
dimnames(m) <- list(c("a", "b"), c("c", "d"))
m


## ------------------------------------------------------------------------------------------------------------------------
colnames(m) <- c("h", "f")
rownames(m) <- c("x", "z")
m


## ------------------------------------------------------------------------------------------------------------------------
options(width = 120)
sessioninfo::session_info()

