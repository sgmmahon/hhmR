#' example_time_series
#'
#' Fake migration dataset used to demonstrate the functionality of the tshm function in the hhmR package.
#' It contains the information on the number of people who have immigrated a series of fictional
#' geographies over the years 2011 to 2015. The geographies themselves have a hierarchical structure,
#' with each county existing within a smaller subset of regions.
#'
#' @docType data
#'
#' @usage data(example_time_series)
#'
#' @format A data frame with 90 rows and 4 variables.
#' \describe{
#'   \item{County     }{The county (lower-level geography) that immigrants move to.    }
#'   \item{Region     }{The region (higher-level geography) that immigrants move to.   }
#'   \item{Year       }{The year during which each wave of immigration occured.        }
#'   \item{Immigration}{The number of immigrants that moved each county in each year.  }
#' }
#'
#' @keywords datasets
#'
#' @examples
#' library(dplyr)
#' library(tidyr)
#'
#' # Define names of fake counties
#' fake_counties = c("Greenridge","Windermoor","Bramblewood","Silverlake",
#'                   "Thornbury","Maplewood","Hawthorne","Pinehurst",
#'                   "Riverton","Meadowbrook","Fairhaven","Oakdale","Stonebridge",
#'                   "Brookfield","Ashford","Glenville","Sunnyvale","Westfield")
#'
#' # Create dataframe of fake migration data
#' set.seed(1234)
#' example_time_series = data.frame(region = c(rep("North",3),rep("Midlands",5),
#'                                             rep("South West",4),rep("South East",6)),
#'                                  county = fake_counties,
#'                                  year_2011 = sample(1:10000,length(fake_counties)),
#'                                  year_2012 = sample(1:10000,length(fake_counties)),
#'                                  year_2013 = sample(1:10000,length(fake_counties)),
#'                                  year_2014 = sample(1:10000,length(fake_counties)),
#'                                  year_2015 = sample(1:10000,length(fake_counties))) %>%
#'   setNames(c("Region","County",2011:2015)) %>%
#'   pivot_longer(cols = `2011`:`2015`,
#'                       names_to = "Year",
#'                       values_to = "Immigration") %>%
#'   mutate(Year = as.numeric(Year))
#' example_time_series[sample(1:(length(fake_counties)*5),5),"Immigration"] = NA
#' @importFrom tidyr pivot_longer
"example_time_series"
