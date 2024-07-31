#' example_migration
#'
#' Fake migration dataset used to demonstrate the functionality of the hhm function in the hhmR package.
#' It contains the information on the number of people who have moved between a series of fictional
#' geographies. The geographies themselves have a hierarchical structure, with each county existing within
#' a smaller subset of regions.
#'
#' @docType data
#'
#' @usage data(example_migration)
#'
#' @format A data frame with 324 rows and 5 variables.
#' \describe{
#'   \item{Origin County     }{The county (lower-level geography) that each migrant began in.    }
#'   \item{Destination County}{The county (lower-level geography) that each migrant ended up in. }
#'   \item{Origin Region     }{The region (higher-level geography) that each migrant began in.   }
#'   \item{Destination Region}{The region (higher-level geography) that each migrant ended up in.}
#'   \item{Migration         }{The number of migrants that moved from each origin county to each destination county.}
#' }
#'
#' @keywords datasets
#'
#' @examples
#' library(dplyr)
#'
#' # Code to create dataset
#'
#' # Define names of fake counties
#' fake_counties = c("Greenridge","Windermoor","Bramblewood","Silverlake",
#'                   "Thornbury","Maplewood","Hawthorne","Pinehurst",
#'                   "Riverton","Meadowbrook","Fairhaven","Oakdale","Stonebridge",
#'                   "Brookfield","Ashford","Glenville","Sunnyvale","Westfield")
#'
#' # Create region county lookup tables
#' rc_lkp = data.frame(region = c(rep("North",3),rep("Midlands",5),
#'                                rep("South West",4),rep("South East",6)),
#'                     county = fake_counties)
#' og_lkp = rc_lkp %>% setNames(c("Origin Region"     ,"Origin County"     ))
#' dn_lkp = rc_lkp %>% setNames(c("Destination Region","Destination County"))
#'
#' # Create dataframe of fake migration data
#' set.seed(1234)
#' example_migration = expand.grid(fake_counties,fake_counties) %>%
#'                     setNames(paste(c("Origin","Destination"),"County",sep=" ")) %>%
#'                     full_join(og_lkp) %>% full_join(dn_lkp) %>%
#'                     mutate(Migration = (1/rgamma(18*18, shape = 17, rate = 0.5)) %>%
#'                                        {. * 1000} %>% round())
#' example_migration[example_migration$`Origin County` ==
#'                   example_migration$`Destination County`,"Migration"] =
#'  example_migration[example_migration$`Origin County` ==
#'                    example_migration$`Destination County`,"Migration"] * 10
"example_migration"
