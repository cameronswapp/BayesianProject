library(tidyverse)
library(tidytuesdayR)

rm(list = ls())

set.seed(432)

gas_data <- tt_load(2025, week = 26)
gas_data <- gas_data$weekly_gas_prices

drop_forms <- c("conventional", "reformulated")
drop_grades <- c("all", "low_sulfur")

gas_clean <- gas_data |>
  filter(!(grade %in% drop_grades) & !(formulation %in% drop_forms) &
           date > as.Date("2007-02-05")) |>
  select(date, grade, price)
