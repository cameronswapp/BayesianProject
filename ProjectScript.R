library(tidyverse)
library(tidytuesdayR)

rm(list = ls())

set.seed(432)

gas_data <- tt_load(2025, week = 26)
gas_data <- gas_data$weekly_gas_prices

regular <- gas_data |>
  filter(fuel == "gasoline", grade == "all", formulation == "all") |>
  select(date, price)
