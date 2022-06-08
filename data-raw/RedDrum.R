## code to prepare RedDrum dataset goes here
RedDrum <- readxl::read_excel("data-raw/Red Drum Data.xlsx") |>
    dplyr::mutate(Collection = factor(Collection),
                  Random = factor(Random),
                  Year = factor(Year),
                  Month = factor(Month),
                  Site = factor(Site),
                  Stratum = factor(Stratum,
                                   levels = c("Port Royal Sound", "ACE Basin",
                                              "Ashley River",
                                              "Charleston Harbor",
                                              "Wando River", "Cape Romain",
                                              "Winyah Bay"),
                                   ordered = TRUE),
                  Area = factor(Area,
                                levels = c("Port Royal Sound", "ACE Basin",
                                           "Charleston Harbor", "Cape Romain",
                                           "Winyah Bay"),
                                ordered = TRUE),
                  Tide = factor(Tide,
                                levels = c("Early Ebb", "Mid Ebb", "Late Ebb"),
                                ordered = TRUE))
usethis::use_data(RedDrum, overwrite = TRUE)
