test_that("Sanity checks work", {
  # Colour-blind friendly colours are returned
  expect_equal(colour_blind_friendly_cols(7),
               c("#009E73", "#D55E00", "#AA4499", "#0072B2",
                 "#F0E442", "#661100", "#332288"))

  # For n > 20, Spectral colour spectrum should be returned with a warning
  expect_message(colour_blind_friendly_cols(21),
                 "There are too many colours")
  expect_equal(colour_blind_friendly_cols(21),
               grDevices::hcl.colors(palette = "Spectral", n = 21))
})
