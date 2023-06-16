library(tidyverse)

# Function to make design file for HTA
# USAGE:
# test <- makeDesign(n.plate = 45, n.row = 8, n.col = 12, assay.name = "LM1", assay.type = "48h", food = "15hHB101_20220727", od = 10)

makeDesign <- function(n.plate, n.row, n.col, assay.name = NULL, assay.type = NULL, food = NULL, od = NULL) {
  
  # setup a plate with correct padding
  Metadata_Plate <- paste0("p", stringr::str_pad(seq(1:n.plate), width = 3, side = "left", pad = "0"))
  row <- LETTERS[1:n.row]
  col <- stringr::str_pad(seq(1:n.col), width = 2, side = "left", pad = "0")
  
  # make a design file of appropriate length
  design <- tidyr::crossing(Metadata_Plate, row, col) %>%
    dplyr::mutate(assay_name = ifelse(is.null(assay.name), NA_character_, assay.name),
                  assay_type = ifelse(is.null(assay.type), NA_character_, assay.type),
                  Metadata_Well = paste0(row, col),
                  plate = as.numeric(stringr::str_replace_all(Metadata_Plate, pattern = "^p00|^p0|^p", replacement = "")),
                  row,
                  col = stringr::str_remove(col, pattern = "^0"),
                  food = ifelse(is.null(food), NA_character_, food),
                  od = ifelse(is.null(od), NA_integer_, od),
                  bleach = NA_integer_,
                  strain = NA_character_,
                  drug = NA_character_,
                  concentration_um = NA_real_,
                  diluent = NA_character_,
                  well_censor = NA_character_,
                  well_censor_reason = NA_character_,
                  notes = NA_character_) %>%
    dplyr::select(assay_name,
                  assay_type,
                  Metadata_Plate,
                  Metadata_Well,
                  plate,
                  row,
                  col,
                  strain,
                  drug,
                  concentration_um,
                  bleach,
                  diluent,
                  food,
                  od,
                  well_censor,
                  well_censor_reason,
                  notes)
  # out
  return(design)
}