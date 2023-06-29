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

# get the genetic distance between two points in the genome
# could make species specific with another parameter.
geneticDistance <- function(left, right, chrom) {
  # message
  message("getting C. elegans genetic map from https://github.com/AndersenLab/post-gatk-nf/raw/main/input_files/annotations/c_elegans/c_elegans_genetic_map.bed.gz")
  # get bed
  bed_data <- data.table::fread("https://github.com/AndersenLab/post-gatk-nf/raw/main/input_files/annotations/c_elegans/c_elegans_genetic_map.bed.gz") %>%
    dplyr::select(CHROM = V1, Start = V2, End = V3, cM = V4)
  
  # get cM distance between our markers
  dist1 <- bed_data %>%
    dplyr::filter(CHROM == chrom) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(bin.median = median(c(Start, End))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(left1.dist = abs(left - bin.median),
                  right1.dist = abs(right - bin.median))
  
  l1cM <- dist1 %>%
    dplyr::arrange(left1.dist) %>%
    dplyr::slice(1)
  
  r1cM <- dist1 %>%
    dplyr::arrange(right1.dist) %>%
    dplyr::slice(1)
  
  gen.dist1 <- r1cM$cM - l1cM$cM
  # return it
  message(glue::glue("The genetic distance between {chrom}:{left} and {chrom}:{right} is {gen.dist1}"))
  return(gen.dist1)
}

# test the probability of observing n recombinants
# n.wells is the number of F2 wells in your cross
# n.recombinants is the number of recombinants observed across n wells
# g.dist is the genetic distance in cM between genetic markers
recombinantProb <- function(n.wells, n.recombinants, g.dist){
  
  # convert to prob from genetic dist
  r.freq <- g.dist/100
  
  # Calculate the binomial probability
  prob <- choose(n.wells, n.recombinants) * r.freq^n.recombinants * (1-r.freq)^(n.wells-n.recombinants)
  
  # Print the result
  message(glue::glue("The probability of getting {n.recombinants} recombinants in {n.wells} wells given a genetic distance between markers of {g.dist} cM is: {prob}"))
  
  # return the result to for assignment
  return(prob)
}
