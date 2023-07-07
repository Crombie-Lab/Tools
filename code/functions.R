library(tidyverse)
library(easyXpress)
library(rebus)
library(data.table)
library(png)
library(cowplot)
library(gridExtra)

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
  message(glue::glue("The genetic distance between {chrom}:{left} and {chrom}:{right} is {gen.dist1} cM
                     The expected frequency of recombination between these positions is {round(gen.dist1/100, digits = 3)}"))
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


#' viewOverlay
#'
#' View and annotate the CellProfiler overlay image for a single well or an array of wells.
#'
#' @param data a single unsummarized dataframe from \code{easyXpress::process}. The
#'   dataframe should be filtered to contain only the wells to be plotted.
#'   Arrange the dataframe in the order that overlays are to be plotted.
#' @param proc.img.dir a variable name in \code{data} that holds the full PATH to the directory holding processed images matching the data.
#' @param well.label a variable name in \code{data} to display as a well label. For example, \code{"Metadata_Well"}.
#' @param obj.label A variable name in \code{data} to label objects by. For example, \code{"model_select"}.
#' @param obj.shape Optional: A variable name in \code{data} to use for the shape of objects. For example, \code{"model_select"}. The default is NULL.
#' @param obj.shape.pal Optional: A shape palette for objects. This is a vector of shape values with names for each unique value in \code{obj.shape}. NULL passes shape values to \code{obj.label}.
#' @param obj.color Optional: A variable name in \code{data} to color objects by. For example, \code{"model_select"}. The default is NULL.
#' @param obj.col.pal Optional: A color palette for obj.labels. This is a vector of colors with names for each unique value in \code{obj.color}. NULL passes color values to \code{obj.label}.
#' @param text.anno Optional: a variable name in \code{data} used to display text annotation at an object center.
#'    by default this value is shifted vertically to avoid plotting over object center.
#' @param file Optional: The full path for saving output plot. Default is NULL which will throw an error asking for a file path.
#' @return A plot showing all the CellProfiler processed well overlays in
#'   \code{data} with objects annotated as desired.
#' @importFrom imager load.image
#' @importFrom dplyr %>%
#' @importFrom gridExtra arrangeGrob
#' @importFrom pals Kelly
#' @export
#'

viewOverlay <- function(data, proc.img.dir, well.label, obj.label,
                           obj.color = NULL, obj.col.pal = NULL,
                           obj.shape = NULL, obj.shape.pal = NULL,
                           text.anno = NULL, file = NULL) {
  
  # flexibility for FileName_RawBF vs. Image_FileName_RawBF in dat
  if("FileName_RawBF" %in% names(data)){
    data <- data %>%
      rename_at(vars(matches("FileName_RawBF")), ~ "Image_FileName_RawBF") 
  }
  
  # rename vars for later, do here to throw errors early
  plot_data_rename <- data %>%
    rename_at(vars(matches(well.label)), ~ "well_label") %>%
    rename_at(vars(matches(obj.label)), ~ "obj_class")
  # rename vars for optional parameters 
  if(!is.null(obj.color)){
    plot_data_rename <- plot_data_rename %>%
      rename_at(vars(matches(obj.color)), ~ "obj_color")
  }
  if(!is.null(obj.shape)){
    plot_data_rename <- plot_data_rename %>%
      rename_at(vars(matches(obj.shape)), ~ "obj_shape")
  }
  if(!is.null(text.anno)){
    plot_data_rename <- plot_data_rename %>%
      rename_at(vars(matches(text.anno)), ~ "text_anno")
  }
  
  #make file list
  file_list_df <- data %>%
    rename_at(vars(matches(proc.img.dir)), ~ "proc.img.dir") %>%
    dplyr::distinct(Image_FileName_RawBF, .keep_all = T) %>%
    dplyr::mutate(file_path = stringr::str_replace(Image_FileName_RawBF, pattern = ".TIF|.tif", replacement = "_overlay.png"),
                  file_path = paste0(proc.img.dir, file_path))
  
  file_list <- file_list_df %>%
    dplyr::pull(file_path)
  
  # Make array dimensions
  n <- length(file_list)
  nrow <- floor(sqrt(n)) # goes in nrow
  ncol <- ceiling(n/nrow) # do we need ceiling?
  max_dim <- max(nrow, ncol)
  
  # send a message that we're making a grob list
  message(glue::glue("Making grob list for {n} overlays"))
  
  # make a grob list
  img_grob_list <- list()
  
  # loop through file list and add to grob list (80sec/36overlays, when local)
  for(i in unique(file_list)) {
    # get the image
    img <- png::readPNG(glue::glue("{i}"))
    grob <- grid::rasterGrob(img) #width=ggplot2::unit(1,"npc"), height=ggplot2::unit(1,"npc"))
    # add to list
    img_grob_list[[i]] <- grob
  }
  
  # get image dimensions from last img
  h<-dim(img)[1] # image height
  w<-dim(img)[2] # image width
  
  # Make gtable with layout
  img_grid <- do.call("arrangeGrob", c(img_grob_list, ncol=ncol, nrow=nrow))
  
  # Get array position df
  array_pos_df <- file_list_df %>%
    dplyr::mutate(array_x_pos = rep(0:(ncol-1), length.out = n()),
                  array_y_pos = rep((nrow-1):0, each = ncol, length.out = n())) %>%
    dplyr::select(Image_FileName_RawBF, array_x_pos, array_y_pos) # just take what's needed
  
  # Get plot df for adding points
  plot_data_df <- plot_data_rename %>%
    dplyr::left_join(array_pos_df, by = "Image_FileName_RawBF") %>%
    dplyr::mutate(obj_center_x = abs(AreaShape_Center_X) * (1/w) + array_x_pos,
                  obj_center_y = abs(AreaShape_Center_Y - h) * (1/h) + array_y_pos) # w = img width, h = img height
  
  # Generate plot message
  message(glue::glue("Making plot for {n} overlays"))
  
  # make color palette for plotting
  if(is.null(obj.col.pal) & is.null(obj.color)){
    # get a typical toxin color palette
    col_pal <- pals::kelly()[3:22] # get a list of colors from kelly
    obj_color_pal <- col_pal[1:length(unique(plot_data_df$obj_class))] 
    names(obj_color_pal) <- unique(plot_data_df$obj_class)
  }
  
  if(is.null(obj.col.pal) & !is.null(obj.color)){
    n_colors <- length(unique(plot_data_df$obj_color))
    obj_col_pal <- seq(1:n_colors)
    names(obj_col_pal) <- unique(plot_data_df$obj_color)
  }
  
  if(!is.null(obj.col.pal) & is.null(obj.color)){
    obj_color_pal <- obj.col.pal
  }
  
  if(!is.null(obj.col.pal) & !is.null(obj.color)){
    obj_color_pal <- obj.col.pal
  }
  
  # make shape palette for plotting
  if(is.null(obj.shape.pal) & is.null(obj.shape)){
    obj_shape_pal <- rep(21, times = length(unique(plot_data_df$obj_class))) # set all to 21 if no shape palette provided
    names(obj_shape_pal) <- unique(plot_data_df$obj_class)
  }
  if(is.null(obj.shape.pal) & !is.null(obj.shape)){
    n_shapes <- length(unique(plot_data_df$obj_shape))
    obj_shape_pal <- seq(21:n_shapes)
    names(obj_shape_pal) <- unique(plot_data_df$obj_shape)
  }
  if(!is.null(obj.shape.pal) & is.null(obj.shape)){
    obj_shape_pal <- obj.shape.pal
    # ADD a catch for mismatch here
  }
  if(!is.null(obj.shape.pal) & !is.null(obj.shape)){
    obj_shape_pal <- obj.shape.pal
    # ADD a catch for mismatch here
  }
  
  # Make plot
  p <- ggplot(plot_data_df) +
    xlim(0, ncol) +
    ylim(-0.25, nrow) +
    {if(is.null(obj.color) & is.null(obj.shape))
      ggplot2::aes(x = obj_center_x, y = obj_center_y, color = obj_class, shape = obj_class)
    } +
    {if(!is.null(obj.color) & is.null(obj.shape))
      ggplot2::aes(x = obj_center_x, y = obj_center_y, color = obj_color, shape = obj_class)
    } +
    {if(is.null(obj.color) & !is.null(obj.shape))
      ggplot2::aes(x = obj_center_x, y = obj_center_y, color = obj_class, shape = obj_shape)
    } +
    {if(!is.null(obj.color) & !is.null(obj.shape))
      ggplot2::aes(x = obj_center_x, y = obj_center_y, color = obj_color, shape = obj_shape)
    } +
    ggplot2::annotation_custom(img_grid, 0, ncol, 0, nrow) +
    ggplot2::scale_shape_manual(name = "Object", values = obj_shape_pal) +
    ggplot2::scale_color_manual(name = "Object", values = obj_color_pal) +
    ggplot2::geom_point(alpha = 0.5, size = 1, stroke = 0.1) + # n = image number (n/12), stroke = 0 should remove outlines
    ggplot2::geom_label(inherit.aes = F,
                        aes(x = array_x_pos + 0.5, y = array_y_pos + 1, label = well_label),
                        label.padding = unit(.05, "line"), size = 1, fill = "white", color = "black",
                        label.size = NA, show.legend = F) + # n = image number (n/6)
    ggplot2::coord_equal() +
    ggplot2::theme_void() +
    theme(legend.justification = 'center',
          legend.direction = "horizontal",
          legend.position = c(0.5, (1/(nrow+1))*.25),
          legend.key.width = unit(1, "line"), #unit(n/12, "line"),
          legend.key.size = unit(1, "line"), #unit(n/12, "line"),
          legend.title = element_text(size = 2), #element_text(size = n/3),
          legend.text = element_text(size = 2)) #element_text(size = n/3))
  
  # add text annotation to objects
  if(!is.null(text.anno)){
    p <- p +
      ggplot2::geom_text(aes(label = text_anno), size = 0.5, nudge_y = 0.03, show.legend = F)
  }
  
  # Save it
  message(glue::glue("Saving plot as: {file}"))
  ggsave(file = file, plot = p, width = ncol, height = nrow + 0.25, dpi = 2048)
  message(glue::glue("DONE"))
  
}

#' viewObjects
#'
#' View one or more CellProfiler objects with annotation.
#'
#' @param data a single unsummarized dataframe from \code{easyXpress::process}. The
#'   dataframe should be filtered to contain only the objects to be plotted.
#' @param proc.img.dir A variable name in \code{data} with the path to the directory holding processed images matching the data.
#' @param out.dir A variable name in \code{data} with the path to the desired output directory.
#' @param file.names A variable name in \code{data} with a unique name for saving the output plots. A .png suffix will be appended for saving files.
#' @param cols A character vector of column names to be added to the plot as a table, e.g. \code{c("modle_select", "AreaShape_Area")}. 
#' @param prompt Logical, ask for permission to process objects? Default is \code{TRUE}.
#' @return A number of plots are saved, one for each object in \code{data}. A .csv file is also saved with file.names corresponding to each plot. This file can be used to help record notes on the objects being plotted.
#' @importFrom imager load.image
#' @importFrom dplyr %>%
#' @importFrom ggfittext geom_fit_text
#' @importFrom pals Kelly
#' @export
#'

viewObjects <- function(data, proc.img.dir, out.dir, file.names, cols, prompt = TRUE){
  
  if(prompt == TRUE){
    # get permission to make directory and object images
    permission = readline(prompt = glue::glue("Do you want to create {data %>%
                                              rename_at(vars(matches(out.dir)), ~ 'out.dir') %>%
                                              dplyr::distinct(out.dir) %>%
                                              dplyr::pull(out.dir)} and save {nrow(data)} object images to it (y/n): "))
  }
  
  if(prompt == FALSE){
    permission = "yes"
  }
  
  if (permission %in% c("Y", "yes", "y")) {
    message(glue::glue("processing {nrow(data)} objects"))
    # flexibility for FileName_RawBF vs. Image_FileName_RawBF in data
    if("FileName_RawBF" %in% names(data)){
      data.n <- data %>%
        rename_at(vars(matches("FileName_RawBF")), ~ "Image_FileName_RawBF") 
    }
    
    # get the columns we care about
    wanted <- data.n %>%
      dplyr::select(dplyr::matches(cols) |
                      dplyr::matches(c(proc.img.dir, out.dir, file.names,
                                       "Image_FileName_RawBF",
                                       "AreaShape_Center_X",
                                       "AreaShape_Center_Y",
                                       "po_AreaShape_BoundingBoxMaximum_X",
                                       "po_AreaShape_BoundingBoxMinimum_X",
                                       "po_AreaShape_BoundingBoxMinimum_Y",
                                       "po_AreaShape_BoundingBoxMaximum_Y")))
    
    # rename vars and add output paths
    d_full <- wanted %>%
      rename_at(vars(matches(proc.img.dir)), ~ "proc.img.dir") %>%
      rename_at(vars(matches(out.dir)), ~ "out.dir") %>%
      rename_at(vars(matches(file.names)), ~ "file.names") %>%
      dplyr::mutate(ol_file = stringr::str_replace(Image_FileName_RawBF, pattern = ".TIF|.tif", replacement = "_overlay.png"),
                    file.names = stringr::str_replace(file.names, pattern = ".png|.jpeg|.TIF|.tif", replacement = ""),
                    ol_path = paste0(proc.img.dir, ol_file),
                    out_path = paste0(out.dir, file.names, ".png")) # ensure .png file
    
    # Create the progress bar
    pb <- progress::progress_bar$new(total = nrow(d_full))
    
    # loop through objects and save
    for(i in 1:nrow(d_full)){
      # get the row and 
      d <- d_full %>%
        dplyr::slice(i) %>%
        dplyr::mutate(raw.y.dim = po_AreaShape_BoundingBoxMaximum_Y - po_AreaShape_BoundingBoxMinimum_Y,
                      raw.x.dim = po_AreaShape_BoundingBoxMaximum_X - po_AreaShape_BoundingBoxMinimum_X,
                      max.dim = max(raw.x.dim, raw.y.dim),
                      exp.fact = case_when(max.dim <= 100 ~ 4,
                                           max.dim > 100 & max.dim <= 200 ~ 2,
                                           max.dim > 200 & max.dim <= 300 ~ 1.5,
                                           max.dim > 300 ~ 1.5), # TAC edited from 1 2023-04-02 testing
                      new.y.min = po_AreaShape_Center_Y - (raw.y.dim * exp.fact) / 2, ####CONTROLS TO AVOID PLOTING OFF THE IMAGE DIM?
                      new.y.max = po_AreaShape_Center_Y + (raw.y.dim * exp.fact) / 2,
                      new.x.min = po_AreaShape_Center_X - (raw.x.dim * exp.fact) / 2,
                      new.x.max = po_AreaShape_Center_X + (raw.x.dim * exp.fact) / 2)
      
      
      # get the overlay image and crop it to object
      #img <- png::readPNG(d$ol_path[1])[d$po_AreaShape_BoundingBoxMinimum_Y[1]:d$po_AreaShape_BoundingBoxMaximum_Y[1], d$po_AreaShape_BoundingBoxMinimum_X[1]:d$po_AreaShape_BoundingBoxMaximum_X[1], ]
      img <- png::readPNG(d$ol_path[1])[d$new.y.min[1]:d$new.y.max[1], d$new.x.min[1]:d$new.x.max[1], ]
      
      # get the dims
      h<-dim(img)[1] # image height
      w<-dim(img)[2] # image width 
      
      # plot it
      obj_img <- ggplot2::ggplot(d) +
        ggplot2::aes(x = AreaShape_Center_X - new.x.min,
                     y = AreaShape_Center_Y - new.y.min) +
        ggplot2::annotation_custom(grid::rasterGrob(img, width = ggplot2::unit(1, "npc"), height = ggplot2::unit(1, "npc")), 0, w, 0, -h) +
        ggplot2::geom_point(shape = 4, alpha = 1, size = 4, stroke = 0.25, color = "black") +
        ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(0, w)) +
        ggplot2::scale_y_reverse(expand = c(0, 0), limits = c(h, 0)) +
        ggplot2::coord_equal() +
        theme_void() +
        theme(legend.position = "none")
      
      # make a nice table
      d_table <- d %>%
        dplyr::select(dplyr::matches("file.names") | dplyr::matches(cols)) %>%
        dplyr::select(file = file.names, everything())
      lt <- tibble::tibble(names = names(d_table),
                           dat = d_table %>%
                             data.table::transpose() %>%
                             dplyr::pull(V1)) %>%
        dplyr::mutate(y = 1:n(),
                      label = paste(names, " = ", dat)) 
      
      # plot table
      tb <- ggplot(lt, aes(x = 1, y = dat, label = label)) +
        geom_tile(fill = "white", colour = "black") +
        ggfittext::geom_fit_text(place = "left", min.size = 1, padding.x = grid::unit(0, "mm")) +
        theme_void()
      
      # add the object and the table together
      full <- cowplot::plot_grid(obj_img, tb, align = "vh", axis = "tblr")
      
      # save it
      cowplot::ggsave2(full, filename = glue::glue("{d$out_path}"), width = 4, height = 2, dpi = case_when(d$exp.fact[1] == 4 ~ 600,
                                                                                                           d$exp.fact[1] == 2 ~ 450,
                                                                                                           d$exp.fact[1] == 1.5 | d$exp.fact[1.5] == 1  ~ 300))
      
      # Update the progress bar
      pb$tick()
    }
    
    # save the .csv file for scoring
    rio::export(d_full %>%
                  dplyr::select(dplyr::matches("file.names") | dplyr::matches(cols)) %>%
                  dplyr::select(file = file.names, everything()),
                file = glue::glue("{unique(d_full$out.dir)}dat_out.csv"))
    
    message(glue::glue("DONE: saved object images and dat_out.csv to {unique(d_full$out.dir)}"))
  }
  
  else{
    message("QUIT viewObject")
  }
  
}
