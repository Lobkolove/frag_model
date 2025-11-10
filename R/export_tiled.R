library(here)

#' Export a Tiled Matrix Image
#'
#' Creates a larger tiled matrix by repeating a given matrix in a grid layout 
#' and exports it as a PNG image with specified dimensions and color palette.
#'
#' @param matrix The input matrix to be tiled and exported as an image.
#' @param rows Number of tile rows (default is 3).
#' @param cols Number of tile columns (default is 3).
#' @param output_filename File path for the output PNG (default includes deparse substitute).
#' @param width Width of the output image in pixels (default 2400).
#' @param height Height of the output image in pixels (default 2400).
#' @param color_palette Color palette to use for the image (default is viridis with 100 colors).
#'
#' @return No return value; the function saves a PNG file and outputs a message.
#'
#' @examples
#' \dontrun{
#' # Create a sample matrix
#' sample_matrix <- matrix(runif(100), 10, 10)
#' # Export tiled image with default 3x3 layout
#' export_tiled(sample_matrix)
#' # Export tiled image with 4x2 layout and custom filename
#' export_tiled(sample_matrix, rows = 4, cols = 2, output_filename = "pics/custom.png")
#' }
#'
#' @export
export_tiled <- function(
    matrix,
    rows = 3,
    cols = 3,
    output_filename = paste0("pics/", deparse(substitute(grf_matrix)), ".png"),
    width = 2400,
    height = 2400,
    color_palette = viridis(100)
) {
  
  # Get dimensions of the single matrix
  ny_single <- nrow(matrix)
  nx_single <- ncol(matrix)
  
  # Create the larger tiled matrix
  # Initialize an empty matrix of the new total size
  tiled_matrix <- matrix(NA, nrow = ny_single * rows, 
                             ncol = nx_single * cols)
  
  # Fill the tiled matrix by repeating the original matrix
  for (i in 0:(rows - 1)) {
    for (j in 0:(cols - 1)) {
      row_start <- i * ny_single + 1
      row_end <- (i + 1) * ny_single
      col_start <- j * nx_single + 1
      col_end <- (j + 1) * nx_single
      tiled_matrix[row_start:row_end, col_start:col_end] <- matrix
    }
  }
  
  # Export the tiled image
  png(
    filename = output_filename,
    width = width,
    height = height,
    res = 300 # You can adjust resolution if needed for print quality
  )
  
  # Set plotting parameters to remove all fluff
  par(mar = c(0, 0, 0, 0), # No margins
      omi = c(0, 0, 0, 0), # No outer margins
      xaxs = "i", yaxs = "i") # "i" for internal axes to prevent extra space
  
  image(
    x = 1:ncol(tiled_matrix),
    y = 1:nrow(tiled_matrix),
    z = tiled_matrix,
    col = color_palette,
    axes = FALSE,      # No axes
    xlab = "",         # No x-label
    ylab = "",         # No y-label
    main = "",         # No main title
    asp = 1            # Maintain 1:1 aspect ratio
  )
  
  dev.off()
  message(paste("Tiled grid saved to:", output_filename))
}
