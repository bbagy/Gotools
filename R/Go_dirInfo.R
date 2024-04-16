#' Set Working Directory and Return Paths for Image Storage
#'
#' This function sets the working directory to a specified path and returns the
#' paths for storing .pdf images and differential abundance test images.
#' If no directory is specified, it displays the current working directory,
#' lists its contents, and describes the usage of other parameters.
#'
#' @param Current_woringing_dir A character string specifying the path to set
#'   as the current working directory. If \code{NA}, the function will output
#'   the current working directory and its contents.
#' @param Image_locations A subdirectory within the current working directory
#'   where .pdf images are stored. This is optional.
#' @param DA_image_location A subdirectory within the current working directory
#'   where .pdf images for differential abundance tests are stored. This is optional.
#'
#' @return A list containing the paths for the current working directory,
#'   image directory, and differential abundance test image directory.
#'   If \code{Corrent_woringing_dir} is \code{NA}, no list is returned and
#'   function outputs directory information directly.
#'
#' @examples
#' \dontrun{
#'   dir_info <- Go_dirInfo("path/to/directory", "images", "DA_images")
#'   print(dir_info)
#'
#'   # If Corrent_woringing_dir is NA
#'   Go_dirInfo()
#' }
#'
#' @export

Go_dirInfo <- function(Current_working_dir=NA,
                       Image_locations=NA,
                       DA_image_location=NA) {
  # Check if all arguments are missing and print options if they are
  if (is.na(Current_working_dir)) {
    cat(
      "Current_working_dir: Add the current working directory. \n",
      sprintf("Current directory is %s. \n\n", getwd()),

      sprintf("Current directory's contents:\n %s \n\n",paste(list.files(), collapse=", ")),

      "image.wd: Add the directory of .pdf image \n\n",
      "da.wd: Add the directory of .pdf differential abundnat test image.\n\n")
    return(invisible())
  }


  # If 'Current_working_dir' is not NA, then set the working directory to it
  if (!is.na(Current_working_dir)) {
    setwd(Current_working_dir)
  }

  # Construct the full paths for image and DA image locations
  image_path <- if (!is.na(Image_locations)) {
    file.path(getwd(), Image_locations)
  } else {
    NA  # Return NA if Image_locations is not specified
  }

  da_image_path <- if (!is.na(DA_image_location)) {
    file.path(getwd(), DA_image_location)
  } else {
    NA  # Return NA if DA_image_location is not specified
  }

  # Return a list
  return(list(
    currentwd = getwd(),
    image.wd = image_path,
    da.wd = da_image_path
  ))
}





