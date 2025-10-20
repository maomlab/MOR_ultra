library(httr2)
library(curl)

#' Check if the MARCUS server is healthy
#'
#' @description Returns TRUE if healthy
#'
#' @param url character url for MARCUS server
#' @export
marcus_is_healthy <- function(url) {
    response <- url |>
        paste0("/health") |>
        httr2::request() |>
        httr2::req_perform()
    response$status_code == 200
}

#' Use DECIMER to segment a given PDF file
#'
#' @description given a pdf filename it will create a directory of segmented
#'   images of compounds in the output it gives the directory path on the
#'   server for the segmented images that can used to retrieve them
#'
#'   Takes about 1 minute per document
#'
#' @param url character url for MARCUS server
#' @param pdf_fname character path to pdf file
#'
#' @return data.frame with columns
#'   * segments_extracted
#'   * segments_already_existed
#'   * segments_count
#'   * segments_directory
#'
#' @export
marcus_decimer_extract_segments <- function(url, pdf_fname) {
    response <- url |>
        paste0("/api/decimer/extract_segments") |>
        httr2::request() |>
        httr2::req_body_multipart(
            collect_all = as.character(TRUE),
            pdf_file = curl::form_file(pdf_fname, type = "application/pdf")) |>
        httr2::req_perform()
    data <- response |>
        httr2::resp_body_json() |>
        as.data.frame()

    if ("segments_directory" %in% names(data)) {
        data <- data |>
            dplyr::mutate(
                segments_directory = segments_directory |>
                    stringr::str_replace("/all_segments$", ""))
    }
    data
}

#' List all the directories of segmented images on the server
#'
#'
#' @param url character url for MARCUS server
#'
#' @return data.frame with columns [directory, segments_count, has_all_segments]
#'
#' @export
marcus_list_segments <- function(url) {
    response <- url |>
        paste0("/api/decimer/list_segments") |>
        httr2::request() |>
        httr2::req_perform()
    response |>
        httr2::resp_body_json() |>
        dplyr::bind_rows()
}

#' List the segment files for a given directory on the server
#'
#'
#' @param url character url for MARCUS server
#' @param directory character directory parameter to access pre-computed
#'    segments on the server for a givne document
#'
#' @return list
#' @export
marcus_list_segment_files <- function(url, directory) {
    response <- url |>
        paste0("/api/decimer/list_directory/", directory, "/all_segments") |>
        httr2::request() |>
        httr2::req_perform()
    z <- response |>
        httr2::resp_body_json()
    data.frame(
        directory = directory,
        segment_file = z$files |> unlist(),
        stringsAsFactors = FALSE)
}

#' Get segmented image given the directory and image path
#'
#'
#' @param url character url for MARCUS server
#' @param directory character path to segments directory on MARCUS server
#' @param image_name character path to image name on MARCUS server
#'
#' @return data.frame with columns
#'   * segments_directory
#'   * image_name
#'   * image_bytes
#'
#' @export
marcus_get_segment_image <- function(url, directory, image_name) {
    response <- url |>
        paste0("/api/decimer/get_segment_image/", directory, "/", image_name) |>
        httr2::request() |>
        httr2::req_perform()

    if (httr2::resp_header(response, "Content-Type") == "image/png") {
        image_bytes <- response |>
            httr2::resp_body_raw()
        return(tibble::tibble(
            segments_directory = directory,
            image_name = image_name,
            image_bytes = image_bytes))
    } else {
        message("The response did not have Content-Type: image/png")
        return()
    }
}
#' Generate smiles from image on MARCUS server
#'
#' @description Given a segments directory and image name, gneerate a smiles
#' with the given engine and indicate if the image is hand-drawn or not.
#'
#' @param url character url for MARCUS server
#' @param directory character path to segments directory on MARCUS server
#' @param image_name character path to image name on MARCUS server
#' @param engine character one of ['deciminer', 'molnextr', or 'molscribe']
#'
#' @return data.frame
#'
#' @export
marcus_generate_smiles <- function(url, directory, image_name, engine, hand_drawn) {
    response <- url |>
        paste0("/api/latest/ocsr/generate_smiles") |>
        httr2::request() |>
        httr2::req_body_multipart(
            image_path = paste0(directory, "/", image_name),
            engine = engine,
            hand_drawn =  as.character(hand_drawn)) |>
        httr2::req_perform()
    response |>
        httr2::resp_body_json() |>
        as.data.frame()
}
