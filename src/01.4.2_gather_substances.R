


library(httr2)
parameters <- yaml::read_yaml("parameters.yaml")


intermediate_path <- "intermediate/01.4.2_gather_substances"
if (!dir.exists(intermediate_path)) {
    cat("Creating intermediate path ", intermediate_path, "\n", sep = "")
    dir.create(intermediate_path, recursive = TRUE)
}


# make sure the server is accessible. E.g. if it is running a remote computer
# you may need to do ssh with port-forwarding make it appear like it is running locally
#
#  ssh -L 8080:localhost:8080 <USERNAME>@<URL_OF_REMOTE_SERVER>
#


source("src/marcus_client.R")

documents <- data.frame(
    pdf_fname = list.files(
        path = "data/KeBOB/documents",
        pattern = "*.pdf",
        full.names = TRUE))


problematic_documents <- data.frame(
    pdf_fname = c(
        # Payload too large
        "data/KeBOB/documents/(Chen, 2022, 10.1021:acs.jcim.2c00962).pdf",
        "data/KeBOB/documents/(Baby, 2024, 10.3389:fphar.2024.1381073).pdf",
        "data/KeBOB/documents/(Gaston, 2020, 10.1172:jci.insight.134174).pdf",
        "data/KeBOB/documents/(Gomes, 2024, 10.1124:molpharm.124.000947).pdf",
        "data/KeBOB/documents/(Lešnik, 2023, 10.1021:acs.jcim.3c00197).pdf",
        "data/KeBOB/documents/(Möller, 2020, 10.1038:s41589-020-0566-1).pdf",
        "data/KeBOB/documents/(O’Brien, 2024, 10.1038:s41586-024-07587-7).pdf",
        "data/KeBOB/documents/(Root-Bernstein, 2022, 10.3390:ph15020214).pdf",
        "data/KeBOB/documents/(Root-Bernstein, 2022, 10.3390:ph15020214).pdf",
        "data/KeBOB/documents/(Sadler, 2023, 10.1038:s41586-023-05789-z).pdf",
        "data/KeBOB/documents/(Scott, 2024, 10.1021:acs.jpcb.4c05214).pdf",
        "data/KeBOB/documents/(Uprety, 2021, 10.7554:eLife.56519).pdf",
        "data/KeBOB/documents/(Viswanadham, 2021, 10.3390:pharmaceutics13070927).pdf",
        "data/KeBOB/documents/(Wang, 2018, 10.1016:j.neuron.2018.03.002).pdf",
        "data/KeBOB/documents/(Wang, 2021, 10.4155:fmc-2020-0308).pdf",
        "data/KeBOB/documents/(Wang, 2023, 10.1016:j.cell.2022.12.026).pdf",
        "data/KeBOB/documents/(Zhang, 2024, 10.1038:s41392-024-01803-6).pdf",
        "data/KeBOB/documents/(Zhuang, 2022, 10.1016:j.cell.2022.09.041).pdf"))
        

documents <- documents |>
    dplyr::anti_join(problematic_documents, by = "pdf_fname")



is_healthy <- marcus_is_healthy(url = parameters$MARCUS_server$url)

document_segments_info <- documents |>
    dplyr::rowwise() |>
    dplyr::do({
        document_info <- .
        cat("Segmenting document '", document_info$pdf_fname[1], "'\n", sep = "")
        tryCatch({
            info <- marcus_decimer_extract_segments(
                url = parameters$MARCUS_server$url,
                pdf_fname = document_info$pdf_fname) |>
                dplyr::mutate(local_pdf_fname = document_info$pdf_fname)

            cat("  Got ", info$segments_count, " segments\n", sep = "")
            info
        }, error = function(err) {
            cat("ERROR: ", err$message, "\n")
            data.frame()
        })
    }) |>
    dplyr::ungroup()



document_segment_files <- document_segments_info |>
    dplyr::filter(segments_count > 0) |>
    dplyr::rowwise() |>
    dplyr::do({
        x <- .
        cat("getting segment files for '", x$segments_directory[1], "'\n", sep = "")
        z <- marcus_list_segment_files(
            url = parameters$MARCUS_server$url,
            directory = x$segments_directory[1])
        z
    }) |>
    dplyr::ungroup()

segment_images <- document_segment_files |>
    dplyr::rowwise() |>
    dplyr::do({
        x <- .
        cat("getting segment image for '", x$segments_directory[1], " ", x$segment_file[1], "'\n", sep = "")
        marcus_get_segment_image(
            url = parameters$MARCUS_server$url,
            directory = x$segments_directory[1],
            image_name = x$segment_file[1])
    })

segment_smiles <- document_segment_files |>
    dplyr::rowwise() |>
    dplyr::do({
        x <- .
        cat("getting segment smiles for '", x$segments_directory[1], " ", x$segment_file[1], "'\n", sep = "")
        tryCatch({
            data <- marcus_generate_smiles(
                url = parameters$MARCUS_server$url,
                directory = x$segments_directory[1],
                image_name = x$segment_file[1],
                engine = "decimer",
                hand_drawn = "false")
        }, error = function(err) {
            cat("ERROR: ", err$message, "\n", sep = "")
            data <- NULL
        })
        data
    })

document_segments <- dplyr::bind_cols(
    document_segment_files,
    segment_smiles) |>
    dplyr::left_join(
        document_segments_info |>
        dplyr::filter(pdf_filename != "_Krotulski__2020__10.1093_jat_bkaa016_.pdf") |>
        dplyr::select(
            segments_directory,
            pdf_filename,
            local_pdf_fname),
        by = "segments_directory")

document_segments |>
    readr::write_tsv(
        file = paste0(intermediate_path, "/document_segments_20251020.tsv"))
        
