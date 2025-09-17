
library(xml2)

# some of the data files have an extra empty columns that seem to cause problems
# when trying to parse it using the pzfx package
pzfx_remove_empty_ycolumns <- function(input_file, output_file) {
  doc <- xml2::read_xml(input_file)

  # Find all YColumn nodes
  ycols <- xml2::xml_find_all(doc, ".//YColumn")

  for (yc in ycols) {
    subcols <- xml2::xml_find_all(yc, "./Subcolumn")

    # Check if all Subcolumns are empty
    all_empty <- all(sapply(subcols, function(sc) {
      # Empty if it has no children and no text
      length(xml2::xml_children(sc)) == 0 && xml2::xml_text(sc) == ""
    }))

    if (all_empty) {
      xml2::xml_remove(yc)
    }
  }

  xml2::write_xml(doc, output_file)
}
