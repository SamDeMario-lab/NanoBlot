#' @name subsetNanoblot
#' @title Subsetting Bam files for Nanoblot generation
#' @description Essentially calls bedtools intersect and samtools ampliconclip with various different parameters
#'
#' @param
#' @param
#' @export
#' @examples
#'
#'

subsetNanoblot <- function(BamFileList,
                           probesFile,
                           targetProbes,
                           targetAntiProbes = NULL) {
  # So basically, we skip the bash script where it deals with plotting files and metadata file
  # locations --> we only really just need a probes bed file

  # Instead of taking in a metadata file --> we can just take in a BamFileList
  # And essentially, there are 3 modes of subsetting, basic mode, RT-PCR mode, and RACE mode
  # There is the additional argument of whether the reads are treated as cDNA or not
  # And there is the option to clear all files from ./temp/ after plot generation

  probesData <- tryCatch({
    read.delim(probesFile, sep = "\t", header = FALSE)},
    error=function(cond) {
      message(paste("probes bed file caused an error:", probesFile))
      message("Here's the original error message:")
      message(cond)
      return(NULL)
    },
    warning=function(cond) {
      message(paste("probes bed caused a warning:", probesFile))
      message("Here's the original warning message:")
      message(cond)
      return(NULL)
    }
    )
  # maybe we can have some sort of bed check for probes file too if we want
  if (!identical(probesData[[4]], unique(probesData[[4]]))) {
    message(paste("Duplicate probes found in ", probesFile, ", please fix and rerun."))
    return(NULL)
  }

  # If temp directory does not exist --> create it
  # Since we are running the markdown folder in the Nanoblot directory, this is okay for now

  previousProbe <- ""
  # For loop through each target probe
  for (probe in targetProbes) {
    if (!(probe %in% probesData[[4]])){
      message(paste("Probe: ", probe, " not found. Check bed file or probe spelling. Exiting"))
      return(NULL)
    }

    probeLine <- dplyr::filter(probesData, probesData[4] == probe)
    print(paste("Probe:", probe, sep = ""))
    print(paste(probeLine[[1]], ":", probeLine[[2]], "-", probeLine[[3]], sep = ""))
    write.table(probeLine, file = "./temp/temp_bed.bed",
                sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

    BamFileListNames = names(BiocGenerics::path(BamFileList))
    if (!isUnique(BiocGenerics::path(BamFileList))) {
      stop("BamFileList paths contain non-unique names. All file paths must be unique.")
    }

    if (!isUnique(BamFileListNames)) {
      stop("BamFileList names contain non-unique names. All names must be unique.")
    }

    for (bamFileIndex in seq_along(BamFileList)) {
      bamFilePath <- BiocGenerics::path(BamFileList[[bamFileIndex]])
      sampleName <- strsplit(names(BiocGenerics::path(BamFileList))[[bamFileIndex]], split = "[.]")[[1]][[1]]
      print(sampleName)
    }

  }

  # Then for loop through each antitarget probe

  # If RT-PCR is true, then go through this loop here
  # Includes the RACE option too


}
