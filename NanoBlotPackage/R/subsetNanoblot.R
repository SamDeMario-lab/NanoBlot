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
                           targetAntiProbes = NULL,
                           viewingWindow = NULL,
                           cDNA = FALSE,
                           RTPCR = FALSE,
                           RACE = FALSE) {
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

  BamFileListNames = names(BiocGenerics::path(BamFileList))
  if (!isUnique(BiocGenerics::path(BamFileList))) {
    stop("BamFileList paths contain non-unique names. All file paths must be unique.") }

  if (!isUnique(BamFileListNames)) {
    stop("BamFileList names contain non-unique names. All names must be unique.") }

  # TO DO
  # If temp directory does not exist --> create it, OR, we can let the user specify what they want as the temp folder
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

    dataLocation <- ""
    tempName <- ""
    for (bamFileIndex in seq_along(BamFileList)) {
      bamFilePath <- BiocGenerics::path(BamFileList[[bamFileIndex]])
      sampleName <- strsplit(names(BiocGenerics::path(BamFileList))[[bamFileIndex]], split = "[.]")[[1]][[1]]

      # checks if its an empty string
      if (previousProbe != "") {
        dataLocation <- paste("./temp/",sampleName,"_",previousProbe,".bam",sep = "")
        tempName <- paste(sampleName,"_",previousProbe,"_",probe,".bam",sep = "")
      }
      else {
        dataLocation <- bamFilePath
        tempName <- paste(sampleName,"_",probe,".bam", sep = "")
      }
      print(paste("Subsetting: ",dataLocation, sep = ""))
      print(paste("Naming Subset: ",tempName, sep = ""))

      outFile <- paste("./temp/", tempName, sep = "")
      if (cDNA == TRUE) {
        system2("bedtools",
                args=c("intersect", "-a", dataLocation,
                       "-b", "./temp/temp_bed.bed", "-wa", "-split", "-nonamecheck"),
                stdout = outFile)
      }
      else
      {
        system2("bedtools",
                args=c("intersect", "-a", dataLocation,
                       "-b", "./temp/temp_bed.bed", "-wa", "-split", "-s", "-nonamecheck"),
                stdout = outFile)
      }
      system2("samtools",
              args=c("index", outFile))
    } #end for loop through each bam file in BamFileList

    if (previousProbe != "")
    { previousProbe <- paste(previousProbe, "_", probe, sep = "") }
    else { previousProbe <- probe }

  } #end for loop through each probe

  # Then for loop through each antitarget probe
  previousAntiProbe <- previousProbe
  if (is.null(targetAntiProbes)) {
    print("No negative probe used")
  }
  else {print("Negative Probe used")}

  for (antiprobe in targetAntiProbes) {
    if (!(antiprobe %in% probesData[[4]])){
      message(paste("Antiprobe: ", antiprobe, " not found. Check bed file or antiprobe spelling. Exiting"))
      return(NULL)
    }

    antiProbeLine <- dplyr::filter(probesData, probesData[4] == antiprobe)
    print(paste("Antiprobe:", antiprobe, sep = ""))
    print(paste(antiProbeLine[[1]], ":", antiProbeLine[[2]], "-", antiProbeLine[[3]], sep = ""))
    write.table(antiProbeLine, file = "./temp/temp_anti_bed.bed",
                sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

    dataLocation <- ""
    tempName <- ""
    for (bamFileIndex in seq_along(BamFileList)) {
      sampleName <- strsplit(names(BiocGenerics::path(BamFileList))[[bamFileIndex]], split = "[.]")[[1]][[1]]

      dataLocation <- paste("./temp/",sampleName,"_",previousAntiProbe,".bam",sep = "")
      tempName <- paste(sampleName,"_",previousAntiProbe,"_anti_",antiprobe,".bam",sep = "")

      print(paste("Subsetting: ",dataLocation, sep = ""))
      print(paste("Naming Subset: ",tempName, sep = ""))
      outFile <- paste("./temp/", tempName, sep = "")
      if (cDNA == TRUE) {
        system2("bedtools",
                args=c("intersect", "-a", dataLocation,
                       "-b", "./temp/temp_anti_bed.bed", "-wa", "-split", "-v", "-nonamecheck"),
                stdout = outFile)
      }
      else
      {
        system2("bedtools",
                args=c("intersect", "-a", dataLocation,
                       "-b", "./temp/temp_anti_bed.bed", "-wa", "-split", "-v", "-s", "-nonamecheck"),
                stdout = outFile)
      }
      system2("samtools",
              args=c("index", outFile))
    } #end for loop through each bam file in BamFileList

    previousAntiProbe <- paste(previousAntiProbe, "_anti_", antiprobe, sep = "")
  }

  # If RT-PCR is true, then go through this loop here
  # Includes the RACE option too
  if (RTPCR == TRUE)
  {
    BUFFER_SIZE <- 5
    PARIS_JAPONICA <- 149000000000
    print("=======")
    print(paste("Running viewing window ", viewingWindow, " subset now for RT-PCR or RACE mode"))

    if (!(viewingWindow %in% probesData[[4]])){
      message(paste("Viewing window: ", viewingWindow, " not found. Check probes file. Exiting"))
      return(NULL)
    }

    viewingWindowLine <- dplyr::filter(probesData, probesData[4] == viewingWindow)
    WINDOW_START <- viewingWindowLine[[2]]
    WINDOW_END <- viewingWindowLine[[3]]
    STRAND <- viewingWindowLine[[6]]

    dataLocation <- ""
    tempName <- ""
    for (bamFileIndex in seq_along(BamFileList)) {
      sampleName <- strsplit(names(BiocGenerics::path(BamFileList))[[bamFileIndex]], split = "[.]")[[1]][[1]]

      dataLocation <- paste("./temp/",sampleName,"_",previousAntiProbe,".bam",sep = "")
      tempDataLocation <- "./temp/RTPCR_temp.bam"
      file.copy(dataLocation, tempDataLocation, overwrite = TRUE)

      if (RACE == TRUE) {
        if (STRAND == "+") {
          # If strand of RACE viewing window is positive, we want inclusive start
          start_frame <- data.frame(viewingWindowLine[[1]],
                                    WINDOW_START - BUFFER_SIZE,
                                    WINDOW_START)
          write.table(start_frame, file = "./temp/temp_start.bed",
                      sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
          system2("bedtools",
                  args=c("intersect", "-a", dataLocation,
                         "-b", "./temp/temp_start.bed", "-wa", "-split", "-nonamecheck"),
                  stdout = tempDataLocation)
          rm <- file.remove("./temp/temp_start.bed")
        }
        else if (STRAND == "-") {
          # If strand of RACE viewing window is negative, we want inclusive end
          end_frame <- data.frame(viewingWindowLine[[1]],
                                    WINDOW_END,
                                    WINDOW_END + BUFFER_SIZE)
          write.table(end_frame, file = "./temp/temp_end.bed",
                      sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
          system2("bedtools",
                  args=c("intersect", "-a", dataLocation,
                         "-b", "./temp/temp_end.bed", "-wa", "-split", "-nonamecheck"),
                  stdout = tempDataLocation)
          rm <- file.remove("./temp/temp_end.bed")
        }
        else {
          message("Correct strand for RACE not found. Check probes file. Exiting")
          return(NULL)
        }
      }
      else {
        start_frame <- data.frame(viewingWindowLine[[1]],
                                  WINDOW_START - BUFFER_SIZE,
                                  WINDOW_START)
        write.table(start_frame, file = "./temp/temp_start.bed",
                    sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
        end_frame <- data.frame(viewingWindowLine[[1]],
                                WINDOW_END,
                                WINDOW_END + BUFFER_SIZE)
        write.table(end_frame, file = "./temp/temp_end.bed",
                    sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
        system2("bedtools",
                args=c("intersect", "-a", tempDataLocation,
                       "-b", "./temp/temp_start.bed", "-wa", "-split", "-nonamecheck"),
                stdout = dataLocation)
        system2("bedtools",
                args=c("intersect", "-a", dataLocation,
                       "-b", "./temp/temp_end.bed", "-wa", "-split", "-nonamecheck"),
                stdout = tempDataLocation)
        rm <- file.remove(c("./temp/temp_start.bed", "./temp/temp_end.bed"))
      }
      #Performing the bedtools complement then running ampliconclip
      print(paste("Clipping ", sampleName, " to ", viewingWindow))
      amplicon_frame <- data.frame(a = c(viewingWindowLine[[1]],viewingWindowLine[[1]]),
                                b = c(0, WINDOW_END),
                                c = c(WINDOW_START, PARIS_JAPONICA))
      write.table(amplicon_frame, file = "./temp/temp.bed",
                  sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
      system2("samtools",
              args=c("ampliconclip", "--hard-clip", "--both-ends",
                     "-b", "./temp/temp.bed", tempDataLocation, "|",
                     "samtools", "sort"),
              stdout = dataLocation)
      system2("samtools",
              args=c("index", dataLocation))
      rm <- file.remove(tempDataLocation)
    } #end for loop through each bam file in BamFileList

  } # end if statement for RT-PCR
print("=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=")
}
