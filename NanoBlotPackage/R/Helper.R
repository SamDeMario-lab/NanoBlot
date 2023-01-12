#' @name isUnique
#' @title Unique Vector Checker
#' @description This function is required to run NanoBlot. It is not intended
#' for use by an end user. This function takes a vector input and checks that each element
#' in that vector is Unique.
#' @param vector Any vector to check uniqueness on. 
#' Returns true if all elements are unique.
#' @export
#' @examples
#' testVector_1 <- c('alpha', 'beta', 'gamma')
#' testVector_2 <- c('alpha', 'beta', 'alpha')
#' 
#' isUnique(testVector_1)
#'   returns TRUE
#' isUnique(testVector_2)
#'   returns FALSE

isUnique <-
	function(vector) {
		return(!any(duplicated(vector)))
	}

#' @name extractNanoblotData
#' @title Extract Nanoblot Data from Rsamtools::scanBam()
#' @description This function is required to run NanoBlot. It is not intended for use by an end user.
#' @param SampleList A list of the outputs of Rsamtools::scanBam() to extract the qname and qwidth from
#' @export
#' @examples
#' extractNanoblotData(SampleList)
#' 
#' returns a dataframe with 2 columns

extractNanoblotData <-
	function(SampleList) {
		NBdf <- data.frame("qname"=SampleList[[1]][['qname']],
											 "qwidth"=SampleList[[1]][['qwidth']])
		return(NBdf)
	}

