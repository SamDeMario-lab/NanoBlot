isUnique <-
	function(vector) {
		return(!any(duplicated(vector)))
	}

extractNanoblotData <-
	function(SampleList, SampleNames, index) {
		NBdf <- data.frame("qname"=SampleList[[index]][[1]][['qname']],
											 "qwidth"=SampleList[[index]][[1]][['qwidth']],
											 SampleID = as.factor(SampleNames[index]))
		return(NBdf)
	}
