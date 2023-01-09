isUnique <-
	function(vector) {
		return(!any(duplicated(vector)))
	}

extractNanoblotData <-
	function(SampleList) {
		NBdf <- data.frame("qname"=SampleList[[1]][['qname']],
											 "qwidth"=SampleList[[1]][['qwidth']])
		return(NBdf)
	}

