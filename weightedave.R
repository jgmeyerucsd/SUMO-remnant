


weightedAve=function(x=c(10,1), weights=c(0.9,0.1)){
	tempsum=0
	for(i in 1:length(weights)){
		tempsum=tempsum+x[i]*(weights[i]/sum(weights))
		}
	print(tempsum)
	}

weightedAve()