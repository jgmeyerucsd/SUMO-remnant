


removeLowLocalization=function(object=mgpl.all,minscore=0.75){
	PTMscorelines<-as.character(object@data[,"ptm_peptide"])
	scores<-list()
	keepthese<-c()
	modposition<-list()
	line<-1
	for(i in 1:length(PTMscorelines)){
		tempscore<-as.numeric(unlist(regmatches(PTMscorelines[i],gregexpr("[[:digit:]]+\\.*[[:digit:]]*",PTMscorelines[i]))))
		if(length(which(tempscore>=minscore)>0)>0){
			#print("isnumeric")
			keepthese<-c(keepthese,i)
			n=which(tempscore>=minscore)
			modposition[[line]]<-unlist(gregexpr("[[:digit:]]+\\.*[[:digit:]]*",PTMscorelines[i]))[which(tempscore>=minscore)]-((n-1)*7+2)
			line=line+1
			}
		}
	### now have index of rows to keep
	### and positions as a list of positions
	### convert the list of peptide positions where length>1 into csv

	### this works but the double mods are form of c(2,3)
	#test<-cbind(as.character(modposition),object@data[keepthese,])
	object@modposition.peptide<-modposition
	object@data<-object@data[keepthese,]
	object@pepvec<-object@pepvec[keepthese]
	return(object)
	}
