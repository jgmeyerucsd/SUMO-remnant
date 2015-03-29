##################################################################################
### this function is for cleaning the peptide sequences output by IDpicker
########3
#########		accepts inputfile and returns unique sequences in character form


cleanPeps=function(object=MGpl.norm[[2]],type="PeptideProphet"){
	object@pepraw<-as.character(object@data$peptide)
	peptempvec<-as.character(object@data$peptide)
	peptempvec.l<-length(peptempvec)
	###loop over all the strings in pepvec
	peps<-c(rep(0,times=peptempvec.l))
	if(type=="PeptideProphet"){
		for(i in 1:peptempvec.l){
			###	remove the first 2 and last 2 characters
			peptempvec[i]<-substr(peptempvec[i],start=3,stop=nchar(peptempvec[i])-2)
			}
		for(i in 1:peptempvec.l){
			###match the string starting with [/[- followed by n integers and "]"
			### and replace with nothing
			peps[i]<-gsub("(\\[|\\[-)([0-9]+)(.)([0-9]+)(]+)", replacement="", peptempvec[i]) 
			}
		for(i in 1:peptempvec.l){
			###remove "n" from peptides with n-term acetyl
			if(unlist(strsplit(peps[i],split=""))[1]=="n"){
				peps[i]<-substr(peps[i],start=2,stop=nchar(peps[i]))
				}
			}
		}	
	#### section to add digly position in peptide
	object@pepraw
	test<-as.list(rep(0,times=peptempvec.l))
	typeof(test)
	test
	for(i in 1:peptempvec.l){
		if(length(unlist(grepRaw(pattern="250",object@pepraw[i],fixed=T,all=T)))>=1){
			test[[i]]<-unlist(grepRaw(pattern="250",object@pepraw[i],fixed=T,all=T))-4
			}
		if(length(unlist(grepRaw(pattern="242",object@pepraw[i],fixed=T,all=T)))>=1){
			test[[i]]<-unlist(grepRaw(pattern="242",object@pepraw[i],fixed=T,all=T))-4	
			}
		}
	unique(object@pepraw)
	as.character(test)[1847]
	as.data.frame(test)[2,]
	#peps<-unique(peps)
	cbind(object@data,test)
	print(length(peps))
	### put peptide positions into their slot
	object@modposition.peptide<-test
	object@pepvec<-peps
	

	return(object)
	
	}


