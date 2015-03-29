


comparemodsum=function(object1=mgpl.single.all.filt,
	object2=mgpl.brp.all.filt)
	{
	modsum1<-object1@modsummary
	modsum2<-object2@modsummary
	modsum1.len<-length(modsum1)
	newmodlist<-list()
	for(x in names(modsum2)){
		#print(x)
		which(names(modsum1)==x)
		modsum1[x]

		modsum2[x]
		is.na(match(modsum2[[x]],modsum1[[x]]))
		newmodlist[[x]]<-c(modsum2[[x]][is.na(match(modsum2[[x]],modsum1[[x]]))],modsum1[[x]])
		}	
	
	length(unlist(newmodlist))
