
compareSites=function(object1=mgpl.sum,object2=nature){
	modsummary1<-object1@modsummary
	modsummary2<-object2@modsummary
	#names(proteins1)<-proteins1
	proteins1<-names(modsummary1)
	proteins2<-names(modsummary2)
	proteins1<-substr(start=4,stop=9,proteins1)
	#proteins2<-substr(start=4,stop=9,proteins2)
	names(modsummary1)<-proteins1
	names(modsummary2)<-proteins2

	overlap.list<-list()
	count=0
	for(x in proteins1){
		print(x)
		temppos<-which(x==proteins2)
		print(temppos)
		if(length(temppos)>0){
			#modsummary1[[x]]
			#length(modsummary1[[x]])
			#length(modsummary2[[x]])
			if(length(na.omit(modsummary2[[x]][match(modsummary1[[x]],modsummary2[[x]])]))>0){
				overlap.list[[x]]<-na.omit(modsummary2[[x]][match(modsummary1[[x]],modsummary2[[x]])])
				}
			#tempmatch<-match(modsummary1[[x]],modsummary2[[x]])
			#x<-"Q13547"
			
			#print(tempmatch)
			count=count+1
			}
		}
	return(overlap.list)
	
	}
