

compareSites=function(object1=mgpl.sum,object2=mgneg.sum){
	modsummary1<-object1@modsummary
	modsummary2<-object2@modsummary
	#names(proteins1)<-proteins1
	proteins1<-names(modsummary1)
	proteins2<-names(modsummary2)
	proteins1<-substr(start=4,stop=9,proteins1)
	proteins2<-substr(start=4,stop=9,proteins2)
	names(modsummary1)<-proteins1
	names(modsummary2)<-proteins2

	overlap.list<-list()
	count=0
	for(x in proteins1){
		print(x)
		temppos<-which(x==proteins2)
		print(temppos)
		if(length(temppos)>0){
			modsummary1[[x]]
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

lapply(overlap.list, write, "test2.txt", append=TRUE)
?con
writeLines(unlist(lapply(overlap.list, paste, collapse=" ")))


library(marray)
data(swirl)
test <- list(A = 1:10, B= maM(swirl)[1:10,], C=list(x=1:10, y=1:4),
             D = summary(maA(swirl[,1])))
write.list(overlap.list, filename="test4.txt")

length(unlist(modsummary1))
length(unlist(modsummary2))
length(unlist(overlap.list))

