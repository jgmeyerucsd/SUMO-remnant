
compareSites=function(object1=mgpl.sig[[2]],object2=nature){
	modsummary1<-object1@modsummary
	modsummary2<-object2@modsummary
	#names(proteins1)<-proteins1
	proteins1<-names(modsummary1)
	proteins2<-names(modsummary2)
	proteins1<-substr(start=4,stop=9,proteins1)
	#proteins2<-substr(start=4,stop=9,proteins2)
	names(modsummary1)<-proteins1
	names(modsummary2)<-proteins2

	length(modsummary2)
	overlap.list<-list()
	nooverlap.list<-list()
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
			if(length(na.omit(modsummary2[[x]][match(modsummary1[[x]],modsummary2[[x]])]))==0){
				nooverlap.list[[x]]<-modsummary1[[x]]
				#names(x)
				}


			#print(tempmatch)
			count=count+1
			}
		}

	newcol<-rep("",times=length(object@data[,1]))
	#### part to add a column with our study and insert pluses
	overlap.l<-length(overlap.list)
	for(j in 1:length(overlap.list)){
		tempprotname<-names(overlap.list[j])
		#unlist(x)
		print(j)
		templines<-which(object2@data[,"Protein"]==tempprotname)
		tempsites<-unlist(overlap.list[j])
		for(i in 1:length(tempsites)){
			newcol[templines[which(object2@data[templines,2]==tempsites[i])]]<-"+"
			print(newcol[templines[which(object2@data[templines,2]==tempsites[i])]])
			}
		}


	for(x in overlap.list){
		tempprotname<-names(x)
		#unlist(x)
		print(x)
		templines<-which(object2@data[,"Protein"]==tempprotname)
		tempsites<-unlist(x)
		for(i in 1:length(tempsites)){
			newcol[templines[which(object2@data[templines,2]==tempsites[i])]]<-"+"
			print(newcol[templines[which(object2@data[templines,2]==tempsites[i])]])
			
			}
		
		}

	return(overlap.list)
	
	}


newcol[newcol=="+"]
length(unlist(overlap.list))
new<-cbind(object2@data,newcol)
write.table(cbind(object2@data,newcol),file="nature.ourMG132.tsv",quote=F,sep="\t",col.names=T,row.names=F)
### this column gives + for Vertegaal
object2@data[,4]
lapply(nooverlap.list, write, "ourunique.txt", append=TRUE, ncolumns=1000)
library(erer)
write.list(nooverlap.list, file="ourunique.txt", t.name = NULL, row.names = TRUE,quote=F)

