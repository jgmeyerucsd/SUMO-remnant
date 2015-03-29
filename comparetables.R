


compare.table=function(object1=mgpl.sig[[2]],
	table=mgpl.sig[[3]],
	natureobject=nature,
	ub=ub.psp,
	sumo=sumo.psp)
	{
	modsummary1<-object1@modsummary
	modsummary2<-natureobject@modsummary
	sumomods<-sumo.psp@modsummary
	ubmods<-ub.psp@modsummary
	methylmods<-methyl.psp@modsummary
	acetylmods<-acetyl.psp@modsummary


	#names(proteins1)<-proteins1
	proteins1<-names(modsummary1)
	proteins2<-names(modsummary2)
	proteins3<-substr(start=4,stop=9,table[,"protein"])
	sumo.proteins<-names(sumomods)
	ub.proteins<-names(ubmods)
	acetyl.proteins<-names(acetylmods)
	methyl.proteins<-names(methylmods)

	proteins1<-substr(start=4,stop=9,proteins1)

	#proteins2<-substr(start=4,stop=9,proteins2)
	names(modsummary1)<-proteins1
	names(modsummary2)<-proteins2
	
	overlap.list<-list()
	nooverlap.list<-list()
	sumo.overlap<-list()
	ub.overlap<-list()
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
		templines<-which(natureobject@data[,"Protein"]==tempprotname)
		tempsites<-unlist(overlap.list[j])
		for(i in 1:length(tempsites)){
			newcol[templines[which(natureobject@data[templines,2]==tempsites[i])]]<-"+"
			print(newcol[templines[which(natureobject@data[templines,2]==tempsites[i])]])
			}
		}

	
	#### part to add a column into our data table for presence in Nature
	tableproteins<-substr(start=4,stop=9,table[,"protein"])
	naturecol<-rep("",times=length(table[,1]))
	overlap.l<-length(overlap.list)
	table[1,]
	for(j in 1:length(overlap.list)){
		tempprotname<-names(overlap.list[j])
		#unlist(x)
		print(j)
		templines<-which(tableproteins==tempprotname)
		tempsites<-unlist(overlap.list[j])
		for(i in 1:length(tempsites)){
			naturecol[templines[which(table[templines,1]==tempsites[i])]]<-"+"
			print(naturecol[templines[which(natureobject@data[templines,2]==tempsites[i])]])
			}
		}
	### find what positions from the table overlap with the sumo psp sites and return list
	count=0
	sumocol<-rep("",times=length(table[,1]))

	for(x in tableproteins){
		print(x)
		pos.in.sumo.modsummary<-which(x==sumo.proteins)
		pos.in.tableproteins<-which(x==tableproteins)
		print(temppos)
		tableproteins.sites<-table[pos.in.tableproteins,1]
		### if the protein accession appears in the SUMO proteins
		### check if the table
		if(length(pos.in.sumo.modsummary)>0){
			#modsummary1[[x]]
			#length(modsummary1[[x]])
			#length(modsummary2[[x]])
			if(length(na.omit(pos.in.tableproteins[match(tableproteins.sites,sumomods[[x]])]))>0){

				### this give the positions that should be + in the 
				sumocol[na.omit(pos.in.tableproteins[match(tableproteins.sites,sumomods[[x]])])]<-"+"
				}
			#tempmatch<-match(modsummary1[[x]],modsummary2[[x]])
			#x<-"Q13547"
			#print(tempmatch)
			count=count+1
			}
		}
	count=0
	ubcol<-rep("",times=length(table[,1]))

	for(x in tableproteins){
		print(x)
		pos.in.ub.modsummary<-which(x==ub.proteins)
		pos.in.tableproteins<-which(x==tableproteins)
		print(temppos)
		tableproteins.sites<-table[pos.in.tableproteins,1]
		### if the protein accession appears in the SUMO proteins
		### check if the table
		if(length(pos.in.ub.modsummary)>0){
			#modsummary1[[x]]
			#length(modsummary1[[x]])
			#length(modsummary2[[x]])
			if(length(na.omit(pos.in.tableproteins[match(tableproteins.sites,ubmods[[x]])]))>0){

				### this give the positions that should be + in the 
				ubcol[na.omit(pos.in.tableproteins[match(tableproteins.sites,ubmods[[x]])])]<-"+"
				}
			#tempmatch<-match(modsummary1[[x]],modsummary2[[x]])
			#x<-"Q13547"
			#print(tempmatch)
			count=count+1
			}
		}

	count=0
	methyl.col<-rep("",times=length(table[,1]))

	for(x in tableproteins){
		print(x)
		pos.in.methyl.modsummary<-which(x==methyl.proteins)
		pos.in.tableproteins<-which(x==tableproteins)
		print(temppos)
		tableproteins.sites<-table[pos.in.tableproteins,1]
		### if the protein accession appears in the SUMO proteins
		### check if the table
		if(length(pos.in.methyl.modsummary)>0){
			#modsummary1[[x]]
			#length(modsummary1[[x]])
			#length(modsummary2[[x]])
			if(length(na.omit(pos.in.tableproteins[match(tableproteins.sites,methylmods[[x]])]))>0){

				### this give the positions that should be + in the 
				methyl.col[na.omit(pos.in.tableproteins[match(tableproteins.sites,methylmods[[x]])])]<-"+"
				}
			#tempmatch<-match(modsummary1[[x]],modsummary2[[x]])
			#x<-"Q13547"
			#print(tempmatch)
			count=count+1
			}
		}

	count=0
	acetyl.col<-rep("",times=length(table[,1]))

	for(x in tableproteins){
		print(x)
		pos.in.acetyl.modsummary<-which(x==acetyl.proteins)
		pos.in.tableproteins<-which(x==tableproteins)
		print(temppos)
		tableproteins.sites<-table[pos.in.tableproteins,1]
		### if the protein accession appears in the SUMO proteins
		### check if the table
		if(length(pos.in.acetyl.modsummary)>0){
			if(length(na.omit(pos.in.tableproteins[match(tableproteins.sites,acetylmods[[x]])]))>0){

				### this give the positions that should be + in the 
				acetyl.col[na.omit(pos.in.tableproteins[match(tableproteins.sites,acetylmods[[x]])])]<-"+"
				}
			count=count+1
			}
		}


	output.table<-cbind(table,in.nature=naturecol,in.PSP.sumo=sumocol,in.PSP.ub=ubcol,in.PSP.methyl=methyl.col,in.PSP.acetyl=acetyl.col)


	output.table[1,]




	#### part to add a column into our data table for presence in sumo PSP db
	tableproteins<-substr(start=4,stop=9,table[,"protein"])
	sumocol<-rep("",times=length(table[,1]))
	overlap.l<-length(sumo.overlap)
	for(j in 1:length(sumo.overlap)){
		tempprotname<-names(overlap.list[j])
		#unlist(x)
		print(j)
		templines<-which(tableproteins==tempprotname)
		tempsites<-unlist(overlap.list[j])
		for(i in 1:length(tempsites)){
			naturecol[templines[which(table[templines,1]==tempsites[i])]]<-"+"
			print(naturecol[templines[which(natureobject@data[templines,2]==tempsites[i])]])
			}
		}


	return(overlap.list)
	
	}


newcol[newcol=="+"]
length(unlist(overlap.list))
new<-cbind(object2@data,newcol)
write.table(output.table,file="MG132plus.unchanged.tsv",quote=F,sep="\t",col.names=T,row.names=F)
### this column gives + for Vertegaal
object2@data[,4]
lapply(nooverlap.list, write, "ourunique.txt", append=TRUE, ncolumns=1000)
library(erer)
write.list(nooverlap.list, file="ourunique.txt", t.name = NULL, row.names = TRUE,quote=F)

	}