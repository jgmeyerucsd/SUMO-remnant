nature@data
	modsum1<-mgpl.brp.uniquelines
	modsum2<-mgpl.single.uniquelines
	getwd()
	fasta2="C:/Users/JgMeyer/Documents/R/SILACquant/ubi.full.fasta"
files<-

files<-list.files()


getMotif=function(object1=mgpl.brp.all.filt,
	fasta1="C:/MSGFplus/database/110712_human.cc.fasta",
	size=15)
	{
	datamat1<-object1@data
	datamat2<-object2@data
	modsum1<-object1@modsummary
	modsum2<-object2@modsummary

	sumopositions<-object1@data[,1]

	sumoproteins<-substr(as.character(object1@data$protein),start=4,stop=9)	
	ubiproteins<-as.character(mg132.ubi.lines[,"Protein.Id"])
	datamat2[1,]
	mg132.ubi.lines<-datamat2[which(datamat2[,"MG132_HPRPbeforeIP_spectral_counts"]>0),]
	mg132.ubi.lines[,1:10]
	ubiproteins<-as.character(mg132.ubi.lines[,"Protein.Id"])
	ubipositions<-mg132.ubi.lines[,5]
	#mylist<-strsplit(as.character(datamat2[,"Protein.Id"]),split="|",fixed=T)
	#ubiproteins<-sapply(mylist,function(x) x[2])
	#ubiproteins<-substr(ubiproteins,start=1,stop=6)
	#unlist(strsplit(as.character(datamat2[,"Protein.Id"]),split="|",fixed=T))[seq(from=2,to=5316,by=3)]

	write(unique(ubiproteins),file="ubiacc.tsv")
	for(x in unique(ubiproteins)){
		object2@modsummary[[x]]<-datamat2[which(ubiproteins==x),5]
		}
	



	require(seqinr)
	require(Biostrings)
	fastaobj1<-read.fasta(fasta1,seqtype="AA",as.string=TRUE)
	fastaobj2<-read.fasta(fasta2,seqtype="AA",as.string=TRUE)
	object1@modsummary
	object2@modsummary

	fastaacc1<-substr(names(fastaobj1), start=4,stop=9)
	#fastaacc2<-substr(names(fastaobj2), start=4, stop=9)
	mylist<-strsplit(names(fastaobj2),split="|",fixed=T)
	fastaacc2<-sapply(mylist,function(x) x[2])
	fastaacc2<-names(fastaobj2)

	sort(fastaacc2)
	####
	### retrieve the previous 15 residues and post 15 residues
	### fix to add blanks for missing n-term or c-term values
	sites.l<-length(unlist(object1@modsummary))
	window.vec<-rep(0,times=sites.l)
	window.list<-list()
	significant.vec<-rep(0,times=sites.l)
	targetprot.names<-names(object1@modsummary)
	#targetprot.names<-substr(start=4,stop=9,targetprot.names)
	targetprot.l<-length(targetprot.names)
	

	window.vec<-rep(0,times=sumo.l)
	window.list<-list()
	### loop through lines in sumopositions and ubipositions
	sumo.l<-length(sumopositions)
	sumopos2<-as.numeric(as.character(sumopositions))
	j=1
	for(i in 1:sumo.l){
		print(i)
		position<-sumopos2[i]+size
		tempprotname<-sumoproteins[i]
		seq<-fastaobj1[[which(tempprotname==fastaacc1)]][1]
		#### add "-" to each end of the sequence based on the max motif size
		seq<-paste(paste(rep("_",times=size),collapse=""),seq,paste(rep("_",times=size),collapse=""),sep="")
		seq.l<-nchar(seq)
		window.vec[j]<-substr(seq,start=position-size,stop=position+size)
		if(nchar(window.vec[j])!=((size*2)+1)){print(i)}
		j=j+1
		#print(x)
		}
	
	object1@sequence<-window.vec

	#### now ubi

	window.list<-list()
	### loop through lines in sumopositions and ubipositions
	ubi.l<-length(ubipositions)
	window.vec<-rep(0,times=ubi.l)
	j=1
	for(i in 1:ubi.l){
		print(i)
		position<-ubipositions[i]+size
		tempprotname<-ubiproteins[i]
		seq<-fastaobj2[[which(tempprotname==fastaacc2)]][1]
		#### add "-" to each end of the sequence based on the max motif size
		seq<-paste(paste(rep("_",times=size),collapse=""),seq,paste(rep("_",times=size),collapse=""),sep="")
		seq.l<-nchar(seq)
		window.vec[j]<-substr(seq,start=position-size,stop=position+size)
		if(nchar(window.vec[j])!=((size*2)+1)){print(i)}
		j=j+1
		#print(x)
		}
	
	object2@sequence<-window.vec
	

	#### repeat the same thing for object 2
	### retrieve the previous 15 residues and post 15 residues
	### fix to add blanks for missing n-term or c-term values
	sites.l<-length(unlist(object1@modsummary))
	window.vec<-rep(0,times=sites.l)
	window.list<-list()
	significant.vec<-rep(0,times=sites.l)
	dbprotnames<-names(fastaobj2)
	targetprot.names<-names(object1@modsummary)
	#targetprot.names<-substr(start=4,stop=9,targetprot.names)
	targetprot.l<-length(targetprot.names)
	
	j=1
	for(i in 1:targetprot.l){
		print(i)
		positions<-object1@modsummary[[i]]+size
		tempprotname<-targetprot.names[i]
		seq<-fastaobj1[[which(tempprotname==fastaacc1)]][1]
		#### add "-" to each end of the sequence based on the max motif size
		seq<-paste(paste(rep("_",times=size),collapse=""),seq,paste(rep("_",times=size),collapse=""),sep="")
		seq.l<-nchar(seq)
		for(x in positions){
			window.vec[j]<-substr(seq,start=x-size,stop=x+size)
			if(nchar(window.vec[j])!=((size*2)+1)){print(i)}
			j=j+1
			#print(x)
			}
		}
	object1@sequence<-window.vec

	match(object1@sequence,table=object2@sequence)
	object1@sequence[624]
	object2@sequence[104]
 	overlap.mg.table<-mg132.ubi.lines[na.omit(match(object1@sequence,object2@sequence)),]
	write.table(file="matched3.mg.ubi.tsv",cbind(object2@sequence[na.omit(match(object1@sequence,object2@sequence))],overlap.mg.table),quote=F,sep="\t",col.names=T,row.names=F)
	write.table(file="MG132sumo.sites.w.windows.tsv",cbind(window.vec,object1@data),quote=F,sep="\t",col.names=T,row.names=F)

	#### repeat the same thing for object 2

overlap.mg.table[,1:8]

	mods.vec<-unlist(object@modsummary)
	names(mods.vec)<-substr(names(mods.vec),start=1,stop=6)
	i=1
	position.in.nature.vec<-rep(0,times=length(mods.vec))
	for(x in object@sequence){
		print(which(nature@data[,3]==x))
		if(length(which(nature@data[,3]==x))>0){
			position.in.nature.vec[i]<-which(nature@data[,3]==x)
			}
		i=i+1
		}
	newcol<-rep("",times=nrow(nature@data))
	newcol[position.in.nature.vec]<-"+"
	write.table(file="matched.with.nature.tsv",cbind(nature@data,meyer=newcol),quote=F,sep="\t",col.names=T,row.names=F)

	
	object@data<-cbind(object@data,sequencewindow=object@sequence)
	filtered<-unique(window.vec)
	filtered.l<-length(window.vec)
	de<-c("D","E")
	motifvec<-rep(0,times=filtered.l)
	change.vec<-rep(0,times=filtered.l)
	for(i in 1:filtered.l){
	### assign 1 for known normal , 2 for inverted, 0 for other
		nextAA<-substr(filtered[i],start=12,stop=12)
		prevAA<-substr(filtered[i],start=10,stop=10)
		next.test<-which(hydrophobic==nextAA)
		prev.test<-which(hydrophobic==prevAA)
		filter
		if(length(next.test)>0){
		### check number 4 for D/E
			if(length(which(substr(filtered[i],start=9,stop=9)==de)>=1)>0){
				print("inverted")
				motifvec[i]<-2
				}
			}

		if(length(prev.test)>0){
			### check number 8 for D/E
			if(length(which(substr(filtered[i],start=13,stop=13)==de)>=1)>0){
				print("normal")
				motifvec[i]<-1
				}
			}

		}
		}
	print("number of normal")
	print(length(which(motifvec==1)))

	print("number of inverted")
	print(length(which(motifvec==2)))

	print("number of new")
	print(length(which(motifvec==0)))


	na.omit(unique(window.vec[which(motifvec==1)]))
	na.omit(unique(window.vec[which(motifvec==2)]))
	na.omit(unique(window.vec[which(motifvec==0)]))

	newmotif<-na.omit(unique(filtered[which(motifvec==0)]))
	write(newmotif,file="newmotiflines.tsv")	
	generalmotif<-na.omit(unique(filtered[which(motifvec==1)]))
	write(generalmotif,file="generalmotiflines.tsv")	
	invertedmotif<-na.omit(unique(filtered[which(motifvec==2)]))
	write(invertedmotif,file="invertedmotiflines.tsv")	
	allmotif<-na.omit(unique(filtered))
	write(allmotif,file="allmotiflines.tsv")	

	### now report those ==0 and print for icelogo
	
	
	which(hydrophobic)

	###

	###
	#### filter to keep unique sequence windows

	###  
	
		}


