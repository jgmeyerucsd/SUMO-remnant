nature@data
	modsum1<-mgpl.brp.uniquelines
	modsum2<-mgpl.single.uniquelines
	

getMotif=function(object=mgpl.brp.all.filt,
	fasta="C:/MSGFplus/database/110712_human.cc.fasta",
	size=15)
	{
	datamat<-object@data
	require(seqinr)	
	require(Biostrings)
	fastaobj<-read.fasta(fasta,seqtype="AA",as.string=TRUE)
	fastaacc<-substr(names(fastaobj),start=4,stop=9)
	#proteins<-substr(as.character(datamat[,"protein"]),start=4,stop=9)
	proteins<-as.character(datamat[,"protein"])
	uniqproteins<-substr(start=4,stop=9,proteins)
	allprotacc<-substr(as.character(datamat[,"protein"]),start=4,stop=9)
	prot.l<-length(proteins)
	hydrophobic<-c("F", "M", "P", "C", "L", "I", "W", "A", "V", "Q", "Y")
	####

	### retrieve the previous 15 residues and post 15 residues

	### fix to add blanks for missing n-term or c-term values
	sites.l<-length(unlist(object@modsummary))
	window.vec<-rep(0,times=sites.l)
	window.list<-list()
	significant.vec<-rep(0,times=sites.l)
	dbprotnames<-names(fastaobj)
	targetprot.names<-names(object@modsummary)
	#targetprot.names<-substr(start=4,stop=9,targetprot.names)
	targetprot.l<-length(targetprot.names)
	
	j=1
	for(i in 1:targetprot.l){
		print(i)
		positions<-object@modsummary[[i]]+size
		tempprotname<-targetprot.names[i]
		seq<-fastaobj[[which(tempprotname==fastaacc)]][1]
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
	object@sequence<-window.vec
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


