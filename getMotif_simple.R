#fasta="C:/MSGFplus/database/110712_human.cc.fasta"

#object<-mgpl.ave

getMotif=function(object=mgpl.brp.all,
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
		seq<-paste(paste(rep("-",times=size),collapse=""),seq,paste(rep("-",times=size),collapse=""),sep="")
		seq.l<-nchar(seq)
		for(x in positions){
			window.vec[j]<-substr(seq,start=x-size,stop=x+size)
			if(nchar(window.vec[j])!=21){print(i)}
			j=j+1
			#print(x)
			}
		}

	filtered<-unique(window.vec)
	filtered.l<-length(window.vec)
	de<-c("D","E")
	motifvec<-rep(0,times=filtered.l)
	change.vec<-rep(0,times=filtered.l)
	for(i in 1:filtered.l){
	### assign 1 for known normal , 2 for inverted, 0 for other
		nextAA<-substr(filtered[i],start=12,stop=12)
		prevAA<-substr(filtered[i],start=10,stop=10)
		### check number 4 for D/E
		if(length(which(substr(filtered[i],start=size+3,stop=size+3)==de)>=1)>0){
			print("forward")
			motifvec[i]<-2
			}
		

		### check number 8 for D/E
		if(length(which(substr(filtered[i],start=size-1,stop=size-1)==de)>=1)>0){
			print("inverted")
			motifvec[i]<-1
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
	write(newmotif,file="20150310neithermotiflines.tsv")	
	generalmotif<-na.omit(unique(filtered[which(motifvec==2)]))
	write(generalmotif,file="20150310_forwardmotiflines.tsv")	
	invertedmotif<-na.omit(unique(filtered[which(motifvec==1)]))
	write(invertedmotif,file="20150310_invertedmotiflines.tsv")	
	allmotif<-na.omit(unique(filtered))
	write(allmotif,file="20150310_allmotiflines.tsv")	

	### now report those ==0 and print for icelogo
	
	
	which(hydrophobic)

	###

	###
	#### filter to keep unique sequence windows

	###  
	
		}

