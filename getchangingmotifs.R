


#### separate motifs by up/down/unchanged


### unchanged lins


nochangelines

#### separate those that are up/down

change.len<-nrow(changelines)

colnames(changelines)
log.heavy.over.light<-changelines$log.ratios*-1
	

	up<-changelines[log.heavy.over.light>0,]
	down<-changelines[log.heavy.over.light<0,]
	
changemotif=function(datamat=nochangelines,size=10){
	require(seqinr)	
	require(Biostrings)
	fastaobj<-read.fasta(fasta,seqtype="AA",as.string=TRUE)
	fastaacc<-substr(names(fastaobj),start=4,stop=9)
	#proteins<-substr(as.character(datamat[,"protein"]),start=4,stop=9)
	proteins<-as.character(datamat[,"protein"])
	uniqproteins<-unique(substr(start=4,stop=9,proteins))
	allprotacc<-substr(as.character(datamat[,"protein"]),start=4,stop=9)
	prot.l<-length(proteins)
	hydrophobic<-c("F", "M", "P", "C", "L", "I", "W", "A", "V", "Q", "Y")
	####

	### retrieve the previous 15 residues and post 15 residues

	### fix to add blanks for missing n-term or c-term values
	sites.l<-nrow(datamat)
	window.vec<-rep(0,times=sites.l)
	window.list<-list()
	significant.vec<-rep(0,times=sites.l)
	dbprotnames<-names(fastaobj)
	targetprot.names<-allprotacc
	#targetprot.names<-substr(start=4,stop=9,targetprot.names)
	targetprot.l<-length(targetprot.names)
	positions<-as.numeric(as.character(datamat[,1]))+size

	j=1
	for(i in 1:sites.l){
		print(i)
		tempprotname<-targetprot.names[i]
		position<-positions[i]
		seq<-fastaobj[[which(tempprotname==fastaacc)]][1]
		#### add "-" to each end of the sequence based on the max motif size
		seq<-paste(paste(rep("-",times=size),collapse=""),seq,paste(rep("-",times=size),collapse=""),sep="")
		seq.l<-nchar(seq)
		window.vec[j]<-substr(seq,start=position-size,stop=position+size)
		#if(nchar(window.vec[j])!=21){print(i)}
		j=j+1
		#print(x)
		}
	
	write(window.vec,file="nochange_motiflines.tsv")	


nrow(changelines[log.heavy.over.light>0,])

changelines$weighted.ratios


for(i in 1:change.len){
	

	}


down<-
up<-

