getwd()
setwd("C:/Users/JgMeyer/Documents/R/SILACquant")

files<-list.files()
files

uv0<-read.PepProph(input=files[124])
uv2<-read.PepProph(input=files[125])
uv8<-read.PepProph(input=files[126])
filterAll(object=uv0)->uv0.filt
filterAll(object=uv2)->uv2.filt
filterAll(object=uv8)->uv8.filt

length(unlist(uv0.filt@modsummary))
length(unlist(uv2.filt@modsummary))
length(unlist(uv8.filt@modsummary))


uv2.filt@data[1,]
mgpl.brp.all@data
object<-mgpl.brp.all
quantlines<-object@data[which(as.numeric(as.character(object@data$xpress))>0),]
quantlines[1:100,]
proteins<-substr(quantlines$protein,start=4,stop=9)
quantlist<-list()

for(i in 1:length(proteins)){
	quantlist[[proteins[i]]]<-quantlines$xpress[i]
	}
quantlines



mgpl.single.all<-read.PepProph(input=files[94])
mgpl.single.all@data$xpress

#### filter by 	(1) localization
#### 			(2) protein position
####			(3) quantification
####			(4) 

object1->mgpl.single.all.filt
mgpl.single.all@modsummary
uniquelines->mgpl.brp.uniquelines
uniquelines->mgpl.single.uniquelines

object@data
filterAll(
filterAll=function(object=mgpl.single.all,
	minscore=0.75,
	fasta="C:/Users/JgMeyer/Documents/R/pepsum/pepsum/110712_human.cc.fasta",
	name="test.tsv")
	{
	PTMscorelines<-as.character(object@data[,"ptm_peptide"])
	object1<-object
	#### object 1 only those with localization scores
	#which(PTMscorelines!="unavailable")
	object1@data<-object@data[PTMscorelines!="[unavailable]",]
	#nrow(object@data[PTMscorelines!="[unavailable]",])
	object1@pepvec<-object@pepvec[PTMscorelines!="[unavailable]"]
	#object1@pepvec[1]
	#object1@data[1,]

	### also filter based on 242\250 in peptide
	peptides<-as.character(object1@data$peptide)
	diglylines<-c(grep("242",peptides),grep("250",peptides))
	object1@data<-object1@data[diglylines,]
	#nrow(object1@data)
	object1@pepvec<-object1@pepvec[diglylines]
	
	object1@data$peptide
	object1@pepvec
	scorelines<-as.character(object1@data[,"ptm_peptide"])

	scores<-list()
	keepthese<-c()
	modposition<-list()
	line<-1
	for(i in 1:length(scorelines)){
		tempscore<-as.numeric(unlist(regmatches(scorelines[i],gregexpr("[[:digit:]]+\\.*[[:digit:]]*",scorelines[i]))))
		if(length(which(tempscore>=minscore)>0)>0){
			#print("isnumeric")
			keepthese<-c(keepthese,i)
			n=which(tempscore>=minscore)
			modposition[[line]]<-unlist(gregexpr("[[:digit:]]+\\.*[[:digit:]]*",scorelines[i]))[which(tempscore>=minscore)]-((n-1)*7+2)
			line=line+1
			}
		}
	### now have index of rows to keep
	### and positions as a list of positions
	### convert the list of peptide positions where length>1 into csv

	object1@data<-object1@data[keepthese,]
	object1@data[,2]
	object1@modposition.peptide<-modposition

	object1@pepvec<-object1@pepvec[keepthese]
	#### now get protein positions and summarize
	### read in fasta file
	require(seqinr)
	require(Biostrings)
	fastaobj<-read.fasta(fasta,seqtype="AA",as.string=TRUE)
	fastaacc<-substr(names(fastaobj),start=4,stop=9)
	datamat<-object1@data
	#proteins<-levels(datamat[,"protein"])
	proteins<-as.character(unique(datamat[,"protein"]))
	uniqproteins<-substr(start=4,stop=9,proteins)
	allprotacc<-substr(as.character(datamat[,"protein"]),start=4,stop=9)
	### replace object@pepvec with cleaned peptides
	peptempvec<-as.character(object1@pepvec)
	peptempvec.l<-length(peptempvec)
	### make empty vector for protein positions
	proteinposition<-rep(0,times=peptempvec.l)
	####  loop through each line in the datamat
	proteinposition<-list()
	#### fix to correctly assign n-terminal position
	for(i in 1:peptempvec.l){
		print(i)
		currentprot<-fastaobj[[which(fastaacc==substr(start=4,stop=9,datamat[i,"protein"]))]][1]
		currentpep<-peptempvec[i]
		#print(currentpro)
		#print(currentpep)
		matched<-matchPattern(currentpep,currentprot)
		proteinposition[[i]]<-matched@ranges[[1]][1]+(unlist(object1@modposition.peptide[i])-1)
		}
	object1@data<-cbind(protein.position=as.character(proteinposition),object1@data)
	#object@data[1,]
	object1@modposition.protein<-proteinposition
	#name="testMGneg.tsv"

	############################
	#### part to summarize and output unique positions

	proteinIDs<-unique(as.character(object1@data[,"protein"]))
	proteinIDs<-substr(proteinIDs,start=4,stop=9)
	prot.pos.list<-list()
	###  gives a list of proteins with their corresponding unique locations
	for(i in 1:length(proteinIDs)){
		#print(i)
		prot.pos.list[[proteinIDs[i]]]<-unique(unlist(object1@modposition.protein[which(substr(start=4,stop=9,object1@data$protein)==proteinIDs[i])]))
		#unique(unlist(object1@modposition.protein[which(substr(start=4,stop=9,object1@data$protein)=="Q96T23")]))
		
		#prot.pos.list[["Q96T23"]]
		#print(proteinIDs[i])
		#print(prot.pos.list[[proteinIDs[i]]])
		}
	length(unlist(prot.pos.list))
	prot.lines.list<-list()
	#### gives the lines in object@data that correspond to each protein
	for(i in 1:length(proteinIDs)){
		#print(i)
		prot.lines.list[[proteinIDs[i]]]<-which(substr(start=4,stop=9,as.character(object1@data$protein))==proteinIDs[i])
		}
	#range(unlist(prot.lines.list))
	position.index.list<-list()

	#### assign each unique position an index 
	### loop through the proteins
	index=1
	index2=-1
	for(i in 1:length(proteinIDs)){
		print(i)
		temp.positions<-prot.pos.list[[which(names(prot.pos.list)==proteinIDs[i])]]
		temp.prot.lines<-which(substr(object1@data$protein,start=4,stop=9)==proteinIDs[i])  ### gives row numbers of of mods in object@data
		protein.position.list<-object1@modposition.protein[temp.prot.lines]  ### gives the values of mod positions as vector
		### set lines that have multiple mods to 0
		for(j in 1:length(protein.position.list)){
			#print(length(protein.position.list[[j]]))
			if(length(protein.position.list[[j]])>1){
				protein.position.list[[j]]<-0
				}
			}
		unique.positions<-unique(unlist(protein.position.list))
		unique.positions.l<-length(unique.positions)

		### loop through the positions and assign those each an index number
		for(x in unique.positions){
			#print(x)
			### if x!=0, do give an index
			if(x!=0){
				temprowlen<-length(temp.prot.lines[which(protein.position.list==x)])
				for(j in 1:temprowlen){
					#print(j)
					position.index.list[[temp.prot.lines[which(protein.position.list==x)][j]]]<-index
					}
				index=index+1
				}
			###  if x==0 (peptide with two mods), assign negative index
			if(x==0){
				temprowlen<-length(temp.prot.lines[which(protein.position.list==x)])
				for(j in 1:temprowlen){
					#print(j)
					position.index.list[[temp.prot.lines[which(protein.position.list==x)][j]]]<-index2
					}
				index2=index2-1
				}
			}
		}


	object1@modsummary<-prot.pos.list
	object1@modindex<-position.index.list
	numuniq<-length(unique(unlist(object1@modindex)))
	uniq.indexes<-unique(unlist(object1@modindex))
	unique.position.singleline<-rep(0,times=numuniq)
	unique.position.singleweight<-rep(0,times=numuniq)
	#targetIDs@data$weighted.ratios
	i=1
	for(x in uniq.indexes){
		print(x)
		unique.position.singleline[i]<-which(object1@modindex==x)[1]
		i=i+1
		#unique.position.singleweight[i]<-log(targetIDs@data[which(object1@modindex==i)[1],"weighted.ratios"],base=2)
		}

	uniquelines<-object1@data[unique.position.singleline,]
	object1@data<-uniquelines



	
	
	if(writetsv==TRUE){
		write.table(uniquelines,file=name,quote=FALSE,sep="\t",row.names=F)
		}



	return(object1)


	}



	for(i in 1:length(proteinIDs)){
		print(i)
		temp.positions<-prot.pos.list[[which(names(prot.pos.list)=="Q96T23")]]
		temp.prot.lines<-which(substr(object@data$protein,start=4,stop=9)=="Q96T23")  ### gives row numbers of of mods in object@data
		protein.position.list<-object@modposition.protein[temp.prot.lines]  ### gives the values of mod positions as vector
		### set lines that have multiple mods to 0
		for(j in 1:length(protein.position.list)){
			#print(length(protein.position.list[[j]]))
			if(length(protein.position.list[[j]])>1){
				protein.position.list[[j]]<-0
				}
			}
		unique.positions<-unique(unlist(protein.position.list))
		unique.positions.l<-length(unique.positions)

