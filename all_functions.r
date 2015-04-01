
#### this file contains all needed functions to go through SILAC quantification analysis
#### select all and run all within R


### create class for peptide identification data

setClass("pepsum",representation(summary="matrix",scans="ANY",specdir="character",
	specfiles="list",data="ANY", totalheat="matrix",residues="character",
	sequence="ANY",score="character",filetype="character",scanvec="ANY",
	modposition.protein="ANY", chargevec="ANY",proteinmodlist="ANY",index="ANY", 
	fraction="character",pepvec="ANY",provec="ANY",filenames="ANY",p1.index="ANY",
	modposition.peptide="ANY",countlist="list",pepraw="ANY",ionCoverage="list",
	ionCovSum="list",modindex="list",modsummary="list"))

	
### function to read peptide prophet results in tab-delim format
	
read.PepProph=function(input=files[2],type="PeptideProphet"){
	object<-new("pepsum")
	object@filetype="PeptideProphet"

	object@data<-read.delim(input,header=T)

	### make scan and charge vectors
	
	scanlist<-strsplit(as.character(object@data$spectrum),split=".",fixed=T)
	scanlen<-length(scanlist)
	scanvec<-as.numeric(unlist(scanlist)[seq(2,length(scanlist)*4,by=4)])
	chargevec<-as.numeric(unlist(scanlist)[seq(4,length(scanlist)*4,by=4)])

	### make peptide vector 
	rawpepvec<-as.character(object@data$peptide)
	
	### clean sequences into format for matchions
	#rawpepvec<-substr(rawpepvec,start=3, stop=nchar(rawpepvec)-2)
	
	peptempvec<-substr(rawpepvec,start=3, stop=nchar(rawpepvec)-2)
	peptempvec.l<-length(rawpepvec)
	peps<-rep(0,times=peptempvec.l)
	if(type=="PeptideProphet"){
		for(i in 1:peptempvec.l){
			###match the string starting with [/[- followed by n integers and "]"
			### and replace with nothing
			peps[i]<-gsub("(\\[|\\[-)([0-9]+)(.)([0-9]+)(]+)", replacement="", peptempvec[i]) 
			}
		for(i in 1:peptempvec.l){
			###remove "n" from peptides with n-term acetyl
			if(unlist(strsplit(peps[i],split=""))[1]=="n"){
				peps[i]<-substr(peps[i],start=2,stop=nchar(peps[i]))
				}
			}
		}

	#### if format is inspect-style
	if(type=="inspect"){
		### start by replacing those with nterminal pyroglu
		rawpepvec<-gsub(rawpepvec,pattern="(n)(\\[)(-16.02)(\\])(Q)",replacement="Q[-17.027]")
		#### pyro glu done, move to the next mod
		rawpepvec<-gsub(rawpepvec,pattern="(n)(\\[)(43.02)(\\])",replacement="[42.011]")
		##### metox next
		rawpepvec<-gsub(rawpepvec,pattern="(\\[)(147.04)(\\])",replacement="[15.995]")
		#### finally replace NEM-modified Cysteine
		rawpepvec<-gsub(rawpepvec,pattern="(\\[)(228.06)(\\])",replacement="[125.048]")
		}
	if(type=="specnets"){
		### start by replacing those with nterminal pyroglu
		rawpepvec<-gsub(rawpepvec,pattern="(n)(\\[)(-16.02)(\\])(Q)",replacement="(Q,-17.027)")
		#### pyro glu done, move to the next mod
		nacetylindex<-grep(rawpepvec,pattern="(n)(\\[)(43.02)(\\])")
		nacetylstartres<-substr(rawpepvec[nacetylindex],start=9, stop=9)
		nacetyl_len<-length(rawpepvec[nacetylindex])
		for(i in 1:nacetyl_len){
			rawpepvec[nacetylindex][i]<-paste("(",nacetylstartres[i],",+42.011)",rawpepvec[nacetylindex][i],sep="",collapse="")
			}
		rawpepvec<-gsub(rawpepvec,pattern="(n)(\\[)(43.02)(\\])([A-Z]){1}",replacement="")
		##### metox next
		rawpepvec<-gsub(rawpepvec,pattern="(M)(\\[)(147.04)(\\])",replacement="(M,+15.995)")
		#### finally replace NEM-modified Cysteine
		rawpepvec<-gsub(rawpepvec,pattern="(C)(\\[)(228.06)(\\])",replacement="(C,+125.048)")
		}
	if(type=="R"){
		### start by replacing those with nterminal pyroglu
		rawpepvec<-gsub(rawpepvec,pattern="(n)(\\[)(-16.02)(\\])(Q)",replacement="-17.027Q")
		#### pyro glu done, move to the next mod
		rawpepvec<-gsub(rawpepvec,pattern="(n)(\\[)(43.02)(\\])",replacement="+42.011")
		##### metox next
		rawpepvec<-gsub(rawpepvec,pattern="(\\[)(147.04)(\\])",replacement="+15.995")
		#### finally replace NEM-modified Cysteine
		rawpepvec<-gsub(rawpepvec,pattern="(\\[)(228.06)(\\])",replacement="+125.05")
		}	
	object@pepvec<-peps
	object@scanvec<-scanvec
	object@chargevec<-chargevec
	return(object)
	}


#### function to filter identifications based on a minimum score

#### usage:
#### newobject<-removeLowLocalizations(object=[your pepsum object], minscore = [your choice of localization cutoff score])
removeLowLocalization=function(object=mgpl.all,minscore=0.75){
	PTMscorelines<-as.character(object@data[,"ptm_peptide"])
	scores<-list()
	keepthese<-c()
	modposition<-list()
	line<-1
	for(i in 1:length(PTMscorelines)){
		tempscore<-as.numeric(unlist(regmatches(PTMscorelines[i],gregexpr("[[:digit:]]+\\.*[[:digit:]]*",PTMscorelines[i]))))
		if(length(which(tempscore>=minscore)>0)>0){
			#print("isnumeric")
			keepthese<-c(keepthese,i)
			n=which(tempscore>=minscore)
			modposition[[line]]<-unlist(gregexpr("[[:digit:]]+\\.*[[:digit:]]*",PTMscorelines[i]))[which(tempscore>=minscore)]-((n-1)*7+2)
			line=line+1
			}
		}
	### now have index of rows to keep
	### and positions as a list of positions
	### convert the list of peptide positions where length>1 into csv

	### this works but the double mods are form of c(2,3)
	#test<-cbind(as.character(modposition),object@data[keepthese,])
	object@modposition.peptide<-modposition
	object@data<-object@data[keepthese,]
	object@pepvec<-object@pepvec[keepthese]
	return(object)
	}

#### function to determine the mod postion within a protein using the fasta DB
####  usage:
### 	#### determine to position of each site ID in their protein

### mgpl.s.pos<-proteinPositions(object=mgpl.s,
###		fasta="C:/MSGFplus/database/110712_human.cc.fasta",
###		writetsv=FALSE,
###		name="mgplus.localized.tsv")
		
		

proteinPositions=function(fasta="C:/Users/JgMeyer/Documents/R/pepsum/pepsum/110712_human.cc.fasta",
	object=mgpl.single.pos,
	writetsv=FALSE,
	name="MG132.all.filtered.tsv"){

	### read in fasta file
	require(seqinr)
	require(Biostrings)
	fastaobj<-read.fasta(fasta,seqtype="AA",as.string=TRUE)
	fastaacc<-substr(names(fastaobj),start=4,stop=9)
	datamat<-object@data
	#proteins<-levels(datamat[,"protein"])
	proteins<-as.character(unique(datamat[,"protein"]))
	uniqproteins<-substr(start=4,stop=9,proteins)
	allprotacc<-substr(as.character(datamat[,"protein"]),start=4,stop=9)
	### replace object@pepvec with cleaned peptides
	peptempvec<-as.character(object@pepvec)
	peptempvec.l<-length(peptempvec)
	### make empty vector for protein positions
	proteinposition<-rep(0,times=peptempvec.l)
	####  loop through each line in the datamat
	proteinposition<-list()
	#### fix to correctly assign n-terminal position
	for(i in 1:peptempvec.l){
		print(i)
		#print(peptempvec[i])
		#print(datamat[i,"protein"])		
		#currentpro<-substr(start=4,stop=9,datamat[i,"protein"])
		currentprot<-fastaobj[[which(fastaacc==substr(start=4,stop=9,datamat[i,"protein"]))]][1]
		currentpep<-peptempvec[i]
		#print(currentpro)
		#print(currentpep)
		matched<-matchPattern(currentpep,currentprot)
		proteinposition[[i]]<-matched@ranges[[1]][1]+(unlist(object@modposition.peptide[i])-1)
		}
	object@data<-cbind(protein.position=as.character(proteinposition),object@data)
	#object@data[1,]
	object@modposition.protein<-proteinposition
	#name="testMGneg.tsv"
	if(writetsv==TRUE){
		write.table(object@data,file=name,quote=FALSE,sep="\t",row.names=F)
		}
	return(object)
	}

#### 	function to summarize and combine the unique protein site IDs
##  	usage:  
### 	mgpl.s.ave<-summarizeProtPositions(object=mgpl.s.pos)


#mgpl.ave<-summarizeProtPositions()
#mgpl.ave@data[1,]
#mgpl.ave@modindex

### works with peptides containing multiple sites, sets their weighted ratio =0


summarizeProtPositions=function(object=mgpl.pos){
	proteinIDs<-unique(as.character(object@data[,"protein"]))
	prot.pos.list<-list()
	###  gives a list of proteins with their corresponding unique locations
	for(i in 1:length(proteinIDs)){
		#print(i)
		prot.pos.list[[proteinIDs[i]]]<-unique(unlist(object@modposition.protein[which(object@data$protein==proteinIDs[i])]))
		#print(proteinIDs[i])
		#print(prot.pos.list[[proteinIDs[i]]])
		}
	prot.lines.list<-list()
	#### gives the lines in object@data that correspond to each protein
	for(i in 1:length(proteinIDs)){
		#print(i)
		prot.lines.list[[proteinIDs[i]]]<-which(object@data$protein==proteinIDs[i])
		}
	#range(unlist(prot.lines.list))
	position.index.list<-list()

	#### assign each unique position an index 
	### loop through the proteins
	index=1
	for(i in 1:length(proteinIDs)){
		#print(i)
		temp.positions<-prot.pos.list[[which(names(prot.pos.list)==proteinIDs[i])]]
		
		temp.prot.lines<-which(object@data$protein==proteinIDs[i])  ### gives row numbers of of mods in object@data
		
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
			###  if x==0 (peptide with two mods), assign index=0
			if(x==0){
				temprowlen<-length(temp.prot.lines[which(protein.position.list==x)])
				for(j in 1:temprowlen){
					#print(j)
					position.index.list[[temp.prot.lines[which(protein.position.list==x)][j]]]<-0
					}
				}
			}
		}
		### which(position.index.list==111)
		#### now loop through those values and take weighted averages of the ones with more than 1 line
		unique.indexes.len<-length(unique(unlist(position.index.list)))
		unique.indexes<-unique(unlist(position.index.list))
		#unique.indexes
		weighted.ratios<-list()
		linesum=0
		for(j in 1:unique.indexes.len){
			#print(j)

			templines<-which(position.index.list==unique.indexes[[j]])
			#print(templines)
			object@data[templines,]
			if(unique.indexes[[j]]>=1){
				### divide by temparea
				temparea=sum(object@data[templines,"light_area"])+sum(object@data[templines,"heavy_area"])
				templines.l<-length(templines)
				tempsum<-0
				linesum=linesum+templines.l
				weights<-rep(0,times=templines.l)
				for(i in 1:templines.l){
					weights[i]<-(object@data[templines[i],"light_area"]+object@data[templines[i],"heavy_area"])/temparea
					}
				for(i in 1:templines.l){
					tempsum=tempsum+object@data[templines[i],"xpress"]*weights[i]
					#print(object@data[templines[i],"light_area"]+object@data[templines[i],"heavy_area"])
					#print(tempsum)
					}
				for(i in 1:templines.l){
					weighted.ratios[[templines[i]]]<-tempsum
					}
				}	
			if(unique.indexes[[j]]==0){
				templines.l<-length(templines)
				for(i in 1:templines.l){
					weighted.ratios[[templines[i]]]<-0
					}
				}
			}
	object@modsummary<-prot.pos.list
	object@modindex<-position.index.list
	object@data<-cbind(object@data,weighted.ratios=unlist(weighted.ratios))
	print("unique SUMO modified protein IDs")
	print(length(proteinIDs))
	print("unique SUMO modification sites from single-site peptides")
	print(length(index))
	length(prot.pos.list)
	### write something to make these text and put them in one column
	#paste(unlist(prot.pos.list[1]))
	return(object)	
	}
	
	
	
#### compute the weighted average of sites identified by with multiple sequences,
#### then determine which sites are outside at least one standard deviation
#### will produce a plot of quantification values as given in figure 3a
###		usage:
###		mgpl.s.sig<-getSignificant(allIDs=mgpl.s.all,targetIDs=mgpl.s.ave)


getSignificant=function(allIDs=mgpl.all,targetIDs=mgpl.ave,
		name="MG132plus.filtered",
		writetsv=F,
		stdev=0.9404181){
	logRatiosAll<-log(allIDs@data$xpress,base=2)
	allIDs@data<-cbind(allIDs@data,logRatiosAll)
	meanAll<-mean(logRatiosAll)
	medianAll<-median(logRatiosAll)
	logRatiosAllnorm<-logRatiosAll-medianAll
	norm.median.all<-median(logRatiosAllnorm)
	### get distribution of weighted averages for unique sites
	#targetIDs@modindex
	numuniq<-length(unique(unlist(targetIDs@modindex)))-1
	unique.position.singleline<-rep(0,times=numuniq)
	unique.position.singleweight<-rep(0,times=numuniq)
	#targetIDs@data$weighted.ratios
	for(i in 1:numuniq){
		unique.position.singleline[i]<-which(targetIDs@modindex==i)[1]
		unique.position.singleweight[i]<-log(targetIDs@data[which(targetIDs@modindex==i)[1],"weighted.ratios"],base=2)
		}
	median(unique.position.singleweight)
	range(unique.position.singleweight)

	norm.uniq.pos.singleweight<-unique.position.singleweight-medianAll
	#2*sd(norm.uniq.pos.singleweight)
	if(length(stdev)==1){
		upperlim<-norm.median.all+2*stdev
		lowerlim<-norm.median.all-2*stdev
		onesdup<-norm.median.all+1*stdev
		onesddown<-norm.median.all-1*stdev
		}
	if(length(stdev)==0){
		upperlim<-norm.median.all+2*sd(logRatiosAllnorm)
		lowerlim<-norm.median.all-2*sd(logRatiosAllnorm)
		onesdup<-norm.median.all+1*sd(logRatiosAllnorm)
		onesddown<-norm.median.all-1*sd(logRatiosAllnorm)
		}
	#?pairlist
	#log(10)
	
	#unique.position.singleline[1]
	#norm.uniq.pos.singleweight[1]
	names(norm.uniq.pos.singleweight)<-unique.position.singleline
	sort.norm.uniq.pos.singleweight<-sort(norm.uniq.pos.singleweight*-1)
	#median(sort.norm.uniq.pos.singleweight)
	
	plot(x=1:numuniq,y=sort.norm.uniq.pos.singleweight)
	abline(h=upperlim,col="red")
	abline(h=lowerlim,col="red")
	###
	print("one standard deviation is")
	print(onesdup)
	abline(h=onesdup,col="red")
	abline(h=onesddown,col="red")
	#### now flag those IDs that are outside 2 std. dev.
	idlen<-length(targetIDs@data$xpress)

	oversigma<-rep(0,times=idlen)
	count=0
	#### loop through and test which are outside 2 std. dev
	for(i in 1:numuniq){
		if(norm.uniq.pos.singleweight[i]>=onesdup | norm.uniq.pos.singleweight[i]<=onesddown){
			#print(norm.uniq.pos.singleweight[i])
			count=count+1
			#print(count)
			oversigma[unique.position.singleline[i]]<-1
			}
		}
	over2sigma<-rep(0,times=idlen)
	count=0
	#### loop through and test which are outside 2 std. dev
	for(i in 1:numuniq){
		if(norm.uniq.pos.singleweight[i]>=upperlim | norm.uniq.pos.singleweight[i]<=lowerlim){
			#print(norm.uniq.pos.singleweight[i])
			count=count+1
			#print(count)
			oversigma[unique.position.singleline[i]]<-2
			}
		}
	# for those missing their weighted average value, put in xpress value
	for(i in 1:length(targetIDs@data[,"weighted.ratios"])){
		if(targetIDs@data[i,"weighted.ratios"]==0){
			targetIDs@data[i,"weighted.ratios"]<-targetIDs@data[i,"xpress"]
			}
		}

	targetIDs@data<-cbind(targetIDs@data,log.ratios=log(targetIDs@data[,"weighted.ratios"],base=2)-medianAll)
	targetIDs@data<-cbind(targetIDs@data,oversigma)
	#targetIDs@data<-cbind(targetIDs@data,over1sigma)

	### part to output only unique significant changers
	#unique.position.singleline
	#targetIDs@data[,"log.ratios"]
	uniquelines<-targetIDs@data[unique.position.singleline,]
	changelines<-uniquelines[uniquelines[,"oversigma"]>=1,]
	nochangelines<-uniquelines[uniquelines[,"oversigma"]==0,]


	### part to output only unique non changers


	### now have new object with binary whether outside 2*sigma
	print(i)
	#### maybe write new table?
	print("unique sites")
	print(numuniq)
	print("over 2 stdev")
	print(length(oversigma[oversigma==2]))
	print("over 1 stdev")
	print(length(oversigma[oversigma==1]))
	if(writetsv){
		write.table(file=paste(name,".all.tsv",sep=""),targetIDs@data,quote=F,sep="\t",col.names=T,row.names=F)
		write.table(file=paste(name,".changing.tsv",sep=""),changelines,quote=F,sep="\t",col.names=T,row.names=F)
		write.table(file=paste(name,".nochange.tsv",sep=""),nochangelines,quote=F,sep="\t",col.names=T,row.names=F)
		}
	return(list(allIDs,targetIDs,nochangelines,changelines))
	}


