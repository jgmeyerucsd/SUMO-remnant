
#mgpl.ave<-summarizeProtPositions()
#mgpl.ave@data[1,]
#mgpl.ave@modindex

### works with peptides containing multiple sites, sets their weighted ratio =0


summarizeProtPositions=function(object=mgpl.pos){
	proteinIDs<-unique(as.character(object@data[,"protein"]))
	proteinIDs<-substr(proteinIDs,start=4,stop=9)
	prot.pos.list<-list()
	###  gives a list of proteins with their corresponding unique locations
	for(i in 1:length(proteinIDs)){
		#print(i)
		prot.pos.list[[proteinIDs[i]]]<-unique(unlist(object@modposition.protein[which(substr(start=4,stop=9,object@data$protein)==proteinIDs[i])]))
		#unique(unlist(object@modposition.protein[which(substr(start=4,stop=9,object@data$protein)=="Q96T23")]))
		
		#prot.pos.list[["Q96T23"]]
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
	index2=-1
	for(i in 1:length(proteinIDs)){
		print(i)
		temp.positions<-prot.pos.list[[which(names(prot.pos.list)==proteinIDs[i])]]
		temp.prot.lines<-which(substr(object@data$protein,start=4,stop=9)==proteinIDs[i])  ### gives row numbers of of mods in object@data
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
					position.index.list[[temp.prot.lines[which(protein.position.list==x)][j]]]<-index2
					}
				index2=index2-1
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
			if(unique.indexes[[j]]<=0){
				templines.l<-length(templines)
				for(i in 1:templines.l){
					weighted.ratios[[templines[i]]]<-0
					}
				}
			}

	object@modsummary<-prot.pos.list
	object@modindex<-position.index.list
	object@data<-cbind(object@data,weighted.ratios=unlist(weighted.ratios),uniquesite.indexes=unlist(position.index.list))
	print("unique SUMO modified protein IDs")
	print(length(proteinIDs))
	print("unique SUMO modification sites from single-site peptides")
	print(index)
	length(unlist(prot.pos.list))
	### write something to make these text and put them in one column
	write.table(object@data,file="testallmods.tsv",quote=FALSE,sep="\t",row.names=F)
		

	#paste(unlist(prot.pos.list[1]))
	return(object)	
	}
