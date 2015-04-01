
### function to read peptide prophet output in tab-delim format
### usage:
###  newobject<-read.PepProph(input=[your file path])

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

	
	