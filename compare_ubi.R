hist(logRatiosAllnorm,seq(from=-12,to=10,by=0.1))
hist(logRatiosAllnorm,seq(from=-12,to=10,by=0.1))


files<-list.files()
files

read.nature=function(input=files[126],type="nature",any=F,studies=T){
	object2<-new("pepsum")
	object@filetype=type

	object2@data<-read.delim(input,header=T)
	object2@data[1,]

	motifs<-as.character(object@data$Motif)
	matched<-match(motifs,object@sequence)
	object@sequence[na.omit(matched)]
	object@data[na.omit(matched),]
	write.table(file="matched.with.ubi.tsv",object@data[na.omit(matched),],quote=F,sep="\t",col.names=T,row.names=F)

	matched[matched!="NA"]
	object@which((matched!=is.na(matched))==TRUE)

	protacc<-levels(object2@data[,2])

	study.lines<-object@data[object@data[,"Count....Studies.SUMO.2.Modified."]>=1,]
	study.lines.l<-length(study.lines[,1])
	protacc<-levels(as.factor(as.character(study.lines[,1])))
	for(x in protacc){
		object@modsummary[[x]]<-object@data[which(object@data[,1]==x),2]
		}
	}

	return(object)
	}




#object<-nature
#read.nature(input=files[90])->nature
nature@data[1,]
