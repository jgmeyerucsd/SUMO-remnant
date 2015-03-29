

read.nature=function(input=files[90],type="nature",any=F,studies=T){
	object<-new("pepsum")
	object@filetype=type

	object@data<-read.delim(input,header=T)

	if(any){
		protacc<-levels(object@data[,1])
		for(x in protacc){
			object@modsummary[[x]]<-object@data[which(object@data[,1]==x),2]
			}
		}
	if(studies){
		protacc<-levels(object@data[,1])
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
