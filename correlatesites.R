
correlate.positions()

correlate.positions=function(object1=mgpl.s.sig,object2=mgpl.sig)
	{
	modsum1<-object1[[2]]@modsummary
	modsum2<-object2[[2]]@modsummary
	names(modsum1)<-substr(names(modsum1),start=4,stop=9)
	names(modsum2)<-substr(names(modsum2),start=4,stop=9)
	i=1
	matchvec<-c()
	for(x in names(modsum1)){	
		print(which(names(modsum2)==x))
		if(length(which(names(modsum2)==x))>0){
			matchvec[i]<-x
			i=i+1
			}
		}
	#### loop through the match names and check their positions
	#### in the larger matrix
	#### take value of quant for each site put into a list of pairs
	quantpairs<-list()
	length(unlist(modsum1[matchvec]))
	modsum1.sub<-modsum1[matchvec]	
	modsum2.sub<-modsum2[matchvec]
	uniprot<-substr(object2[[2]]@data[,"protein"],start=4,stop=9)

	object2[[2]]@data<-cbind(object2[[2]]@data,uniprot=uniprot)
	uniprot<-substr(object1[[2]]@data[,"protein"],start=4,stop=9)

	object1[[2]]@data<-cbind(object1[[2]]@data,uniprot=uniprot)
	hprp.values<-c()
	single.values<-c()
	for(x in names(modsum1.sub)){
		matched.positions<-modsum2.sub[[x]][na.omit(match(modsum1.sub[[x]],modsum2.sub[[x]]))]
		for(y in matched.positions){
			obj2.protlines<-object2[[2]]@data[which(as.character(object2[[2]]@data[,"uniprot"])==x),]
			hprp.value<-obj2.protlines[,"weighted.ratios"][obj2.protlines[,1]==y][1]
			obj1.protlines<-object1[[2]]@data[which(as.character(object1[[2]]@data[,"uniprot"])==x),]
			single.value<-obj1.protlines[,"weighted.ratios"][obj1.protlines[,1]==y][1]
			if(is.na(single.value)==is.na(hprp.value)){
				quantpairs[[paste(x,"_",y,sep="")]]<-c(single.value,hprp.value)
				hprp.values<-c(hprp.values,hprp.value)
				single.values<-c(single.values,single.value)
				}
			}				
		}
	
	par(cex.axis=1.5,cex=1.25)
	plot(log(hprp.values,base=2),log(single.values,base=2),ylim=c(),ylab="standard",xlab="hprp",main="correlation of SILAC quant values",pch=20)
	y=log(hprp.values,base=2)
	x=log(single.values,base=2)
	lm(x ~ y)
	abline(lm(x~y+0),col="red",lwd=2)
	abline(h=0)
	abline(v=0)
	print(summary(lm(x~y+0)))
	}
	
	