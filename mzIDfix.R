
#### set current working directory to the location of your mzid file
setwd("C:/Users/JgMeyer/Documents/R/pepsum/pepsum")

### make an R object containing file names that have mzid suffix
files<-list.files(pattern="mzid")
#### print the files in the above object
files
### load required external libraries
library(mzID)
library(XML)
#test<-xmlTreeParse(files, useInternal = TRUE, , asText=TRUE)

### prints usage for parse command
?xmlParse

### parse and find peptide elements
### set the number in brackets to the correct file number in the "files" object
doc <- xmlParse(files[3])
r = xmlRoot(doc) #gives content of root

#r[[3]]["Peptide"]  ### gives me a list of the peptide nodes

#r[[3]]["Peptide"][[1]]  ### gives me the first peptide element 
### output from the above line
#<Peptide id="Pep1">
#  <PeptideSequence>GGSRNMGGPYGGGNYGPGGSGGSGGYGGRSRY</PeptideSequence>
#</Peptide> 

#typeof(r[[3]]["Peptide"][[1]])  ### reports type "externalptr"

#typeof(r[[3]]["Peptide"][[1]])

### get the peptide tags and store in tags object
tags <- xmlElementsByTagName(r[[3]], "Peptide")

#### function to fix all tags, run all lines in order
manipulate <- function(tag) {
    ## get 'Modification' node set
    dMod <- tag["Modification"]
    ## get 'location' numbers
    loc <- sapply(dMod, xmlGetAttr, "location")
    ## get the sum of 'monoisotopicMassDelta' 
    lapply(unique(loc), function(i) {
        if(length(d <- dMod[loc == i]) > 1) {
            nm <- "monoisotopicMassDelta"
            s <- sapply(d, xmlGetAttr, nm)
            xmlAttrs(d[[1]])[nm] <- sum(as.numeric(s))
        }
    })
    ## remove duplicated location nodes
    removeNodes(dMod[duplicated(loc)])
    ## return the adjusted tag
    tag
}

#### this command takes about 20 minutes, R appears to freeze but it does not actually, just let it work until you can interact with the GUI again
lapply(tags, manipulate)   ### applies manipulate to the XML structure to fix tags
?saveXML

#saveXML(r)
#### saves new XML file
saveXML(r, file="dh_fix.mzid", indent=TRUE)

### now have mzid file with tags fixed