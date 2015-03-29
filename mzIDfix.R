setwd("C:/Users/JgMeyer/Documents/R/pepsum/pepsum")

files<-list.files(pattern="mzid")
files
library(mzID)
library(XML)
test<-xmlTreeParse(files, useInternal = TRUE, , asText=TRUE)
?xmlParse

### parse and gind peptide elements
doc <- xmlParse(files[3])
r = xmlRoot(doc) #gives content of root
r[[3]]["Peptide"]  ### gives me a list of the peptide nodes

r[[3]]["Peptide"][[1]]  ### gives me the first peptide element 
### output from the above line
#<Peptide id="Pep1">
#  <PeptideSequence>GGSRNMGGPYGGGNYGPGGSGGSGGYGGRSRY</PeptideSequence>
#</Peptide> 

typeof(r[[3]]["Peptide"][[1]])  ### reports type "externalptr"

typeof(r[[3]]["Peptide"][[1]])
tags <- xmlElementsByTagName(r[[3]], "Peptide")

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

lapply(tags, manipulate)
?saveXML

saveXML(r)
saveXML(r, file="dh_fix.mzid", indent=TRUE)
xmlSize(r)
xmlName(r[[1]])
xmlName(r[[2]])
xmlName(r[[3]]["Peptide"])
xmlName(r[[1]])
xmlName(r[[2]])
xmlName(r[[3]])



xmltop2 = xmlRoot(test)
pepelem<-xmlElementsByTagName(xmltop[[3]],name="Peptide")
pepelem2<-xmlElementsByTagName(xmltop2[[3]],name="Peptide")

pepelem[[1]]
typeof(pepelem[[1]])



pepelem2[[1239]]

typeof(pepelem2[[1]])


testchildren<-xmlChildren(xmltop[[3]])
typeof(testchildren)
testchildren$peptide
names(testchildren)=="Peptide"

xmlElementSummary(xmlfile)
?readLines

## Parse the XML file
doc <- xmlTreeParse("sample.xml", useInternal = TRUE)
## Select the nodes we want to update
nodes <- getNodeSet(test, "//peptide")
## For each node, apply gsub on the content of the node
lapply(nodes, function(n) {
  xmlValue(n) <- gsub("ABC","foobar",xmlValue(n))
})


nct_url <- "http://clinicaltrials.gov/ct2/show/NCT00112281?resultsxml=true"
xml_doc <- xmlParse(nct_url, useInternalNode=TRUE)

mod_path <- "/MzIdentML/peptide/modification" 
mod_text <- xpathSApply(xmlfile, mod_path, xmlValue)


# read in to a tree:
x = xmlParse("test.xml")

# this returns a *list* of text nodes under sequence
# and NOT the text nodes under taxon
nodeSet = xpathApply(x,"//sequence/text()")
nodeSet = xpathApply(xmlfile,"//SequenceCollection/text()")

# now we loop over the list returned, and get and modify the node value:
sapply(nodeSet,function(G){
  text = paste("CCCCC",xmlValue(G),"GGGGGGG",sep="")
  text = gsub("[^A-Z]","",text)
  xmlValue(G) = text
})


