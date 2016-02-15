### Open input ###

vars.tmp <- commandArgs()
vars <- vars.tmp[length(vars.tmp)]
split.vars <- unlist(strsplit(vars, ","))
IdVec <- split.vars[1]
PathToGold <- split.vars[2]
WorkingPath <- split.vars[3]
Gold <- read.table(PathToGold,sep="\t",header=F,quote = "\"",fill=T)

### Header Parsing ###

StringHeader<-paste(" head -200 ",PathToGold," | grep '#'",sep="")
TotalHeader<-system(StringHeader,intern = TRUE)
NLine<-length(TotalHeader)
HeaderNames<-TotalHeader[NLine]
Header<-TotalHeader[-NLine]
HeaderSplit<-unlist(strsplit(HeaderNames,"\t"))

### IdVec Parsing ###

IdVecSplit<-unlist(strsplit(IdVec, "_"))
IdNumber <- length(IdVecSplit)
ControlSamp <- IdVecSplit[1]

MatchingTab<-c()
for (s in 2:IdNumber)
{
	TumSamp <- IdVecSplit[s]
	ControlInGold <- match(ControlSamp, HeaderSplit)
	TumInGold <- match(TumSamp, HeaderSplit)
	MatchingTabIdx <- sort(c(which(substr(Gold[,ControlInGold], 1, 3)=="0/0" & substr(Gold[,TumInGold], 1, 3)=="0/1"),which(substr(Gold[,ControlInGold], 1, 3)=="0/0" & substr(Gold[,TumInGold], 1, 3)=="1/1")))
	MatchingTab <- rbind(MatchingTab,cbind(Gold[MatchingTabIdx,c(1,2,4,5)], substr(Gold[MatchingTabIdx,ControlInGold],start=1,stop=3), substr(Gold[MatchingTabIdx,TumInGold],start=1,stop=3),rep(TumSamp,length(MatchingTabIdx))))
	TabName <- file.path(WorkingPath,paste(TumSamp, ".somatic.txt", sep=""))
	write(t(MatchingTab),file=TabName,sep = "\t",ncolumns=ncol(MatchingTab),append=TRUE)
}
