####################
# Define functions #
####################

### Interval Generator ###

Randomizer <- function(ChooseChr, EventSize, BedFile, BedGz, Chrs, i, Intervals){
  repeat {
    f <- as.numeric(1e-20)
    Start1 <- sample(as.numeric(SizeCentromerFile[ChooseChr,2]):as.numeric(SizeCentromerFile[ChooseChr,3]), 1, replace=F)
    Start2 <- sample(as.numeric(SizeCentromerFile[ChooseChr,4]):as.numeric(SizeCentromerFile[ChooseChr,5]), 1, replace=F)
    StartVec <- c(Start1, Start2)
    RandomStart <- sample(StartVec, 1)
    Stop <- as.numeric(RandomStart) + as.numeric(EventSize)
    
    ## Non overlap control ##
    
    SameIdx <- which(as.numeric(SizeCentromerFile[ChooseChr,1]) == as.numeric(Chrs))
    PreviousLines <- SameIdx[which(SameIdx < i)]
    OverlapControl <- 1
    if (length(PreviousLines) != 0)
    {
     OverlapVec <- c()
     for(k in 1:length(PreviousLines))
     {
       RecOverlap <- ReciprocalOverlap(c(RandomStart, Stop), Intervals[PreviousLines[k],], f)
       OverlapVec <- c(OverlapVec, RecOverlap)
     }
     if(sum(OverlapVec)!=0)
       {
       OverlapControl <- 0
       }
    }

    ## Bed control ##    

    BedControl <- 1
    if(BedFile == "Yes")
    {
      TabixString <- paste("tabix ", BedGz, " ", ChooseChr, ":", RandomStart, "-", Stop, " | wc -l", sep="")
      ExonsInRegion <- system(TabixString, intern=TRUE)
      if(as.numeric(ExonsInRegion) < 0)
      {
        BedControl <- 0
      }
    }
    if((Stop<as.numeric(SizeCentromerFile[ChooseChr,3]) | RandomStart>as.numeric(SizeCentromerFile[ChooseChr,4])) & (BedControl == 1) & (OverlapControl == 1)) break
  }
  return(c(RandomStart, Stop))
}

### Reciprocal Overlap Control ###

ReciprocalOverlap <- function(IntervalA,IntervalB,f)
{
  SA <- IntervalA[1]
  EA <- IntervalA[2]
  SB <- IntervalB[1]
  EB <- IntervalB[2]
  LA <- EA-SA
  LB <- EB-SB
  LOverlap <- 0
  if (SA<=SB & EA>=SB & EA<=EB)
  {
    LComm <- EA-SB
    LAComm <- LComm/LA
    LBComm <- LComm/LB
    if (LAComm>f & LBComm>f)
    {
      LOverlap<-1
    }
  }
  if (SA<=SB & EA>=EB)
  {
    LComm <- EB-SB
    LAComm <- LComm/LA
    LBComm <- LComm/LB
    if (LAComm>f & LBComm>f)
    {
      LOverlap <- 1
    }
  }
  if (SA>=SB & SA<=EB & EA>=EB)
  {
    LComm <- EB-SA
    LAComm <- LComm/LA
    LBComm <- LComm/LB
    if (LAComm>f & LBComm>f)
    {
      LOverlap <- 1
    }
  }
  if (SA>=SB & SA<=EB & EA<=EB)
  {
    LComm <- EA-SA
    LAComm <- LComm/LA
    LBComm <- LComm/LB
    if (LAComm>f & LBComm>f)
    {
      LOverlap <- 1
    }
  }
  Results <- LOverlap
  Results
}

### Open input ###

vars.tmp <- commandArgs()
vars <- vars.tmp[length(vars.tmp)]
split.vars <- unlist(strsplit(vars, ","))
FinalRData <- split.vars[1]
BedFile <- ""

load(FinalRData)

if(EventSize < 1000)
{
  EventSize <- 1000
}

load(file.path(FilesFolder,"xomeblender/scripts","ChrSizCentr.Rdata"))

### Bed Indexing ###

BedGz <- ""

if(Bed != '')
{
  BedGz <- paste(Bed, ".gz", sep="")
  if(file.exists(BedGz)==FALSE)
  {
    SortingString <- paste("sort -k1,1 -V -s", Bed, "| sed '/^$/d' >", paste(Bed, '.srt', sep =""), sep = " ")
    system(SortingString, wait = TRUE)
    BgzipString <- paste("bgzip -c",  paste(Bed, '.srt', sep =""), ">", BedGz, sep=" ")
    system(BgzipString, wait = TRUE)
    IndexingString <- paste("tabix -p bed", BedGz, sep=" ")
    system(IndexingString, wait = TRUE)
    rmline <-  paste("rm", Bed, sep = " ")
    mvline <- paste("mv", paste(Bed, '.srt', sep =""), Bed, sep = " ")
    system(rmline, wait = TRUE)
    system(mvline, wait = TRUE)
  }
  BedFile <- "Yes"
}

### Id Assignment ### 

#SampVec <- unlist(strsplit(IDs, " "))[-1]
SampVec <- IDs[-1]
SampCol <- c()

for(s in 1:EventsNumb)
{
  MatrixComboSample <- combn(SampVec, round(runif(1,1,length(SampVec))))
  indSel <- round(runif(1,1,ncol(MatrixComboSample)))
  SampCol <- rbind(SampCol, paste(MatrixComboSample[,indSel],collapse="-"))
}

### Chr Generation ###

Chrs <- sample(rep(c(1:22),c(22:1)),EventsNumb)
# Chrs <- round(runif(EventsNumb, min=1, max=22))

### Intervals Generation ###

Intervals <- c()
for(i in 1:length(Chrs))
{
  ChooseChr <- which(Chrs[i] == as.numeric(SizeCentromerFile[,1]))
  RandomResult <- Randomizer(ChooseChr, EventSize, BedFile, BedGz, Chrs, i, Intervals)
  Intervals <- rbind(Intervals, RandomResult)
}
ChrIntervals <- cbind(Chrs, Intervals)

### Outputting ###

## Chr Check ##

if(grepl("chr", ChrStyle)==TRUE)
{
  ChrVec <- rep("chr", length(Chrs))
  Chrs <- paste(ChrVec,ChrIntervals[,1], sep = "")
  ChrIntervals <- cbind(Chrs, ChrIntervals[,2:3])
}

## Del Dup Assignment ##

if(EventsNumb==2)
{
  DelDupVec <- c("Dup", "Del")
}else{
  DelDupVec <- sample(c("Dup","Del"),EventsNumb,replace=TRUE)
}

RefVec <- rep(Ref, length(Chrs))

ChrIntervals <- cbind(ChrIntervals, DelDupVec, SampCol, RefVec)
TabOut <- file.path(PathIn, paste(Label, "_CNV.txt", sep=""))
write.table(ChrIntervals, file=TabOut, sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)

