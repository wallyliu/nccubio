####################################
#         Initialization           #
####################################
# read parameters
args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0)
{
	stop( "USAGE: Rscript pro2_104753035.R --input test.fasta --score pam250.txt --aln global --gap_open -10 --gap_extend -2 --output result.fasta --gap open", call.=FALSE)
}

# parse argument in command line
input_flag <- pmatch( "--input", args)
output_flag <- pmatch( "--output", args)
score_flag <- pmatch( "--score", args)
aln_flag <- pmatch( "--aln", args)
gop_flag <- pmatch( "--gap_open", args)
gep_flag <- pmatch( "--gap_extend", args)
gap_flag <- pmatch( "--gap", args)

input<-c(args[(input_flag+1)])
output<-c(args[output_flag+1])
score<-c(args[score_flag+1])
aln<-c(args[aln_flag+1])
gop<-as.numeric(args[gop_flag+1])
gep<-as.numeric(args[gep_flag+1])
gap<-c(args[gap_flag+1])


####################################
#       READ DATA FROM FILE        #
####################################
# VARABLE:                         #
#   sequence, seq_name, seq_length # 
#   aln_length, scoreMatrix        #
####################################

# (1) Read and Parse Fasta File
# VARIABLE: sequence, seq_name, aln_length
fastaFile <- read.table(input)
seq_name <- c(substring(fastaFile[1,1],2,), substring(fastaFile[3,1],2,))
sequence <- c(substring(fastaFile[2,1],1,), substring(fastaFile[4,1],1,))
seq_length <- c(nchar(sequence[1]), nchar(sequence[2]))
aln_length <- nchar(sequence[1])

# (2) Read Score File (PAM100 or 250)
scoreMatrix<-read.table(score)
scoreMatrix<-as.matrix(scoreMatrix)

####################################
#    CALCULATE GLOBAL ALIGNMENT    #
####################################
alignmentFunc <- function(sequence, seq_length, aln_length, scoreMatrix, gop, gap, aln)
{
    # constuct ALLIGNMENT SCORE TABLE
    seqTitle1 <- c( c(strsplit(sequence[1], "")[[1]]) )
    seqTitle2 <- c( c(strsplit(sequence[2], "")[[1]]) )
    globalTable <- matrix(data=0, nrow=seq_length[1]+1, ncol=seq_length[2]+1, dimname=list(c("-", seqTitle1), c("-", seqTitle2)))

    # Initialize globalTable which is be used to store ALIGNMENT SCORE.
    if( aln == "global")
    {
        tmp_gep <- ifelse(gep<0, -gep, gep)
        t_v1 <- seq(0, seq_length[2]*tmp_gep, tmp_gep)
        t_v2 <- seq(0, seq_length[1]*tmp_gep, tmp_gep)
        if( gep<0 ){
            globalTable[1,] <- -t_v1
            globalTable[,1] <- -t_v2
        }else{
            globalTable[1,] <- t_v1
            globalTable[,1] <- t_v2
        }
    }

    # calculate score word by word
    for(i in 1:length(seqTitle1))
    {
        for(j in 1:length(seqTitle2))
        {
            a<-substring(sequence[1], i, i)
            b<-substring(sequence[2], j, j)
 
            openPenalty <- 0
            if( gap =="open"){
                if( i==1 && j==1 )
                {
                    openPenalty <- gop
                }else if( !is.na( which( c( tracebackMatrix[i-1][j-1], tracebackMatrix[i][j-1], tracebackMatrix[i-1][j]) =="substitution")[1] ))
                {
                    openPenalty <- gop
                }
            }
            if((a != "-")&&(b != "-"))
            {
                # substitution, insertion, deletion
                if( aln == "global")
                {
                    value <- c( globalTable[i,j]+scoreMatrix[a,b], globalTable[i+1,j]+gep+openPenalty, globalTable[i,j+1]+gep+openPenalty)
                    globalTable[i+1,j+1] <- max(value)
                }else{
                    value <- c( globalTable[i,j]+scoreMatrix[a,b], globalTable[i+1,j]+gep+openPenalty, globalTable[i,j+1]+gep+openPenalty, 0)
                    globalTable[i+1,j+1] <- max(value)
                }
            }else{
                if( aln == "global")
                {
                    value <- max( globalTable[i+1,j]+gep+openPenalty, globalTable[i,j+1]+gep+openPenalty)
                    globalTable[i+1,j+1] <- max(value)
                }else{
                    value <- max( globalTable[i+1,j]+gep+openPenalty, globalTable[i,j+1]+gep+openPenalty, 0)
                    globalTable[i+1,j+1] <- max(value)
                }
            }

            # constuct TRACEBACKMATRIX used to traceback sequence
            direction <- which( value==max(value))[1]
            if( direction == 1 ){
                tracebackMatrix[i,j] <<- "substitution"
            }else if( direction == 2){
                tracebackMatrix[i,j] <<- "insertion"
            }else{
                tracebackMatrix[i,j] <<- "deletion"
            }
        }
    }
    return(globalTable)
}

####################################
#     TRACKBACK TWO SEQUENCES      #
####################################
tracebackString <- function( finalMatrix, sequence, aln)
{
    finalString1 <- ""
    finalString2 <- ""
    
    if( aln == "global")
    {
        i <- nrow(finalMatrix)
        j <- ncol(finalMatrix)
    
        # substitution, insertion, deletion
        while( i>1 && j>1 )
        {
            if( tracebackMatrix[i-1,j-1] == "substitution")
            {   # substitution
                finalString1 <- paste( substring(sequence[1], i-1, i-1), finalString1)
                finalString2 <- paste( substring(sequence[2], j-1, j-1), finalString2)
                i <- i-1
                j <- j-1
            }else if( tracebackMatrix[i-1,j-1] == "insertion")
            {   # insertion
                finalString1 <- paste( "-", finalString1)
                finalString2 <- paste( substring(sequence[2], j-1, j-1), finalString2)
                i <- i
                j <- j-1
            }else
            {   # deletion
                finalString1 <- paste( substring(sequence[1], i-1, i-1), finalString1)
                finalString2 <- paste( "-", finalString2)
                i <- i-1
                j <- j
            }
        }
    }else{
        index <- which( finalMatrix == max(finalMatrix), arr.ind = TRUE)
        i <- index[1]
        j <- index[2]
        
        # substitution, insertion, deletion
        while( i>1 && j>1 )
        {
            if( tracebackMatrix[i-1,j-1] == "substitution")
            {   # substitution
                finalString1 <- paste( substring(sequence[1], i-1, i-1), finalString1)
                finalString2 <- paste( substring(sequence[2], j-1, j-1), finalString2)
                i <- i-1
                j <- j-1
            }else if( tracebackMatrix[i-1,j-1] == "insertion")
            {   # insertion
                finalString1 <- paste( "-", finalString1)
                finalString2 <- paste( substring(sequence[2], j-1, j-1), finalString2)
                i <- i
                j <- j-1
            }else
            {   # deletion
                finalString1 <- paste( substring(sequence[1], i-1, i-1), finalString1)
                finalString2 <- paste( "-", finalString2)
                i <- i-1
                j <- j
            }
        }
    }
    return(c(finalString1, finalString2))
}

####################################
#          MAIN FUNCTION           #
####################################
tracebackMatrix <<- matrix(data="", nrow=nchar(sequence[1]), ncol=nchar(sequence[2]))
finalMatrix <- alignmentFunc(sequence, seq_length, aln_length, scoreMatrix, gop, gap, aln)
final <- tracebackString(finalMatrix, sequence, aln)
finalData <- c( paste(">",seq_name[1]), final[1], paste(">",seq_name[2]), final[2])
finalData <- gsub( " ", "", finalData)

# write data to output file
write(finalData, file = output)
