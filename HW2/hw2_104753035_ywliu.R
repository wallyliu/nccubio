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

input<-c(args[(input_flag+1)])
output<-c(args[output_flag+1])
score<-c(args[score_flag+1])
aln<-c(args[aln_flag+1])
gop<-as.numeric(args[gop_flag+1])
gep<-as.numeric(args[gep_flag+1])

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
    # Constuct ALLIGNMENT SCORE TABLE (globalTable)
    seqTitle1 <- c( c(strsplit(sequence[1], "")[[1]]) )
    seqTitle2 <- c( c(strsplit(sequence[2], "")[[1]]) )
    globalTable <- matrix(data=0, nrow=seq_length[1]+1, ncol=seq_length[2]+1, dimname=list(c("-", seqTitle1), c("-", seqTitle2)))

    # Initialize globalTable used to store ALIGNMENT SCORE.
    if( aln == "global")
    {
        tmp_gep <- ifelse(gep<0, -gep, gep)
        tmp_gop <- ifelse(gep<0, -gop, gop)

        t_v1 <- seq( tmp_gop-tmp_gep, tmp_gop+(seq_length[2]-1)*tmp_gep, tmp_gep)
        t_v2 <- seq( tmp_gop-tmp_gep, tmp_gop+(seq_length[1]-1)*tmp_gep, tmp_gep)
        
        if( gep<0 ){
            globalTable[1,] <- -t_v1
            globalTable[,1] <- -t_v2
        }else{
            globalTable[1,] <- t_v1
            globalTable[,1] <- t_v2
        }
        globalTable[1,1] <- 0
    }
    # Initialize tracebackMatrix used to stroe TRACKBACK INFORMATION
    if( aln == "global")
    {
        tracebackMatrix[1,] <<- "insertion"
        tracebackMatrix[,1] <<- "deletion"
        tracebackMatrix[1,1] <<- "substitution"
    }else
    {
        tracebackMatrix[1,] <<- "no"
        tracebackMatrix[,1] <<- "no"
    }

    # Calculate score word by word
    for(i in 1:length(seqTitle1))
    {
        for(j in 1:length(seqTitle2))
        {
            # decide insertion penalty
            insertion_penalty <- gop
            if( tracebackMatrix[i+1, j] == "substitution")
            {
                insertion_penalty <- gop
            }else if( tracebackMatrix[i+1, j] == "insertion" || tracebackMatrix[i+1, j] == "deletion" )
            {
                insertion_penalty <- gep
            }
            
            # decide deletion penalty
            deletion_penalty <- gop
            if( tracebackMatrix[i, j+1] == "substitution")
            {
                deletion_penalty <- gop
            }else if( tracebackMatrix[i, j+1] == "insertion" || tracebackMatrix[i, j+1] == "deletion" )
            {
                deletion_penalty <- gep
            }

            # get comparison character
            a<-substring(sequence[1], i, i)
            b<-substring(sequence[2], j, j)
 
            # calculate socre of each column
            if((a != "-")&&(b != "-"))
            {
                # Substitution, Insertion, Deletion
                if( aln == "global")
                {   # global alignment
                    value <- c( globalTable[i,j]+scoreMatrix[a,b], globalTable[i+1,j]+insertion_penalty, globalTable[i,j+1]+deletion_penalty)
                    globalTable[i+1,j+1] <- max(value)
                }else
                {   # local alignment
                    value <- c( globalTable[i,j]+scoreMatrix[a,b], globalTable[i+1,j]+insertion_penalty, globalTable[i,j+1]+deletion_penalty, 0)
                    globalTable[i+1,j+1] <- max(value)
                }
            }else{
                if( aln == "global")
                {   # global alignment
                    value <- max( globalTable[i+1,j]+insertion_penalty, globalTable[i,j+1]+deletion_penalty)
                    globalTable[i+1,j+1] <- max(value)
                }else
                {   # local alignment
                    value <- max( globalTable[i+1,j]+insertion_penalty, globalTable[i,j+1]+deletion_penalty, 0)
                    globalTable[i+1,j+1] <- max(value)
                }
            }

            # Constuct TRACEBACKMATRIX used to traceback sequence
            direction <- which( value==max(value))[1]
            if( direction == 1 ){
                tracebackMatrix[i+1,j+1] <<- "substitution"
            }else if( direction == 2){
                tracebackMatrix[i+1,j+1] <<- "insertion"
            }else if( direction == 3){
                tracebackMatrix[i+1,j+1] <<- "deletion"
            }else{
                tracebackMatrix[i+1,j+1] <<- "no"
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
    
    # Determine initial index i and j
    if( aln == "global")
    {   # global alignment
        i <- nrow(finalMatrix)
        j <- ncol(finalMatrix)
    }else
    {   # local alignment
        index <- which( finalMatrix == max(finalMatrix), arr.ind = TRUE)
        i <- index[nrow(index),][1]
        j <- index[nrow(index),][2]    
    }

    # Substitution, Insertion, Deletion
    while( i>1 && j>1 )
    {
        if( tracebackMatrix[i,j] == "substitution")
        {   # Substitution
            finalString1 <- paste( substring(sequence[1], i-1, i-1), finalString1)
            finalString2 <- paste( substring(sequence[2], j-1, j-1), finalString2)
            i <- i-1
            j <- j-1
        }else if( tracebackMatrix[i,j] == "insertion")
        {   # Insertion
            finalString1 <- paste( "-", finalString1)
            finalString2 <- paste( substring(sequence[2], j-1, j-1), finalString2)
            i <- i
            j <- j-1
        }else if( tracebackMatrix[i,j] == "deletion")
        {   # Deletion
            finalString1 <- paste( substring(sequence[1], i-1, i-1), finalString1)
            finalString2 <- paste( "-", finalString2)
            i <- i-1
            j <- j
        }else
        {   # Do Nothing -> local alignment
            finalString1 <- paste( substring(sequence[1], i-1, i-1), finalString1)
            finalString2 <- paste( substring(sequence[2], j-1, j-1), finalString2)
            return(c(finalString1, finalString2))
        }
    }
    return(c(finalString1, finalString2))
}

####################################
#          MAIN FUNCTION           #
####################################
tracebackMatrix <<- matrix(data="", nrow=nchar(sequence[1])+1, ncol=nchar(sequence[2])+1)
finalMatrix <- alignmentFunc(sequence, seq_length, aln_length, scoreMatrix, gop, gap, aln)
final <- tracebackString(finalMatrix, sequence, aln)
finalData <- c( paste(">",seq_name[1]), final[1], paste(">",seq_name[2]), final[2])
finalData <- gsub( " ", "", finalData)

# write data to output file
write(finalData, file = output)
