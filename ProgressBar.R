PB<-function(Methods,Rep)
{
###################################################################
## USAGE
 # setting the progress bar, used to shown the progress of the program
## ARGUMENTS
 # Methods: the name showed by the progress bar
 # Rep: number of loops
## VALUES
 # pb: progress_bar object
	library(progress)
	pb <- progress_bar$new(
         format = paste(Methods,":completed [:bar] :percent, Execute time::elapsed",sep=""),
         total = Rep, clear = FALSE, width= 60)
	return(pb)
}