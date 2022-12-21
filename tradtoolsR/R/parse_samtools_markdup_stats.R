
  # input_file  <- '/Users/dk432/Desktop/multiqc/SUV39H1_200J_UVB_1_markdup_stats.txt'
  # output_file <- '/Users/dk432/Desktop/multiqc/SUV39H1_200J_UVB_1_trad_markdups.csv'
  # parse_samtools_markdup_stats(input_file,output_file)

  #' @export
  parse_samtools_markdup_stats <- function(input_file,output_file) {
    ID <- tools:::file_path_sans_ext(basename(input_file))
    x <- readLines(input_file)
    x <- lapply(x,function(i){
      dat <- trimws(unlist(strsplit(i,':')))
      if(length(dat)>2) NULL else {
        c(dat[1],dat[2])
      }
    })
    x <- do.call(rbind,x)
    x <- as.data.frame(x,row.names=x[,1])
    duplicate_pair_optical <- as.numeric(x['DUPLICATE PAIR OPTICAL',2])
    duplicate_pair_other   <- as.numeric(x['DUPLICATE PAIR',2]) - duplicate_pair_optical
    unique_pair            <- as.numeric(x['PAIRED',2]) - (duplicate_pair_optical+duplicate_pair_other)
    unpaired_excluded      <- as.numeric(x['READ',2]) - (unique_pair+duplicate_pair_optical+duplicate_pair_other)
    Line1 <- c('','unique_pair','unpaired_excluded','duplicate_pair_optical','duplicate_pair_other')
    Line1 <- paste(Line1,collapse=',')
    Line2 <- c(ID, unique_pair,  unpaired_excluded,  duplicate_pair_optical,  duplicate_pair_other)
    Line2 <- paste(Line2,collapse=',')
    #Line1 <- c("",x[,1])
    #Line1 <- paste(Line1,collapse=',')
    #Line2 <- c(ID,as.numeric(x[,2]))
    #Line2 <- paste(Line2,collapse=',')
    writeLines(c(Line1,Line2),con=output_file)
  }  
  
  
