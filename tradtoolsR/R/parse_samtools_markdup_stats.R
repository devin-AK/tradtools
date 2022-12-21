
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
    Line1 <- c("",x[,1])
    Line1 <- paste(Line1,collapse=',')
    Line2 <- c(ID,as.numeric(x[,2]))
    Line2 <- paste(Line2,collapse=',')
    writeLines(c(Line1,Line2),con=output_file)
  }  
  
  
