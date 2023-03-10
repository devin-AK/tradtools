% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build_signal_track.R
\name{build_signal_track}
\alias{build_signal_track}
\title{TRADtools build_signal_track}
\usage{
\code{build_signal_track(...)}
}
\arguments{
\item{input_bam}{PATH to the input BAM file. BAM file should be pre-processed and sorted}

\item{input_blacklist}{PATH to the blacklist BED file to be used for filtering (e.g. "/path/to/hg19-blacklist.v2.bed.gz"), can be NULL (or "NULL" or "null")}

\item{output_bed}{PATH to the output BED file (will not be overwritten if already exists, unless overwrite=TRUE)}

\item{output_bw}{PATH to the output BW file (will not be overwritten if already exists, unless overwrite =T)}

\item{reference_fasta}{PATH to the reference FASTA file that was used for alignment. This provides the sequences that will be retrieved after applying ranges_operation() to the BAM reads}

\item{ranges_operation}{can be "default" (two bases immediately 5' to R1 read), "identity" (no transformation), or a user-defined function}

\item{mapQ}{Integer. Reads with mapQ below this number will be discarded}

\item{standard_chromosomes}{TRUE or FALSE. If TRUE, retain only the standard chromosomes (by using \code{GenomeInfoDb::keepStandardChromosomes()})}

\item{overwrite}{TRUE or FALSE. If TRUE, output BED file will be overwritten even if it already exists}

\item{verbose}{TRUE or FALSE. Whether to report helpful messages and summary statistics}

\item{export_stats}{TRUE or FALSE. Whether to write out stats files from the signal track processing (to the current directory)}
}
\value{
Nothing returned. Files are exported.
}
\description{
This function processes a BAM file into a TRAD-seq signal track file in BED format
}
\details{
The exported BED file contains reads after applying mapQ filter and transformation defined by ranges_operation. The exported BW file reflects the "final" signal after further filtering out non-dipyrimidine records and records that overlap a black-listed region.
}
\examples{
\code{build_signal_track(gen1.bam, bl.bed.gz, gen1.bed, gen1.bw, gen1.fasta)}
}
