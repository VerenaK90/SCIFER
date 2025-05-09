% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Extract_info_from_vcf.R
\name{Extract.info.from.vcf}
\alias{Extract.info.from.vcf}
\title{Extracts information from a vcf file.}
\usage{
Extract.info.from.vcf(
  vcf,
  info = "readcounts",
  type = "snvs",
  mutationcaller = "Strelka",
  tumor.col.mutect = 10,
  normal.col.mutect = 11,
  sample.col.mpileup = NA,
  tumor.id = NULL
)
}
\arguments{
\item{vcf}{Mutation information in VCF format represented as a list (as returned by read.vcf from package bedR).}

\item{info}{Variant information to be retrieved. Possible values are
\describe{
\item{\code{readcounts}}{ returns the number of reference and variant reads for each variant position.}
\item{\code{varCounts}}{ the number of variant reads (if mutationcaller is 'Strelka' or 'Manta').}
\item{\code{depth}}{ returns the sequencing depths for each variant position.}
\item{\code{VAF}}{ returns the variant allele frequency at each variant position.}
\item{\code{VAF.control}}{ returns the variant allele frequency in the germline control at each variant position.}
\item{\code{depth.control}}{ returns the sequencing depth in the germline control at each variant position}
\item{\code{AA_change}}{ returns the amino acid change at each variant position (only works if vcf-file had been annotated with annovar).}
\item{\code{cDNA_change}}{ returns the cDNA change at each variant position (only works if vcf-file had been annotated with annovar).}
\item{\code{Exon}}{ returns the exon targeted by each variant (only works if vcf-file had been annotated with annovar).}
\item{\code{Gene}}{ returns the gene targeted by each variant (only works if vcf-file had been annotated with annovar).}
\item{\code{annovar_function}}{ returns the consequence of the variant (e.g. "intronic"; only works if vcf-file had been annotated with annovar).}
\item{\code{exonic_function}}{ returns the exonic consequence of each variant (e.g. 'nonsynonymous SNV'; only works if vcf-file had been annotated with annovar).}
\item{\code{svtype}}{ type of structural variant at each variant position; only works on output generated with Manta.}
\item{\code{somaticScore}}{ the somatic score computed by Manta for each variant; only works on output generated with Manta.}
\item{\code{End}}{ the end position of a structural variant; only works with output generated by Manta.}
}}

\item{type}{Specify the variant type (\code{"snvs"} or \code{"indel"}, defaults to \code{"snvs"}). Ignored if \code{mutationcaller = "Mutect2"} or  \code{mutationcaller = "Manta"}.}

\item{mutationcaller}{The mutation caller that generated the vcf-files. Must be either \code{"Strelka"}, \code{"Mutect2"}, \code{"mpileup"} or \code{"Manta"}. Defaults to \code{"Strelka"}.}

\item{tumor.col.mutect}{Column index or name of the tumor information when \code{mutationcaller="Mutect2"}.}

\item{normal.col.mutect}{Column index or name of the germline control information when \code{mutationcaller="Mutect2"}.}

\item{tumor.id}{The name of the tumor sample. Only used if vcf-file was generated with Manta.}
}
\value{
The requested information for each variant position.
}
\description{
Extracts information from a vcf file.
}
