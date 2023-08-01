#' SNV data from individual A1
#'
#' Exemplary SNV data from individual A1 of the study Körber et al., Detecting and quantifying clonal selection in somatic mosaicism. The dataset is a list object, containing variant information in vcf format. 
#'
#' @format ## `snvs`
#' A list containing a data frame with 447 rows and 45 columns:
#' \describe{
#'   \item{Chr}{Chromosome}
#'   \item{Start, End}{Start and end position of the variant}
#'   \item{Ref, Alt}{Reference and alternative base par}
#'   \item{VAF, Depth, varCounts}{Variant allele frequency, read depth and number of variant reads}
#'   \item{VAF.control}{Variant allele frequency in the control data set}
#'   \item{Func.refGene}{Annovar annotation of the functional change (e.g., exonic, intergentic)}
#'   \item{GeneDetail.refGene}{Annovar annotation of the gene ID.}
#'   \item{ExonicFunc.refGene}{Annovar annotation of the exonic change in the gene (e.g. nonsynonymous)}
#'   \item{Gene.refGene}{Annovar annotation of the gene symbol}
#'   \item{AAChange.refGene}{Annovar annotation of the amino acid substitution}
#'   \item{avsnp150}{dbSNP identifier}
#'   \item{ExAC_ALL, ExAC_AFR, ExAC_AMR, ExAC_EAS, ExAC_FIN, ExAC_NFE, ExAC_OTH, ExAC_SAS}{exome aggregation consortium information}
#'   \item{AF, AF_popmax, AF_male, AF_female, AF_raw, AF_afr, AF_sas, AF_amr, AF_eas, AF_nfe, AF_fin, AF_asj, AF_oth, non_topmed_AF_popmax, non_neuro_AF_popmax, non_cancer_AF_popmax, controls_AF_popmax}{GnomAD annotated population-wide allele frequencies}
#'   \item{CLINALLELEID, CLNDN, CLNDISDB, CLNREVSTAT, CLNSIG}{Clinvar annotation}
#' }
"snvs"