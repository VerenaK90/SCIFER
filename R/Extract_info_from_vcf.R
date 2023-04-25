#' Extracts information from a vcf file
#'
#' @param vcf Mutation file in VCF representation as a list (as returned by read.vcf from package bedR).
#' @param info Which information is to be retrieved? Possible values are "readcounts", "depth", "VAF", "AA_change", "cDNA_change", "Exon", "VAF.control", "depth.control", "Gene", "annovar_function", "exonic_function". Defaults to "readcounts".
#' @param type Is the VCF file reporting snvs or indels? Possible values are "snvs" and "indel". Defaults to "snvs".
#' @param mutationcaller which mutation caller was used to generate the vcf-files? Defaults to "Strelka".
#' @return purity and ploidy
#' @export


Extract.info.from.vcf <- function(vcf, info="readcounts", type="snvs", mutationcaller="Strelka", tumor.col.mutect=10,
                                  normal.col.mutect=11, sample.col.mpileup=NA){
  
  if(length(vcf$vcf)==0){
    vcf <- list(vcf=vcf)
  }
  
  if(mutationcaller=="Strelka"){
    
    ref.column <- which(colnames(vcf$vcf)=="REF")
    alt.column <-  which(colnames(vcf$vcf)=="ALT")
    normal.column <- which(colnames(vcf$vcf)=="NORMAL")
    tumor.column <- which(colnames(vcf$vcf)=="TUMOR")
    if("FORMAT" %in% colnames(vcf$vcf)){
      SUBDP <- which(strsplit(vcf$vcf$FORMAT[1], split=":")[[1]]=="SUBDP")
    }
    
    if(type=="snvs"){
      if(info=="readcounts"){
        readcounts <- t(apply(vcf$vcf, 1, function(x){
          
          ref <- which(c("A", "C", "G", "T")==as.character(x[ref.column]))
          alt <- which(c("A", "C", "G", "T")==as.character(x[alt.column]))
          x <- strsplit(as.character(x[11]), split=":")[[1]]
          ref <- strsplit(x[SUBDP+ref], split=",")[[1]][1]
          alt <- strsplit(x[SUBDP+alt], split=",")[[1]][1]
          return(c(as.numeric(ref),as.numeric(alt)))
        }))
        colnames(readcounts) <- c("REF", "ALT")
        
        return(readcounts)
      }
      if(info=="varCounts"){
        alt <- t(apply(vcf$vcf, 1, function(x){
          
          alt <- which(c("A", "C", "G", "T")==as.character(x[alt.column]))
          x <- strsplit(as.character(x[11]), split=":")[[1]]
          alt <- strsplit(x[SUBDP+alt], split=",")[[1]][1]
          return(as.numeric(alt))
        }))

        return(alt)
      }
      if(info=="depth"){
        depth <- apply(vcf$vcf, 1, function(x){
          
          ref <- which(c("A", "C", "G", "T")==as.character(x[ref.column]))
          alt <- which(c("A", "C", "G", "T")==as.character(x[alt.column]))
          x <- strsplit(as.character(x[tumor.column]), split=":")[[1]]
          ref <- strsplit(x[SUBDP+ref], split=",")[[1]][1]
          alt <- strsplit(x[SUBDP+alt], split=",")[[1]][1]
          return(as.numeric(ref) +as.numeric(alt))
        })
        
        return(depth)
      }
      if(info=="VAF"){
        vaf <- apply(vcf$vcf, 1, function(x){
          
          ref <- which(c("A", "C", "G", "T")==as.character(x[ref.column]))
          alt <- which(c("A", "C", "G", "T")==as.character(x[alt.column]))
          x <- strsplit(as.character(x[tumor.column]), split=":")[[1]]
          ref <- strsplit(x[SUBDP+ref], split=",")[[1]][1]
          alt <- strsplit(x[SUBDP+alt], split=",")[[1]][1]
          return(as.numeric(alt)/(as.numeric(ref) +as.numeric(alt)))
        })
        return(vaf)
      }
      if(info=="AA_change"){
        protein_change <- sapply(vcf$vcf$INFO, function(x){
          x <- strsplit(x, split="AAChange")[[1]][2]
          x <- strsplit(x, split=":")[[1]][5]
          x <- strsplit(x, split=";")[[1]][1]
          x <- strsplit(x, split=",")[[1]][1]
        })
        aa_change <- sapply(protein_change, function(x){
          x <- strsplit(x, split="[.]")[[1]][2]
          #x <- strsplit(x, split="")[[1]]
          #paste0(x[c(1,length(x))], collapse = "")
        })
        return(aa_change)
      }
      if(info=="cDNA_change"){
        cDNA_change <- sapply(vcf$vcf$INFO, function(x){
          x <- strsplit(x, split="AAChange")[[1]][2]
          x <- strsplit(x, split=":")[[1]][4]
        })
        return(cDNA_change)
      }
      if(info=="Exon"){
        cDNA_change <- sapply(vcf$vcf$INFO, function(x){
          x <- strsplit(x, split="AAChange")[[1]][2]
          x <- strsplit(x, split=":")[[1]][3]
        })
        return(cDNA_change)
      }
      if(info=="VAF.control"){
        VAF.control <- apply(vcf$vcf, 1, function(x){
          ref <- which(c("A", "C", "G", "T")==as.character(x[ref.column]))
          alt <- which(c("A", "C", "G", "T")==as.character(x[alt.column]))
          x <- strsplit(as.character(x[normal.column]), split=":")[[1]]
          ref <- strsplit(x[SUBDP+ref], split=",")[[1]][1]
          alt <- strsplit(x[SUBDP+alt], split=",")[[1]][1]
          return(as.numeric(alt)/(as.numeric(alt) + as.numeric(ref)))
        })
        return(VAF.control)
      }
      if(info=="depth.control"){
        depth.control <- apply(vcf$vcf, 1, function(x){
          ref <- which(c("A", "C", "G", "T")==as.character(x[ref.column]))
          alt <- which(c("A", "C", "G", "T")==as.character(x[alt.column]))
          x <- strsplit(as.character(x[normal.column]), split=":")[[1]]
          ref <- strsplit(x[SUBDP+ref], split=",")[[1]][1]
          alt <- strsplit(x[SUBDP+alt], split=",")[[1]][1]
          return(as.numeric(alt) + as.numeric(ref))
        })
        return(depth.control)
      }
      if(info=="Gene"){
        gene <- sapply(vcf$vcf$INFO, function(x){
          x <- strsplit(x, split="Gene.refGene=")[[1]][2]
          x <- strsplit(x, split=";")[[1]][1]
          x
        })
        return(unname(gene))
      }
      if(info=="annovar_function"){
        annovar_function <- sapply(vcf$vcf$INFO, function(x){
          x <- strsplit(x, split="Func.refGene=")[[1]][2]
          x <- strsplit(x, split=";")[[1]][1]
          x
        })
        return(unname(annovar_function))
      }
      if(info=="exonic_function"){
        exonic_function <- sapply(vcf$vcf$INFO, function(x){
          x <- strsplit(x, split="ExonicFunc.refGene=")[[1]][2]
          x <- strsplit(x, split=";")[[1]][1]
          x
        })
        return(unname(exonic_function))
      }
    }else if(type=="indel"){
      
      if(info=="readcounts"){
        readcounts <- t(sapply(vcf$vcf$TUMOR, function(x){
          ref <- strsplit(x, split=":")[[1]][ref.column]
          alt <- strsplit(x, split=":")[[1]][alt.column]
          ref <- strsplit(ref, split=",")[[1]][1]
          alt <- strsplit(alt, split=",")[[1]][1]
          alt <- unname(alt)
          ref <- unname(ref)
          return(c(as.numeric(ref),as.numeric(alt)))
        }))
        colnames(readcounts) <- c("REF", "ALT")
        return(readcounts)
      }
      if(info=="varCounts"){
        varCounts <- sapply(vcf$vcf$TUMOR, function(x){
          ref <- strsplit(x, split=":")[[1]][ref.column]
          alt <- strsplit(x, split=":")[[1]][alt.column]
          ref <- strsplit(ref, split=",")[[1]][1]
          alt <- strsplit(alt, split=",")[[1]][1]
          alt <- unname(alt)
          ref <- unname(ref)
          return(as.numeric(alt))
        })
        return(varCounts)
      }
      if(info=="depth"){
        depth <- sapply(vcf$vcf$TUMOR, function(x){
          ref <- strsplit(x, split=":")[[1]][ref.column]
          alt <- strsplit(x, split=":")[[1]][alt.column]
          ref <- strsplit(ref, split=",")[[1]][1]
          alt <- strsplit(alt, split=",")[[1]][1]
          return(as.numeric(ref)+as.numeric(alt))
        })
        return(unname(depth))
      }
      if(info=="VAF"){
        vaf <- sapply(vcf$vcf$TUMOR, function(x){
          ref <- strsplit(x, split=":")[[1]][ref.column]
          alt <- strsplit(x, split=":")[[1]][alt.column]
          ref <- strsplit(ref, split=",")[[1]][1]
          alt <- strsplit(alt, split=",")[[1]][1]
          return(as.numeric(alt)/(as.numeric(ref)+as.numeric(alt)))
        })
        return(unname(vaf))
      }
      if(info=="AA_change"){
        protein_change <- sapply(vcf$vcf$INFO, function(x){
          x <- strsplit(x, split="AAChange")[[1]][2]
          x <- strsplit(x, split=":")[[1]][5]
          x <- strsplit(x, split=";")[[1]][1]
          x <- strsplit(x, split=",")[[1]][1]
          
        })
        aa_change <- sapply(protein_change, function(x){
          x <- strsplit(x, split="[.]")[[1]][2]
          #x <- strsplit(x, split="")[[1]]
          #paste0(x[c(1,length(x))], collapse = "")
        })
        return(aa_change)
      }
      if(info=="cDNA_change"){
        cDNA_change <- sapply(vcf$vcf$INFO, function(x){
          x <- strsplit(x, split="AAChange")[[1]][2]
          x <- strsplit(x, split=":")[[1]][4]
        })
        return(cDNA_change)
      }
      if(info=="Exon"){
        cDNA_change <- sapply(vcf$vcf$INFO, function(x){
          x <- strsplit(x, split="AAChange")[[1]][2]
          x <- strsplit(x, split=":")[[1]][3]
        })
        return(cDNA_change)
      }
      if(info=="VAF.control"){
        VAF.control <- sapply(vcf$vcf$NORMAL, function(x){
          ref <- strsplit(x, split=":")[[1]][ref.column]
          alt <- strsplit(x, split=":")[[1]][alt.column]
          ref <- strsplit(ref, split=",")[[1]][1]
          alt <- strsplit(alt, split=",")[[1]][1]
          return(as.numeric(alt)/(as.numeric(alt) + as.numeric(ref)))
        })
        return(unname(VAF.control))
      }
      if(info=="Gene"){
        gene <- sapply(vcf$vcf$INFO, function(x){
          x <- strsplit(x, split=";")[[1]][17]
          strsplit(x, split="=")[[1]][2]
        })
        return(unname(gene))
      }
      if(info=="annovar_function"){
        annovar_function <-  sapply(vcf$vcf$INFO, function(x){
          x <- strsplit(x, split="Func.refGene=")[[1]][2]
          x <- strsplit(x, split=";")[[1]][1]
          x
        })
        return(unname(annovar_function))
      }
      if(info=="exonic_function"){
        exonic_function <- sapply(vcf$vcf$INFO, function(x){
            x <- strsplit(x, split="ExonicFunc.refGene=")[[1]][2]
            x <- strsplit(x, split=";")[[1]][1]
            x
        })
        return(unname(exonic_function))
      }
      
    }
  }else if(mutationcaller=="Mutect2"){
    
      if(info=="readcounts"){
        readcounts <- t(apply(vcf$vcf, 1, function(x){
          
          ref.alt <- strsplit(x[9], split=":")[[1]]
          ref.alt <- which(ref.alt=="AD")

          ref.alt <- strsplit(x[tumor.col.mutect], split=":")[[1]][ref.alt]
          ref <- strsplit(ref.alt, split=",")[[1]][1]
          alt <- strsplit(ref.alt, split=",")[[1]][2]

          return(c(as.numeric(ref),as.numeric(alt)))
        }))
        colnames(readcounts) <- c("REF", "ALT")
        
        return(readcounts)
      }
      if(info=="depth"){
        depth <- apply(vcf$vcf, 1, function(x){
          
          ref.alt <- strsplit(x[9], split=":")[[1]]
          ref.alt <- which(ref.alt=="AD")
          
          ref.alt <- strsplit(x[tumor.col.mutect], split=":")[[1]][ref.alt]
          ref <- strsplit(ref.alt, split=",")[[1]][1]
          alt <- strsplit(ref.alt, split=",")[[1]][2]
          

          return(as.numeric(ref) + as.numeric(alt))
        })
        
        return(depth)
      }
      if(info=="VAF"){
        vaf <- apply(vcf$vcf, 1, function(x){
 
          ref.alt <- strsplit(x[9], split=":")[[1]]
          ref.alt <- which(ref.alt=="AD")
          
          ref.alt <- strsplit(x[tumor.col.mutect], split=":")[[1]][ref.alt]
          ref <- strsplit(ref.alt, split=",")[[1]][1]
          alt <- strsplit(ref.alt, split=",")[[1]][2]
          
          return(as.numeric(alt)/(as.numeric(ref) + as.numeric(alt)))
        })
        return(vaf)
      }
      
      if(info=="VAF.control"){
        VAF.control <- apply(vcf$vcf, 1, function(x){
          ref.alt <- strsplit(x[9], split=":")[[1]]
          ref.alt <- which(ref.alt=="AD")
          
          ref.alt <- strsplit(x[normal.col.mutect], split=":")[[1]][ref.alt]
          ref <- strsplit(ref.alt, split=",")[[1]][1]
          alt <- strsplit(ref.alt, split=",")[[1]][2]
          
          return(as.numeric(alt)/(as.numeric(ref) + as.numeric(alt)))
        })
        return(VAF.control)
      }
      
    
  }else if(mutationcaller=="mpileup"){
    vcf$vcf <- as.data.frame(vcf$vcf)
    if(type=="snvs"){
      if(info=="readcounts"){
        readcounts <- t(sapply(vcf$vcf[,sample.col.mpileup], function(x){
          x <- strsplit(x, split=":")[[1]]
          x <- x[length(x)]
          x <- as.numeric(strsplit(x, split=",")[[1]][c(1,2)])
          c(x[c(2)], x[1])
        }))
        rownames(readcounts) <- paste(vcf$vcf$GENE, vcf$vcf$POS, sep=".")
        colnames(readcounts) <- c("REF", "ALT")
        
        return(readcounts)
      }
      if(info=="depth"){
        depth <- sapply(vcf$vcf[,sample.col.mpileup], function(x){
          x <- strsplit(x, split=":")[[1]]
          x <- x[length(x)]
          x <- as.numeric(strsplit(x, split=",")[[1]][c(1,2)])
          sum(x)
        })
        
        return(depth)
      }
      if(info=="VAF"){
        vaf <- sapply(vcf$vcf[,sample.col.mpileup], function(x){
          x <- strsplit(x, split=":")[[1]]
          x <- x[length(x)]
          x <- as.numeric(strsplit(x, split=",")[[1]][c(1,2)])
          x[2]/sum(x)
        })
        vaf <- unname(vaf)
        return(vaf)
      }
      if(info=="Gene"){
        gene <- sapply(vcf$vcf$INFO, function(x){
          x <- strsplit(x, split="Gene.refGene=")[[1]][2]
          x <- strsplit(x, split=";")[[1]][1]
        })
        return(gene)
      }
      if(info=="annovar_function"){
        annovarfunc <- sapply(vcf$vcf$INFO, function(x){
          x <- strsplit(x, split="Func.refGene=")[[1]][2]
          x <- strsplit(x, split=";")[[1]][1]
        })
        return(annovarfunc)
      }
      if(info=="exonic_function"){
        exonicfunc <- sapply(vcf$vcf$INFO, function(x){
          x <- strsplit(x, split="ExonicFunc.refGene=")[[1]][2]
          x <- strsplit(x, split=";")[[1]][1]
        })
        return(exonicfunc)
      }
      if(info=="AA_change"){
        aa_change <- sapply(vcf$vcf$INFO, function(x){
          x <- strsplit(x, split="AAChange.refGene=")[[1]][2]
          x <- strsplit(x, split=";")[[1]][1]
          x <- strsplit(x, split="p.")[[1]][2]
          x <- strsplit(x, split=",")[[1]][1]
        })
        return(aa_change)
      }
    }
  }
  
  
}