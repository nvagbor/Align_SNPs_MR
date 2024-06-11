
# Function to help align SNPs for Mendelian randomisation analysis 
# Function dependencies: dplyr and rlang 
# Function arguments: 
#     dataset: Dataset containing the exposure (a dataframe or tibble)
#     var_beta: Name of column containing the betas (Character vector of length 1)
#     var_effect_allele: Name of column containing the effect allele (Character vector of length 1)
#     var_other_allele: Name of column containing the non-effect allele (Character vector of length 1)
#     var_effect_allele_freq: Name of column containing the minor allele frequency (Character vector of length 1)
#     positive_beta: Logical vector that determines the sign of the new beta. 
#                    If TRUE, all negative betas will be flipped to have negative signs 
#                    OR negative sign if FALSE. 
#     The function will preserve the old variables using the suffix ".old" 


AlignSNPs <- function(dataset, 
                      var_beta, 
                      var_effect_allele, 
                      var_other_allele, 
                      var_effect_allele_freq,
                      positive_beta = TRUE  # Will make all beta positive by default
                      ){
  
  # Rename old beta.value and preserve them for verification 
  inputVars <- c(var_beta, var_effect_allele, var_other_allele, var_effect_allele_freq)
  
  for (k in inputVars) { dataset[ , paste0(k, ".old")] <- dataset[[k]] }
  
  if(positive_beta){
    old_beta <- paste0(var_beta, ".old")
    old_EA <- paste0(var_effect_allele, ".old")
    old_AA <- paste0(var_other_allele, ".old")
    
    df_temp <- dplyr::mutate(
      dataset, 
      pos_beta=positive_beta,
      !!var_beta:=case_when(.data[[old_beta]]<0~.data[[old_beta]]*-1, TRUE~.data[[old_beta]]), 
      !!var_effect_allele:=case_when(.data[[old_beta]]<0~.data[[old_AA]], TRUE~.data[[old_EA]]),
      !!var_other_allele:=case_when(.data[[old_beta]]<0~.data[[old_EA]], TRUE~.data[[old_AA]]),
      !!var_effect_allele_freq:=case_when(.data[[old_beta]]<0~1 - .data[[var_effect_allele_freq]], TRUE~.data[[var_effect_allele_freq]]),
    )
    
  }else{
    old_beta <- paste0(var_beta, ".old")
    old_EA <- paste0(var_effect_allele, ".old")
    old_AA <- paste0(var_other_allele, ".old")
    
    df_temp <- dplyr::mutate(
      dataset, 
      pos_beta=positive_beta,
      !!var_beta:=case_when(.data[[var_beta]]>0~.data[[var_beta]]*-1, TRUE~.data[[var_beta]]), 
      !!var_effect_allele:=case_when(.data[[var_beta]]>0~.data[[old_AA]], TRUE~.data[[old_EA]]),
      !!var_other_allele:=case_when(.data[[var_beta]]>0~.data[[old_EA]], TRUE~.data[[old_AA]]),
      !!var_effect_allele_freq:=case_when(.data[[var_beta]]>0~1 - .data[[var_effect_allele_freq]], TRUE~.data[[var_effect_allele_freq]]),
    )
   
  }
  
  return(df_temp)
  
}
