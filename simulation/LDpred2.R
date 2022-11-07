rm(list=ls())
#setwd('/dcl01/chatterj/data/jin/prs/logfile/')
setwd('/fastscratch/myscratch/jjin/')
library(data.table)
#library(GENESIS)
library(Rcpp)
library(readr)
library(MASS) # for mvrnorm
library(reshape) # for melt
library(parallel)
#library(mvtnorm) # for dmvnorm
library(RcppTN)
library(devtools)
library(Rcpp)
library(RcppArmadillo)
library(inline)
library(mvnfast)
library(genio) # for read_plink
library(dplyr)
library(stringr)
library(gdata)
library(R.utils) # for gzip
library(pROC)
library(bigsnpr)
temp <- commandArgs(TRUE)
RACES = c('EUR', 'AFR', 'AMR', 'EAS', 'SAS')
race = RACES[as.numeric(temp[1])]
chr =  as.numeric(temp[2])
rho = as.numeric(temp[3])
size = as.numeric(temp[4])
GA = as.numeric(temp[5])
rep = as.numeric(temp[6])

tem = paste0('/dcs04/nilanjan/data/jjin/prs/sim/results/ldpred2/ldpred2prs-mega-',race,'-rho=',rho,
             '-size=',size,'-chr',chr,'-rep',rep,'-GA',GA,'.txt')
if (file.exists(tem)) print('Completed; file exists.')

if (!file.exists(tem))
{
  ldr = 3/1000
  sum.raw = bigreadr::fread2(paste0('/dcl01/chatterj/data/jin/prs/simulation/',race,'/sumdata/megasum-rho',rho,'-size',size,'-rep',rep,'-GA',GA,'-chr',chr,'.txt'))
  sum.raw = as.data.frame(sum.raw)
  
  # save phenotype data into the same folder with the same file name
  iid = read.table(paste0('/dcl01/chatterj/data/jin/prs/simulation/',race,'/geno/test_size',size,'.id.txt'),header=F)
  iid=as.character(iid[,1])
  which.test = rep
  pheno = read.table(paste0('/dcl01/chatterj/data/hzhang1/multi_ethnic/LD_simulation_new/',race,'/pheno_summary_out_GA/phenotypes_rho',rho,'_',GA,'.phen'),header=F,
                     colClasses = c(rep('character',2),rep('numeric',which.test),rep('NULL',100-which.test)))
  rownames(pheno) = pheno[,1]
  pheno = pheno[iid,2+which.test]
  names(pheno) = iid
  phenomat = pheno
  # ------------------------ Run LDpred2
  temfile = paste0('/dcl01/chatterj/data/jin/prs/simulation/',race,'/geno/mega/ref_chr',chr,'.bk')
  system(paste0('rm -rf ',temfile))
  temfile = paste0('/dcl01/chatterj/data/jin/prs/simulation/',race,'/geno/mega/ref_chr',chr,'.rds')
  system(paste0('rm -rf ',temfile))
  temfile = paste0('/dcl01/chatterj/data/jin/prs/simulation/',race,'/geno/mega/ref_chr',chr,'.bk')
  if (!file.exists(temfile)){
    snp_readBed(paste0('/dcl01/chatterj/data/jin/prs/simulation/',race,'/geno/mega/ref_chr',chr,'.bed'))
  }
  obj.bigSNP <- snp_attach(paste0('/dcl01/chatterj/data/jin/prs/simulation/',race,'/geno/mega/ref_chr',chr,'.rds'))
  
  G   <- obj.bigSNP$genotypes
  CHR <- obj.bigSNP$map$chromosome
  POS <- obj.bigSNP$map$physical.pos
  y   <-  pheno # obj.bigSNP$fam$affection - 1 # outcome
  #NCORES <- nb_cores()
  # Read external summary statistics
  sumstats = sum.raw[sum.raw$CHR == chr,c('CHR', 'SNP_ID', 'POS', 'REF', 'ALT', 'BETA', 'SE', 'PVAL', 'N')]
  #str(sumstats)
  
  set.seed(2020)
  #ind.val <- sample(1:nrow(G), nrow(G))
  #ind.test <- setdiff(rows_along(G), ind.val)
  #sumstats$n_eff <- sumstats$N
  #sumstats$n_case <- sumstats$n_control <- NULL
  names(sumstats) <- c("chr", "rsid", "pos", "a0", "a1", "beta", "beta_se", "p", "n_eff")
  map <- obj.bigSNP$map[-(2:3)]
  names(map) <- c("chr", "pos", "a0", "a1")
  info_snp <- snp_match(sumstats, map, strand_flip = T)
  rownames(info_snp) = info_snp$rsid
  ## compute correlation
  POS2 <- snp_asGeneticPos(CHR, POS, dir = paste0('/dcs04/nilanjan/data/jjin/prs/sim/',race,'/sumdata/'), ncores = 2)
  ## indices in info_snp
  ind.chr <- which(info_snp$chr == chr)
  df_beta <- info_snp[ind.chr, c("beta", "beta_se", "n_eff")]
  ## indices in G
  ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
  corr0 <- snp_cor(G, ind.col = ind.chr2, ncores = 4, #size = ldradius)
                   infos.pos = POS2[ind.chr2], size = ldr) # default
  corr <- bigsparser::as_SFBM(as(corr0, "dgCMatrix"))
  
  # Automatic model
  ldsc <- snp_ldsc2(corr0, df_beta)
  h2_est <- abs(ldsc[["h2"]])
  
  # takes a few minutes
  if (T == F){
    multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est,
                                   vec_p_init = seq_log(1e-4, 0.9, length.out = NCORES),
                                   ncores = NCORES)
    str(multi_auto)
    
    beta_auto <- sapply(multi_auto, function(auto) auto$beta_est)
    pred_auto <- big_prodMat(G, beta_auto, ind.row = ind.test, ind.col = ind.chr2)
    #apply(pred_auto, 2, mad)
    final_pred_auto <- rowMeans(pred_auto)
    cor(final_pred_auto, pheno1[ind.test])
    
    # compared with infinitesimal model
    beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = h2_est)
    pred_inf <- big_prodVec(G, beta_inf, ind.row = ind.test, ind.col = ind.chr2)
    cor(pred_inf, pheno1[ind.test])
  }
  
  
  
  # grid of models:
  H2_seq = signif(abs(h2_est) * c(0.7, 1, 1.4), 3)
  h2_seq <- abs(h2_est) #signif(abs(h2_est) * c(0.7, 1, 1.4), 3)
  p_seq <- signif(seq_log(1e-4, 0.5, length.out = 15), 2)[-c(1:3)]
  #p_seq = p_seq[p_seq<0.4]
  params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE))
  
  beta_grid <- snp_ldpred2_grid(corr, df_beta, params, burn_in = 50, num_iter = 200, ncores = min(24,nrow(params)))
  beta_grid = as.data.frame(beta_grid)
  rownames(beta_grid) = info_snp$rsid
  beta_grid = cbind(info_snp$rsid, info_snp$a0, info_snp$a1, beta_grid)
  colnames(beta_grid) = c(c('marker.ID', 'a0', 'a1'),paste0('e',1:nrow(params)))
  
  
  bigreadr::fwrite2(beta_grid, paste0('/dcs04/nilanjan/data/jjin/prs/sim/results/ldpred2/ldpred2effect-mega-',race,'-rho=',rho,
                                      '-size=',size,'-chr',chr,'-rep',rep,'-GA',GA,'.txt'))
  
  temfile = paste0('/dcl01/chatterj/data/jin/prs/simulation/',race,'/geno/mega/chr',chr,'.bk')
  system(paste0('rm -rf ',temfile))
  temfile = paste0('/dcl01/chatterj/data/jin/prs/simulation/',race,'/geno/mega/chr',chr,'.rds')
  system(paste0('rm -rf ',temfile))
  temfile = paste0('/dcl01/chatterj/data/jin/prs/simulation/',race,'/geno/mega/chr',chr,'.bk')
  if (!file.exists(temfile)){
    snp_readBed(paste0('/dcl01/chatterj/data/jin/prs/simulation/',race,'/geno/mega/chr',chr,'.bed'))
  }
  obj.bigSNP <- snp_attach(paste0('/dcl01/chatterj/data/jin/prs/simulation/',race,'/geno/mega/chr',chr,'.rds'))
  combined = merge(beta_grid, obj.bigSNP$map[,c('marker.ID', 'allele1', 'allele2')], by='marker.ID')
  snps = combined$marker.ID
  flipped = ((combined$a0 == combined$allele1)&(combined$a1 == combined$allele2))
  combined[flipped,paste0('e',1:nrow(params))] = - combined[flipped,paste0('e',1:nrow(params))]
  
  #ind.rows = sapply(1:length(snps),function(x){which(obj.bigSNP$fam$sample.ID == iid[x])})
  n.training = c(15000, 45000, 80000, 100000)
  n.total = 1.2e5
  split = list()
  # extract test + validation
  for (ethnicity in RACES){
    split[[ethnicity]] = list()
    for (siz in 1:4){
      split[[ethnicity]][[siz]] = list()
      n.test = n.val = (n.total-n.training[siz])/2#n.training[siz]*0.1
      sample.test = (n.training[siz]+1):(n.training[siz]+n.val)
      sample.val = (n.training[siz]+n.val+1):(n.total)
      split[[ethnicity]][[siz]] = list(1:n.training[siz], sample.test, sample.val)
    }
  }
  ind.rows = unlist(split[[race]][[size]][2:3])
  ind.cols = sapply(1:length(snps),function(x){which(obj.bigSNP$map$marker.ID == snps[x])})
  pred_grid <- big_prodMat(obj.bigSNP$genotypes, as.matrix(combined[,paste0('e',1:nrow(params))]), ind.row = ind.rows, ind.col = ind.cols)
  pred_grid = as.data.frame(pred_grid)
  rownames(pred_grid) = ind.rows
  
  pathname = '/dcs04/nilanjan/data/jjin/prs/sim/results/ldpred2/'
  if (!dir.exists(pathname)){dir.create(pathname)}
  bigreadr::fwrite2(pred_grid, paste0('/dcs04/nilanjan/data/jjin/prs/sim/results/ldpred2/ldpred2prs-mega-',race,'-rho=',rho,
                                      '-size=',size,'-chr',chr,'-rep',rep,'-GA',GA,'.txt'))
  
  print('Complete')
}
