


require(VineCopula)
require(foreach)
require(doParallel)
require(parallel)
require(pracma)
require(lmom)
require(BSDA)
require(pls)
require(segmented)
require(abind)

##
###################lmf calculation###############################


fpall = '/GPUFS/ygo_fwf_1/00junli/01next_step/GCM_all_data_his_126_245_370_585/all_gcm_da'
get_lmf_all = function(fpall, scenario = '_ssp585_', modelnai, threshold=0.9, window_yr=40){
  t1 = proc.time()
  get_lmf_ij_sd = function(tas__bgc, mrso_bgc, i, j, threshold = 0.9, window_yr=30){#i=200;j=80
    
    sca=3
    yr=(2100-1850+1)
    tas__ij  = rowMeans(embed(c( tas__bgc[i,j,],rep(NA, sca-1)), sca), na.rm=T)
    loca_bgc = which.max(rowMeans(matrix(tas__ij, ncol=yr), na.rm=T))
    locaij = sort(embed(c(1:12,1,2),3)[loca_bgc,])
    len = yr-window_yr+1
    loca_sel = embed(1:yr, window_yr)
    
    if(loca_bgc==12){
      tas__ijk = cbind(rep(NA, 3), matrix( tas__bgc[i,j,12:3011], ncol=(yr-1))[c(1,2,3),])
      mrso_ijk = cbind(rep(NA, 3), matrix(-mrso_bgc[i,j,12:3011], ncol=(yr-1))[c(1,2,3),])
    }else if(loca_bgc==11){
      tas__ijk = cbind(rep(NA, 3), matrix( tas__bgc[i,j,11:3010], ncol=(yr-1))[c(1,2,3),])
      mrso_ijk = cbind(rep(NA, 3), matrix(-mrso_bgc[i,j,11:3010], ncol=(yr-1))[c(1,2,3),])
    }else if(loca_bgc==1){
      tas__ijk = cbind(rep(NA, 3), matrix( tas__bgc[i,j,13:3012], ncol=(yr-1))[c(1,2,3),])
      mrso_ijk = cbind(rep(NA, 3), matrix(-mrso_bgc[i,j,13:3012], ncol=(yr-1))[c(1,2,3),])
    }else{
      tas__ijk = matrix( tas__bgc[i,j,], ncol=yr)[locaij,]
      mrso_ijk = matrix(-mrso_bgc[i,j,], ncol=yr)[locaij,] }
    
    lmf_sm_mean   = c()
    lmf_ta_mean   = c()
    lmf_sm_ta_sd  = c()
    lmf_dep       = c()
    lmf_total     = c()
    sm_mean_c_raw = c()
    ta_mean_c_raw = c()
    sm_ta_sd_c_raw= c()
    cor_c_raw     = c()
    for (kk in 2:len) {#kk=202
      mrso_2 = mrso_ijk[ ,loca_sel[kk,]]
      tas_2  = tas__ijk[ ,loca_sel[kk,]]
      mrso_1 = mrso_ijk[ ,loca_sel[01,]]      
      tas_1  = tas__ijk[ ,loca_sel[01,]]      
      
      sm_tas_2 = na.omit(cbind(matrix(mrso_2, ncol=1), matrix(tas_2, ncol=1)))
      sm_tas_1 = na.omit(cbind(matrix(mrso_1, ncol=1), matrix(tas_1, ncol=1)))
      
      mrso_2 = sm_tas_2[,1]
      tas_2  = sm_tas_2[,2]
      mrso_1 = sm_tas_1[,1]
      tas_1  = sm_tas_1[,2]
      if(!is.na(mrso_2[1]) & round(sd(mrso_2, na.rm=T),2)!=0 & 
         !is.na(mrso_1[1]) & round(sd(mrso_1, na.rm=T),2)!=0 & 
         !is.na(tas_2 [1]) & round(sd(tas_2 , na.rm=T),2)!=0 & 
         !is.na(tas_1 [1]) & round(sd(tas_1 , na.rm=T),2)!=0 ){
        ######Mean sd#####
        m_bgc_2   = cbind(pobs(mrso_2), pobs(tas_2))
        m_bgc_1   = cbind(pobs(mrso_1), pobs(tas_1))
        pa_bgc_2  = BiCopEst(m_bgc_2[,1], m_bgc_2[,2], family = 5)
        pa_bgc_1  = BiCopEst(m_bgc_1[,1], m_bgc_1[,2], family = 5)
        cop_bgc_2 = BiCop(family=5, par=pa_bgc_2$par, par2=pa_bgc_2$par2)
        cop_bgc_1 = BiCop(family=5, par=pa_bgc_1$par, par2=pa_bgc_1$par2)
        
        ta__thre = cdfnor(quanor(threshold, para = c(mean(tas_1 , na.rm=T), sd(tas_1 , na.rm=T))), para = c(mean(tas_2 , na.rm=T), sd(tas_2 , na.rm=T)))
        sm__thre = cdfnor(quanor(threshold, para = c(mean(mrso_1, na.rm=T), sd(mrso_1, na.rm=T))), para = c(mean(mrso_2, na.rm=T), sd(mrso_2, na.rm=T)))
        
        lmf2_dep_sm_sd_mean = (1 - sm__thre  - threshold + BiCopCDF(sm__thre, threshold, cop_bgc_2))
        lmf2_dep_ta_sd_mean = (1 - threshold - ta__thre  + BiCopCDF(threshold, ta__thre, cop_bgc_2))
        
        ta__thre_detrend = cdfnor(quanor(threshold, para=c(mean(tas_1 , na.rm=T), sd(tas_1 , na.rm=T))), para=c(mean(tas_1 , na.rm=T), sd(tas_2 , na.rm=T)))
        sm__thre_detrend = cdfnor(quanor(threshold, para=c(mean(mrso_1, na.rm=T), sd(mrso_1, na.rm=T))), para=c(mean(mrso_1, na.rm=T), sd(mrso_2, na.rm=T)))
        
        lmf2_dep_sm_sd    = (1 - sm__thre_detrend  - threshold         + BiCopCDF(sm__thre_detrend, threshold       , cop_bgc_2))
        lmf2_dep_ta_sd    = (1 - threshold         - ta__thre_detrend  + BiCopCDF(threshold       , ta__thre_detrend, cop_bgc_2))
        lmf2_dep          = (1 - threshold         - threshold         + BiCopCDF(threshold       , threshold       , cop_bgc_2))
        lmf2_dep_sm_ta_sd = (1 - sm__thre_detrend  - ta__thre_detrend  + BiCopCDF(sm__thre_detrend, ta__thre_detrend, cop_bgc_2))
        
        lmf_sm_mean = c(lmf_sm_mean , log(lmf2_dep_sm_sd_mean/lmf2_dep_sm_sd))
        lmf_ta_mean = c(lmf_ta_mean , log(lmf2_dep_ta_sd_mean/lmf2_dep_ta_sd))
        lmf_sm_ta_sd= c(lmf_sm_ta_sd, log(lmf2_dep_sm_ta_sd/lmf2_dep         ))
        
        sm_mean_c_raw = c(sm_mean_c_raw , mean(mrso_2) - mean(mrso_1))
        ta_mean_c_raw = c(ta_mean_c_raw , mean(tas_2 ) - mean(tas_1 ))
        cor_c_raw     = c(cor_c_raw     , cov(mrso_2, tas_2) - cov(mrso_1, tas_1))
        
        sm_ta_sd_c_raw= c(sm_ta_sd_c_raw, sd(mrso_2)*sd(tas_2) - sd(mrso_1)*sd(tas_1))
        ######Dep####
        m_bgc_2   = cbind(pobs(mrso_2), pobs(tas_2))
        pa_bgc_2  = BiCopEst(m_bgc_2[,1], m_bgc_2[,2], family = 5)
        cop_bgc_2 = BiCop(family=5, par=pa_bgc_2$par, par2=pa_bgc_2$par2)
        
        m_bgc_1   = cbind(pobs(mrso_1), pobs(tas_1))
        pa_bgc_1  = BiCopEst(m_bgc_1[,1], m_bgc_1[,2], family = 5)
        cop_bgc_1 = BiCop(family=5, par=pa_bgc_1$par, par2=pa_bgc_1$par2)
        
        lmf2_dep = 1 - threshold - threshold + BiCopCDF(threshold, threshold, cop_bgc_2)
        lmf1_dep = 1 - threshold - threshold + BiCopCDF(threshold, threshold, cop_bgc_1)
        
        lmf_dep = c(lmf_dep, log(lmf2_dep/lmf1_dep))
        ######Total#####
        m_bgc_2   = cbind(pobs(mrso_2), pobs(tas_2))
        pa_bgc_2  = BiCopEst(m_bgc_2[,1], m_bgc_2[,2], family = 5)
        cop_bgc_2 = BiCop(family=5, par=pa_bgc_2$par, par2=pa_bgc_2$par2)
        
        m_bgc_1   = cbind(pobs(mrso_1), pobs(tas_1))
        pa_bgc_1  = BiCopEst(m_bgc_1[,1], m_bgc_1[,2], family = 5)
        cop_bgc_1 = BiCop(family=5, par=pa_bgc_1$par, par2=pa_bgc_1$par2)
        
        ta__thre = cdfnor(quanor(threshold, para = c(mean(tas_1 , na.rm=T), sd(tas_1 , na.rm=T))), para = c(mean(tas_2 , na.rm=T), sd(tas_2 , na.rm=T)))
        sm__thre = cdfnor(quanor(threshold, para = c(mean(mrso_1, na.rm=T), sd(mrso_1, na.rm=T))), para = c(mean(mrso_2, na.rm=T), sd(mrso_2, na.rm=T)))
        
        lmf2_dep_sm_ta_sd_mean = (1 - sm__thre  - ta__thre  + BiCopCDF(sm__thre , ta__thre , cop_bgc_2))
        lmf1_dep               = (1 - threshold - threshold + BiCopCDF(threshold, threshold, cop_bgc_1))
        
        lmf_total = c(lmf_total, log(lmf2_dep_sm_ta_sd_mean/lmf1_dep))
      }else{
        lmf_sm_mean    = c(lmf_sm_mean   , NA)
        lmf_ta_mean    = c(lmf_ta_mean   , NA)
        lmf_sm_ta_sd   = c(lmf_sm_ta_sd  , NA)
        lmf_dep        = c(lmf_dep       , NA)
        lmf_total      = c(lmf_total     , NA)
        sm_mean_c_raw  = c(sm_mean_c_raw , NA)
        ta_mean_c_raw  = c(ta_mean_c_raw , NA)
        sm_ta_sd_c_raw = c(sm_ta_sd_c_raw, NA)
        cor_c_raw      = c(cor_c_raw     , NA)
      }
    }
    
    out_da = c(lmf_sm_mean   ,lmf_ta_mean   ,lmf_sm_ta_sd  ,lmf_dep       ,lmf_total     ,
               sm_mean_c_raw ,ta_mean_c_raw ,sm_ta_sd_c_raw,cor_c_raw     ) 
    return(out_da)
  }
  
  lucc = readRDS('/GPUFS/ygo_fwf_1/00junli/01next_step/lucc2016_1.5degree_vegetated_area.rds')
  #fphis=fpall; modelnai=modelna[01]; window_yr=30;scenario='_ssp585_';threshold = 0.9
  setwd(fpall)
  fp_file = list.files(fpall)
  for (k in 1:length(fp_file)) {#k=10
    con_tas     = grepl(c('tas'  ), fp_file[k])
    con_mrso_   = grepl(c('mrso_'), fp_file[k])
    con_his     = grepl(c('_historical_'), fp_file[k])
    con_ssp585  = grepl(c(scenario      ), fp_file[k])
    con_modelna = grepl(c(modelnai      ), fp_file[k])
    
    if(con_tas   & con_his & con_modelna) {tas__his_bgcna = fp_file[k]}
    if(con_mrso_ & con_his & con_modelna) {mrso_his_bgcna = fp_file[k]}
    
    if(con_tas   & con_ssp585 & con_modelna) {tas__ssp_bgcna = fp_file[k]}
    if(con_mrso_ & con_ssp585 & con_modelna) {mrso_ssp_bgcna = fp_file[k]}
  }
  tas__his = readRDS(tas__his_bgcna)
  mrso_his = readRDS(mrso_his_bgcna)
  tas__ssp = readRDS(tas__ssp_bgcna)
  mrso_ssp = readRDS(mrso_ssp_bgcna)
  
  #if(modelnai == "TaiESM1" & dim(tas__his)[3] != 1980){ tas__his = abind(matrix(NA, nrow=240, ncol=93), tas__his)}
  #if(modelnai == "TaiESM1" & dim(mrso_his)[3] != 1980){ mrso_his = abind(matrix(NA, nrow=240, ncol=93), mrso_his)}
  #if(modelnai == "TaiESM1" & dim(tas__ssp)[3] != 1032){ tas__ssp = abind(matrix(NA, nrow=240, ncol=93), tas__ssp)}
  #if(modelnai == "TaiESM1" & dim(mrso_ssp)[3] != 1032){ mrso_ssp = abind(matrix(NA, nrow=240, ncol=93), mrso_ssp)}
  
  tas__bgc = abind(tas__his, tas__ssp); rm(tas__his); rm(tas__ssp)
  mrso_bgc = abind(mrso_his, mrso_ssp); rm(mrso_his); rm(mrso_ssp)
  gc()
  #print('data in')
  
  yr = 2100-1850+1
  len = yr - window_yr + 1
  da_len = (len - 1)*9
  core=36
  cl = makeCluster(core)
  registerDoParallel(cl)
  lmf_bgc1 = foreach (k1 = 1:dim(tas__bgc)[1], .packages=c('VineCopula','lmom','pracma','abind'))  %:% 
    foreach (k2 = 1:dim(tas__bgc)[2],.combine='rbind') %dopar% {  #k1=50;k2=80; 
      if(is.na(lucc[k1,k2])         |
         is.na(tas__bgc[k1,k2,100]) | 
         is.na(mrso_bgc[k1,k2,100])){
        out_i = rep(NA, da_len)        
      }else{ 
        out_i = get_lmf_ij_sd(tas__bgc, mrso_bgc, i=k1, j=k2, threshold, window_yr) 
      }
      out_i
    }
  stopCluster(cl)
  
  #plot(out_i[((len - 1)*(06-1)+1):((len - 1)*06)], out_i[((len - 1)*(1-1)+1):((len - 1)*1)])
  #plot(out_i[((len - 1)*(07-1)+1):((len - 1)*07)], out_i[((len - 1)*(2-1)+1):((len - 1)*2)])
  #plot(out_i[((len - 1)*(08-1)+1):((len - 1)*08)], out_i[((len - 1)*(3-1)+1):((len - 1)*3)])
  #plot(out_i[((len - 1)*(09-1)+1):((len - 1)*09)], out_i[((len - 1)*(4-1)+1):((len - 1)*4)])
  #plot(out_i[((len - 1)*(1-1)+1):((len - 1)*1)]+
  #     out_i[((len - 1)*(2-1)+1):((len - 1)*2)]+
  #     out_i[((len - 1)*(3-1)+1):((len - 1)*3)]+
  #     out_i[((len - 1)*(4-1)+1):((len - 1)*4)],
  #     out_i[((len - 1)*(5-1)+1):((len - 1)*5)])
  
  lmf_bgc_timing = array(NA, c(dim(tas__bgc)[1], dim(tas__bgc)[2], da_len))
  for (k3 in 1:length(lmf_bgc1)) {#k3=1
    lmf_bgc_timing[k3,,] = lmf_bgc1[[k3]]
  }
  
  out_na = 'lmf-sm-ta-mean-sm-ta-var-sm-ta-dep-total_'
  setwd('/GPUFS/ygo_fwf_1/00junli/01next_step/7lmf_all_and_ensemble')
  saveRDS(lmf_bgc_timing, file = paste(out_na, scenario, modelnai, '_', threshold, '_', window_yr, sep=''))
} #contain raw values

#for (kk in 1:23) {
#  tryCatch({get_lmf_all(fpall, '_ssp585_', modelna[kk], 0.90, 30);print(kk)}, error=function(e){print(paste('Error:', kk))})}
#for (kk in 1:23) {
#  tryCatch({get_lmf_all(fpall, '_ssp126_', modelna[kk], 0.90, 30);print(kk)}, error=function(e){print(paste('Error:', kk))})}
#for (kk in 1:23) {
#  tryCatch({get_lmf_all(fpall, '_ssp245_', modelna[kk], 0.90, 30);print(kk)}, error=function(e){print(paste('Error:', kk))})}
#for (kk in 1:23) {
#  tryCatch({get_lmf_all(fpall, '_ssp370_', modelna[kk], 0.90, 30);print(kk)}, error=function(e){print(paste('Error:', kk))})}
#
#
#for (kk in 1:23) {
#  tryCatch({get_lmf_all(fpall, '_ssp585_', modelna[kk], 0.90, 20);print(kk)}, error=function(e){print(paste('Error:', kk))})}
#for (kk in 1:23) {
#  tryCatch({get_lmf_all(fpall, '_ssp585_', modelna[kk], 0.90, 40);print(kk)}, error=function(e){print(paste('Error:', kk))})}
#for (kk in 1:23) {
#  tryCatch({get_lmf_all(fpall, '_ssp585_', modelna[kk], 0.95, 30);print(kk)}, error=function(e){print(paste('Error:', kk))})}
#for (kk in 1:23) {
#  tryCatch({get_lmf_all(fpall, '_ssp585_', modelna[kk], 0.80, 30);print(kk)}, error=function(e){print(paste('Error:', kk))})}
#for (kk in 1:23) {
#  tryCatch({get_lmf_all(fpall, '_ssp585_', modelna[kk], 0.85, 30);print(kk)}, error=function(e){print(paste('Error:', kk))})}

get_lmf_all = function(modelnai, experiment){
  threshold=0.9; window_yr=30
  get_lmf_ij_sd = function(tas__bgc, mrso_bgc, i, j, threshold = 0.9, window_yr=30){#i=200;j=80
    
    sca=3
    yr=(2100-1850+1)
    tas__ij  = rowMeans(embed(c( tas__bgc[i,j,],rep(NA, sca-1)), sca), na.rm=T)
    loca_bgc = which.max(rowMeans(matrix(tas__ij, ncol=yr), na.rm=T))
    locaij = sort(embed(c(1:12,1,2),3)[loca_bgc,])
    len = yr-window_yr+1
    loca_sel = embed(1:yr, window_yr)
    
    if(loca_bgc==12){
      tas__ijk = cbind(rep(NA, 3), matrix( tas__bgc[i,j,12:3011], ncol=(yr-1))[c(1,2,3),])
      mrso_ijk = cbind(rep(NA, 3), matrix(-mrso_bgc[i,j,12:3011], ncol=(yr-1))[c(1,2,3),])
    }else if(loca_bgc==11){
      tas__ijk = cbind(rep(NA, 3), matrix( tas__bgc[i,j,11:3010], ncol=(yr-1))[c(1,2,3),])
      mrso_ijk = cbind(rep(NA, 3), matrix(-mrso_bgc[i,j,11:3010], ncol=(yr-1))[c(1,2,3),])
    }else if(loca_bgc==1){
      tas__ijk = cbind(rep(NA, 3), matrix( tas__bgc[i,j,13:3012], ncol=(yr-1))[c(1,2,3),])
      mrso_ijk = cbind(rep(NA, 3), matrix(-mrso_bgc[i,j,13:3012], ncol=(yr-1))[c(1,2,3),])
    }else{
      tas__ijk = matrix( tas__bgc[i,j,], ncol=yr)[locaij,]
      mrso_ijk = matrix(-mrso_bgc[i,j,], ncol=yr)[locaij,] }
    
    lmf_sm_mean   = c()
    lmf_ta_mean   = c()
    lmf_sm_ta_sd  = c()
    lmf_dep       = c()
    lmf_total     = c()
    sm_mean_c_raw = c()
    ta_mean_c_raw = c()
    sm_ta_sd_c_raw= c()
    cor_c_raw     = c()
    for (kk in 2:len) {#kk=202
      mrso_2 = mrso_ijk[ ,loca_sel[kk,]]
      tas_2  = tas__ijk[ ,loca_sel[kk,]]
      mrso_1 = mrso_ijk[ ,loca_sel[01,]]      
      tas_1  = tas__ijk[ ,loca_sel[01,]]      
      
      sm_tas_2 = na.omit(cbind(matrix(mrso_2, ncol=1), matrix(tas_2, ncol=1)))
      sm_tas_1 = na.omit(cbind(matrix(mrso_1, ncol=1), matrix(tas_1, ncol=1)))
      
      mrso_2 = sm_tas_2[,1]
      tas_2  = sm_tas_2[,2]
      mrso_1 = sm_tas_1[,1]
      tas_1  = sm_tas_1[,2]
      if(!is.na(mrso_2[1]) & round(sd(mrso_2, na.rm=T),2)!=0 & 
         !is.na(mrso_1[1]) & round(sd(mrso_1, na.rm=T),2)!=0 & 
         !is.na(tas_2 [1]) & round(sd(tas_2 , na.rm=T),2)!=0 & 
         !is.na(tas_1 [1]) & round(sd(tas_1 , na.rm=T),2)!=0 ){
        ######Mean sd#####
        m_bgc_2   = cbind(pobs(mrso_2), pobs(tas_2))
        m_bgc_1   = cbind(pobs(mrso_1), pobs(tas_1))
        pa_bgc_2  = BiCopEst(m_bgc_2[,1], m_bgc_2[,2], family = 5)
        pa_bgc_1  = BiCopEst(m_bgc_1[,1], m_bgc_1[,2], family = 5)
        cop_bgc_2 = BiCop(family=5, par=pa_bgc_2$par, par2=pa_bgc_2$par2)
        cop_bgc_1 = BiCop(family=5, par=pa_bgc_1$par, par2=pa_bgc_1$par2)
        
        ta__thre = cdfnor(quanor(threshold, para = c(mean(tas_1 , na.rm=T), sd(tas_1 , na.rm=T))), para = c(mean(tas_2 , na.rm=T), sd(tas_2 , na.rm=T)))
        sm__thre = cdfnor(quanor(threshold, para = c(mean(mrso_1, na.rm=T), sd(mrso_1, na.rm=T))), para = c(mean(mrso_2, na.rm=T), sd(mrso_2, na.rm=T)))
        
        lmf2_dep_sm_sd_mean = (1 - sm__thre  - threshold + BiCopCDF(sm__thre, threshold, cop_bgc_2))
        lmf2_dep_ta_sd_mean = (1 - threshold - ta__thre  + BiCopCDF(threshold, ta__thre, cop_bgc_2))
        
        ta__thre_detrend = cdfnor(quanor(threshold, para=c(mean(tas_1 , na.rm=T), sd(tas_1 , na.rm=T))), para=c(mean(tas_1 , na.rm=T), sd(tas_2 , na.rm=T)))
        sm__thre_detrend = cdfnor(quanor(threshold, para=c(mean(mrso_1, na.rm=T), sd(mrso_1, na.rm=T))), para=c(mean(mrso_1, na.rm=T), sd(mrso_2, na.rm=T)))
        
        lmf2_dep_sm_sd    = (1 - sm__thre_detrend  - threshold         + BiCopCDF(sm__thre_detrend, threshold       , cop_bgc_2))
        lmf2_dep_ta_sd    = (1 - threshold         - ta__thre_detrend  + BiCopCDF(threshold       , ta__thre_detrend, cop_bgc_2))
        lmf2_dep          = (1 - threshold         - threshold         + BiCopCDF(threshold       , threshold       , cop_bgc_2))
        lmf2_dep_sm_ta_sd = (1 - sm__thre_detrend  - ta__thre_detrend  + BiCopCDF(sm__thre_detrend, ta__thre_detrend, cop_bgc_2))
        
        lmf_sm_mean = c(lmf_sm_mean , log(lmf2_dep_sm_sd_mean/lmf2_dep_sm_sd))
        lmf_ta_mean = c(lmf_ta_mean , log(lmf2_dep_ta_sd_mean/lmf2_dep_ta_sd))
        lmf_sm_ta_sd= c(lmf_sm_ta_sd, log(lmf2_dep_sm_ta_sd/lmf2_dep        ))
        
        sm_mean_c_raw = c(sm_mean_c_raw , mean(mrso_2) - mean(mrso_1))
        ta_mean_c_raw = c(ta_mean_c_raw , mean(tas_2 ) - mean(tas_1 ))
        cor_c_raw     = c(cor_c_raw     , cov(mrso_2, tas_2) - cov(mrso_1, tas_1))
        
        sm_ta_sd_c_raw= c(sm_ta_sd_c_raw, sd(mrso_2)*sd(tas_2) - sd(mrso_1)*sd(tas_1))
        ######Dep####
        m_bgc_2   = cbind(pobs(mrso_2), pobs(tas_2))
        pa_bgc_2  = BiCopEst(m_bgc_2[,1], m_bgc_2[,2], family = 5)
        cop_bgc_2 = BiCop(family=5, par=pa_bgc_2$par, par2=pa_bgc_2$par2)
        
        m_bgc_1   = cbind(pobs(mrso_1), pobs(tas_1))
        pa_bgc_1  = BiCopEst(m_bgc_1[,1], m_bgc_1[,2], family = 5)
        cop_bgc_1 = BiCop(family=5, par=pa_bgc_1$par, par2=pa_bgc_1$par2)
        
        lmf2_dep = 1 - threshold - threshold + BiCopCDF(threshold, threshold, cop_bgc_2)
        lmf1_dep = 1 - threshold - threshold + BiCopCDF(threshold, threshold, cop_bgc_1)
        
        lmf_dep = c(lmf_dep, log(lmf2_dep/lmf1_dep))
        ######Total#####
        m_bgc_2   = cbind(pobs(mrso_2), pobs(tas_2))
        pa_bgc_2  = BiCopEst(m_bgc_2[,1], m_bgc_2[,2], family = 5)
        cop_bgc_2 = BiCop(family=5, par=pa_bgc_2$par, par2=pa_bgc_2$par2)
        
        m_bgc_1   = cbind(pobs(mrso_1), pobs(tas_1))
        pa_bgc_1  = BiCopEst(m_bgc_1[,1], m_bgc_1[,2], family = 5)
        cop_bgc_1 = BiCop(family=5, par=pa_bgc_1$par, par2=pa_bgc_1$par2)
        
        ta__thre = cdfnor(quanor(threshold, para = c(mean(tas_1 , na.rm=T), sd(tas_1 , na.rm=T))), para = c(mean(tas_2 , na.rm=T), sd(tas_2 , na.rm=T)))
        sm__thre = cdfnor(quanor(threshold, para = c(mean(mrso_1, na.rm=T), sd(mrso_1, na.rm=T))), para = c(mean(mrso_2, na.rm=T), sd(mrso_2, na.rm=T)))
        
        lmf2_dep_sm_ta_sd_mean = (1 - sm__thre  - ta__thre  + BiCopCDF(sm__thre , ta__thre , cop_bgc_2))
        lmf1_dep               = (1 - threshold - threshold + BiCopCDF(threshold, threshold, cop_bgc_1))
        
        lmf_total = c(lmf_total, log(lmf2_dep_sm_ta_sd_mean/lmf1_dep))
      }else{
        lmf_sm_mean    = c(lmf_sm_mean   , NA)
        lmf_ta_mean    = c(lmf_ta_mean   , NA)
        lmf_sm_ta_sd   = c(lmf_sm_ta_sd  , NA)
        lmf_dep        = c(lmf_dep       , NA)
        lmf_total      = c(lmf_total     , NA)
        sm_mean_c_raw  = c(sm_mean_c_raw , NA)
        ta_mean_c_raw  = c(ta_mean_c_raw , NA)
        sm_ta_sd_c_raw = c(sm_ta_sd_c_raw, NA)
        cor_c_raw      = c(cor_c_raw     , NA)
      }
    }
    
    out_da = c(lmf_sm_mean   ,lmf_ta_mean   ,lmf_sm_ta_sd  ,lmf_dep       ,lmf_total     ,
               sm_mean_c_raw ,ta_mean_c_raw ,sm_ta_sd_c_raw,cor_c_raw     ) 
    return(out_da)
  }
  
  lucc = readRDS('/GPUFS/ygo_fwf_1/00junli/01next_step/lucc2016_1.5degree_vegetated_area.rds')
  setwd('/GPUFS/ygo_fwf_1/00junli/01next_step/GCM_ensemble')
  fp_file = list.files('/GPUFS/ygo_fwf_1/00junli/01next_step/GCM_ensemble')
  for (k in 1:length(fp_file)) {#k=10
    con_tas     = grepl(c('tas'         ), fp_file[k])
    con_mrso_   = grepl(c('mrso_'       ), fp_file[k])
    con_his     = grepl(c('_historical_'), fp_file[k])
    con_ssp585  = grepl(c('_ssp585_'    ), fp_file[k])
    con_modelna = grepl(c(modelnai      ), fp_file[k])
    con_expna   = grepl(c(experiment    ), fp_file[k])
    
    if(con_tas   & con_his    & con_modelna & con_expna) {tas__his_bgcna = fp_file[k]}
    if(con_mrso_ & con_his    & con_modelna & con_expna) {mrso_his_bgcna = fp_file[k]}
    if(con_tas   & con_ssp585 & con_modelna & con_expna) {tas__ssp_bgcna = fp_file[k]}
    if(con_mrso_ & con_ssp585 & con_modelna & con_expna) {mrso_ssp_bgcna = fp_file[k]}
  }
  tas__his = readRDS(tas__his_bgcna)
  mrso_his = readRDS(mrso_his_bgcna)
  tas__ssp = readRDS(tas__ssp_bgcna)
  mrso_ssp = readRDS(mrso_ssp_bgcna)
  tas__bgc = abind(tas__his, tas__ssp); rm(tas__his); rm(tas__ssp)
  mrso_bgc = abind(mrso_his, mrso_ssp); rm(mrso_his); rm(mrso_ssp)
  gc()
  
  yr = 2100-1850+1
  len = yr - window_yr + 1
  da_len = (len - 1)*9
  core=36
  cl = makeCluster(core)
  registerDoParallel(cl)
  lmf_bgc1 = foreach (k1 = 1:dim(tas__bgc)[1], .packages=c('VineCopula','lmom','pracma','abind'))  %:% 
    foreach (k2 = 1:dim(tas__bgc)[2],.combine='rbind') %dopar% {  #k1=200;k2=80; 
      if(is.na(lucc[k1,k2])         |
         is.na(tas__bgc[k1,k2,100]) | 
         is.na(mrso_bgc[k1,k2,100])){
        out_i = rep(NA, da_len)        
      }else{ 
        out_i = get_lmf_ij_sd(tas__bgc, mrso_bgc, i=k1, j=k2, threshold=threshold, window_yr=window_yr) 
      }
      out_i }
  stopCluster(cl)
  
  #print('data fi')
  #plot(out_i[((len - 1)*(06-1)+1):((len - 1)*06)], out_i[((len - 1)*(1-1)+1):((len - 1)*1)])
  #plot(out_i[((len - 1)*(07-1)+1):((len - 1)*07)], out_i[((len - 1)*(2-1)+1):((len - 1)*2)])
  #plot(out_i[((len - 1)*(08-1)+1):((len - 1)*08)], out_i[((len - 1)*(3-1)+1):((len - 1)*3)])
  #plot(out_i[((len - 1)*(09-1)+1):((len - 1)*09)], out_i[((len - 1)*(4-1)+1):((len - 1)*4)])
  lmf_bgc_timing = array(NA, c(dim(tas__bgc)[1], dim(tas__bgc)[2], da_len))
  for (k3 in 1:length(lmf_bgc1)) {#k3=1
    lmf_bgc_timing[k3,,] = lmf_bgc1[[k3]]
  }
  
  scenario = '_ssp585_'
  out_na = 'lmf_ensemble' 
  setwd('/GPUFS/ygo_fwf_1/00junli/01next_step/7lmf_all_and_ensemble')
  saveRDS(lmf_bgc_timing, file = paste(out_na, scenario, modelnai, '_', experiment, '_', threshold, '_', window_yr, sep=''))
} #contain raw values

#for (kk in 1:50) {#kk=1
#  tryCatch({get_lmf_all(modelna[1], canesm5[kk]);print(kk)}, error=function(e){print(paste('Error:', kk))})}
#
#for (kk in 1:20) {#kk=1
#  tryCatch({get_lmf_all(modelna[2], canesm51[kk]);print(kk)}, error=function(e){print(paste('Error:', kk))})}
#
#for (kk in 1:20) {#kk=1
#  tryCatch({get_lmf_all(modelna[3], mpi[kk]);print(kk)}, error=function(e){print(paste('Error:', kk))})}

###################tran-induced change##########################################################

decompose_own = function(x, type = c("additive", "multiplicative"), filter = NULL, window_yr){
  #x = x_mrso_ij
  type <- match.arg(type)
  
  if (anyNA(x)) {
    x <- na.interp(x, lambda = NULL) }
  l <- length(x)
  f <- frequency(x)
  if (f <= 1 || length(na.omit(x)) < 2 * f) 
    stop("time series has no or less than 2 periods")
  if (is.null(filter)) 
    filter <- if (!f%%2) 
      c(0.5, rep_len(1, f - 1), 0.5)/f
  else rep_len(1, f)/f
  trend <- filter(x, filter)
  season <- if (type == "additive") 
    x - trend
  else x/trend
  
  periods <- l%/%f
  index <- (seq.int(1L, l, by = f) - 1L)[1:window_yr]
  #index <- seq.int(1L, 349, by = f) - 1L
  figure <- numeric(f)
  for (i in 1L:f) figure[i] <- mean(season[index + i], na.rm = TRUE)
  figure <- if (type == "additive") 
    figure - mean(figure)
  else figure/mean(figure)
  seasonal <- ts(rep(figure, periods + 1)[seq_len(l)], start = start(x), 
                 frequency = f)
  structure(list(x = x, seasonal = seasonal, trend = trend, 
                 random = if (type == "additive") x - seasonal - trend else x/seasonal/trend, 
                 figure = figure, type = type), class = "decomposed.ts")
}
cal_tran_induced_change = function(fpall, modelna, scenario, core1, window_yr){
  
  
  get_re_sm_mean_g = function(fpall, modelnai, scenario, window_yr=window_yr){
    
    #window_yr=30; modelnai = modelna[01]; scenario='_ssp585_'
    setwd(fpall)
    fp_file = list.files(fpall)
    for (k in 1:length(fp_file)) {#k=10
      con_tas    = grepl(c('tas'         ), fp_file[k])
      con_mrso_  = grepl(c('mrso_'       ), fp_file[k])
      con_mrro_  = grepl(c('mrro_'       ), fp_file[k])
      con__lai   = grepl(c('lai_L'       ), fp_file[k])
      con_et     = grepl(c('evspsbl_A'   ), fp_file[k])#
      con_tran   = grepl(c('tran_L'      ), fp_file[k])#evspsbl_A
      con_rsds   = grepl(c('rsds'        ), fp_file[k])
      con_rsus   = grepl(c('rsus'        ), fp_file[k])
      con_ts     = grepl(c('ts_A'        ), fp_file[k])
      con_hfss   = grepl(c('hfss_A'      ), fp_file[k])
      con_wind   = grepl(c('sfcWind_A'   ), fp_file[k])
      con_pr     = grepl(c('pr_A'        ), fp_file[k])
      con_his    = grepl(c('_historical_'), fp_file[k])
      con_ssp585 = grepl(c(scenario      ), fp_file[k])
      con_modelna= grepl(c(modelnai      ), fp_file[k])
      
      if(con_tas   & con_his & con_modelna) {tas__his_bgcna = fp_file[k]}
      if(con_mrso_ & con_his & con_modelna) {mrso_his_bgcna = fp_file[k]}
      if(con_mrro_ & con_his & con_modelna) {mrro_his_bgcna = fp_file[k]}
      if(con__lai  & con_his & con_modelna) {lai__his_bgcna = fp_file[k]}
      if(con_et    & con_his & con_modelna) {et_his_bgcna   = fp_file[k]}
      if(con_tran  & con_his & con_modelna) {tranhis_bgcna  = fp_file[k]}
      if(con_rsds  & con_his & con_modelna) {rsds_his_bgcna = fp_file[k]}
      if(con_rsus  & con_his & con_modelna) {rsus_his_bgcna = fp_file[k]}
      if(con_ts    & con_his & con_modelna) {ts_his_bgcna   = fp_file[k]}
      if(con_hfss  & con_his & con_modelna) {hfss_his_bgcna = fp_file[k]}
      if(con_wind  & con_his & con_modelna) {wind_his_bgcna = fp_file[k]}
      if(con_pr    & con_his & con_modelna) {pr_his_bgcna   = fp_file[k]}
      
      if(con_tas   & con_ssp585 & con_modelna) {tas__ssp_bgcna = fp_file[k]}
      if(con_mrso_ & con_ssp585 & con_modelna) {mrso_ssp_bgcna = fp_file[k]}
      if(con_mrro_ & con_ssp585 & con_modelna) {mrro_ssp_bgcna = fp_file[k]}
      if(con__lai  & con_ssp585 & con_modelna) {lai__ssp_bgcna = fp_file[k]}
      if(con_et    & con_ssp585 & con_modelna) {et_ssp_bgcna   = fp_file[k]}
      if(con_tran  & con_ssp585 & con_modelna) {transsp_bgcna  = fp_file[k]}
      if(con_rsds  & con_ssp585 & con_modelna) {rsds_ssp_bgcna = fp_file[k]}
      if(con_rsus  & con_ssp585 & con_modelna) {rsus_ssp_bgcna = fp_file[k]}
      if(con_ts    & con_ssp585 & con_modelna) {ts_ssp_bgcna   = fp_file[k]}
      if(con_hfss  & con_ssp585 & con_modelna) {hfss_ssp_bgcna = fp_file[k]}
      if(con_wind  & con_ssp585 & con_modelna) {wind_ssp_bgcna = fp_file[k]}
      if(con_pr    & con_ssp585 & con_modelna) {pr_ssp_bgcna   = fp_file[k]}
    }
    tas__his = readRDS(tas__his_bgcna)
    tran_his = readRDS(tranhis_bgcna )
    mrso_his = readRDS(mrso_his_bgcna)
    mrro_his = readRDS(mrro_his_bgcna)
    rsds_his = readRDS(rsds_his_bgcna)
    ts___his = readRDS(ts_his_bgcna  )
    pr___his = readRDS(pr_his_bgcna  )#rsus_his = readRDS(rsus_his_bgcna)
   #wind_his = readRDS(wind_his_bgcna)#lai__his = readRDS(lai__his_bgcna)
    
    tas__ssp = readRDS(tas__ssp_bgcna)
    tran_ssp = readRDS(transsp_bgcna )
    mrso_ssp = readRDS(mrso_ssp_bgcna)
    mrro_ssp = readRDS(mrro_ssp_bgcna)
    rsds_ssp = readRDS(rsds_ssp_bgcna)
    ts___ssp = readRDS(ts_ssp_bgcna  )
    pr___ssp = readRDS(pr_ssp_bgcna  )#lai__ssp = readRDS(lai__ssp_bgcna)
   #wind_ssp = readRDS(wind_ssp_bgcna)#rsus_ssp = readRDS(rsus_ssp_bgcna)
    
    tas__bgc = abind(tas__his, tas__ssp); rm(tas__his); rm(tas__ssp); print(dim(tas__bgc))
    mrso_bgc = abind(mrso_his, mrso_ssp); rm(mrso_his); rm(mrso_ssp); print(dim(mrso_bgc))
    mrro_bgc = abind(mrro_his, mrro_ssp); rm(mrro_his); rm(mrro_ssp); print(dim(mrro_bgc))
    tran_bgc = abind(tran_his, tran_ssp); rm(tran_his); rm(tran_ssp); print(dim(tran_bgc))
    rsds_bgc = abind(rsds_his, rsds_ssp); rm(rsds_his); rm(rsds_ssp); print(dim(rsds_bgc))
    ts___bgc = abind(ts___his, ts___ssp); rm(ts___his); rm(ts___ssp); print(dim(ts___bgc))
    pr___bgc = abind(pr___his, pr___ssp); rm(pr___his); rm(pr___ssp); print(dim(pr___bgc))#rsus_bgc = abind(rsus_his, rsus_ssp); rm(rsus_his); rm(rsus_ssp); print(dim(rsus_bgc))
   #wind_bgc = abind(wind_his, wind_ssp); rm(wind_his); rm(wind_ssp); print(dim(wind_bgc))#lai__bgc = abind(lai__his, lai__ssp); rm(lai__his); rm(lai__ssp); print(dim(lai__bgc))
    
    rm(wind_bgc); gc()
    
    sca=3
    yr=2100-1850+1
    len = yr-window_yr+1
    nrow1=240
    ncol1=93
    da_tran_sm         = array(NA, c(nrow1,ncol1,3  ))     
    lai_tas_mean1_c    = array(NA, c(nrow1,ncol1,len))     
    lai_tas_sd1_c      = array(NA, c(nrow1,ncol1,len))   
    lai_sm_mean0_c     = array(NA, c(nrow1,ncol1,len))    
    lai_sm_mean1_c     = array(NA, c(nrow1,ncol1,len))    
    lai_sm_sd0_c       = array(NA, c(nrow1,ncol1,len))  
    lai_sm_sd1_c       = array(NA, c(nrow1,ncol1,len))  
    lai_cov_c0         = array(NA, c(nrow1,ncol1,len))
    lai_cov_c1         = array(NA, c(nrow1,ncol1,len))
    cov_tran_c01       = array(NA, c(nrow1,ncol1,len))
    lai_sm_lai_tas_cor = matrix(NA, nrow=nrow1, ncol=ncol1)
    lucc = readRDS('/GPUFS/ygo_fwf_1/00junli/01next_step/lucc2016_1.5degree_vegetated_area.rds');col=1
    for (i in 1:nrow(tas__bgc)) {
      for (j in 1:ncol(tas__bgc)) {
        if(is.na(lucc[i,j]) | is.na(tas__bgc[i,j,1]) | 
           is.na(mrso_bgc[i,j,1]) | is.na(pr___bgc[i,j,1]) | is.na(tran_bgc[i,j,1])) next
        #i=200;j=80;image.plot(1:240,1:93,tas__bgc[,,8]);abline(v=i,h=j); gc()
        tas__ij  = rowMeans(embed(c(tas__bgc[i,j,],rep(NA, sca-1)), sca), na.rm=T)
        loca_bgc = which.max(rowMeans(matrix(tas__ij, ncol=yr), na.rm=T))
        locaij = sort(embed(c(1:12,1,2),3)[loca_bgc,])
        loca_sel = embed(1:yr, window_yr)
        
        x_mrso_ij     = ts(mrso_bgc[i,j,]         , frequency=12)
        x_mrro_ij     = ts(mrro_bgc[i,j,]*86400*30, frequency=12)
        x_tran_ij_sm  = ts(tran_bgc[i,j,]*86400*30, frequency=12)
        x_pr___ij_sm  = ts(pr___bgc[i,j,]*86400*30, frequency=12)
        x_tran_ij_tas = ts(tran_bgc[i,j,]         , frequency=12)
        
        mrso_ij    = try(decompose_own(x_mrso_ij    , type = "additive", window_yr = window_yr), silent = T)
        mrro_ij    = try(decompose_own(x_mrro_ij    , type = "additive", window_yr = window_yr), silent = T)
        tran_ij_sm = try(decompose_own(x_tran_ij_sm , type = "additive", window_yr = window_yr), silent = T)
        pr___ij_sm = try(decompose_own(x_pr___ij_sm , type = "additive", window_yr = window_yr), silent = T)
        tran_ij_tas= try(decompose_own(x_tran_ij_tas, type = "additive", window_yr = window_yr), silent = T)
        
        tranij_raw  = tran_bgc[i,j,]
        tasij_raw   = tas__bgc[i,j,]
        #ra_ij_raw   = ra___bgc[i,j,]
        #ts_ij_raw   = ts___bgc[i,j,]
        #sw_ij_raw   = rsds_bgc[i,j,]
        if(class(mrso_ij    )[1] != "try-error" &
           class(mrro_ij    )[1] != "try-error" &
           class(tran_ij_sm )[1] != "try-error" &
           class(pr___ij_sm )[1] != "try-error" &
           class(tran_ij_tas)[1] != "try-error"){
          
          getre = function(mrso_ij  ,mrro_ij, tran_ij_sm,pr___ij_sm,tran_ij_tas, tranij_raw,tasij_raw,id1=09:3008, id2=c(1:6), yr=251){
            return(list(mrso_ijk     = matrix(mrso_ij    [id1], ncol=yr)[id2,],
                        mrro_ijk     = matrix(mrro_ij    [id1], ncol=yr)[id2,],
                        tran_ijk_sm  = matrix(tran_ij_sm [id1], ncol=yr)[id2,],
                        pr___ijk_sm  = matrix(pr___ij_sm [id1], ncol=yr)[id2,],
                        tran_ijk_tas = matrix(tran_ij_tas[id1], ncol=yr)[id2,],
                        tranij_raw   = matrix(tranij_raw [id1], ncol=yr)[id2,],
                        tasij_raw    = matrix(tasij_raw  [id1], ncol=yr)[id2,])) 
            }
          
          if(loca_bgc==12){
            id1 = 09:3008; id2 = c(1:6); yr1 = 2100-1850 }else if(loca_bgc==11){
              id1 = 08:3007; id2 = c(1:6); yr1 = 2100-1850 }else if(loca_bgc==1){
                id1 = 10:3009; id2 = c(1:6); yr1 = 2100-1850 }else if(loca_bgc==2){
                  id1 = 11:3010; id2 = c(1:6); yr1 = 2100-1850 }else if(loca_bgc==3){
                    id1 = 12:3011; id2 = c(1:6); yr1 = 2100-1850 }else{
                      id1 = 1:3012; id2 = c((min(locaij)-3):max(locaij)); yr1 = 2100-1850+1
                    }
          
          mrso_ij_1     = mrso_ij    $random# - mrso_ij    [,2] - rep(rowMeans(matrix(mrso_ij    [,3], ncol=(2100-1850+1)), na.rm=T), (2100-1850+1))
          mrro_ij_1     = mrro_ij    $random# - mrso_ij    [,2] - rep(rowMeans(matrix(mrso_ij    [,3], ncol=(2100-1850+1)), na.rm=T), (2100-1850+1))
          tran_ij_sm_1  = tran_ij_sm $random# - tran_ij_sm [,2] - rep(rowMeans(matrix(tran_ij_sm [,3], ncol=(2100-1850+1)), na.rm=T), (2100-1850+1))
          pr___ij_sm_1  = pr___ij_sm $random# - pr___ij_sm [,2] - rep(rowMeans(matrix(pr___ij_sm [,3], ncol=(2100-1850+1)), na.rm=T), (2100-1850+1))
          tran_ij_tas_1 = tran_ij_tas$random# - tran_ij_tas[,2] - rep(rowMeans(matrix(tran_ij_tas[,3], ncol=(2100-1850+1)), na.rm=T), (2100-1850+1))
          tranij_raw_1  = tran_ij_tas$trend
          
          dare1 = getre(mrso_ij_1,mrro_ij_1,tran_ij_sm_1,pr___ij_sm_1,tran_ij_tas_1, tranij_raw_1,
                        tasij_raw,id1, id2, yr1)
          
          mrso_ijk        = dare1$mrso_ijk    
          mrro_ijk        = dare1$mrro_ijk    
          tran_ijk        = dare1$tran_ijk_sm 
          pr___ijk        = dare1$pr___ijk_sm 
          tran_ijk_tas_sd = dare1$tran_ijk_tas
          tran_raw_ijk    = dare1$tranij_raw    
          tas_raw_ijk     = dare1$tasij_raw 

          tran_raw_ijk_sm0 = rbind(colMeans(tran_raw_ijk[1:3,], na.rm=T),
                                   colMeans(tran_raw_ijk[2:4,], na.rm=T),
                                   colMeans(tran_raw_ijk[3:5,], na.rm=T))
          tran_raw_ijk_sm1 = rbind(tran_raw_ijk[4,],
                                   tran_raw_ijk[5,],
                                   tran_raw_ijk[6,])
          tran_dtr_dsea_ijk_sm0 = rbind(colMeans(tran_ijk_tas_sd[1:3,]),
                                        colMeans(tran_ijk_tas_sd[2:4,]),
                                        colMeans(tran_ijk_tas_sd[3:5,]))
          tran_dtr_dsea_ijk_sm1 = rbind(tran_ijk_tas_sd[4,],
                                        tran_ijk_tas_sd[5,],
                                        tran_ijk_tas_sd[6,])
          
          smt1 = try(lm(mrso_ijk[4,] ~ colMeans(mrso_ijk[1:3,], na.rm=T) + colMeans(pr___ijk[1:3,], na.rm=T) + colMeans(tran_ijk[1:3,], na.rm=T) + colMeans(mrro_ijk[1:3,], na.rm=T)), silent = T)
          smt2 = try(lm(mrso_ijk[5,] ~ colMeans(mrso_ijk[2:4,], na.rm=T) + colMeans(pr___ijk[2:4,], na.rm=T) + colMeans(tran_ijk[2:4,], na.rm=T) + colMeans(mrro_ijk[2:4,], na.rm=T)), silent = T)
          smt3 = try(lm(mrso_ijk[6,] ~ colMeans(mrso_ijk[3:5,], na.rm=T) + colMeans(pr___ijk[3:5,], na.rm=T) + colMeans(tran_ijk[3:5,], na.rm=T) + colMeans(mrro_ijk[3:5,], na.rm=T)), silent = T)
          
          if(class(smt1) != "try-error" &
             class(smt2) != "try-error" &
             class(smt3) != "try-error"){
            
            da_tran_sm[i,j,] = c(coef(smt1)[4], coef(smt2)[4],  coef(smt3)[4])
            deta_tran_sm     = mean(c(coef(smt1)[4], coef(smt2)[4],  coef(smt3)[4]), na.rm=T)
            ####sm mean####
            
            lai_sm0 = (tran_raw_ijk_sm0*86400*30)*deta_tran_sm
            lai_sm1 = (tran_raw_ijk_sm1*86400*30)*deta_tran_sm
            
            if(ncol(lai_sm0)!=(2100-1850+1)){lai_sm0 = cbind(rep(NA, 3), lai_sm0)}
            if(ncol(lai_sm1)!=(2100-1850+1)){lai_sm1 = cbind(rep(NA, 3), lai_sm1)}
            
            lai_sm_mean0_c[i,j,] = rowMeans(embed(colMeans(lai_sm0), window_yr), na.rm=T)
            lai_sm_mean1_c[i,j,] = rowMeans(embed(colMeans(lai_sm1), window_yr), na.rm=T)
            ####sm sd####
            lai_sm_sd0 = (tran_dtr_dsea_ijk_sm0*86400*30)*deta_tran_sm 
            lai_sm_sd1 = (tran_dtr_dsea_ijk_sm1*86400*30)*deta_tran_sm 
            if(ncol(lai_sm_sd0)!=(2100-1850+1)){lai_sm_sd0 = cbind(rep(NA,3), lai_sm_sd0)}
            if(ncol(lai_sm_sd1)!=(2100-1850+1)){lai_sm_sd1 = cbind(rep(NA,3), lai_sm_sd1)}
            
            lai_c_sm_sd0 = c()
            lai_c_sm_sd1 = c()
            for (ij in 1:nrow(loca_sel)) {#ij=1
              lai_c_sm_sd0 = c(lai_c_sm_sd0, sd(lai_sm_sd0[,loca_sel[ij,]], na.rm=T))
              lai_c_sm_sd1 = c(lai_c_sm_sd1, sd(lai_sm_sd1[,loca_sel[ij,]], na.rm=T))}
            
            lai_sm_sd0_c[i,j,] = lai_c_sm_sd0
            lai_sm_sd1_c[i,j,] = lai_c_sm_sd1
            ####tas mean####
            tas_ctr = mean(tas_raw_ijk[4:6,1:window_yr],na.rm=T)
            f = 1/(4*0.95*5.67*10^-8*tas_ctr^3)
            if(ncol(tran_raw_ijk)!=(2100-1850+1)){ tran_raw_ijk = cbind(rep(NA,6), tran_raw_ijk) }
            
            tran_raw_ijk_yr = tran_raw_ijk[4:6,] 
            lai_tas = tran_raw_ijk_yr*f*(-2.54*10^6)
            
            lai_tas_mean1_c[i,j,] = rowMeans(embed(colMeans(lai_tas), window_yr), na.rm=T)
            ####tas sd####
            lai_tas_sd1 = tran_dtr_dsea_ijk_sm1*f*(-2.54*10^6) 
            if(ncol(lai_tas_sd1)!=(2100-1850+1)){lai_tas_sd1 = cbind(rep(NA,3), lai_tas_sd1)}
            
            lai_c_sd1 = c()
            for (ij in 1:nrow(loca_sel)) {lai_c_sd1 = c(lai_c_sd1, sd(lai_tas_sd1[,loca_sel[ij,]], na.rm=T))  }
            
            lai_tas_sd1_c[i,j,] = lai_c_sd1
            ####cor####
            lai_c_sm_tas0 = c()
            lai_c_sm_tas1 = c()
            for (ijk in 1:nrow(loca_sel)) {#ijk=1
              lai_c_sm_tas1 = c(lai_c_sm_tas1, cov(na.omit(cbind(matrix(lai_sm_sd1 [,loca_sel[ijk,]], ncol=1), 
                                                                 matrix(lai_tas_sd1[,loca_sel[ijk,]], ncol=1))))[1,2])
              lai_c_sm_tas0 = c(lai_c_sm_tas0, cov(na.omit(cbind(matrix(lai_sm_sd0 [,loca_sel[ijk,]], ncol=1), 
                                                                 matrix(lai_tas_sd1[,loca_sel[ijk,]], ncol=1))))[1,2])}
            
            lai_cov_c0[i,j,] = lai_c_sm_tas0
            lai_cov_c1[i,j,] = lai_c_sm_tas1
            lai_sm_lai_tas_cor[i,j ] = cor(na.omit(cbind(matrix(lai_sm_sd0 , ncol=1),
                                                         matrix(lai_tas_sd1, ncol=1))))[1,2]
            ####tran_change#######
            
            if(ncol(tran_dtr_dsea_ijk_sm0)!=251){tran_dtr_dsea_ijk_sm0 = cbind(rep(NA,3), tran_dtr_dsea_ijk_sm0)}
            if(ncol(tran_dtr_dsea_ijk_sm1)!=251){tran_dtr_dsea_ijk_sm1 = cbind(rep(NA,3), tran_dtr_dsea_ijk_sm1)}
            tran_dtr_dsea_ijk_sm0 = tran_dtr_dsea_ijk_sm0*86400*30
            tran_dtr_dsea_ijk_sm1 = tran_dtr_dsea_ijk_sm1*86400*30
            
            cov_tran_01_ij = c()
            for (ijk in 1:nrow(loca_sel)) {#ijk=1
              cov_tran_01_ij = c(cov_tran_01_ij, 
                                 cov(na.omit(cbind(matrix(tran_dtr_dsea_ijk_sm0[,loca_sel[ijk,]], ncol=1), 
                                                   matrix(tran_dtr_dsea_ijk_sm1[,loca_sel[ijk,]], ncol=1))))[1,2])}
            
            cov_tran_c01[i,j,] = cov_tran_01_ij
          }
        }
      }
    }
    
    out = list(da_tran_sm           = da_tran_sm        ,
               tran_tas_mean1_c     = lai_tas_mean1_c   ,
               tran_tas_sd1_c       = lai_tas_sd1_c     ,
               tran_sm_mean0_c      = lai_sm_mean0_c    ,
               tran_sm_mean1_c      = lai_sm_mean1_c    ,
               tran_sm_sd0_c        = lai_sm_sd0_c      ,
               tran_sm_sd1_c        = lai_sm_sd1_c      ,
               tran_cov_c0          = lai_cov_c0        ,
               tran_cov_c1          = lai_cov_c1        ,
               cov_tran_c01         = cov_tran_c01      ,
               tran_sm_tran_tas_cor = lai_sm_lai_tas_cor
               )
    return(out)
  }
  
  t1 = proc.time()
  core=core1
  cl = makeCluster(core)
  registerDoParallel(cl)
  lmf_contribution = foreach (kk = 1:17, .packages=c('forecast','pracma','abind')) %dopar% {#kk=1
    out_re = get_re_sm_mean_g(fpall, modelna[kk], scenario, window_yr=window_yr)
    out_re
  }
  stopCluster(cl)
  t2 = proc.time()
  print((t2 - t1)[3]/60)
  
  outna = '3Month_no_ra_sm1_sm0_pr0_tran0_mrro0_ano_3mon_sen_mean_decompose_own_tas_by_sm_first30yr_win_'
  setwd('/GPUFS/ygo_fwf_1/00junli/01next_step/10results')
  saveRDS(lmf_contribution[[01]], paste(outna, scenario, '_',modelna[01], window_yr,sep=''))
  saveRDS(lmf_contribution[[02]], paste(outna, scenario, '_',modelna[02], window_yr,sep=''))
  saveRDS(lmf_contribution[[03]], paste(outna, scenario, '_',modelna[03], window_yr,sep=''))
  saveRDS(lmf_contribution[[04]], paste(outna, scenario, '_',modelna[04], window_yr,sep=''))
  saveRDS(lmf_contribution[[05]], paste(outna, scenario, '_',modelna[05], window_yr,sep=''))
  saveRDS(lmf_contribution[[06]], paste(outna, scenario, '_',modelna[06], window_yr,sep=''))
  saveRDS(lmf_contribution[[07]], paste(outna, scenario, '_',modelna[07], window_yr,sep=''))
  saveRDS(lmf_contribution[[08]], paste(outna, scenario, '_',modelna[08], window_yr,sep=''))
  saveRDS(lmf_contribution[[09]], paste(outna, scenario, '_',modelna[09], window_yr,sep=''))
  saveRDS(lmf_contribution[[10]], paste(outna, scenario, '_',modelna[10], window_yr,sep=''))
  saveRDS(lmf_contribution[[11]], paste(outna, scenario, '_',modelna[11], window_yr,sep=''))
  saveRDS(lmf_contribution[[12]], paste(outna, scenario, '_',modelna[12], window_yr,sep=''))
  saveRDS(lmf_contribution[[13]], paste(outna, scenario, '_',modelna[13], window_yr,sep=''))
  saveRDS(lmf_contribution[[14]], paste(outna, scenario, '_',modelna[14], window_yr,sep=''))
  saveRDS(lmf_contribution[[15]], paste(outna, scenario, '_',modelna[15], window_yr,sep=''))
  saveRDS(lmf_contribution[[16]], paste(outna, scenario, '_',modelna[16], window_yr,sep=''))
  saveRDS(lmf_contribution[[17]], paste(outna, scenario, '_',modelna[17], window_yr,sep=''))

  rm(lmf_contribution)
  gc()
}
#cal_tran_induced_change(fpall, modelna, '_ssp585_', 17, window_yr=30)
#cal_tran_induced_change(fpall, modelna, '_ssp370_', 17, window_yr=30)
#cal_tran_induced_change(fpall, modelna, '_ssp245_', 17, window_yr=30)
#cal_tran_induced_change(fpall, modelna, '_ssp126_', 17, window_yr=30)
#cal_tran_induced_change(fpall, modelna, '_ssp585_', 17, window_yr=20)
#cal_tran_induced_change(fpall, modelna, '_ssp585_', 17, window_yr=40)

###################Albedo-induced change##########################################################

#modelna
#fpall = '/GPUFS/ygo_fwf_1/00junli/01next_step/GCM_all_data_his_126_245_370_585/all_gcm_da'
cal_albedo_induced_change = function(fpall, modelna, scenario, core1, window_yr){
  
  get_re_sm_mean_g = function(fpall, modelnai, scenario, window_yr=window_yr){
    
    #window_yr=20; modelnai = modelna[01]; scenario='_ssp585_'
    setwd(fpall)
    fp_file = list.files(fpall)
    for (k in 1:length(fp_file)) {#k=10
      con_tas    = grepl(c('tas'         ), fp_file[k])
      con_mrso_  = grepl(c('mrso_'       ), fp_file[k])
      con_mrro_  = grepl(c('mrro_'       ), fp_file[k])
      con__lai   = grepl(c('lai_L'       ), fp_file[k])
      con_et     = grepl(c('evspsbl_A'   ), fp_file[k])#
      con_tran   = grepl(c('tran_L'      ), fp_file[k])#evspsbl_A
      con_rsds   = grepl(c('rsds'        ), fp_file[k])
      con_rsus   = grepl(c('rsus'        ), fp_file[k])
      con_ts     = grepl(c('ts_A'        ), fp_file[k])
      con_hfss   = grepl(c('hfss_A'      ), fp_file[k])
      con_wind   = grepl(c('sfcWind_A'   ), fp_file[k])
      con_pr     = grepl(c('pr_A'        ), fp_file[k])
      con_his    = grepl(c('_historical_'), fp_file[k])
      con_ssp585 = grepl(c(scenario      ), fp_file[k])
      con_modelna= grepl(c(modelnai      ), fp_file[k])
      
      if(con_tas   & con_his & con_modelna) {tas__his_bgcna = fp_file[k]}
      if(con_mrso_ & con_his & con_modelna) {mrso_his_bgcna = fp_file[k]}
      if(con_mrro_ & con_his & con_modelna) {mrro_his_bgcna = fp_file[k]}
      if(con__lai  & con_his & con_modelna) {lai__his_bgcna = fp_file[k]}
      if(con_et    & con_his & con_modelna) {et_his_bgcna   = fp_file[k]}
      if(con_tran  & con_his & con_modelna) {tranhis_bgcna  = fp_file[k]}
      if(con_rsds  & con_his & con_modelna) {rsds_his_bgcna = fp_file[k]}
      if(con_rsus  & con_his & con_modelna) {rsus_his_bgcna = fp_file[k]}
      if(con_ts    & con_his & con_modelna) {ts_his_bgcna   = fp_file[k]}
      if(con_hfss  & con_his & con_modelna) {hfss_his_bgcna = fp_file[k]}
      if(con_wind  & con_his & con_modelna) {wind_his_bgcna = fp_file[k]}
      if(con_pr    & con_his & con_modelna) {pr_his_bgcna   = fp_file[k]}
      
      if(con_tas   & con_ssp585 & con_modelna) {tas__ssp_bgcna = fp_file[k]}
      if(con_mrso_ & con_ssp585 & con_modelna) {mrso_ssp_bgcna = fp_file[k]}
      if(con_mrro_ & con_ssp585 & con_modelna) {mrro_ssp_bgcna = fp_file[k]}
      if(con__lai  & con_ssp585 & con_modelna) {lai__ssp_bgcna = fp_file[k]}
      if(con_et    & con_ssp585 & con_modelna) {et_ssp_bgcna   = fp_file[k]}
      if(con_tran  & con_ssp585 & con_modelna) {transsp_bgcna  = fp_file[k]}
      if(con_rsds  & con_ssp585 & con_modelna) {rsds_ssp_bgcna = fp_file[k]}
      if(con_rsus  & con_ssp585 & con_modelna) {rsus_ssp_bgcna = fp_file[k]}
      if(con_ts    & con_ssp585 & con_modelna) {ts_ssp_bgcna   = fp_file[k]}
      if(con_hfss  & con_ssp585 & con_modelna) {hfss_ssp_bgcna = fp_file[k]}
      if(con_wind  & con_ssp585 & con_modelna) {wind_ssp_bgcna = fp_file[k]}
      if(con_pr    & con_ssp585 & con_modelna) {pr_ssp_bgcna   = fp_file[k]}
    }
    tas__his = readRDS(tas__his_bgcna)
    tran_his = readRDS(tranhis_bgcna )
    mrso_his = readRDS(mrso_his_bgcna)
    mrro_his = readRDS(mrro_his_bgcna)
    rsds_his = readRDS(rsds_his_bgcna)
    pr___his = readRDS(pr_his_bgcna  )
    rsus_his = readRDS(rsus_his_bgcna)
    #ts___his = readRDS(ts_his_bgcna  )
    wind_his = readRDS(wind_his_bgcna)#lai__his = readRDS(lai__his_bgcna)
    
    tas__ssp = readRDS(tas__ssp_bgcna)
    tran_ssp = readRDS(transsp_bgcna )
    mrso_ssp = readRDS(mrso_ssp_bgcna)
    mrro_ssp = readRDS(mrro_ssp_bgcna)
    rsds_ssp = readRDS(rsds_ssp_bgcna)
    pr___ssp = readRDS(pr_ssp_bgcna  )#lai__ssp = readRDS(lai__ssp_bgcna)
    rsus_ssp = readRDS(rsus_ssp_bgcna)
    #ts___ssp = readRDS(ts_ssp_bgcna  )
    wind_ssp = readRDS(wind_ssp_bgcna)
    
    tas__bgc = abind(tas__his, tas__ssp); rm(tas__his); rm(tas__ssp); print(dim(tas__bgc))
    mrso_bgc = abind(mrso_his, mrso_ssp); rm(mrso_his); rm(mrso_ssp); print(dim(mrso_bgc))
    mrro_bgc = abind(mrro_his, mrro_ssp); rm(mrro_his); rm(mrro_ssp); print(dim(mrro_bgc))
    tran_bgc = abind(tran_his, tran_ssp); rm(tran_his); rm(tran_ssp); print(dim(tran_bgc))
    rsds_bgc = abind(rsds_his, rsds_ssp); rm(rsds_his); rm(rsds_ssp); print(dim(rsds_bgc))
    pr___bgc = abind(pr___his, pr___ssp); rm(pr___his); rm(pr___ssp); print(dim(pr___bgc))
    rsus_bgc = abind(rsus_his, rsus_ssp); rm(rsus_his); rm(rsus_ssp); print(dim(rsus_bgc))
    #ts___bgc = abind(ts___his, ts___ssp); rm(ts___his); rm(ts___ssp); print(dim(ts___bgc))
    wind_bgc = abind(wind_his, wind_ssp); rm(wind_his); rm(wind_ssp); print(dim(wind_bgc))
    #lai__bgc = abind(lai__his, lai__ssp); rm(lai__his); rm(lai__ssp); print(dim(lai__bgc))
    
    albedo_bgc = rsus_bgc/rsds_bgc
    albedo_bgc[is.infinite(albedo_bgc)] = NA
    ra___bgc = 208/wind_bgc
    #rm(wind_bgc); gc()     #ra___bgc   = 1.2*1013*(ts___bgc - tas__bgc)/hfss_bgc
    rm(rsus_bgc); gc() 
    
    sca=3
    yr=2100-1850+1
    len = yr-window_yr+1
    nrow1=240
    ncol1=93
    da_tran_sm         = array(NA, c(nrow1,ncol1,3  ))     
    lai_tas_mean1_c    = array(NA, c(nrow1,ncol1,len))     
    lai_tas_sd1_c      = array(NA, c(nrow1,ncol1,len))   
    lai_sm_mean0_c     = array(NA, c(nrow1,ncol1,len))    
    lai_sm_mean1_c     = array(NA, c(nrow1,ncol1,len))    
    lai_sm_sd0_c       = array(NA, c(nrow1,ncol1,len))  
    lai_sm_sd1_c       = array(NA, c(nrow1,ncol1,len))  
    lai_cov_c0         = array(NA, c(nrow1,ncol1,len))
    lai_cov_c1         = array(NA, c(nrow1,ncol1,len))
    cov_tran_c01       = array(NA, c(nrow1,ncol1,len))
    lai_sm_lai_tas_cor = matrix(NA, nrow=nrow1, ncol=ncol1)
    lucc = readRDS('/GPUFS/ygo_fwf_1/00junli/01next_step/lucc2016_1.5degree_vegetated_area.rds');col=1
    for (i in 1:nrow(tas__bgc)) {
      for (j in 1:ncol(tas__bgc)) {
        if(is.na(lucc[i,j]) | is.na(tas__bgc[i,j,1]) | 
           is.na(mrso_bgc[i,j,1]) | is.na(pr___bgc[i,j,1]) | is.na(tran_bgc[i,j,1])) next
        #i=200;j=80;image.plot(1:240,1:93,tas__bgc[,,8]);abline(v=i,h=j); gc()
        tas__ij  = rowMeans(embed(c(tas__bgc[i,j,],rep(NA, sca-1)), sca), na.rm=T)
        loca_bgc = which.max(rowMeans(matrix(tas__ij, ncol=yr), na.rm=T))
        locaij = sort(embed(c(1:12,1,2),3)[loca_bgc,])
        loca_sel = embed(1:yr, window_yr)
        
        mrso_ij    = try(decompose_own(ts(mrso_bgc  [i,j,]         , frequency=12), type = "additive", window_yr = window_yr), silent = T)
        mrro_ij    = try(decompose_own(ts(mrro_bgc  [i,j,]*86400*30, frequency=12), type = "additive", window_yr = window_yr), silent = T)
        tran_ij_sm = try(decompose_own(ts(tran_bgc  [i,j,]*86400*30, frequency=12), type = "additive", window_yr = window_yr), silent = T)
        pr___ij_sm = try(decompose_own(ts(pr___bgc  [i,j,]*86400*30, frequency=12), type = "additive", window_yr = window_yr), silent = T)
        tran_ij_tas= try(decompose_own(ts(tran_bgc  [i,j,]         , frequency=12), type = "additive", window_yr = window_yr), silent = T)
        albe_ij_tas= try(decompose_own(ts(albedo_bgc[i,j,]         , frequency=12), type = "additive", window_yr = window_yr), silent = T)
        
        tranij_raw  = tran_bgc[i,j,]
        tasij_raw   = tas__bgc[i,j,]
        sw_ij_raw   = rsds_bgc[i,j,]
        if(class(mrso_ij    )[1] != "try-error" &
           class(mrro_ij    )[1] != "try-error" &
           class(tran_ij_sm )[1] != "try-error" &
           class(pr___ij_sm )[1] != "try-error" &
           class(albe_ij_tas)[1] != "try-error" &
           class(tran_ij_tas)[1] != "try-error"){
          
          getre = function(mrso_ij   ,mrro_ij  , tran_ij_sm,pr___ij_sm ,tran_ij_tas,albe_ij_tas, 
                           tranij_raw,tasij_raw, sw_ij_raw ,albe_ij_raw, id1=09:3008, id2=c(1:6), yr=251){
            return(list(mrso_ijk     = matrix(mrso_ij    [id1], ncol=yr)[id2,],
                        mrro_ijk     = matrix(mrro_ij    [id1], ncol=yr)[id2,],
                        tran_ijk_sm  = matrix(tran_ij_sm [id1], ncol=yr)[id2,],
                        pr___ijk_sm  = matrix(pr___ij_sm [id1], ncol=yr)[id2,],
                        tran_ijk_tas = matrix(tran_ij_tas[id1], ncol=yr)[id2,],
                        albe_ij_tas  = matrix(albe_ij_tas[id1], ncol=yr)[id2,],
                        tranij_raw   = matrix(tranij_raw [id1], ncol=yr)[id2,],
                        tasij_raw    = matrix(tasij_raw  [id1], ncol=yr)[id2,],
                        sw_ij_raw    = matrix(sw_ij_raw  [id1], ncol=yr)[id2,],
                        albe_ij_raw  = matrix(albe_ij_raw[id1], ncol=yr)[id2,])) }
          
          if(loca_bgc==12){
            id1 = 09:3008; id2 = c(1:6); yr1 = 2100-1850 }else if(loca_bgc==11){
              id1 = 08:3007; id2 = c(1:6); yr1 = 2100-1850 }else if(loca_bgc==1){
                id1 = 10:3009; id2 = c(1:6); yr1 = 2100-1850 }else if(loca_bgc==2){
                  id1 = 11:3010; id2 = c(1:6); yr1 = 2100-1850 }else if(loca_bgc==3){
                    id1 = 12:3011; id2 = c(1:6); yr1 = 2100-1850 }else{
                      id1 = 1:3012; id2 = c((min(locaij)-3):max(locaij)); yr1 = 2100-1850+1
                    }
          
          mrso_ij_1     = mrso_ij    $random
          mrro_ij_1     = mrro_ij    $random
          tran_ij_sm_1  = tran_ij_sm $random
          pr___ij_sm_1  = pr___ij_sm $random
          tran_ij_tas_1 = tran_ij_tas$random
          albe_ij_tas_1 = albe_ij_tas$random
          tranij_raw_1  = tran_ij_tas$trend
          albe_ij_raw   = albe_ij_tas$trend
          
          dare1 = getre(mrso_ij_1   ,mrro_ij_1,tran_ij_sm_1,pr___ij_sm_1,tran_ij_tas_1, albe_ij_tas_1,
                        tranij_raw_1,tasij_raw, sw_ij_raw  ,albe_ij_raw , id1, id2, yr1)
          
          mrso_ijk        = dare1$mrso_ijk    
          mrro_ijk        = dare1$mrro_ijk    
          tran_ijk        = dare1$tran_ijk_sm 
          pr___ijk        = dare1$pr___ijk_sm 
          tran_ijk_tas_sd = dare1$tran_ijk_tas
          albe_ijk_tas_sd = dare1$albe_ij_tas
          tran_raw_ijk    = dare1$tranij_raw    
          tas_raw_ijk     = dare1$tasij_raw 
          sw__raw_ijk     = dare1$sw_ij_raw 
          albe_ij_raw     = dare1$albe_ij_raw 
          
          tran_raw_ijk_sm0 = rbind(colMeans(tran_raw_ijk[1:3,], na.rm=T),
                                   colMeans(tran_raw_ijk[2:4,], na.rm=T),
                                   colMeans(tran_raw_ijk[3:5,], na.rm=T))
          tran_raw_ijk_sm1 = rbind(tran_raw_ijk[4,],
                                   tran_raw_ijk[5,],
                                   tran_raw_ijk[6,])
          tran_dtr_dsea_ijk_sm0 = rbind(colMeans(tran_ijk_tas_sd[1:3,]),
                                        colMeans(tran_ijk_tas_sd[2:4,]),
                                        colMeans(tran_ijk_tas_sd[3:5,]))
          tran_dtr_dsea_ijk_sm1 = rbind(tran_ijk_tas_sd[4,],
                                        tran_ijk_tas_sd[5,],
                                        tran_ijk_tas_sd[6,])
          
          smt1 = try(lm(mrso_ijk[4,] ~ colMeans(mrso_ijk[1:3,], na.rm=T) + colMeans(pr___ijk[1:3,], na.rm=T) + colMeans(tran_ijk[1:3,], na.rm=T) + colMeans(mrro_ijk[1:3,], na.rm=T)), silent = T)
          smt2 = try(lm(mrso_ijk[5,] ~ colMeans(mrso_ijk[2:4,], na.rm=T) + colMeans(pr___ijk[2:4,], na.rm=T) + colMeans(tran_ijk[2:4,], na.rm=T) + colMeans(mrro_ijk[2:4,], na.rm=T)), silent = T)
          smt3 = try(lm(mrso_ijk[6,] ~ colMeans(mrso_ijk[3:5,], na.rm=T) + colMeans(pr___ijk[3:5,], na.rm=T) + colMeans(tran_ijk[3:5,], na.rm=T) + colMeans(mrro_ijk[3:5,], na.rm=T)), silent = T)
          
          if(class(smt1) != "try-error" &
             class(smt2) != "try-error" &
             class(smt3) != "try-error"){
            
            #da_tran_sm[i,j,] = c(coef(smt1)[4], coef(smt2)[4],  coef(smt3)[4])
            deta_tran_sm = mean(c(coef(smt1)[4], coef(smt2)[4],  coef(smt3)[4]), na.rm=T)
            ####sm mean####
            
            #lai_sm0 = (tran_raw_ijk_sm0*86400*30)*deta_tran_sm
            #lai_sm1 = (tran_raw_ijk_sm1*86400*30)*deta_tran_sm
            #
            #if(ncol(lai_sm0)!=(2100-1850+1)){lai_sm0 = cbind(rep(NA, 3), lai_sm0)}
            #if(ncol(lai_sm1)!=(2100-1850+1)){lai_sm1 = cbind(rep(NA, 3), lai_sm1)}
            #
            #lai_sm_mean0_c[i,j,] = rowMeans(embed(colMeans(lai_sm0), window_yr), na.rm=T)
            #lai_sm_mean1_c[i,j,] = rowMeans(embed(colMeans(lai_sm1), window_yr), na.rm=T)
            ####sm sd####
            lai_sm_sd0 = (tran_dtr_dsea_ijk_sm0*86400*30)*deta_tran_sm 
            lai_sm_sd1 = (tran_dtr_dsea_ijk_sm1*86400*30)*deta_tran_sm 
            if(ncol(lai_sm_sd0)!=(2100-1850+1)){lai_sm_sd0 = cbind(rep(NA,3), lai_sm_sd0)}
            if(ncol(lai_sm_sd1)!=(2100-1850+1)){lai_sm_sd1 = cbind(rep(NA,3), lai_sm_sd1)}
            
            lai_c_sm_sd0 = c()
            lai_c_sm_sd1 = c()
            for (ij in 1:nrow(loca_sel)) {#ij=1
              lai_c_sm_sd0 = c(lai_c_sm_sd0, sd(lai_sm_sd0[,loca_sel[ij,]], na.rm=T))
              lai_c_sm_sd1 = c(lai_c_sm_sd1, sd(lai_sm_sd1[,loca_sel[ij,]], na.rm=T))}
            
            lai_sm_sd0_c[i,j,] = lai_c_sm_sd0
            lai_sm_sd1_c[i,j,] = lai_c_sm_sd1
            ####tas mean####
            rsdsctr = mean(sw__raw_ijk[4:6,1:window_yr],na.rm=T)
            tas_ctr = mean(tas_raw_ijk[4:6,1:window_yr],na.rm=T)
            
            f = 1/(4*0.95*5.67*10^-8*tas_ctr^3)
            if(ncol(albe_ij_raw)!=(2100-1850+1)){ albe_ij_raw = cbind(rep(NA,6), albe_ij_raw) }
            
            lai_tas = albe_ij_raw[4:6,]*f*(-rsdsctr)
            
            lai_tas_mean1_c[i,j,] = rowMeans(embed(colMeans(lai_tas), window_yr), na.rm=T)
            ####tas sd####
            albe_dtr_dsea_ijk_sm1 = rbind(albe_ijk_tas_sd[4,],
                                          albe_ijk_tas_sd[5,],
                                          albe_ijk_tas_sd[6,])#albe_ijk_tas_sd
            
            lai_tas_sd1 = albe_dtr_dsea_ijk_sm1*f*(-rsdsctr)
            if(ncol(lai_tas_sd1)!=(2100-1850+1)){lai_tas_sd1 = cbind(rep(NA,3), lai_tas_sd1)}
            
            lai_c_sd1 = c()
            for (ij in 1:nrow(loca_sel)) {lai_c_sd1 = c(lai_c_sd1, sd(lai_tas_sd1[,loca_sel[ij,]], na.rm=T))  }
            
            lai_tas_sd1_c[i,j,] = lai_c_sd1
            ####cor####
            lai_c_sm_tas0 = c()
            lai_c_sm_tas1 = c()
            for (ijk in 1:nrow(loca_sel)) {#ijk=1
              lai_c_sm_tas1 = c(lai_c_sm_tas1, cov(na.omit(cbind(matrix(lai_sm_sd1 [,loca_sel[1,]], ncol=1), 
                                                                 matrix(lai_tas_sd1[,loca_sel[ijk,]], ncol=1))))[1,2])
              lai_c_sm_tas0 = c(lai_c_sm_tas0, cov(na.omit(cbind(matrix(lai_sm_sd0 [,loca_sel[1,]], ncol=1), 
                                                                 matrix(lai_tas_sd1[,loca_sel[ijk,]], ncol=1))))[1,2])}
            
            lai_cov_c0[i,j,] = lai_c_sm_tas0
            lai_cov_c1[i,j,] = lai_c_sm_tas1
            ##
            ####tran_change#######
            
            #if(ncol(tran_dtr_dsea_ijk_sm0)!=251){tran_dtr_dsea_ijk_sm0 = cbind(rep(NA,3), tran_dtr_dsea_ijk_sm0)}
            #if(ncol(tran_dtr_dsea_ijk_sm1)!=251){tran_dtr_dsea_ijk_sm1 = cbind(rep(NA,3), tran_dtr_dsea_ijk_sm1)}
            #tran_dtr_dsea_ijk_sm0 = tran_dtr_dsea_ijk_sm0*86400*30
            #tran_dtr_dsea_ijk_sm1 = tran_dtr_dsea_ijk_sm1*86400*30
            #
            #cov_tran_01_ij = c()
            #for (ijk in 1:nrow(loca_sel)) {#ijk=1
            #  cov_tran_01_ij = c(cov_tran_01_ij, 
            #                     cov(na.omit(cbind(matrix(tran_dtr_dsea_ijk_sm0[,loca_sel[ijk,]], ncol=1), 
            #                                       matrix(tran_dtr_dsea_ijk_sm1[,loca_sel[ijk,]], ncol=1))))[1,2])}
            #
            #cov_tran_c01[i,j,] = cov_tran_01_ij
            #col=col+1; print(col)
          }
        }
      }
    }
    
    out = list(alb_tas_mean1_c = lai_tas_mean1_c,
               alb_tas_sd1_c   = lai_tas_sd1_c  ,
               alb_sm_mean0_c  = lai_sm_mean0_c ,
               alb_sm_mean1_c  = lai_sm_mean1_c ,
               alb_sm_sd0_c    = lai_sm_sd0_c   ,
               alb_sm_sd1_c    = lai_sm_sd1_c   ,
               alb_cov_c0      = lai_cov_c0     ,
               alb_cov_c1      = lai_cov_c1     )
    return(out)
  }
  
  t1 = proc.time()
  core=core1
  cl = makeCluster(core)
  registerDoParallel(cl)
  lmf_contribution = foreach (kk = 1:17, .packages=c('forecast','pracma','abind')) %dopar% {#kk=1
    out_re = get_re_sm_mean_g(fpall, modelna[kk], scenario, window_yr=window_yr)
    out_re
  }
  stopCluster(cl)
  t2 = proc.time()
  print((t2 - t1)[3]/60)
  
  outna = 'Albedo_induce_change_3mon_'
  setwd('/GPUFS/ygo_fwf_1/00junli/01next_step/10results_Albedo_induce_change')
  saveRDS(lmf_contribution[[01]], paste(outna, scenario, '_',modelna[01], window_yr,sep=''))
  saveRDS(lmf_contribution[[02]], paste(outna, scenario, '_',modelna[02], window_yr,sep=''))
  saveRDS(lmf_contribution[[03]], paste(outna, scenario, '_',modelna[03], window_yr,sep=''))
  saveRDS(lmf_contribution[[04]], paste(outna, scenario, '_',modelna[04], window_yr,sep=''))
  saveRDS(lmf_contribution[[05]], paste(outna, scenario, '_',modelna[05], window_yr,sep=''))
  saveRDS(lmf_contribution[[06]], paste(outna, scenario, '_',modelna[06], window_yr,sep=''))
  saveRDS(lmf_contribution[[07]], paste(outna, scenario, '_',modelna[07], window_yr,sep=''))
  saveRDS(lmf_contribution[[08]], paste(outna, scenario, '_',modelna[08], window_yr,sep=''))
  saveRDS(lmf_contribution[[09]], paste(outna, scenario, '_',modelna[09], window_yr,sep=''))
  saveRDS(lmf_contribution[[10]], paste(outna, scenario, '_',modelna[10], window_yr,sep=''))
  saveRDS(lmf_contribution[[11]], paste(outna, scenario, '_',modelna[11], window_yr,sep=''))
  saveRDS(lmf_contribution[[12]], paste(outna, scenario, '_',modelna[12], window_yr,sep=''))
  saveRDS(lmf_contribution[[13]], paste(outna, scenario, '_',modelna[13], window_yr,sep=''))
  saveRDS(lmf_contribution[[14]], paste(outna, scenario, '_',modelna[14], window_yr,sep=''))
  saveRDS(lmf_contribution[[15]], paste(outna, scenario, '_',modelna[15], window_yr,sep=''))
  saveRDS(lmf_contribution[[16]], paste(outna, scenario, '_',modelna[16], window_yr,sep=''))
  saveRDS(lmf_contribution[[17]], paste(outna, scenario, '_',modelna[17], window_yr,sep=''))

  rm(lmf_contribution)
  gc()
}
#cal_albedo_induced_change(fpall, modelna, '_ssp585_', 17, window_yr=30)
#cal_albedo_induced_change(fpall, modelna, '_ssp370_', 17, window_yr=30)
#cal_albedo_induced_change(fpall, modelna, '_ssp245_', 17, window_yr=30)
#cal_albedo_induced_change(fpall, modelna, '_ssp126_', 17, window_yr=30)
#cal_albedo_induced_change(fpall, modelna, '_ssp585_', 17, window_yr=20)
#cal_albedo_induced_change(fpall, modelna, '_ssp585_', 17, window_yr=40)

###################lmf Prediction####################################################

lmf_pre1 = function(rena, lai_gs_indc, modelna, threshold=0.9, window_yr=50, out_na){
  t1 = proc.time()
  get_da = function(lai_gs_ind,lmf_sm_mean ,lmf_ta_mean ,lmf_sm_ta_sd,lmf_dep,
                    lmf_total ,sm_mean_c   ,ta_mean_c   ,sm_ta_sdc   ,cor_c  ,i,j, id=len, da_len=7*len){#i=70;j=78
    
    #id = 221;i=200;j=70
    lmf_sm_mean_ij  = lmf_sm_mean [i,j,]
    lmf_ta_mean_ij  = lmf_ta_mean [i,j,]
    lmf_sm_ta_sd_ij = lmf_sm_ta_sd[i,j,]
    lmf_dep_ij      = lmf_dep     [i,j,]
    #lmf_total_ij    = lmf_total   [i,j,]
    sm_mean_c_ij    = sm_mean_c   [i,j,]
    ta_mean_c_ij    = ta_mean_c   [i,j,]
    sm_ta_sdc_ij    = sm_ta_sdc   [i,j,]
    cor_c_ij        = cor_c       [i,j,]
    
    dasm_mean_ij = cbind(lmf_sm_mean_ij  = lmf_sm_mean_ij , sm_mean_c_ij = sm_mean_c_ij)
    data_mean_ij = cbind(lmf_ta_mean_ij  = lmf_ta_mean_ij , ta_mean_c_ij = ta_mean_c_ij)
    dasm_ta_v_ij = cbind(lmf_sm_ta_sd_ij = lmf_sm_ta_sd_ij, sm_ta_sdc_ij = sm_ta_sdc_ij)
    dadep_ij     = cbind(lmf_dep_ij      = lmf_dep_ij     , cor_c_ij     = cor_c_ij    )
    
    dasm_mean_ij[is.infinite(dasm_mean_ij)] = NA
    data_mean_ij[is.infinite(data_mean_ij)] = NA
    dasm_ta_v_ij[is.infinite(dasm_ta_v_ij)] = NA
    dadep_ij    [is.infinite(dadep_ij    )] = NA
    
    dasm_mean_ij = as.data.frame(na.omit(dasm_mean_ij))
    data_mean_ij = as.data.frame(na.omit(data_mean_ij))
    dasm_ta_v_ij = as.data.frame(na.omit(dasm_ta_v_ij))
    dadep_ij     = as.data.frame(na.omit(dadep_ij    ))
    
    mod_sm_mean_ij = try(segmented(lm(lmf_sm_mean_ij  ~ sm_mean_c_ij, data=dasm_mean_ij), seg.Z = ~sm_mean_c_ij), silent = T)
    mod_ta_mean_ij = try(segmented(lm(lmf_ta_mean_ij  ~ ta_mean_c_ij, data=data_mean_ij), seg.Z = ~ta_mean_c_ij), silent = T)
    mod_sm_ta_v_ij = try(segmented(lm(lmf_sm_ta_sd_ij ~ sm_ta_sdc_ij, data=dasm_ta_v_ij), seg.Z = ~sm_ta_sdc_ij), silent = T)
    mod_dep_ij     = try(segmented(lm(lmf_dep_ij      ~ cor_c_ij    , data=dadep_ij    ), seg.Z = ~cor_c_ij    ), silent = T)
    
    #plot(dasm_mean_ij[,2], dasm_mean_ij[,1], col='orange'); points(dasm_mean_ij[,2], predict(mod_sm_mean_ij));abline(h=0,v=0,lty=2, col='blue')
    #plot(data_mean_ij[,2], data_mean_ij[,1], col='orange'); points(data_mean_ij[,2], predict(mod_ta_mean_ij));abline(h=0,v=0,lty=2, col='blue')
    #plot(dasm_ta_v_ij[,2], dasm_ta_v_ij[,1], col='orange'); points(dasm_ta_v_ij[,2], predict(mod_sm_ta_v_ij));abline(h=0,v=0,lty=2, col='blue')
    #plot(dadep_ij    [,2], dadep_ij    [,1], col='orange'); points(dadep_ij    [,2], predict(mod_dep_ij    ));abline(h=0,v=0,lty=2, col='blue')
    con1 = sum(!is.na(sm_mean_c_ij))==id & sum(!is.na(lmf_sm_mean_ij )) ==id
    con2 = sum(!is.na(ta_mean_c_ij))==id & sum(!is.na(lmf_ta_mean_ij )) ==id
    con3 = sum(!is.na(sm_ta_sdc_ij))==id & sum(!is.na(lmf_sm_ta_sd_ij)) ==id
    con4 = sum(!is.na(cor_c_ij    ))==id & sum(!is.na(lmf_dep_ij     )) ==id
    
    lai_tas_mean_all1 = lai_gs_ind$lai_tas_mean1_c                               [i,j,]    
    lai_sm_mean_all1  = lai_gs_ind$lai_sm_mean1_c                                [i,j,]    
    lai_sm_mean_all0  = lai_gs_ind$lai_sm_mean0_c                                [i,j,]    
    lai_var_all1      = (lai_gs_ind$lai_sm_sd1_c[i,j,])*(lai_gs_ind$lai_tas_sd1_c[i,j,])                          
    lai_cov_all1      = lai_gs_ind$lai_cov_c1                                    [i,j,]
    lai_var_all0      = (lai_gs_ind$lai_sm_sd0_c[i,j,])*(lai_gs_ind$lai_tas_sd1_c[i,j,])                             
    lai_cov_all0      = lai_gs_ind$lai_cov_c0                                    [i,j,]                            
    
    lai_tas_mean_c1 = lai_tas_mean_all1[2:(id+1)] - lai_tas_mean_all1[1]
    lai_sm_mean_c1  = lai_sm_mean_all1 [2:(id+1)] - lai_sm_mean_all1 [1]
    lai_sm_mean_c0  = lai_sm_mean_all0 [2:(id+1)] - lai_sm_mean_all0 [1]
    lai_var_c1      = lai_var_all1     [2:(id+1)] - lai_var_all1     [1]
    lai_cov_c1      = lai_cov_all1     [2:(id+1)] - lai_cov_all1     [1]
    lai_var_c0      = lai_var_all0     [2:(id+1)] - lai_var_all0     [1]
    lai_cov_c0      = lai_cov_all0     [2:(id+1)] - lai_cov_all0     [1]
    
    if(class(mod_sm_mean_ij)[1] != "try-error" & con1 &  #& class(mod_sm_mean_ij)[2] != "try-error"   
       class(mod_ta_mean_ij)[1] != "try-error" & con2 &  #& class(mod_ta_mean_ij)[2] != "try-error"   
       class(mod_sm_ta_v_ij)[1] != "try-error" & con3 &  #& class(mod_ta_mean_ij)[2] != "try-error"   
       class(mod_dep_ij    )[1] != "try-error" & con4){  #& class(mod_dep_ij    )[2] != "try-error"
      
      if(is.na(mean(lai_tas_mean_c1, na.rm=T))){
        out_da = rep(NA, da_len)
      }else{
        ctr_ta_mean = predict(mod_ta_mean_ij, data.frame(ta_mean_c_ij = 0))
        ctr_sm_mean = predict(mod_sm_mean_ij, data.frame(sm_mean_c_ij = 0))
        ctr_sd      = predict(mod_sm_ta_v_ij, data.frame(sm_ta_sdc_ij = 0))
        ctr_dep     = predict(mod_dep_ij    , data.frame(cor_c_ij     = 0))
        
        out_da = c(predict(mod_ta_mean_ij, data.frame(ta_mean_c_ij  =  lai_tas_mean_c1)) - ctr_ta_mean,
                   predict(mod_sm_mean_ij, data.frame(sm_mean_c_ij  = -lai_sm_mean_c1 )) + ctr_sm_mean,
                   predict(mod_sm_mean_ij, data.frame(sm_mean_c_ij  = -lai_sm_mean_c0 )) + ctr_sm_mean,
                  (predict(mod_sm_ta_v_ij, data.frame(sm_ta_sdc_ij  =  lai_var_c1     )) - ctr_sd ),
                  (predict(mod_dep_ij    , data.frame(cor_c_ij      = -lai_cov_c1     )) - ctr_dep),
                  (predict(mod_sm_ta_v_ij, data.frame(sm_ta_sdc_ij  =  lai_var_c0     )) - ctr_sd ),
                  (predict(mod_dep_ij    , data.frame(cor_c_ij      = -lai_cov_c0     )) - ctr_dep))}
    }else{
      out_da = rep(NA, da_len) }
    return(out_da)
  }
  cal_re = function(rena, lai_gs_indc, id=len, da_len=7*len){
    #id = 221;len=221;rena = rena585_0.90_30[01]; lai_gs_indc = tran_indc585_30yr[01];  id=len; da_len=6*len
    lucc = readRDS('/GPUFS/ygo_fwf_1/00junli/01next_step/lucc2016_1.5degree_vegetated_area.rds')
    setwd('/GPUFS/ygo_fwf_1/00junli/01next_step/7lmf_all_and_ensemble')
    dare  = readRDS(rena)
    
    setwd('/GPUFS/ygo_fwf_1/00junli/01next_step/10results')
    lai_gs_ind = readRDS(lai_gs_indc)
    
    lmf_sm_mean = dare[,,1        : (id*01)]
    lmf_ta_mean = dare[,,(id*01+1): (id*02)]
    lmf_sm_ta_sd= dare[,,(id*02+1): (id*03)]
    lmf_dep     = dare[,,(id*03+1): (id*04)]
    lmf_total   = dare[,,(id*04+1): (id*05)]
    sm_mean_c   = dare[,,(id*05+1): (id*06)]
    ta_mean_c   = dare[,,(id*06+1): (id*07)]
    sm_ta_sdc   = dare[,,(id*07+1): (id*08)]
    cor_c       = dare[,,(id*08+1): (id*09)]#
    rm(dare)
    gc()
    
    lmf_predict = array(NA, c(dim(sm_mean_c)[1],dim(sm_mean_c)[2],da_len));col=1
    for(k1 in 1:dim(sm_mean_c)[1]) { 
      for(k2 in 1:dim(sm_mean_c)[2]) {  #k1=200;k2=80; 
        if(is.na(lucc[k1,k2]) | is.na(sm_mean_c[k1,k2,50]) | is.na(lmf_sm_mean[k1,k2,50])) next
        
        lmf_predict[k1,k2,] = get_da(lai_gs_ind,lmf_sm_mean ,lmf_ta_mean ,lmf_sm_ta_sd,lmf_dep,
                                     lmf_total ,sm_mean_c   ,ta_mean_c   ,sm_ta_sdc   ,cor_c  ,      
                                     k1, k2, id, da_len) 
        #col=col+1
        #print(col)
      }
    }
    return(lmf_predict) }
  
  #threshold=0.90; window_yr=20;
  sca = 3
  yr = 2100-1850+1
  len = yr - window_yr
  
  core=36
  cl = makeCluster(core)
  registerDoParallel(cl)
  lmf_contribution = foreach (kk = 1:17, .packages=c('forecast','pracma','segmented')) %dopar% {#kk=1
    out_re = cal_re(rena[kk], lai_gs_indc[kk], id = len, da_len = 7*len)
    out_re
  }
  stopCluster(cl)
  
  setwd('/GPUFS/ygo_fwf_1/00junli/01next_step/10test_has_ant_Ec_LMF')
  saveRDS(lmf_contribution[[01]], paste(out_na, '_',modelna[01], '_', threshold,'_',window_yr,sep=''))
  saveRDS(lmf_contribution[[02]], paste(out_na, '_',modelna[02], '_', threshold,'_',window_yr,sep=''))
  saveRDS(lmf_contribution[[03]], paste(out_na, '_',modelna[03], '_', threshold,'_',window_yr,sep=''))
  saveRDS(lmf_contribution[[04]], paste(out_na, '_',modelna[04], '_', threshold,'_',window_yr,sep=''))
  saveRDS(lmf_contribution[[05]], paste(out_na, '_',modelna[05], '_', threshold,'_',window_yr,sep=''))
  saveRDS(lmf_contribution[[06]], paste(out_na, '_',modelna[06], '_', threshold,'_',window_yr,sep=''))
  saveRDS(lmf_contribution[[07]], paste(out_na, '_',modelna[07], '_', threshold,'_',window_yr,sep=''))
  saveRDS(lmf_contribution[[08]], paste(out_na, '_',modelna[08], '_', threshold,'_',window_yr,sep=''))
  saveRDS(lmf_contribution[[09]], paste(out_na, '_',modelna[09], '_', threshold,'_',window_yr,sep=''))
  saveRDS(lmf_contribution[[10]], paste(out_na, '_',modelna[10], '_', threshold,'_',window_yr,sep=''))
  saveRDS(lmf_contribution[[11]], paste(out_na, '_',modelna[11], '_', threshold,'_',window_yr,sep=''))
  saveRDS(lmf_contribution[[12]], paste(out_na, '_',modelna[12], '_', threshold,'_',window_yr,sep=''))
  saveRDS(lmf_contribution[[13]], paste(out_na, '_',modelna[13], '_', threshold,'_',window_yr,sep=''))
  saveRDS(lmf_contribution[[14]], paste(out_na, '_',modelna[14], '_', threshold,'_',window_yr,sep=''))
  saveRDS(lmf_contribution[[15]], paste(out_na, '_',modelna[15], '_', threshold,'_',window_yr,sep=''))
  saveRDS(lmf_contribution[[16]], paste(out_na, '_',modelna[16], '_', threshold,'_',window_yr,sep=''))
  saveRDS(lmf_contribution[[17]], paste(out_na, '_',modelna[17], '_', threshold,'_',window_yr,sep=''))

  rm(lmf_contribution)
  gc()
  t2 = proc.time()
  print((t2 - t1)[3]/60)
}

#out_na = 'lmf-sm-ta-mean-sm-ta-var-sm-ta-dep-total_'
#threshold = 0.9
#window_yr = 30
#scenario= '_ssp585_'; rena585_0.90_30 = c(paste(out_na, scenario, modelna[01], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[02], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[03], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[04], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[05], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[06], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[07], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[08], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[09], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[10], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[11], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[12], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[13], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[14], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[15], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[16], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[17], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[18], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[19], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[20], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[21], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[22], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[23], '_', threshold, '_', window_yr, sep=''))
#scenario= '_ssp370_'; rena370_0.90_30 = c(paste(out_na, scenario, modelna[01], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[02], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[03], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[04], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[05], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[06], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[07], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[08], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[09], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[10], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[11], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[12], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[13], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[14], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[15], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[16], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[17], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[18], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[19], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[20], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[21], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[22], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[23], '_', threshold, '_', window_yr, sep=''))
#scenario= '_ssp245_'; rena245_0.90_30 = c(paste(out_na, scenario, modelna[01], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[02], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[03], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[04], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[05], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[06], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[07], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[08], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[09], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[10], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[11], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[12], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[13], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[14], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[15], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[16], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[17], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[18], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[19], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[20], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[21], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[22], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[23], '_', threshold, '_', window_yr, sep=''))
#scenario= '_ssp126_'; rena126_0.90_30 = c(paste(out_na, scenario, modelna[01], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[02], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[03], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[04], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[05], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[06], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[07], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[08], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[09], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[10], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[11], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[12], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[13], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[14], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[15], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[16], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[17], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[18], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[19], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[20], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[21], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[22], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[23], '_', threshold, '_', window_yr, sep=''))
#
#scenario= '_ssp585_'; 
#rena585_0.95_30 = c(paste(out_na, scenario, modelna[01], '_', 0.95, '_', 30, sep=''),paste(out_na, scenario, modelna[02], '_', 0.95, '_', 30, sep=''),paste(out_na, scenario, modelna[03], '_', 0.95, '_', 30, sep=''),paste(out_na, scenario, modelna[04], '_', 0.95, '_', 30, sep=''),paste(out_na, scenario, modelna[05], '_', 0.95, '_', 30, sep=''),paste(out_na, scenario, modelna[06], '_', 0.95, '_', 30, sep=''),paste(out_na, scenario, modelna[07], '_', 0.95, '_', 30, sep=''),paste(out_na, scenario, modelna[08], '_', 0.95, '_', 30, sep=''),paste(out_na, scenario, modelna[09], '_', 0.95, '_', 30, sep=''),paste(out_na, scenario, modelna[10], '_', 0.95, '_', 30, sep=''),paste(out_na, scenario, modelna[11], '_', 0.95, '_', 30, sep=''),paste(out_na, scenario, modelna[12], '_', 0.95, '_', 30, sep=''),paste(out_na, scenario, modelna[13], '_', 0.95, '_', 30, sep=''),paste(out_na, scenario, modelna[14], '_', 0.95, '_', 30, sep=''),paste(out_na, scenario, modelna[15], '_', 0.95, '_', 30, sep=''),paste(out_na, scenario, modelna[16], '_', 0.95, '_', 30, sep=''),paste(out_na, scenario, modelna[17], '_', 0.95, '_', 30, sep=''),paste(out_na, scenario, modelna[18], '_', 0.95, '_', 30, sep=''),paste(out_na, scenario, modelna[19], '_', 0.95, '_', 30, sep=''),paste(out_na, scenario, modelna[20], '_', 0.95, '_', 30, sep=''),paste(out_na, scenario, modelna[21], '_', 0.95, '_', 30, sep=''),paste(out_na, scenario, modelna[22], '_', 0.95, '_', 30, sep=''),paste(out_na, scenario, modelna[23], '_', 0.95, '_', 30, sep=''))
#rena585_0.80_30 = c(paste(out_na, scenario, modelna[01], '_', 0.80, '_', 30, sep=''),paste(out_na, scenario, modelna[02], '_', 0.80, '_', 30, sep=''),paste(out_na, scenario, modelna[03], '_', 0.80, '_', 30, sep=''),paste(out_na, scenario, modelna[04], '_', 0.80, '_', 30, sep=''),paste(out_na, scenario, modelna[05], '_', 0.80, '_', 30, sep=''),paste(out_na, scenario, modelna[06], '_', 0.80, '_', 30, sep=''),paste(out_na, scenario, modelna[07], '_', 0.80, '_', 30, sep=''),paste(out_na, scenario, modelna[08], '_', 0.80, '_', 30, sep=''),paste(out_na, scenario, modelna[09], '_', 0.80, '_', 30, sep=''),paste(out_na, scenario, modelna[10], '_', 0.80, '_', 30, sep=''),paste(out_na, scenario, modelna[11], '_', 0.80, '_', 30, sep=''),paste(out_na, scenario, modelna[12], '_', 0.80, '_', 30, sep=''),paste(out_na, scenario, modelna[13], '_', 0.80, '_', 30, sep=''),paste(out_na, scenario, modelna[14], '_', 0.80, '_', 30, sep=''),paste(out_na, scenario, modelna[15], '_', 0.80, '_', 30, sep=''),paste(out_na, scenario, modelna[16], '_', 0.80, '_', 30, sep=''),paste(out_na, scenario, modelna[17], '_', 0.80, '_', 30, sep=''),paste(out_na, scenario, modelna[18], '_', 0.80, '_', 30, sep=''),paste(out_na, scenario, modelna[19], '_', 0.80, '_', 30, sep=''),paste(out_na, scenario, modelna[20], '_', 0.80, '_', 30, sep=''),paste(out_na, scenario, modelna[21], '_', 0.80, '_', 30, sep=''),paste(out_na, scenario, modelna[22], '_', 0.80, '_', 30, sep=''),paste(out_na, scenario, modelna[23], '_', 0.80, '_', 30, sep=''))
#rena585_0.85_30 = c(paste(out_na, scenario, modelna[01], '_', 0.85, '_', 30, sep=''),paste(out_na, scenario, modelna[02], '_', 0.85, '_', 30, sep=''),paste(out_na, scenario, modelna[03], '_', 0.85, '_', 30, sep=''),paste(out_na, scenario, modelna[04], '_', 0.85, '_', 30, sep=''),paste(out_na, scenario, modelna[05], '_', 0.85, '_', 30, sep=''),paste(out_na, scenario, modelna[06], '_', 0.85, '_', 30, sep=''),paste(out_na, scenario, modelna[07], '_', 0.85, '_', 30, sep=''),paste(out_na, scenario, modelna[08], '_', 0.85, '_', 30, sep=''),paste(out_na, scenario, modelna[09], '_', 0.85, '_', 30, sep=''),paste(out_na, scenario, modelna[10], '_', 0.85, '_', 30, sep=''),paste(out_na, scenario, modelna[11], '_', 0.85, '_', 30, sep=''),paste(out_na, scenario, modelna[12], '_', 0.85, '_', 30, sep=''),paste(out_na, scenario, modelna[13], '_', 0.85, '_', 30, sep=''),paste(out_na, scenario, modelna[14], '_', 0.85, '_', 30, sep=''),paste(out_na, scenario, modelna[15], '_', 0.85, '_', 30, sep=''),paste(out_na, scenario, modelna[16], '_', 0.85, '_', 30, sep=''),paste(out_na, scenario, modelna[17], '_', 0.85, '_', 30, sep=''),paste(out_na, scenario, modelna[18], '_', 0.85, '_', 30, sep=''),paste(out_na, scenario, modelna[19], '_', 0.85, '_', 30, sep=''),paste(out_na, scenario, modelna[20], '_', 0.85, '_', 30, sep=''),paste(out_na, scenario, modelna[21], '_', 0.85, '_', 30, sep=''),paste(out_na, scenario, modelna[22], '_', 0.85, '_', 30, sep=''),paste(out_na, scenario, modelna[23], '_', 0.85, '_', 30, sep=''))
#rena585_0.90_20 = c(paste(out_na, scenario, modelna[01], '_', 0.90, '_', 20, sep=''),paste(out_na, scenario, modelna[02], '_', 0.90, '_', 20, sep=''),paste(out_na, scenario, modelna[03], '_', 0.90, '_', 20, sep=''),paste(out_na, scenario, modelna[04], '_', 0.90, '_', 20, sep=''),paste(out_na, scenario, modelna[05], '_', 0.90, '_', 20, sep=''),paste(out_na, scenario, modelna[06], '_', 0.90, '_', 20, sep=''),paste(out_na, scenario, modelna[07], '_', 0.90, '_', 20, sep=''),paste(out_na, scenario, modelna[08], '_', 0.90, '_', 20, sep=''),paste(out_na, scenario, modelna[09], '_', 0.90, '_', 20, sep=''),paste(out_na, scenario, modelna[10], '_', 0.90, '_', 20, sep=''),paste(out_na, scenario, modelna[11], '_', 0.90, '_', 20, sep=''),paste(out_na, scenario, modelna[12], '_', 0.90, '_', 20, sep=''),paste(out_na, scenario, modelna[13], '_', 0.90, '_', 20, sep=''),paste(out_na, scenario, modelna[14], '_', 0.90, '_', 20, sep=''),paste(out_na, scenario, modelna[15], '_', 0.90, '_', 20, sep=''),paste(out_na, scenario, modelna[16], '_', 0.90, '_', 20, sep=''),paste(out_na, scenario, modelna[17], '_', 0.90, '_', 20, sep=''),paste(out_na, scenario, modelna[18], '_', 0.90, '_', 20, sep=''),paste(out_na, scenario, modelna[19], '_', 0.90, '_', 20, sep=''),paste(out_na, scenario, modelna[20], '_', 0.90, '_', 20, sep=''),paste(out_na, scenario, modelna[21], '_', 0.90, '_', 20, sep=''),paste(out_na, scenario, modelna[22], '_', 0.90, '_', 20, sep=''),paste(out_na, scenario, modelna[23], '_', 0.90, '_', 20, sep=''))
#rena585_0.90_40 = c(paste(out_na, scenario, modelna[01], '_', 0.90, '_', 40, sep=''),paste(out_na, scenario, modelna[02], '_', 0.90, '_', 40, sep=''),paste(out_na, scenario, modelna[03], '_', 0.90, '_', 40, sep=''),paste(out_na, scenario, modelna[04], '_', 0.90, '_', 40, sep=''),paste(out_na, scenario, modelna[05], '_', 0.90, '_', 40, sep=''),paste(out_na, scenario, modelna[06], '_', 0.90, '_', 40, sep=''),paste(out_na, scenario, modelna[07], '_', 0.90, '_', 40, sep=''),paste(out_na, scenario, modelna[08], '_', 0.90, '_', 40, sep=''),paste(out_na, scenario, modelna[09], '_', 0.90, '_', 40, sep=''),paste(out_na, scenario, modelna[10], '_', 0.90, '_', 40, sep=''),paste(out_na, scenario, modelna[11], '_', 0.90, '_', 40, sep=''),paste(out_na, scenario, modelna[12], '_', 0.90, '_', 40, sep=''),paste(out_na, scenario, modelna[13], '_', 0.90, '_', 40, sep=''),paste(out_na, scenario, modelna[14], '_', 0.90, '_', 40, sep=''),paste(out_na, scenario, modelna[15], '_', 0.90, '_', 40, sep=''),paste(out_na, scenario, modelna[16], '_', 0.90, '_', 40, sep=''),paste(out_na, scenario, modelna[17], '_', 0.90, '_', 40, sep=''),paste(out_na, scenario, modelna[18], '_', 0.90, '_', 40, sep=''),paste(out_na, scenario, modelna[19], '_', 0.90, '_', 40, sep=''),paste(out_na, scenario, modelna[20], '_', 0.90, '_', 40, sep=''),paste(out_na, scenario, modelna[21], '_', 0.90, '_', 40, sep=''),paste(out_na, scenario, modelna[22], '_', 0.90, '_', 40, sep=''),paste(out_na, scenario, modelna[23], '_', 0.90, '_', 40, sep=''))
#
#outna = '3Month_no_ra_sm1_sm0_pr0_tran0_mrro0_ano_3mon_sen_mean_decompose_own_tas_by_sm_first30yr_win_'
#window_yr = 30
#scenario ='_ssp585_'; tran_indc585_30yr = c(paste(outna, scenario, '_',modelna[01], window_yr,sep=''),paste(outna, scenario, '_',modelna[02], window_yr,sep=''),paste(outna, scenario, '_',modelna[03], window_yr,sep=''),paste(outna, scenario, '_',modelna[04], window_yr,sep=''),paste(outna, scenario, '_',modelna[05], window_yr,sep=''),paste(outna, scenario, '_',modelna[06], window_yr,sep=''),paste(outna, scenario, '_',modelna[07], window_yr,sep=''),paste(outna, scenario, '_',modelna[08], window_yr,sep=''),paste(outna, scenario, '_',modelna[09], window_yr,sep=''),paste(outna, scenario, '_',modelna[10], window_yr,sep=''),paste(outna, scenario, '_',modelna[11], window_yr,sep=''),paste(outna, scenario, '_',modelna[12], window_yr,sep=''),paste(outna, scenario, '_',modelna[13], window_yr,sep=''),paste(outna, scenario, '_',modelna[14], window_yr,sep=''),paste(outna, scenario, '_',modelna[15], window_yr,sep=''),paste(outna, scenario, '_',modelna[16], window_yr,sep=''),paste(outna, scenario, '_',modelna[17], window_yr,sep=''),paste(outna, scenario, '_',modelna[18], window_yr,sep=''),paste(outna, scenario, '_',modelna[19], window_yr,sep=''),paste(outna, scenario, '_',modelna[20], window_yr,sep=''),paste(outna, scenario, '_',modelna[21], window_yr,sep=''),paste(outna, scenario, '_',modelna[22], window_yr,sep=''),paste(outna, scenario, '_',modelna[23], window_yr,sep=''))
#scenario ='_ssp370_'; tran_indc370_30yr = c(paste(outna, scenario, '_',modelna[01], window_yr,sep=''),paste(outna, scenario, '_',modelna[02], window_yr,sep=''),paste(outna, scenario, '_',modelna[03], window_yr,sep=''),paste(outna, scenario, '_',modelna[04], window_yr,sep=''),paste(outna, scenario, '_',modelna[05], window_yr,sep=''),paste(outna, scenario, '_',modelna[06], window_yr,sep=''),paste(outna, scenario, '_',modelna[07], window_yr,sep=''),paste(outna, scenario, '_',modelna[08], window_yr,sep=''),paste(outna, scenario, '_',modelna[09], window_yr,sep=''),paste(outna, scenario, '_',modelna[10], window_yr,sep=''),paste(outna, scenario, '_',modelna[11], window_yr,sep=''),paste(outna, scenario, '_',modelna[12], window_yr,sep=''),paste(outna, scenario, '_',modelna[13], window_yr,sep=''),paste(outna, scenario, '_',modelna[14], window_yr,sep=''),paste(outna, scenario, '_',modelna[15], window_yr,sep=''),paste(outna, scenario, '_',modelna[16], window_yr,sep=''),paste(outna, scenario, '_',modelna[17], window_yr,sep=''),paste(outna, scenario, '_',modelna[18], window_yr,sep=''),paste(outna, scenario, '_',modelna[19], window_yr,sep=''),paste(outna, scenario, '_',modelna[20], window_yr,sep=''),paste(outna, scenario, '_',modelna[21], window_yr,sep=''),paste(outna, scenario, '_',modelna[22], window_yr,sep=''),paste(outna, scenario, '_',modelna[23], window_yr,sep=''))
#scenario ='_ssp245_'; tran_indc245_30yr = c(paste(outna, scenario, '_',modelna[01], window_yr,sep=''),paste(outna, scenario, '_',modelna[02], window_yr,sep=''),paste(outna, scenario, '_',modelna[03], window_yr,sep=''),paste(outna, scenario, '_',modelna[04], window_yr,sep=''),paste(outna, scenario, '_',modelna[05], window_yr,sep=''),paste(outna, scenario, '_',modelna[06], window_yr,sep=''),paste(outna, scenario, '_',modelna[07], window_yr,sep=''),paste(outna, scenario, '_',modelna[08], window_yr,sep=''),paste(outna, scenario, '_',modelna[09], window_yr,sep=''),paste(outna, scenario, '_',modelna[10], window_yr,sep=''),paste(outna, scenario, '_',modelna[11], window_yr,sep=''),paste(outna, scenario, '_',modelna[12], window_yr,sep=''),paste(outna, scenario, '_',modelna[13], window_yr,sep=''),paste(outna, scenario, '_',modelna[14], window_yr,sep=''),paste(outna, scenario, '_',modelna[15], window_yr,sep=''),paste(outna, scenario, '_',modelna[16], window_yr,sep=''),paste(outna, scenario, '_',modelna[17], window_yr,sep=''),paste(outna, scenario, '_',modelna[18], window_yr,sep=''),paste(outna, scenario, '_',modelna[19], window_yr,sep=''),paste(outna, scenario, '_',modelna[20], window_yr,sep=''),paste(outna, scenario, '_',modelna[21], window_yr,sep=''),paste(outna, scenario, '_',modelna[22], window_yr,sep=''),paste(outna, scenario, '_',modelna[23], window_yr,sep=''))
#scenario ='_ssp126_'; tran_indc126_30yr = c(paste(outna, scenario, '_',modelna[01], window_yr,sep=''),paste(outna, scenario, '_',modelna[02], window_yr,sep=''),paste(outna, scenario, '_',modelna[03], window_yr,sep=''),paste(outna, scenario, '_',modelna[04], window_yr,sep=''),paste(outna, scenario, '_',modelna[05], window_yr,sep=''),paste(outna, scenario, '_',modelna[06], window_yr,sep=''),paste(outna, scenario, '_',modelna[07], window_yr,sep=''),paste(outna, scenario, '_',modelna[08], window_yr,sep=''),paste(outna, scenario, '_',modelna[09], window_yr,sep=''),paste(outna, scenario, '_',modelna[10], window_yr,sep=''),paste(outna, scenario, '_',modelna[11], window_yr,sep=''),paste(outna, scenario, '_',modelna[12], window_yr,sep=''),paste(outna, scenario, '_',modelna[13], window_yr,sep=''),paste(outna, scenario, '_',modelna[14], window_yr,sep=''),paste(outna, scenario, '_',modelna[15], window_yr,sep=''),paste(outna, scenario, '_',modelna[16], window_yr,sep=''),paste(outna, scenario, '_',modelna[17], window_yr,sep=''),paste(outna, scenario, '_',modelna[18], window_yr,sep=''),paste(outna, scenario, '_',modelna[19], window_yr,sep=''),paste(outna, scenario, '_',modelna[20], window_yr,sep=''),paste(outna, scenario, '_',modelna[21], window_yr,sep=''),paste(outna, scenario, '_',modelna[22], window_yr,sep=''),paste(outna, scenario, '_',modelna[23], window_yr,sep=''))
#
#window_yr = 20
#scenario ='_ssp585_'; tran_indc585_20yr = c(paste(outna, scenario, '_',modelna[01], window_yr,sep=''),paste(outna, scenario, '_',modelna[02], window_yr,sep=''),paste(outna, scenario, '_',modelna[03], window_yr,sep=''),paste(outna, scenario, '_',modelna[04], window_yr,sep=''),paste(outna, scenario, '_',modelna[05], window_yr,sep=''),paste(outna, scenario, '_',modelna[06], window_yr,sep=''),paste(outna, scenario, '_',modelna[07], window_yr,sep=''),paste(outna, scenario, '_',modelna[08], window_yr,sep=''),paste(outna, scenario, '_',modelna[09], window_yr,sep=''),paste(outna, scenario, '_',modelna[10], window_yr,sep=''),paste(outna, scenario, '_',modelna[11], window_yr,sep=''),paste(outna, scenario, '_',modelna[12], window_yr,sep=''),paste(outna, scenario, '_',modelna[13], window_yr,sep=''),paste(outna, scenario, '_',modelna[14], window_yr,sep=''),paste(outna, scenario, '_',modelna[15], window_yr,sep=''),paste(outna, scenario, '_',modelna[16], window_yr,sep=''),paste(outna, scenario, '_',modelna[17], window_yr,sep=''),paste(outna, scenario, '_',modelna[18], window_yr,sep=''),paste(outna, scenario, '_',modelna[19], window_yr,sep=''),paste(outna, scenario, '_',modelna[20], window_yr,sep=''),paste(outna, scenario, '_',modelna[21], window_yr,sep=''),paste(outna, scenario, '_',modelna[22], window_yr,sep=''),paste(outna, scenario, '_',modelna[23], window_yr,sep=''))
#window_yr = 40
#scenario ='_ssp585_'; tran_indc585_40yr = c(paste(outna, scenario, '_',modelna[01], window_yr,sep=''),paste(outna, scenario, '_',modelna[02], window_yr,sep=''),paste(outna, scenario, '_',modelna[03], window_yr,sep=''),paste(outna, scenario, '_',modelna[04], window_yr,sep=''),paste(outna, scenario, '_',modelna[05], window_yr,sep=''),paste(outna, scenario, '_',modelna[06], window_yr,sep=''),paste(outna, scenario, '_',modelna[07], window_yr,sep=''),paste(outna, scenario, '_',modelna[08], window_yr,sep=''),paste(outna, scenario, '_',modelna[09], window_yr,sep=''),paste(outna, scenario, '_',modelna[10], window_yr,sep=''),paste(outna, scenario, '_',modelna[11], window_yr,sep=''),paste(outna, scenario, '_',modelna[12], window_yr,sep=''),paste(outna, scenario, '_',modelna[13], window_yr,sep=''),paste(outna, scenario, '_',modelna[14], window_yr,sep=''),paste(outna, scenario, '_',modelna[15], window_yr,sep=''),paste(outna, scenario, '_',modelna[16], window_yr,sep=''),paste(outna, scenario, '_',modelna[17], window_yr,sep=''),paste(outna, scenario, '_',modelna[18], window_yr,sep=''),paste(outna, scenario, '_',modelna[19], window_yr,sep=''),paste(outna, scenario, '_',modelna[20], window_yr,sep=''),paste(outna, scenario, '_',modelna[21], window_yr,sep=''),paste(outna, scenario, '_',modelna[22], window_yr,sep=''),paste(outna, scenario, '_',modelna[23], window_yr,sep=''))

#lmf_pre1(rena585_0.90_30, tran_indc585_30yr, modelna, threshold=0.90, window_yr=30, 'Pre_sm_ta_sd_together_no_ra_sen_mean_each_mon_ssp585_')
#lmf_pre1(rena370_0.90_30, tran_indc370_30yr, modelna, threshold=0.90, window_yr=30, 'Pre_sm_ta_sd_together_no_ra_sen_mean_each_mon_ssp370_')
#lmf_pre1(rena245_0.90_30, tran_indc245_30yr, modelna, threshold=0.90, window_yr=30, 'Pre_sm_ta_sd_together_no_ra_sen_mean_each_mon_ssp245_')
#lmf_pre1(rena126_0.90_30, tran_indc126_30yr, modelna, threshold=0.90, window_yr=30, 'Pre_sm_ta_sd_together_no_ra_sen_mean_each_mon_ssp126_')
#
#lmf_pre1(rena585_0.90_20, tran_indc585_20yr, modelna, threshold=0.90, window_yr=20, 'Pre_sm_ta_sd_together_no_ra_sen_mean_each_mon_ssp585_')
#lmf_pre1(rena585_0.90_40, tran_indc585_40yr, modelna, threshold=0.90, window_yr=40, 'Pre_sm_ta_sd_together_no_ra_sen_mean_each_mon_ssp585_')
#
#lmf_pre1(rena585_0.95_30, tran_indc585_30yr, modelna, threshold=0.95, window_yr=30, 'Pre_sm_ta_sd_together_no_ra_sen_mean_each_mon_ssp585_')
#lmf_pre1(rena585_0.80_30, tran_indc585_30yr, modelna, threshold=0.80, window_yr=30, 'Pre_sm_ta_sd_together_no_ra_sen_mean_each_mon_ssp585_')
#lmf_pre1(rena585_0.85_30, tran_indc585_30yr, modelna, threshold=0.85, window_yr=30, 'Pre_sm_ta_sd_together_no_ra_sen_mean_each_mon_ssp585_')

lmf_pre2 = function(rena, lai_gs_indc, modelna, threshold=0.9, window_yr=50, out_na){
  t1 = proc.time()
  get_da = function(lai_gs_ind,lmf_sm_mean ,lmf_ta_mean ,lmf_sm_ta_sd,lmf_dep,
                    lmf_total ,sm_mean_c   ,ta_mean_c   ,sm_ta_sdc   ,cor_c  ,i,j, id=len, da_len=6*len){#i=70;j=78
    
    #id = 221;i=200;j=80
    lmf_sm_mean_ij  = lmf_sm_mean [i,j,]
    lmf_ta_mean_ij  = lmf_ta_mean [i,j,]
    lmf_sm_ta_sd_ij = lmf_sm_ta_sd[i,j,]
    lmf_dep_ij      = lmf_dep     [i,j,]
    sm_mean_c_ij    = sm_mean_c   [i,j,]
    ta_mean_c_ij    = ta_mean_c   [i,j,]
    sm_ta_sdc_ij    = sm_ta_sdc   [i,j,]
    cor_c_ij        = cor_c       [i,j,]
    
    #dasm_mean_ij = cbind(lmf_sm_mean_ij  = lmf_sm_mean_ij , sm_mean_c_ij = sm_mean_c_ij)
    data_mean_ij = cbind(lmf_ta_mean_ij  = lmf_ta_mean_ij , ta_mean_c_ij = ta_mean_c_ij)
    dasm_ta_v_ij = cbind(lmf_sm_ta_sd_ij = lmf_sm_ta_sd_ij, sm_ta_sdc_ij = sm_ta_sdc_ij)
    dadep_ij     = cbind(lmf_dep_ij      = lmf_dep_ij     , cor_c_ij     = cor_c_ij    )
    
    #dasm_mean_ij[is.infinite(dasm_mean_ij)] = NA
    data_mean_ij[is.infinite(data_mean_ij)] = NA
    dasm_ta_v_ij[is.infinite(dasm_ta_v_ij)] = NA
    dadep_ij    [is.infinite(dadep_ij    )] = NA
    
    #dasm_mean_ij = as.data.frame(na.omit(dasm_mean_ij))
    data_mean_ij = as.data.frame(na.omit(data_mean_ij))
    dasm_ta_v_ij = as.data.frame(na.omit(dasm_ta_v_ij))
    dadep_ij     = as.data.frame(na.omit(dadep_ij    ))
    
    #mod_sm_mean_ij = try(segmented(lm(lmf_sm_mean_ij  ~ sm_mean_c_ij, data=dasm_mean_ij), seg.Z = ~sm_mean_c_ij), silent = T)
    mod_ta_mean_ij = try(segmented(lm(lmf_ta_mean_ij  ~ ta_mean_c_ij, data=data_mean_ij), seg.Z = ~ta_mean_c_ij), silent = T)
    mod_sm_ta_v_ij = try(segmented(lm(lmf_sm_ta_sd_ij ~ sm_ta_sdc_ij, data=dasm_ta_v_ij), seg.Z = ~sm_ta_sdc_ij), silent = T)
    mod_dep_ij     = try(segmented(lm(lmf_dep_ij      ~ cor_c_ij    , data=dadep_ij    ), seg.Z = ~cor_c_ij    ), silent = T)
    
    #plot(dasm_mean_ij[,2], dasm_mean_ij[,1], col='orange'); points(dasm_mean_ij[,2], predict(mod_sm_mean_ij));abline(h=0,v=0,lty=2, col='blue')
    #plot(data_mean_ij[,2], data_mean_ij[,1], col='orange'); points(data_mean_ij[,2], predict(mod_ta_mean_ij));abline(h=0,v=0,lty=2, col='blue')
    #plot(dasm_ta_v_ij[,2], dasm_ta_v_ij[,1], col='orange'); points(dasm_ta_v_ij[,2], predict(mod_sm_ta_v_ij));abline(h=0,v=0,lty=2, col='blue')
    #plot(dadep_ij    [,2], dadep_ij    [,1], col='orange'); points(dadep_ij    [,2], predict(mod_dep_ij    ));abline(h=0,v=0,lty=2, col='blue')
    #con1 = sum(!is.na(sm_mean_c_ij))==id & sum(!is.na(lmf_sm_mean_ij )) ==id
    con2 = sum(!is.na(ta_mean_c_ij))==id & sum(!is.na(lmf_ta_mean_ij )) ==id
    con3 = sum(!is.na(sm_ta_sdc_ij))==id & sum(!is.na(lmf_sm_ta_sd_ij)) ==id
    con4 = sum(!is.na(cor_c_ij    ))==id & sum(!is.na(lmf_dep_ij     )) ==id
    
    lai_tas_mean_all1 = lai_gs_ind$lai_tas_mean1_c                                [i,j,]    
    #lai_sm_mean_all1  = lai_gs_ind$lai_sm_mean1_c                                [i,j,]    
    lai_var_all1      = (lai_gs_ind$lai_sm_sd1_c[i,j,1])*(lai_gs_ind$lai_tas_sd1_c[i,j,])                          
    lai_cov_all1      = lai_gs_ind$lai_cov_c1                                     [i,j,]
    lai_var_all0      = (lai_gs_ind$lai_sm_sd0_c[i,j,1])*(lai_gs_ind$lai_tas_sd1_c[i,j,])                             
    lai_cov_all0      = lai_gs_ind$lai_cov_c0                                     [i,j,]                            
    
    lai_tas_mean_c1 = lai_tas_mean_all1[2:(id+1)] - lai_tas_mean_all1[1]
    #lai_sm_mean_c1  = lai_sm_mean_all1 [2:(id+1)] - lai_sm_mean_all1 [1]
    lai_var_c1      = lai_var_all1     [2:(id+1)] - lai_var_all1     [1]
    lai_cov_c1      = lai_cov_all1     [2:(id+1)] - lai_cov_all1     [1]
    lai_var_c0      = lai_var_all0     [2:(id+1)] - lai_var_all0     [1]
    lai_cov_c0      = lai_cov_all0     [2:(id+1)] - lai_cov_all0     [1]
    
    if(#class(mod_sm_mean_ij)[1] != "try-error" & con1 &  #& class(mod_sm_mean_ij)[2] != "try-error"   
      class(mod_ta_mean_ij)[1] != "try-error" & con2 &  #& class(mod_ta_mean_ij)[2] != "try-error"   
      class(mod_sm_ta_v_ij)[1] != "try-error" & con3 &  #& class(mod_ta_mean_ij)[2] != "try-error"   
      class(mod_dep_ij    )[1] != "try-error" & con4){  #& class(mod_dep_ij    )[2] != "try-error"
      
      if(is.na(mean(lai_tas_mean_c1, na.rm=T))){
        out_da = rep(NA, da_len)
      }else{
        ctr_ta_mean = predict(mod_ta_mean_ij, data.frame(ta_mean_c_ij = 0))
        #ctr_sm_mean= predict(mod_sm_mean_ij, data.frame(sm_mean_c_ij = 0))
        ctr_sd      = predict(mod_sm_ta_v_ij, data.frame(sm_ta_sdc_ij = 0))
        ctr_dep     = predict(mod_dep_ij    , data.frame(cor_c_ij     = 0))
        
        out_da = c(predict(mod_ta_mean_ij, data.frame(ta_mean_c_ij  =  lai_tas_mean_c1)) - ctr_ta_mean,
                   rep(NA, da_len/6), #-predict(mod_sm_mean_ij, data.frame(sm_mean_c_ij  =  lai_sm_mean_c1 )) + ctr_sm_mean,
                   (predict(mod_sm_ta_v_ij, data.frame(sm_ta_sdc_ij  =  lai_var_c1     )) - ctr_sd ),
                   (predict(mod_dep_ij    , data.frame(cor_c_ij      = -lai_cov_c1     )) - ctr_dep)#,
                   #(predict(mod_sm_ta_v_ij, data.frame(sm_ta_sdc_ij  =  lai_var_c0     )) - ctr_sd ),
                   #(predict(mod_dep_ij    , data.frame(cor_c_ij      = -lai_cov_c0     )) - ctr_dep)
                   )}
    }else{
      out_da = rep(NA, da_len) }
    return(out_da)
  }
  cal_re = function(rena, lai_gs_indc, id=len, da_len=6*len){
    #id = 221;len=221;rena = rena585_0.90_30[01]; lai_gs_indc = tran_indc585_1mon[01];  id=len; da_len=6*len
    lucc = readRDS('/GPUFS/ygo_fwf_1/00junli/01next_step/lucc2016_1.5degree_vegetated_area.rds')
    setwd('/GPUFS/ygo_fwf_1/00junli/01next_step/7lmf_all_and_ensemble')
    dare  = readRDS(rena)
    
    setwd('/GPUFS/ygo_fwf_1/00junli/01next_step/10results_Albedo_induce_change')
    lai_gs_ind = readRDS(lai_gs_indc)
    
    lmf_sm_mean = dare[,,1        : (id*01)]
    lmf_ta_mean = dare[,,(id*01+1): (id*02)]
    lmf_sm_ta_sd= dare[,,(id*02+1): (id*03)]
    lmf_dep     = dare[,,(id*03+1): (id*04)]
    lmf_total   = dare[,,(id*04+1): (id*05)]
    sm_mean_c   = dare[,,(id*05+1): (id*06)]
    ta_mean_c   = dare[,,(id*06+1): (id*07)]
    sm_ta_sdc   = dare[,,(id*07+1): (id*08)]
    cor_c       = dare[,,(id*08+1): (id*09)]#
    rm(dare)
    gc()
    
    lmf_predict = array(NA, c(dim(sm_mean_c)[1],dim(sm_mean_c)[2],da_len));col=1
    for(k1 in 1:dim(sm_mean_c)[1]) { 
      for(k2 in 1:dim(sm_mean_c)[2]) {  #k1=200;k2=80; 
        if(is.na(lucc[k1,k2]) | is.na(sm_mean_c[k1,k2,50]) | is.na(lmf_sm_mean[k1,k2,50])) next
        
        lmf_predict[k1,k2,] = get_da(lai_gs_ind,lmf_sm_mean ,lmf_ta_mean ,lmf_sm_ta_sd,lmf_dep,
                                     lmf_total ,sm_mean_c   ,ta_mean_c   ,sm_ta_sdc   ,cor_c  ,      
                                     k1, k2, id, da_len) 
        #col=col+1
        #print(col)
      }
    }
    return(lmf_predict) }
  
  #threshold=0.90; window_yr=20;
  sca = 3
  yr = 2100-1850+1
  len = yr - window_yr
  
  core=36
  cl = makeCluster(core)
  registerDoParallel(cl)
  lmf_contribution = foreach (kk = 1:17, .packages=c('forecast','pracma','segmented')) %dopar% {#kk=1
    out_re = cal_re(rena[kk], lai_gs_indc[kk], id = len, da_len = 6*len)
    out_re
  }
  stopCluster(cl)
  
  setwd('/GPUFS/ygo_fwf_1/00junli/01next_step/10results_albedo_prediction')
  saveRDS(lmf_contribution[[01]], paste(out_na, '_',modelna[01], '_', threshold,'_',window_yr,sep=''))
  saveRDS(lmf_contribution[[02]], paste(out_na, '_',modelna[02], '_', threshold,'_',window_yr,sep=''))
  saveRDS(lmf_contribution[[03]], paste(out_na, '_',modelna[03], '_', threshold,'_',window_yr,sep=''))
  saveRDS(lmf_contribution[[04]], paste(out_na, '_',modelna[04], '_', threshold,'_',window_yr,sep=''))
  saveRDS(lmf_contribution[[05]], paste(out_na, '_',modelna[05], '_', threshold,'_',window_yr,sep=''))
  saveRDS(lmf_contribution[[06]], paste(out_na, '_',modelna[06], '_', threshold,'_',window_yr,sep=''))
  saveRDS(lmf_contribution[[07]], paste(out_na, '_',modelna[07], '_', threshold,'_',window_yr,sep=''))
  saveRDS(lmf_contribution[[08]], paste(out_na, '_',modelna[08], '_', threshold,'_',window_yr,sep=''))
  saveRDS(lmf_contribution[[09]], paste(out_na, '_',modelna[09], '_', threshold,'_',window_yr,sep=''))
  saveRDS(lmf_contribution[[10]], paste(out_na, '_',modelna[10], '_', threshold,'_',window_yr,sep=''))
  saveRDS(lmf_contribution[[11]], paste(out_na, '_',modelna[11], '_', threshold,'_',window_yr,sep=''))
  saveRDS(lmf_contribution[[12]], paste(out_na, '_',modelna[12], '_', threshold,'_',window_yr,sep=''))
  saveRDS(lmf_contribution[[13]], paste(out_na, '_',modelna[13], '_', threshold,'_',window_yr,sep=''))
  saveRDS(lmf_contribution[[14]], paste(out_na, '_',modelna[14], '_', threshold,'_',window_yr,sep=''))
  saveRDS(lmf_contribution[[15]], paste(out_na, '_',modelna[15], '_', threshold,'_',window_yr,sep=''))
  saveRDS(lmf_contribution[[16]], paste(out_na, '_',modelna[16], '_', threshold,'_',window_yr,sep=''))
  saveRDS(lmf_contribution[[17]], paste(out_na, '_',modelna[17], '_', threshold,'_',window_yr,sep=''))

  rm(lmf_contribution)
  gc()
  t2 = proc.time()
  print((t2 - t1)[3]/60)
}

#modelna 
#out_na = 'lmf-sm-ta-mean-sm-ta-var-sm-ta-dep-total_'
#threshold = 0.9
#window_yr = 30
#scenario= '_ssp585_'; rena585_0.90_30 = c(paste(out_na, scenario, modelna[01], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[02], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[03], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[04], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[05], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[06], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[07], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[08], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[09], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[10], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[11], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[12], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[13], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[14], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[15], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[16], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[17], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[18], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[19], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[20], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[21], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[22], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[23], '_', threshold, '_', window_yr, sep=''))
#scenario= '_ssp370_'; rena370_0.90_30 = c(paste(out_na, scenario, modelna[01], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[02], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[03], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[04], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[05], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[06], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[07], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[08], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[09], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[10], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[11], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[12], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[13], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[14], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[15], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[16], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[17], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[18], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[19], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[20], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[21], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[22], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[23], '_', threshold, '_', window_yr, sep=''))
#scenario= '_ssp245_'; rena245_0.90_30 = c(paste(out_na, scenario, modelna[01], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[02], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[03], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[04], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[05], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[06], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[07], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[08], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[09], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[10], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[11], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[12], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[13], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[14], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[15], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[16], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[17], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[18], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[19], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[20], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[21], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[22], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[23], '_', threshold, '_', window_yr, sep=''))
#scenario= '_ssp126_'; rena126_0.90_30 = c(paste(out_na, scenario, modelna[01], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[02], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[03], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[04], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[05], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[06], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[07], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[08], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[09], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[10], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[11], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[12], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[13], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[14], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[15], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[16], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[17], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[18], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[19], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[20], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[21], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[22], '_', threshold, '_', window_yr, sep=''),paste(out_na, scenario, modelna[23], '_', threshold, '_', window_yr, sep=''))
#
#scenario= '_ssp585_'; 
#rena585_0.95_30 = c(paste(out_na, scenario, modelna[01], '_', 0.95, '_', 30, sep=''),paste(out_na, scenario, modelna[02], '_', 0.95, '_', 30, sep=''),paste(out_na, scenario, modelna[03], '_', 0.95, '_', 30, sep=''),paste(out_na, scenario, modelna[04], '_', 0.95, '_', 30, sep=''),paste(out_na, scenario, modelna[05], '_', 0.95, '_', 30, sep=''),paste(out_na, scenario, modelna[06], '_', 0.95, '_', 30, sep=''),paste(out_na, scenario, modelna[07], '_', 0.95, '_', 30, sep=''),paste(out_na, scenario, modelna[08], '_', 0.95, '_', 30, sep=''),paste(out_na, scenario, modelna[09], '_', 0.95, '_', 30, sep=''),paste(out_na, scenario, modelna[10], '_', 0.95, '_', 30, sep=''),paste(out_na, scenario, modelna[11], '_', 0.95, '_', 30, sep=''),paste(out_na, scenario, modelna[12], '_', 0.95, '_', 30, sep=''),paste(out_na, scenario, modelna[13], '_', 0.95, '_', 30, sep=''),paste(out_na, scenario, modelna[14], '_', 0.95, '_', 30, sep=''),paste(out_na, scenario, modelna[15], '_', 0.95, '_', 30, sep=''),paste(out_na, scenario, modelna[16], '_', 0.95, '_', 30, sep=''),paste(out_na, scenario, modelna[17], '_', 0.95, '_', 30, sep=''),paste(out_na, scenario, modelna[18], '_', 0.95, '_', 30, sep=''),paste(out_na, scenario, modelna[19], '_', 0.95, '_', 30, sep=''),paste(out_na, scenario, modelna[20], '_', 0.95, '_', 30, sep=''),paste(out_na, scenario, modelna[21], '_', 0.95, '_', 30, sep=''),paste(out_na, scenario, modelna[22], '_', 0.95, '_', 30, sep=''),paste(out_na, scenario, modelna[23], '_', 0.95, '_', 30, sep=''))
#rena585_0.85_30 = c(paste(out_na, scenario, modelna[01], '_', 0.85, '_', 30, sep=''),paste(out_na, scenario, modelna[02], '_', 0.85, '_', 30, sep=''),paste(out_na, scenario, modelna[03], '_', 0.85, '_', 30, sep=''),paste(out_na, scenario, modelna[04], '_', 0.85, '_', 30, sep=''),paste(out_na, scenario, modelna[05], '_', 0.85, '_', 30, sep=''),paste(out_na, scenario, modelna[06], '_', 0.85, '_', 30, sep=''),paste(out_na, scenario, modelna[07], '_', 0.85, '_', 30, sep=''),paste(out_na, scenario, modelna[08], '_', 0.85, '_', 30, sep=''),paste(out_na, scenario, modelna[09], '_', 0.85, '_', 30, sep=''),paste(out_na, scenario, modelna[10], '_', 0.85, '_', 30, sep=''),paste(out_na, scenario, modelna[11], '_', 0.85, '_', 30, sep=''),paste(out_na, scenario, modelna[12], '_', 0.85, '_', 30, sep=''),paste(out_na, scenario, modelna[13], '_', 0.85, '_', 30, sep=''),paste(out_na, scenario, modelna[14], '_', 0.85, '_', 30, sep=''),paste(out_na, scenario, modelna[15], '_', 0.85, '_', 30, sep=''),paste(out_na, scenario, modelna[16], '_', 0.85, '_', 30, sep=''),paste(out_na, scenario, modelna[17], '_', 0.85, '_', 30, sep=''),paste(out_na, scenario, modelna[18], '_', 0.85, '_', 30, sep=''),paste(out_na, scenario, modelna[19], '_', 0.85, '_', 30, sep=''),paste(out_na, scenario, modelna[20], '_', 0.85, '_', 30, sep=''),paste(out_na, scenario, modelna[21], '_', 0.85, '_', 30, sep=''),paste(out_na, scenario, modelna[22], '_', 0.85, '_', 30, sep=''),paste(out_na, scenario, modelna[23], '_', 0.85, '_', 30, sep=''))
#rena585_0.90_20 = c(paste(out_na, scenario, modelna[01], '_', 0.90, '_', 20, sep=''),paste(out_na, scenario, modelna[02], '_', 0.90, '_', 20, sep=''),paste(out_na, scenario, modelna[03], '_', 0.90, '_', 20, sep=''),paste(out_na, scenario, modelna[04], '_', 0.90, '_', 20, sep=''),paste(out_na, scenario, modelna[05], '_', 0.90, '_', 20, sep=''),paste(out_na, scenario, modelna[06], '_', 0.90, '_', 20, sep=''),paste(out_na, scenario, modelna[07], '_', 0.90, '_', 20, sep=''),paste(out_na, scenario, modelna[08], '_', 0.90, '_', 20, sep=''),paste(out_na, scenario, modelna[09], '_', 0.90, '_', 20, sep=''),paste(out_na, scenario, modelna[10], '_', 0.90, '_', 20, sep=''),paste(out_na, scenario, modelna[11], '_', 0.90, '_', 20, sep=''),paste(out_na, scenario, modelna[12], '_', 0.90, '_', 20, sep=''),paste(out_na, scenario, modelna[13], '_', 0.90, '_', 20, sep=''),paste(out_na, scenario, modelna[14], '_', 0.90, '_', 20, sep=''),paste(out_na, scenario, modelna[15], '_', 0.90, '_', 20, sep=''),paste(out_na, scenario, modelna[16], '_', 0.90, '_', 20, sep=''),paste(out_na, scenario, modelna[17], '_', 0.90, '_', 20, sep=''),paste(out_na, scenario, modelna[18], '_', 0.90, '_', 20, sep=''),paste(out_na, scenario, modelna[19], '_', 0.90, '_', 20, sep=''),paste(out_na, scenario, modelna[20], '_', 0.90, '_', 20, sep=''),paste(out_na, scenario, modelna[21], '_', 0.90, '_', 20, sep=''),paste(out_na, scenario, modelna[22], '_', 0.90, '_', 20, sep=''),paste(out_na, scenario, modelna[23], '_', 0.90, '_', 20, sep=''))
#rena585_0.90_40 = c(paste(out_na, scenario, modelna[01], '_', 0.90, '_', 40, sep=''),paste(out_na, scenario, modelna[02], '_', 0.90, '_', 40, sep=''),paste(out_na, scenario, modelna[03], '_', 0.90, '_', 40, sep=''),paste(out_na, scenario, modelna[04], '_', 0.90, '_', 40, sep=''),paste(out_na, scenario, modelna[05], '_', 0.90, '_', 40, sep=''),paste(out_na, scenario, modelna[06], '_', 0.90, '_', 40, sep=''),paste(out_na, scenario, modelna[07], '_', 0.90, '_', 40, sep=''),paste(out_na, scenario, modelna[08], '_', 0.90, '_', 40, sep=''),paste(out_na, scenario, modelna[09], '_', 0.90, '_', 40, sep=''),paste(out_na, scenario, modelna[10], '_', 0.90, '_', 40, sep=''),paste(out_na, scenario, modelna[11], '_', 0.90, '_', 40, sep=''),paste(out_na, scenario, modelna[12], '_', 0.90, '_', 40, sep=''),paste(out_na, scenario, modelna[13], '_', 0.90, '_', 40, sep=''),paste(out_na, scenario, modelna[14], '_', 0.90, '_', 40, sep=''),paste(out_na, scenario, modelna[15], '_', 0.90, '_', 40, sep=''),paste(out_na, scenario, modelna[16], '_', 0.90, '_', 40, sep=''),paste(out_na, scenario, modelna[17], '_', 0.90, '_', 40, sep=''),paste(out_na, scenario, modelna[18], '_', 0.90, '_', 40, sep=''),paste(out_na, scenario, modelna[19], '_', 0.90, '_', 40, sep=''),paste(out_na, scenario, modelna[20], '_', 0.90, '_', 40, sep=''),paste(out_na, scenario, modelna[21], '_', 0.90, '_', 40, sep=''),paste(out_na, scenario, modelna[22], '_', 0.90, '_', 40, sep=''),paste(out_na, scenario, modelna[23], '_', 0.90, '_', 40, sep=''))
#
#
#outna = 'Albedo_induce_change_3mon_'
#window_yr = 30
#scenario ='_ssp585_'; tran_indc585_30yr = c(paste(outna, scenario, '_',modelna[01], window_yr,sep=''),paste(outna, scenario, '_',modelna[02], window_yr,sep=''),paste(outna, scenario, '_',modelna[03], window_yr,sep=''),paste(outna, scenario, '_',modelna[04], window_yr,sep=''),paste(outna, scenario, '_',modelna[05], window_yr,sep=''),paste(outna, scenario, '_',modelna[06], window_yr,sep=''),paste(outna, scenario, '_',modelna[07], window_yr,sep=''),paste(outna, scenario, '_',modelna[08], window_yr,sep=''),paste(outna, scenario, '_',modelna[09], window_yr,sep=''),paste(outna, scenario, '_',modelna[10], window_yr,sep=''),paste(outna, scenario, '_',modelna[11], window_yr,sep=''),paste(outna, scenario, '_',modelna[12], window_yr,sep=''),paste(outna, scenario, '_',modelna[13], window_yr,sep=''),paste(outna, scenario, '_',modelna[14], window_yr,sep=''),paste(outna, scenario, '_',modelna[15], window_yr,sep=''),paste(outna, scenario, '_',modelna[16], window_yr,sep=''),paste(outna, scenario, '_',modelna[17], window_yr,sep=''),paste(outna, scenario, '_',modelna[18], window_yr,sep=''),paste(outna, scenario, '_',modelna[19], window_yr,sep=''),paste(outna, scenario, '_',modelna[20], window_yr,sep=''),paste(outna, scenario, '_',modelna[21], window_yr,sep=''),paste(outna, scenario, '_',modelna[22], window_yr,sep=''),paste(outna, scenario, '_',modelna[23], window_yr,sep=''))
#scenario ='_ssp370_'; tran_indc370_30yr = c(paste(outna, scenario, '_',modelna[01], window_yr,sep=''),paste(outna, scenario, '_',modelna[02], window_yr,sep=''),paste(outna, scenario, '_',modelna[03], window_yr,sep=''),paste(outna, scenario, '_',modelna[04], window_yr,sep=''),paste(outna, scenario, '_',modelna[05], window_yr,sep=''),paste(outna, scenario, '_',modelna[06], window_yr,sep=''),paste(outna, scenario, '_',modelna[07], window_yr,sep=''),paste(outna, scenario, '_',modelna[08], window_yr,sep=''),paste(outna, scenario, '_',modelna[09], window_yr,sep=''),paste(outna, scenario, '_',modelna[10], window_yr,sep=''),paste(outna, scenario, '_',modelna[11], window_yr,sep=''),paste(outna, scenario, '_',modelna[12], window_yr,sep=''),paste(outna, scenario, '_',modelna[13], window_yr,sep=''),paste(outna, scenario, '_',modelna[14], window_yr,sep=''),paste(outna, scenario, '_',modelna[15], window_yr,sep=''),paste(outna, scenario, '_',modelna[16], window_yr,sep=''),paste(outna, scenario, '_',modelna[17], window_yr,sep=''),paste(outna, scenario, '_',modelna[18], window_yr,sep=''),paste(outna, scenario, '_',modelna[19], window_yr,sep=''),paste(outna, scenario, '_',modelna[20], window_yr,sep=''),paste(outna, scenario, '_',modelna[21], window_yr,sep=''),paste(outna, scenario, '_',modelna[22], window_yr,sep=''),paste(outna, scenario, '_',modelna[23], window_yr,sep=''))
#scenario ='_ssp245_'; tran_indc245_30yr = c(paste(outna, scenario, '_',modelna[01], window_yr,sep=''),paste(outna, scenario, '_',modelna[02], window_yr,sep=''),paste(outna, scenario, '_',modelna[03], window_yr,sep=''),paste(outna, scenario, '_',modelna[04], window_yr,sep=''),paste(outna, scenario, '_',modelna[05], window_yr,sep=''),paste(outna, scenario, '_',modelna[06], window_yr,sep=''),paste(outna, scenario, '_',modelna[07], window_yr,sep=''),paste(outna, scenario, '_',modelna[08], window_yr,sep=''),paste(outna, scenario, '_',modelna[09], window_yr,sep=''),paste(outna, scenario, '_',modelna[10], window_yr,sep=''),paste(outna, scenario, '_',modelna[11], window_yr,sep=''),paste(outna, scenario, '_',modelna[12], window_yr,sep=''),paste(outna, scenario, '_',modelna[13], window_yr,sep=''),paste(outna, scenario, '_',modelna[14], window_yr,sep=''),paste(outna, scenario, '_',modelna[15], window_yr,sep=''),paste(outna, scenario, '_',modelna[16], window_yr,sep=''),paste(outna, scenario, '_',modelna[17], window_yr,sep=''),paste(outna, scenario, '_',modelna[18], window_yr,sep=''),paste(outna, scenario, '_',modelna[19], window_yr,sep=''),paste(outna, scenario, '_',modelna[20], window_yr,sep=''),paste(outna, scenario, '_',modelna[21], window_yr,sep=''),paste(outna, scenario, '_',modelna[22], window_yr,sep=''),paste(outna, scenario, '_',modelna[23], window_yr,sep=''))
#scenario ='_ssp126_'; tran_indc126_30yr = c(paste(outna, scenario, '_',modelna[01], window_yr,sep=''),paste(outna, scenario, '_',modelna[02], window_yr,sep=''),paste(outna, scenario, '_',modelna[03], window_yr,sep=''),paste(outna, scenario, '_',modelna[04], window_yr,sep=''),paste(outna, scenario, '_',modelna[05], window_yr,sep=''),paste(outna, scenario, '_',modelna[06], window_yr,sep=''),paste(outna, scenario, '_',modelna[07], window_yr,sep=''),paste(outna, scenario, '_',modelna[08], window_yr,sep=''),paste(outna, scenario, '_',modelna[09], window_yr,sep=''),paste(outna, scenario, '_',modelna[10], window_yr,sep=''),paste(outna, scenario, '_',modelna[11], window_yr,sep=''),paste(outna, scenario, '_',modelna[12], window_yr,sep=''),paste(outna, scenario, '_',modelna[13], window_yr,sep=''),paste(outna, scenario, '_',modelna[14], window_yr,sep=''),paste(outna, scenario, '_',modelna[15], window_yr,sep=''),paste(outna, scenario, '_',modelna[16], window_yr,sep=''),paste(outna, scenario, '_',modelna[17], window_yr,sep=''),paste(outna, scenario, '_',modelna[18], window_yr,sep=''),paste(outna, scenario, '_',modelna[19], window_yr,sep=''),paste(outna, scenario, '_',modelna[20], window_yr,sep=''),paste(outna, scenario, '_',modelna[21], window_yr,sep=''),paste(outna, scenario, '_',modelna[22], window_yr,sep=''),paste(outna, scenario, '_',modelna[23], window_yr,sep=''))
#
#window_yr = 20
#scenario ='_ssp585_'; tran_indc585_20yr = c(paste(outna, scenario, '_',modelna[01], window_yr,sep=''),paste(outna, scenario, '_',modelna[02], window_yr,sep=''),paste(outna, scenario, '_',modelna[03], window_yr,sep=''),paste(outna, scenario, '_',modelna[04], window_yr,sep=''),paste(outna, scenario, '_',modelna[05], window_yr,sep=''),paste(outna, scenario, '_',modelna[06], window_yr,sep=''),paste(outna, scenario, '_',modelna[07], window_yr,sep=''),paste(outna, scenario, '_',modelna[08], window_yr,sep=''),paste(outna, scenario, '_',modelna[09], window_yr,sep=''),paste(outna, scenario, '_',modelna[10], window_yr,sep=''),paste(outna, scenario, '_',modelna[11], window_yr,sep=''),paste(outna, scenario, '_',modelna[12], window_yr,sep=''),paste(outna, scenario, '_',modelna[13], window_yr,sep=''),paste(outna, scenario, '_',modelna[14], window_yr,sep=''),paste(outna, scenario, '_',modelna[15], window_yr,sep=''),paste(outna, scenario, '_',modelna[16], window_yr,sep=''),paste(outna, scenario, '_',modelna[17], window_yr,sep=''),paste(outna, scenario, '_',modelna[18], window_yr,sep=''),paste(outna, scenario, '_',modelna[19], window_yr,sep=''),paste(outna, scenario, '_',modelna[20], window_yr,sep=''),paste(outna, scenario, '_',modelna[21], window_yr,sep=''),paste(outna, scenario, '_',modelna[22], window_yr,sep=''),paste(outna, scenario, '_',modelna[23], window_yr,sep=''))
#window_yr = 40
#scenario ='_ssp585_'; tran_indc585_40yr = c(paste(outna, scenario, '_',modelna[01], window_yr,sep=''),paste(outna, scenario, '_',modelna[02], window_yr,sep=''),paste(outna, scenario, '_',modelna[03], window_yr,sep=''),paste(outna, scenario, '_',modelna[04], window_yr,sep=''),paste(outna, scenario, '_',modelna[05], window_yr,sep=''),paste(outna, scenario, '_',modelna[06], window_yr,sep=''),paste(outna, scenario, '_',modelna[07], window_yr,sep=''),paste(outna, scenario, '_',modelna[08], window_yr,sep=''),paste(outna, scenario, '_',modelna[09], window_yr,sep=''),paste(outna, scenario, '_',modelna[10], window_yr,sep=''),paste(outna, scenario, '_',modelna[11], window_yr,sep=''),paste(outna, scenario, '_',modelna[12], window_yr,sep=''),paste(outna, scenario, '_',modelna[13], window_yr,sep=''),paste(outna, scenario, '_',modelna[14], window_yr,sep=''),paste(outna, scenario, '_',modelna[15], window_yr,sep=''),paste(outna, scenario, '_',modelna[16], window_yr,sep=''),paste(outna, scenario, '_',modelna[17], window_yr,sep=''),paste(outna, scenario, '_',modelna[18], window_yr,sep=''),paste(outna, scenario, '_',modelna[19], window_yr,sep=''),paste(outna, scenario, '_',modelna[20], window_yr,sep=''),paste(outna, scenario, '_',modelna[21], window_yr,sep=''),paste(outna, scenario, '_',modelna[22], window_yr,sep=''),paste(outna, scenario, '_',modelna[23], window_yr,sep=''))


#lmf_pre2(rena585_0.90_30, tran_indc585_30yr, modelna, threshold=0.90, window_yr=30, 'Pre_albedo_no_ra_ssp585_')
#lmf_pre2(rena370_0.90_30, tran_indc370_30yr, modelna, threshold=0.90, window_yr=30, 'Pre_albedo_no_ra_ssp370_')
#lmf_pre2(rena245_0.90_30, tran_indc245_30yr, modelna, threshold=0.90, window_yr=30, 'Pre_albedo_no_ra_ssp245_')
#lmf_pre2(rena126_0.90_30, tran_indc126_30yr, modelna, threshold=0.90, window_yr=30, 'Pre_albedo_no_ra_ssp126_')
#lmf_pre2(rena585_0.90_20, tran_indc585_20yr, modelna, threshold=0.90, window_yr=20, 'Pre_albedo_no_ra_ssp585_')
#lmf_pre2(rena585_0.90_40, tran_indc585_40yr, modelna, threshold=0.90, window_yr=40, 'Pre_albedo_no_ra_ssp585_')

#lmf_pre2(rena585_0.95_30, tran_indc585_30yr, modelna, threshold=0.95, window_yr=30, 'Pre_albedo_no_ra_ssp585_')
#lmf_pre2(rena585_0.85_30, tran_indc585_30yr, modelna, threshold=0.85, window_yr=30, 'Pre_albedo_no_ra_ssp585_')

##
###################Tran alb variability and cov ############################

#use
Tran_ano_plsr = function(scenario, window_yr, core1){
  
  #modelna
  #fpall = '/GPUFS/ygo_fwf_1/00junli/01next_step/GCM_all_data_his_126_245_370_585/all_gcm_da'
  
  get_re_sm_mean_g = function(fpall, modelnai, scenario, window_yr=30){
    
    #window_yr=30; modelnai = modelna[01]; scenario='_ssp585_'
    setwd(fpall)
    fp_file = list.files(fpall)
    for (k in 1:length(fp_file)) {#k=10
      con_tas    = grepl(c('tas'         ), fp_file[k])
      con_mrso_  = grepl(c('mrso_'       ), fp_file[k])
      con_mrro_  = grepl(c('mrro_'       ), fp_file[k])
      con__lai   = grepl(c('lai_L'       ), fp_file[k])
      con_et     = grepl(c('evspsbl_A'   ), fp_file[k])#
      con_tran   = grepl(c('tran_L'      ), fp_file[k])#evspsbl_A
      con_rsds   = grepl(c('rsds'        ), fp_file[k])
      con_rsus   = grepl(c('rsus'        ), fp_file[k])
      con_ts     = grepl(c('ts_A'        ), fp_file[k])
      con_hfss   = grepl(c('hfss_A'      ), fp_file[k])
      con_hurs   = grepl(c('hurs_A'      ), fp_file[k])
      con_wind   = grepl(c('sfcWind_A'   ), fp_file[k])
      con_pr     = grepl(c('pr_A'        ), fp_file[k])
      con_his    = grepl(c('_historical_'), fp_file[k])
      con_ssp585 = grepl(c(scenario      ), fp_file[k])
      con_modelna= grepl(c(modelnai      ), fp_file[k])
      
      if(con_tas   & con_his & con_modelna) {tas__his_bgcna = fp_file[k]}
      if(con_hurs  & con_his & con_modelna) {hurs_hisna     = fp_file[k]}
      if(con_mrso_ & con_his & con_modelna) {mrso_his_bgcna = fp_file[k]}
      if(con_mrro_ & con_his & con_modelna) {mrro_his_bgcna = fp_file[k]}
      if(con__lai  & con_his & con_modelna) {lai__his_bgcna = fp_file[k]}
      if(con_et    & con_his & con_modelna) {et_his_bgcna   = fp_file[k]}
      if(con_tran  & con_his & con_modelna) {tranhis_bgcna  = fp_file[k]}
      if(con_rsds  & con_his & con_modelna) {rsds_his_bgcna = fp_file[k]}
      if(con_rsus  & con_his & con_modelna) {rsus_his_bgcna = fp_file[k]}
      if(con_ts    & con_his & con_modelna) {ts_his_bgcna   = fp_file[k]}
      if(con_hfss  & con_his & con_modelna) {hfss_his_bgcna = fp_file[k]}
      if(con_wind  & con_his & con_modelna) {wind_his_bgcna = fp_file[k]}
      if(con_pr    & con_his & con_modelna) {pr_his_bgcna   = fp_file[k]}
      
      if(con_tas   & con_ssp585 & con_modelna) {tas__ssp_bgcna = fp_file[k]}
      if(con_hurs  & con_ssp585 & con_modelna) {hurs_sspna     = fp_file[k]}
      if(con_mrso_ & con_ssp585 & con_modelna) {mrso_ssp_bgcna = fp_file[k]}
      if(con_mrro_ & con_ssp585 & con_modelna) {mrro_ssp_bgcna = fp_file[k]}
      if(con__lai  & con_ssp585 & con_modelna) {lai__ssp_bgcna = fp_file[k]}
      if(con_et    & con_ssp585 & con_modelna) {et_ssp_bgcna   = fp_file[k]}
      if(con_tran  & con_ssp585 & con_modelna) {transsp_bgcna  = fp_file[k]}
      if(con_rsds  & con_ssp585 & con_modelna) {rsds_ssp_bgcna = fp_file[k]}
      if(con_rsus  & con_ssp585 & con_modelna) {rsus_ssp_bgcna = fp_file[k]}
      if(con_ts    & con_ssp585 & con_modelna) {ts_ssp_bgcna   = fp_file[k]}
      if(con_hfss  & con_ssp585 & con_modelna) {hfss_ssp_bgcna = fp_file[k]}
      if(con_wind  & con_ssp585 & con_modelna) {wind_ssp_bgcna = fp_file[k]}
      if(con_pr    & con_ssp585 & con_modelna) {pr_ssp_bgcna   = fp_file[k]}
    }
    #mrso_his = readRDS(mrso_his_bgcna)
    tas__his = readRDS(tas__his_bgcna)
    tran_his = readRDS(tranhis_bgcna )
    wind_his = readRDS(wind_his_bgcna)
    et___his = readRDS(et_his_bgcna  )
    pr___his = readRDS(pr_his_bgcna  )
    hurs_his = readRDS(hurs_hisna    )
    nrad_his = 2.54*1000000*et___his + readRDS(hfss_his_bgcna) #net radiation (w m-2)
    lai__his = readRDS(lai__his_bgcna)
    
    #mrso_ssp = readRDS(mrso_ssp_bgcna)
    tas__ssp = readRDS(tas__ssp_bgcna)
    tran_ssp = readRDS(transsp_bgcna )
    wind_ssp = readRDS(wind_ssp_bgcna)
    et___ssp = readRDS(et_ssp_bgcna  )
    pr___ssp = readRDS(pr_ssp_bgcna  )
    hurs_ssp = readRDS(hurs_sspna    )
    nrad_ssp = 2.54*1000000*et___ssp + readRDS(hfss_ssp_bgcna) #net radiation (w m-2)
    lai__ssp = readRDS(lai__ssp_bgcna)
    
    #mrso_bgc = abind(mrso_his, mrso_ssp); rm(mrso_his); rm(mrso_ssp); print(dim(mrso_bgc))
    tas__bgc = abind(tas__his, tas__ssp); rm(tas__his); rm(tas__ssp); print(dim(tas__bgc))
    tran_bgc = abind(tran_his, tran_ssp); rm(tran_his); rm(tran_ssp); print(dim(tran_bgc))
    wind_bgc = abind(wind_his, wind_ssp); rm(wind_his); rm(wind_ssp); print(dim(wind_bgc))
    et___bgc = abind(et___his, et___ssp); rm(et___his); rm(et___ssp); print(dim(et___bgc))
    pr___bgc = abind(pr___his, pr___ssp); rm(pr___his); rm(pr___ssp); print(dim(pr___bgc))
    hurs_bgc = abind(hurs_his, hurs_ssp); rm(hurs_his); rm(hurs_ssp); print(dim(hurs_bgc))
    nrad_bgc = abind(nrad_his, nrad_ssp); rm(nrad_his); rm(nrad_ssp); print(dim(nrad_bgc))
    lai__bgc = abind(lai__his, lai__ssp); rm(lai__his); rm(lai__ssp); print(dim(lai__bgc))
    #lai__bgc = round(lai__bgc, 2); lai__bgc[lai__bgc == 0.00] = NA
    
    gs___bgc = tran_bgc/lai__bgc
    gs___bgc[is.infinite(gs___bgc)] = NA
    hurs_bgc[hurs_bgc > 100] = 100 
    vpd__bgc = (610.8*exp(17.29*(tas__bgc-273.15)/((tas__bgc-273.15)+237.3))/1000)*(1-hurs_bgc/100)# (Pa)
    rm(hurs_bgc); gc()
    rm(et___bgc)
    rm(wind_bgc)
    
    sca=3
    yr=2100-1850+1
    len = yr-window_yr+1
    nrow1=240
    ncol1=93
    sen_lai_tran     = array(NA, c(nrow1,ncol1,7))   
    sen_gs_tran      = array(NA, c(nrow1,ncol1,7))   
    sen_vpd_tran     = array(NA, c(nrow1,ncol1,7))   
    sen_pr_tran      = array(NA, c(nrow1,ncol1,7))   
    sen_nrad_tran    = array(NA, c(nrow1,ncol1,7)) 
    model_r2         = array(NA, c(nrow1,ncol1,7)) 
    tran_sd          = array(NA, c(nrow1,ncol1,len))   
    lai__tran_sd     = array(NA, c(nrow1,ncol1,len))   
    gs___tran_sd     = array(NA, c(nrow1,ncol1,len))   
    vpd__tran_sd     = array(NA, c(nrow1,ncol1,len))   
    pr___tran_sd     = array(NA, c(nrow1,ncol1,len))   
    nrad_tran_sd     = array(NA, c(nrow1,ncol1,len))   
    cov_tran_01      = array(NA, c(nrow1,ncol1,len))   
    cov_lai__tran_01 = array(NA, c(nrow1,ncol1,len)) 
    cov_gs___tran_01 = array(NA, c(nrow1,ncol1,len)) 
    cov_vpd__tran_01 = array(NA, c(nrow1,ncol1,len)) 
    cov_pr___tran_01 = array(NA, c(nrow1,ncol1,len)) 
    cov_nrad_tran_01 = array(NA, c(nrow1,ncol1,len)) 
    lucc = readRDS('/GPUFS/ygo_fwf_1/00junli/01next_step/lucc2016_1.5degree_vegetated_area.rds'); col=1
    for (i in 1:nrow(tas__bgc)) {
      for (j in 1:ncol(tas__bgc)) {
        if(is.na(lucc[i,j])  | is.na(tran_bgc[i,j,10]) | is.na(lai__bgc[i,j,10]) |
           is.na(gs___bgc[i,j,10]) | is.na(vpd__bgc[i,j,10]) | is.na(pr___bgc[i,j,10]) |
           is.na(nrad_bgc[i,j,10])) next#i=200;j=80;image.plot(1:240,1:93,tas__bgc[,,8]);abline(v=i,h=j)
        
        tas__ij  = rowMeans(embed(c(tas__bgc[i,j,],rep(NA, sca-1)), sca), na.rm=T)
        loca_bgc = which.max(rowMeans(matrix(tas__ij, ncol=yr), na.rm=T))
        locaij   = sort(embed(c(1:12,1,2),3)[loca_bgc,])
        loca_sel = embed(1:yr, window_yr)
        
        #tran_bgc_ij = try(decompose_own(ts(et___bgc[i,j,]*86400*30, frequency=12), type = "additive"), silent = T)
        tran_bgc_ij = try(decompose_own(ts(tran_bgc[i,j,]*86400*30, frequency=12), type = "additive"), silent = T)
        lai__bgc_ij = try(decompose_own(ts(lai__bgc[i,j,]         , frequency=12), type = "additive"), silent = T)
        gs___bgc_ij = try(decompose_own(ts(gs___bgc[i,j,]         , frequency=12), type = "additive"), silent = T)
        vpd__bgc_ij = try(decompose_own(ts(vpd__bgc[i,j,]         , frequency=12), type = "additive"), silent = T)
        pr___bgc_ij = try(decompose_own(ts(pr___bgc[i,j,]*86400*30, frequency=12), type = "additive"), silent = T)
        nrad_bgc_ij = try(decompose_own(ts(nrad_bgc[i,j,]         , frequency=12), type = "additive"), silent = T)
        
        if(class(tran_bgc_ij)[1] != "try-error" &
           class(lai__bgc_ij)[1] != "try-error" &
           class(gs___bgc_ij)[1] != "try-error" &
           class(vpd__bgc_ij)[1] != "try-error" &
           class(pr___bgc_ij)[1] != "try-error" &
           class(nrad_bgc_ij)[1] != "try-error"){
          
          getre = function(tran_bgc_ij,lai__bgc_ij,gs___bgc_ij,vpd__bgc_ij,pr___bgc_ij,nrad_bgc_ij, id1=09:3008, id2=c(1:6), yr=251){
            return(list(tran_bgc_ijk = matrix(tran_bgc_ij[id1], ncol=yr)[id2,],
                        lai__bgc_ijk = matrix(lai__bgc_ij[id1], ncol=yr)[id2,],
                        gs___bgc_ijk = matrix(gs___bgc_ij[id1], ncol=yr)[id2,],
                        vpd__bgc_ijk = matrix(vpd__bgc_ij[id1], ncol=yr)[id2,],
                        pr___bgc_ijk = matrix(pr___bgc_ij[id1], ncol=yr)[id2,],
                        nrad_bgc_ijk = matrix(nrad_bgc_ij[id1], ncol=yr)[id2,])) }
          
          if(loca_bgc==12){
            id1 = 09:3008; id2 = c(1:6); yr1 = 2100-1850 }else if(loca_bgc==11){
              id1 = 08:3007; id2 = c(1:6); yr1 = 2100-1850 }else if(loca_bgc==1){
                id1 = 10:3009; id2 = c(1:6); yr1 = 2100-1850 }else if(loca_bgc==2){
                  id1 = 11:3010; id2 = c(1:6); yr1 = 2100-1850 }else if(loca_bgc==3){
                    id1 = 12:3011; id2 = c(1:6); yr1 = 2100-1850 }else{
                      id1 = 1:3012; id2 = c((min(locaij)-3):max(locaij)); yr1 = 2100-1850+1
                    }
          
          dare1 = getre(tran_bgc_ij$random,
                        lai__bgc_ij$random,
                        gs___bgc_ij$random,
                        vpd__bgc_ij$random,
                        pr___bgc_ij$random,
                        nrad_bgc_ij$random,id1, id2, yr1)
          dare2 = getre(tran_bgc[i,j,],
                        lai__bgc[i,j,],
                        gs___bgc[i,j,],
                        vpd__bgc[i,j,],
                        pr___bgc[i,j,],
                        nrad_bgc[i,j,], id1, id2, yr1)
          
          
          tran_bgc_ijk = dare1$tran_bgc_ijk
          lai__bgc_ijk = dare1$lai__bgc_ijk
          gs___bgc_ijk = dare1$gs___bgc_ijk
          vpd__bgc_ijk = dare1$vpd__bgc_ijk
          pr___bgc_ijk = dare1$pr___bgc_ijk  
          nrad_bgc_ijk = dare1$nrad_bgc_ijk
          
          da1 = na.omit(data.frame(tran=tran_bgc_ijk[1,],lai=lai__bgc_ijk[1,],gs=gs___bgc_ijk[1,],vpd=vpd__bgc_ijk[1,],pr=pr___bgc_ijk[1,],nrad=nrad_bgc_ijk[1,]))
          da2 = na.omit(data.frame(tran=tran_bgc_ijk[2,],lai=lai__bgc_ijk[2,],gs=gs___bgc_ijk[2,],vpd=vpd__bgc_ijk[2,],pr=pr___bgc_ijk[2,],nrad=nrad_bgc_ijk[2,]))
          da3 = na.omit(data.frame(tran=tran_bgc_ijk[3,],lai=lai__bgc_ijk[3,],gs=gs___bgc_ijk[3,],vpd=vpd__bgc_ijk[3,],pr=pr___bgc_ijk[3,],nrad=nrad_bgc_ijk[3,]))
          da4 = na.omit(data.frame(tran=tran_bgc_ijk[4,],lai=lai__bgc_ijk[4,],gs=gs___bgc_ijk[4,],vpd=vpd__bgc_ijk[4,],pr=pr___bgc_ijk[4,],nrad=nrad_bgc_ijk[4,]))
          da5 = na.omit(data.frame(tran=tran_bgc_ijk[5,],lai=lai__bgc_ijk[5,],gs=gs___bgc_ijk[5,],vpd=vpd__bgc_ijk[5,],pr=pr___bgc_ijk[5,],nrad=nrad_bgc_ijk[5,]))
          da6 = na.omit(data.frame(tran=tran_bgc_ijk[6,],lai=lai__bgc_ijk[6,],gs=gs___bgc_ijk[6,],vpd=vpd__bgc_ijk[6,],pr=pr___bgc_ijk[6,],nrad=nrad_bgc_ijk[6,]))
          daall = na.omit(data.frame(tran = matrix(dare1$tran_bgc_ijk, ncol=1),
                                     lai  = matrix(dare1$lai__bgc_ijk, ncol=1),
                                     gs   = matrix(dare1$gs___bgc_ijk, ncol=1),
                                     vpd  = matrix(dare1$vpd__bgc_ijk, ncol=1),
                                     pr   = matrix(dare1$pr___bgc_ijk, ncol=1),
                                     nrad = matrix(dare1$nrad_bgc_ijk, ncol=1)))
          
          rc1   = try(plsr(tran ~ lai + gs + vpd + pr + nrad, data=da1  ), silent = T)
          rc2   = try(plsr(tran ~ lai + gs + vpd + pr + nrad, data=da2  ), silent = T)
          rc3   = try(plsr(tran ~ lai + gs + vpd + pr + nrad, data=da3  ), silent = T)
          rc4   = try(plsr(tran ~ lai + gs + vpd + pr + nrad, data=da4  ), silent = T)
          rc5   = try(plsr(tran ~ lai + gs + vpd + pr + nrad, data=da5  ), silent = T)
          rc6   = try(plsr(tran ~ lai + gs + vpd + pr + nrad, data=da6  ), silent = T)
          rcall = try(plsr(tran ~ lai + gs + vpd + pr + nrad, data=daall), silent = T)
          
          #npsel_rc1   = try(selectNcomp(rc1  , "onesigma", plot = F), silent = T)
          #npsel_rc2   = try(selectNcomp(rc2  , "onesigma", plot = F), silent = T)
          #npsel_rc3   = try(selectNcomp(rc3  , "onesigma", plot = F), silent = T)
          #npsel_rc4   = try(selectNcomp(rc4  , "onesigma", plot = F), silent = T)
          #npsel_rc5   = try(selectNcomp(rc5  , "onesigma", plot = F), silent = T)
          #npsel_rc6   = try(selectNcomp(rc6  , "onesigma", plot = F), silent = T)
          #npsel_rcall = try(selectNcomp(rcall, "onesigma", plot = F), silent = T)
          
          if(class(rc1  ) != "try-error" &
             class(rc2  ) != "try-error" &
             class(rc3  ) != "try-error" &
             class(rc4  ) != "try-error" &
             class(rc5  ) != "try-error" &
             class(rc6  ) != "try-error" &
             class(rcall) != "try-error" ){
            
            #if(npsel_rc1==0){
            #  rc1 = plsr(tran ~ lai + gs + vpd + pr + nrad, data=da1, ncomp=5        )  }else{
            #    rc1 = plsr(tran ~ lai + gs + vpd + pr + nrad, data=da1, ncomp=npsel_rc1)}
            #if(npsel_rc2==0){
            #  rc2 = plsr(tran ~ lai + gs + vpd + pr + nrad, data=da2, ncomp=5        )  }else{
            #    rc2 = plsr(tran ~ lai + gs + vpd + pr + nrad, data=da2, ncomp=npsel_rc2)}
            #if(npsel_rc3==0){
            #  rc3 = plsr(tran ~ lai + gs + vpd + pr + nrad, data=da3, ncomp=5        )  }else{
            #    rc3 = plsr(tran ~ lai + gs + vpd + pr + nrad, data=da3, ncomp=npsel_rc3)}
            #if(npsel_rc4==0){
            #  rc4 = plsr(tran ~ lai + gs + vpd + pr + nrad, data=da4, ncomp=5        )  }else{
            #    rc4 = plsr(tran ~ lai + gs + vpd + pr + nrad, data=da4, ncomp=npsel_rc4)}
            #if(npsel_rc5==0){
            #  rc5 = plsr(tran ~ lai + gs + vpd + pr + nrad, data=da5, ncomp=5        )  }else{
            #    rc5 = plsr(tran ~ lai + gs + vpd + pr + nrad, data=da5, ncomp=npsel_rc5)}
            #if(npsel_rc6==0){
            #  rc6 = plsr(tran ~ lai + gs + vpd + pr + nrad, data=da6, ncomp=5        )  }else{
            #    rc6 = plsr(tran ~ lai + gs + vpd + pr + nrad, data=da6, ncomp=npsel_rc6)}
            #if(npsel_rcall==0){
            #  rcall = plsr(tran ~ lai + gs + vpd + pr + nrad, data=daall, ncomp=5          )  }else{
            #    rcall = plsr(tran ~ lai + gs + vpd + pr + nrad, data=daall, ncomp=npsel_rcall)}
            
            sen_lai_tran [i,j,] = c(coef(rc1)[1],coef(rc2)[1],coef(rc3)[1],coef(rc4)[1],coef(rc5)[1],coef(rc6)[1],coef(rcall)[1])
            sen_gs_tran  [i,j,] = c(coef(rc1)[2],coef(rc2)[2],coef(rc3)[2],coef(rc4)[2],coef(rc5)[2],coef(rc6)[2],coef(rcall)[2])
            sen_vpd_tran [i,j,] = c(coef(rc1)[3],coef(rc2)[3],coef(rc3)[3],coef(rc4)[3],coef(rc5)[3],coef(rc6)[3],coef(rcall)[3])
            sen_pr_tran  [i,j,] = c(coef(rc1)[4],coef(rc2)[4],coef(rc3)[4],coef(rc4)[4],coef(rc5)[4],coef(rc6)[4],coef(rcall)[4])
            sen_nrad_tran[i,j,] = c(coef(rc1)[5],coef(rc2)[5],coef(rc3)[5],coef(rc4)[5],coef(rc5)[5],coef(rc6)[5],coef(rcall)[5])
            
            
            #rc1___sum = summary(rc1  )
            #rc2___sum = summary(rc2  )
            #rc3___sum = summary(rc3  )
            #rc4___sum = summary(rc4  )
            #rc5___sum = summary(rc5  )
            #rc6___sum = summary(rc6  )
            #rcall_sum = summary(rcall)
            #model_r2[i,j,] #= c(rc1___sum$adj.r.squared,rc2___sum$adj.r.squared,rc3___sum$adj.r.squared,
            #    rc4___sum$adj.r.squared,rc5___sum$adj.r.squared,rc6___sum$adj.r.squared,
            #    rcall_sum$adj.r.squared)
            
            if(ncol(tran_bgc_ijk) != yr){ tran_bgc_ijk = cbind(rep(NA, 6), tran_bgc_ijk) }
            if(ncol(lai__bgc_ijk) != yr){ lai__bgc_ijk = cbind(rep(NA, 6), lai__bgc_ijk) }
            if(ncol(gs___bgc_ijk) != yr){ gs___bgc_ijk = cbind(rep(NA, 6), gs___bgc_ijk) }
            if(ncol(vpd__bgc_ijk) != yr){ vpd__bgc_ijk = cbind(rep(NA, 6), vpd__bgc_ijk) }
            if(ncol(pr___bgc_ijk) != yr){ pr___bgc_ijk = cbind(rep(NA, 6), pr___bgc_ijk) }
            if(ncol(nrad_bgc_ijk) != yr){ nrad_bgc_ijk = cbind(rep(NA, 6), nrad_bgc_ijk) }
            
            lai__tran_ijk = lai__bgc_ijk*c(coef(rc1)[1],coef(rc2)[1],coef(rc3)[1],coef(rc4)[1],coef(rc5)[1],coef(rc6)[1]) 
            gs___tran_ijk = gs___bgc_ijk*c(coef(rc1)[2],coef(rc2)[2],coef(rc3)[2],coef(rc4)[2],coef(rc5)[2],coef(rc6)[2]) 
            vpd__tran_ijk = vpd__bgc_ijk*c(coef(rc1)[3],coef(rc2)[3],coef(rc3)[3],coef(rc4)[3],coef(rc5)[3],coef(rc6)[3]) 
            pr___tran_ijk = pr___bgc_ijk*c(coef(rc1)[4],coef(rc2)[4],coef(rc3)[4],coef(rc4)[4],coef(rc5)[4],coef(rc6)[4]) 
            nrad_tran_ijk = nrad_bgc_ijk*c(coef(rc1)[5],coef(rc2)[5],coef(rc3)[5],coef(rc4)[5],coef(rc5)[5],coef(rc6)[5]) 
            
            tran_sd_ijk      = c()
            lai__tran_sd_ijk = c()
            gs___tran_sd_ijk = c()
            vpd__tran_sd_ijk = c()
            pr___tran_sd_ijk = c()
            nrad_tran_sd_ijk = c()
            for (ij in 1:nrow(loca_sel)) {#ij=222
              tran_sd_ijk      = c(tran_sd_ijk     , sd(tran_bgc_ijk [4:6,loca_sel[ij,]], na.rm=T))
              lai__tran_sd_ijk = c(lai__tran_sd_ijk, sd(lai__tran_ijk[4:6,loca_sel[ij,]], na.rm=T))
              gs___tran_sd_ijk = c(gs___tran_sd_ijk, sd(gs___tran_ijk[4:6,loca_sel[ij,]], na.rm=T))
              vpd__tran_sd_ijk = c(vpd__tran_sd_ijk, sd(vpd__tran_ijk[4:6,loca_sel[ij,]], na.rm=T))
              pr___tran_sd_ijk = c(pr___tran_sd_ijk, sd(pr___tran_ijk[4:6,loca_sel[ij,]], na.rm=T))
              nrad_tran_sd_ijk = c(nrad_tran_sd_ijk, sd(nrad_tran_ijk[4:6,loca_sel[ij,]], na.rm=T))}
            
            tran_bgc_0  = rbind(colMeans(tran_bgc_ijk [1:3,]),colMeans(tran_bgc_ijk [2:4,]),colMeans(tran_bgc_ijk [3:5,]))
            lai__tran_0 = rbind(colMeans(lai__tran_ijk[1:3,]),colMeans(lai__tran_ijk[2:4,]),colMeans(lai__tran_ijk[3:5,]))
            gs___tran_0 = rbind(colMeans(gs___tran_ijk[1:3,]),colMeans(gs___tran_ijk[2:4,]),colMeans(gs___tran_ijk[3:5,]))
            vpd__tran_0 = rbind(colMeans(vpd__tran_ijk[1:3,]),colMeans(vpd__tran_ijk[2:4,]),colMeans(vpd__tran_ijk[3:5,]))
            pr___tran_0 = rbind(colMeans(pr___tran_ijk[1:3,]),colMeans(pr___tran_ijk[2:4,]),colMeans(pr___tran_ijk[3:5,]))
            nrad_tran_0 = rbind(colMeans(nrad_tran_ijk[1:3,]),colMeans(nrad_tran_ijk[2:4,]),colMeans(nrad_tran_ijk[3:5,]))
            
            tran_bgc_1  = rbind(tran_bgc_ijk [4,],tran_bgc_ijk [5,],tran_bgc_ijk [6,])
            lai__tran_1 = rbind(lai__tran_ijk[4,],lai__tran_ijk[5,],lai__tran_ijk[6,])
            gs___tran_1 = rbind(gs___tran_ijk[4,],gs___tran_ijk[5,],gs___tran_ijk[6,])
            vpd__tran_1 = rbind(vpd__tran_ijk[4,],vpd__tran_ijk[5,],vpd__tran_ijk[6,])
            pr___tran_1 = rbind(pr___tran_ijk[4,],pr___tran_ijk[5,],pr___tran_ijk[6,])
            nrad_tran_1 = rbind(nrad_tran_ijk[4,],nrad_tran_ijk[5,],nrad_tran_ijk[6,])
            
            cov_tran_bgc_01_ijk  = c()
            cov_lai__tran_01_ijk = c()
            cov_gs___tran_01_ijk = c()
            cov_vpd__tran_01_ijk = c()
            cov_pr___tran_01_ijk = c()
            cov_nrad_tran_01_ijk = c()
            for (ijk in 1:nrow(loca_sel)) {#ijk=1
              cov_tran_bgc_01_ijk  = c(cov_tran_bgc_01_ijk , cov(na.omit(cbind(matrix(tran_bgc_0 [,loca_sel[ijk,]], ncol=1), matrix(tran_bgc_1 [,loca_sel[ijk,]], ncol=1))))[1,2])
              cov_lai__tran_01_ijk = c(cov_lai__tran_01_ijk, cov(na.omit(cbind(matrix(lai__tran_0[,loca_sel[ijk,]], ncol=1), matrix(lai__tran_1[,loca_sel[ijk,]], ncol=1))))[1,2])
              cov_gs___tran_01_ijk = c(cov_gs___tran_01_ijk, cov(na.omit(cbind(matrix(gs___tran_0[,loca_sel[ijk,]], ncol=1), matrix(gs___tran_1[,loca_sel[ijk,]], ncol=1))))[1,2])
              cov_vpd__tran_01_ijk = c(cov_vpd__tran_01_ijk, cov(na.omit(cbind(matrix(vpd__tran_0[,loca_sel[ijk,]], ncol=1), matrix(vpd__tran_1[,loca_sel[ijk,]], ncol=1))))[1,2])
              cov_pr___tran_01_ijk = c(cov_pr___tran_01_ijk, cov(na.omit(cbind(matrix(pr___tran_0[,loca_sel[ijk,]], ncol=1), matrix(pr___tran_1[,loca_sel[ijk,]], ncol=1))))[1,2])
              cov_nrad_tran_01_ijk = c(cov_nrad_tran_01_ijk, cov(na.omit(cbind(matrix(nrad_tran_0[,loca_sel[ijk,]], ncol=1), matrix(nrad_tran_1[,loca_sel[ijk,]], ncol=1))))[1,2])
            }
            
            #par(mfrow=c(6,2))
            #plot(tran_sd_ijk     , type='o', main='tran_sd_ijk     ', ylab=''); plot(cov_tran_bgc_01_ijk , type='o', main='cov_tran_bgc_01_ijk ')
            #plot(lai__tran_sd_ijk, type='o', main='lai__tran_sd_ijk', ylab=''); plot(cov_lai__tran_01_ijk, type='o', main='cov_lai__tran_01_ijk')
            #plot(gs___tran_sd_ijk, type='o', main='gs___tran_sd_ijk', ylab=''); plot(cov_gs___tran_01_ijk, type='o', main='cov_gs___tran_01_ijk')
            #plot(vpd__tran_sd_ijk, type='o', main='vpd__tran_sd_ijk', ylab=''); plot(cov_vpd__tran_01_ijk, type='o', main='cov_vpd__tran_01_ijk')
            #plot(pr___tran_sd_ijk, type='o', main='pr___tran_sd_ijk', ylab=''); plot(cov_pr___tran_01_ijk, type='o', main='cov_pr___tran_01_ijk')
            #plot(nrad_tran_sd_ijk, type='o', main='nrad_tran_sd_ijk', ylab=''); plot(cov_nrad_tran_01_ijk, type='o', main='cov_nrad_tran_01_ijk')
            
            tran_sd         [i,j,] = tran_sd_ijk
            lai__tran_sd    [i,j,] = lai__tran_sd_ijk 
            gs___tran_sd    [i,j,] = gs___tran_sd_ijk 
            vpd__tran_sd    [i,j,] = vpd__tran_sd_ijk 
            pr___tran_sd    [i,j,] = pr___tran_sd_ijk 
            nrad_tran_sd    [i,j,] = nrad_tran_sd_ijk 
            cov_tran_01     [i,j,] = cov_tran_bgc_01_ijk   
            cov_lai__tran_01[i,j,] = cov_lai__tran_01_ijk    
            cov_gs___tran_01[i,j,] = cov_gs___tran_01_ijk    
            cov_vpd__tran_01[i,j,] = cov_vpd__tran_01_ijk    
            cov_pr___tran_01[i,j,] = cov_pr___tran_01_ijk    
            cov_nrad_tran_01[i,j,] = cov_nrad_tran_01_ijk    
          }
          #col=col+1;print(col)
        }
      }
    }
    
    out = list(sen_lai_tran     = sen_lai_tran    , 
               sen_gs_tran      = sen_gs_tran     ,
               sen_vpd_tran     = sen_vpd_tran    ,
               sen_pr_tran      = sen_pr_tran     ,
               sen_nrad_tran    = sen_nrad_tran   ,
               model_r2         = model_r2        ,
               tran_sd          = tran_sd         ,
               lai__tran_sd     = lai__tran_sd    ,
               gs___tran_sd     = gs___tran_sd    ,
               vpd__tran_sd     = vpd__tran_sd    ,
               pr___tran_sd     = pr___tran_sd    ,
               nrad_tran_sd     = nrad_tran_sd    ,
               cov_tran_01      = cov_tran_01     ,
               cov_lai__tran_01 = cov_lai__tran_01,
               cov_gs___tran_01 = cov_gs___tran_01,
               cov_vpd__tran_01 = cov_vpd__tran_01,
               cov_pr___tran_01 = cov_pr___tran_01,
               cov_nrad_tran_01 = cov_nrad_tran_01)
    return(out)
  }
  
  t1 = proc.time()
  core=core1
  cl = makeCluster(core)
  registerDoParallel(cl)
  lmf_contribution = foreach (kk = 1:17, .packages=c('forecast','pracma','abind','pls')) %dopar% {#kk=1
    out_re = get_re_sm_mean_g(fpall, modelna[kk], scenario, window_yr=30)
    out_re
  }
  stopCluster(cl)
  t2 = proc.time()
  print((t2 - t1)[3]/60)
  
  outna = 'Tran_ano_cov_by_lai_gs_vpd_pr_nrad_gs_tran_lai_each_sen_plsr_'
  setwd('/GPUFS/ygo_fwf_1/00junli/01next_step/10results/tran_ano_explain')
  saveRDS(lmf_contribution[[01]], paste(outna, '_', scenario, '_', window_yr, modelna[01],sep=''))
  saveRDS(lmf_contribution[[02]], paste(outna, '_', scenario, '_', window_yr, modelna[02],sep=''))
  saveRDS(lmf_contribution[[03]], paste(outna, '_', scenario, '_', window_yr, modelna[03],sep=''))
  saveRDS(lmf_contribution[[04]], paste(outna, '_', scenario, '_', window_yr, modelna[04],sep=''))
  saveRDS(lmf_contribution[[05]], paste(outna, '_', scenario, '_', window_yr, modelna[05],sep=''))
  saveRDS(lmf_contribution[[06]], paste(outna, '_', scenario, '_', window_yr, modelna[06],sep=''))
  saveRDS(lmf_contribution[[07]], paste(outna, '_', scenario, '_', window_yr, modelna[07],sep=''))
  saveRDS(lmf_contribution[[08]], paste(outna, '_', scenario, '_', window_yr, modelna[08],sep=''))
  saveRDS(lmf_contribution[[09]], paste(outna, '_', scenario, '_', window_yr, modelna[09],sep=''))
  saveRDS(lmf_contribution[[10]], paste(outna, '_', scenario, '_', window_yr, modelna[10],sep=''))
  saveRDS(lmf_contribution[[11]], paste(outna, '_', scenario, '_', window_yr, modelna[11],sep=''))
  saveRDS(lmf_contribution[[12]], paste(outna, '_', scenario, '_', window_yr, modelna[12],sep=''))
  saveRDS(lmf_contribution[[13]], paste(outna, '_', scenario, '_', window_yr, modelna[13],sep=''))
  saveRDS(lmf_contribution[[14]], paste(outna, '_', scenario, '_', window_yr, modelna[14],sep=''))
  saveRDS(lmf_contribution[[15]], paste(outna, '_', scenario, '_', window_yr, modelna[15],sep=''))
  saveRDS(lmf_contribution[[16]], paste(outna, '_', scenario, '_', window_yr, modelna[16],sep=''))
  saveRDS(lmf_contribution[[17]], paste(outna, '_', scenario, '_', window_yr, modelna[17],sep=''))

  rm(lmf_contribution)
  gc()
}
#Tran_ano_plsr('_ssp585_', window_yr=30, 23)
#Tran_ano_plsr('_ssp370_', window_yr=30, 23)
#Tran_ano_plsr('_ssp245_', window_yr=30, 23)
#Tran_ano_plsr('_ssp126_', window_yr=30, 23)

alb_ano_plsr = function(fpall, modelna, scenario, core1, window_yr){
  
  
 # modelna  need input
  fpall = '/GPUFS/ygo_fwf_1/00junli/01next_step/GCM_all_data_his_126_245_370_585/all_gcm_da'
  
  get_re_sm_mean_g = function(fpall, modelnai, scenario, window_yr=window_yr){
    
    #window_yr=30; modelnai = modelna[01]; scenario='_ssp585_'
    setwd(fpall)
    fp_file = list.files(fpall)
    for (k in 1:length(fp_file)) {#k=10
      con_tas    = grepl(c('tas'         ), fp_file[k])
      con_mrso_  = grepl(c('mrso_'       ), fp_file[k])
      con_mrro_  = grepl(c('mrro_'       ), fp_file[k])
      con__lai   = grepl(c('lai_L'       ), fp_file[k])
      con_et     = grepl(c('evspsbl_A'   ), fp_file[k])#
      con_tran   = grepl(c('tran_L'      ), fp_file[k])#evspsbl_A
      con_rsds   = grepl(c('rsds'        ), fp_file[k])
      con_rsus   = grepl(c('rsus'        ), fp_file[k])
      con_ts     = grepl(c('ts_A'        ), fp_file[k])
      con_hfss   = grepl(c('hfss_A'      ), fp_file[k])
      #con_wind   = grepl(c('sfcWind_A'   ), fp_file[k])
      con_wind   = grepl(c('clt_A'   ), fp_file[k])
      con_pr     = grepl(c('pr_A'        ), fp_file[k])
      con_his    = grepl(c('_historical_'), fp_file[k])
      con_ssp585 = grepl(c(scenario      ), fp_file[k])
      con_modelna= grepl(c(modelnai      ), fp_file[k])
      
      if(con_tas   & con_his & con_modelna) {tas__his_bgcna = fp_file[k]}
      if(con_mrso_ & con_his & con_modelna) {mrso_his_bgcna = fp_file[k]}
      if(con_mrro_ & con_his & con_modelna) {mrro_his_bgcna = fp_file[k]}
      if(con__lai  & con_his & con_modelna) {lai__his_bgcna = fp_file[k]}
      if(con_et    & con_his & con_modelna) {et_his_bgcna   = fp_file[k]}
      if(con_tran  & con_his & con_modelna) {tranhis_bgcna  = fp_file[k]}
      if(con_rsds  & con_his & con_modelna) {rsds_his_bgcna = fp_file[k]}
      if(con_rsus  & con_his & con_modelna) {rsus_his_bgcna = fp_file[k]}
      if(con_ts    & con_his & con_modelna) {ts_his_bgcna   = fp_file[k]}
      if(con_hfss  & con_his & con_modelna) {hfss_his_bgcna = fp_file[k]}
      if(con_wind  & con_his & con_modelna) {wind_his_bgcna = fp_file[k]}
      if(con_pr    & con_his & con_modelna) {pr_his_bgcna   = fp_file[k]}
      
      if(con_tas   & con_ssp585 & con_modelna) {tas__ssp_bgcna = fp_file[k]}
      if(con_mrso_ & con_ssp585 & con_modelna) {mrso_ssp_bgcna = fp_file[k]}
      if(con_mrro_ & con_ssp585 & con_modelna) {mrro_ssp_bgcna = fp_file[k]}
      if(con__lai  & con_ssp585 & con_modelna) {lai__ssp_bgcna = fp_file[k]}
      if(con_et    & con_ssp585 & con_modelna) {et_ssp_bgcna   = fp_file[k]}
      if(con_tran  & con_ssp585 & con_modelna) {transsp_bgcna  = fp_file[k]}
      if(con_rsds  & con_ssp585 & con_modelna) {rsds_ssp_bgcna = fp_file[k]}
      if(con_rsus  & con_ssp585 & con_modelna) {rsus_ssp_bgcna = fp_file[k]}
      if(con_ts    & con_ssp585 & con_modelna) {ts_ssp_bgcna   = fp_file[k]}
      if(con_hfss  & con_ssp585 & con_modelna) {hfss_ssp_bgcna = fp_file[k]}
      if(con_wind  & con_ssp585 & con_modelna) {wind_ssp_bgcna = fp_file[k]}
      if(con_pr    & con_ssp585 & con_modelna) {pr_ssp_bgcna   = fp_file[k]}
    }
    tas__his = readRDS(tas__his_bgcna)
    tran_his = readRDS(tranhis_bgcna )
    #mrso_his = readRDS(mrso_his_bgcna)
    #mrro_his = readRDS(mrro_his_bgcna)
    rsds_his = readRDS(rsds_his_bgcna)
    pr___his = readRDS(pr_his_bgcna  )
    rsus_his = readRDS(rsus_his_bgcna)
    #ts___his = readRDS(ts_his_bgcna  )
    wind_his = readRDS(wind_his_bgcna)
    lai__his = readRDS(lai__his_bgcna)
    
    tas__ssp = readRDS(tas__ssp_bgcna)
    tran_ssp = readRDS(transsp_bgcna )
    #mrso_ssp = readRDS(mrso_ssp_bgcna)
    #mrro_ssp = readRDS(mrro_ssp_bgcna)
    rsds_ssp = readRDS(rsds_ssp_bgcna)
    pr___ssp = readRDS(pr_ssp_bgcna  )
    lai__ssp = readRDS(lai__ssp_bgcna)
    rsus_ssp = readRDS(rsus_ssp_bgcna)
    #ts___ssp = readRDS(ts_ssp_bgcna  )
    wind_ssp = readRDS(wind_ssp_bgcna)
    
    tas__bgc = abind(tas__his, tas__ssp); rm(tas__his); rm(tas__ssp); print(dim(tas__bgc))
    #mrso_bgc = abind(mrso_his, mrso_ssp); rm(mrso_his); rm(mrso_ssp); print(dim(mrso_bgc))
    #mrro_bgc = abind(mrro_his, mrro_ssp); rm(mrro_his); rm(mrro_ssp); print(dim(mrro_bgc))
    tran_bgc = abind(tran_his, tran_ssp); rm(tran_his); rm(tran_ssp); print(dim(tran_bgc))
    rsds_bgc = abind(rsds_his, rsds_ssp); rm(rsds_his); rm(rsds_ssp); print(dim(rsds_bgc))
    pr___bgc = abind(pr___his, pr___ssp); rm(pr___his); rm(pr___ssp); print(dim(pr___bgc))
    rsus_bgc = abind(rsus_his, rsus_ssp); rm(rsus_his); rm(rsus_ssp); print(dim(rsus_bgc))
    #ts___bgc = abind(ts___his, ts___ssp); rm(ts___his); rm(ts___ssp); print(dim(ts___bgc))
    wind_bgc = abind(wind_his, wind_ssp); rm(wind_his); rm(wind_ssp); print(dim(wind_bgc))
    lai__bgc = abind(lai__his, lai__ssp); rm(lai__his); rm(lai__ssp); print(dim(lai__bgc))
    
    albedo_bgc = rsus_bgc/rsds_bgc
    albedo_bgc[is.infinite(albedo_bgc)] = NA
    #ra___bgc = 208/wind_bgc
    #rm(wind_bgc); gc()     #ra___bgc   = 1.2*1013*(ts___bgc - tas__bgc)/hfss_bgc
    rm(rsus_bgc); gc() 
    
    sca=3
    yr=2100-1850+1
    len = yr-window_yr+1
    nrow1=240
    ncol1=93
    alb_tas_mean_c     = array(NA, c(nrow1,ncol1,len))     
    lai_alb_tas_mean_c = array(NA, c(nrow1,ncol1,len))     
    pr__alb_tas_mean_c = array(NA, c(nrow1,ncol1,len))     
    win_alb_tas_mean_c = array(NA, c(nrow1,ncol1,len))     
    lucc = readRDS('/GPUFS/ygo_fwf_1/00junli/01next_step/lucc2016_1.5degree_vegetated_area.rds');col=1
    for (i in 1:nrow(tas__bgc)) {
      for (j in 1:ncol(tas__bgc)) {
        if(is.na(lucc[i,j]) | is.na(tas__bgc[i,j,1]) | 
           is.na(pr___bgc[i,j,1]) ) next
        #i=200;j=80;image.plot(1:240,1:93,tas__bgc[,,8]);abline(v=i,h=j); gc()
        tas__ij  = rowMeans(embed(c(tas__bgc[i,j,],rep(NA, sca-1)), sca), na.rm=T)
        loca_bgc = which.max(rowMeans(matrix(tas__ij, ncol=yr), na.rm=T))
        locaij = sort(embed(c(1:12,1,2),3)[loca_bgc,])
        loca_sel = embed(1:yr, window_yr)
        
        #mrso_ij    = try(decompose_own(ts(mrso_bgc  [i,j,]         , frequency=12), type = "additive", window_yr = window_yr), silent = T)
        #mrro_ij    = try(decompose_own(ts(mrro_bgc  [i,j,]*86400*30, frequency=12), type = "additive", window_yr = window_yr), silent = T)
        #tran_ij_sm = try(decompose_own(ts(tran_bgc  [i,j,]*86400*30, frequency=12), type = "additive", window_yr = window_yr), silent = T)
        #pr___ij_sm = try(decompose_own(ts(pr___bgc  [i,j,]*86400*30, frequency=12), type = "additive", window_yr = window_yr), silent = T)
        #tran_ij_tas= try(decompose_own(ts(tran_bgc  [i,j,]         , frequency=12), type = "additive", window_yr = window_yr), silent = T)
        albe_ij_tas= try(decompose_own(ts(albedo_bgc[i,j,]         , frequency=12), type = "additive", window_yr = window_yr), silent = T)
        
        tranij_raw  = tran_bgc[i,j,]
        tasij_raw   = tas__bgc[i,j,]
        sw_ij_raw   = rsds_bgc[i,j,]
        
        daalbe = matrix(albedo_bgc[i,j,], ncol=yr)
        dalai_ = matrix(lai__bgc  [i,j,], ncol=yr)
        dapr__ = matrix(pr___bgc  [i,j,], ncol=yr)
        dawind = matrix(wind_bgc  [i,j,], ncol=yr)
        
      
        cal_att = function(daalbe,dalai_,dapr__,dawind){#k=7
          
          alblai_ = matrix(NA, nrow=12, ncol=251)
          albpr__ = matrix(NA, nrow=12, ncol=251)
          albwind = matrix(NA, nrow=12, ncol=251)
          albfit  = matrix(NA, nrow=12, ncol=251)
          for (k in 1:12) {
            rc=try(plsr(daalbe[k,]~
                        dalai_[k,]+
                        dapr__[k,]+
                        dawind[k,]), silent = T)
            if(class(rc)!= "try-error"){
              if(length(rc$fitted.values) == 251){
                #print(k)
                alblai_[k,] = dalai_[k,]*coef(rc)[1]
                albpr__[k,] = dapr__[k,]*coef(rc)[2]
                albwind[k,] = dawind[k,]*coef(rc)[3]
                albfit [k,] = rc$fitted.values
              }
            }
          }
          #length(rc$fitted.values)
          #as.character(rc$fitted.values)
          return(list(alblai_ = alblai_,
                      albpr__ = albpr__,
                      albwind = albwind,
                      albfit  = albfit))
        }
        re = cal_att(daalbe,dalai_,dapr__,dawind)
        
        da_albedo  = try(decompose_own(ts(albedo_bgc          [i,j,], frequency = 12), type = "additive", window_yr=30), silent = T)
        da_alblai_ = try(decompose_own(ts(matrix(re$alblai_, ncol=1), frequency = 12), type = "additive", window_yr=30), silent = T)
        da_albpr__ = try(decompose_own(ts(matrix(re$albpr__, ncol=1), frequency = 12), type = "additive", window_yr=30), silent = T)
        da_albwind = try(decompose_own(ts(matrix(re$albwind, ncol=1), frequency = 12), type = "additive", window_yr=30), silent = T)
        
        if(class(da_albedo  )[1] != "try-error" &
           class(da_alblai_ )[1] != "try-error" &
           class(da_albpr__ )[1] != "try-error" &
           class(da_albwind )[1] != "try-error" &
           class(albe_ij_tas)[1] != "try-error"){
          
          getre = function(albe_ij_tas,tasij_raw, sw_ij_raw, albe_ij_raw, id1=09:3008, id2=c(1:6), yr=251){
            return(list(albe_ij_tas  = matrix(albe_ij_tas[id1], ncol=yr)[id2,],
                        #tranij_raw   = matrix(tranij_raw [id1], ncol=yr)[id2,],
                        tasij_raw    = matrix(tasij_raw  [id1], ncol=yr)[id2,],
                        sw_ij_raw    = matrix(sw_ij_raw  [id1], ncol=yr)[id2,],
                        albe_ij_raw  = matrix(albe_ij_raw[id1], ncol=yr)[id2,])) }
          
          if(loca_bgc==12){
            id1 = 09:3008; id2 = c(1:6); yr1 = 2100-1850 }else if(loca_bgc==11){
              id1 = 08:3007; id2 = c(1:6); yr1 = 2100-1850 }else if(loca_bgc==1){
                id1 = 10:3009; id2 = c(1:6); yr1 = 2100-1850 }else if(loca_bgc==2){
                  id1 = 11:3010; id2 = c(1:6); yr1 = 2100-1850 }else if(loca_bgc==3){
                    id1 = 12:3011; id2 = c(1:6); yr1 = 2100-1850 }else{
                      id1 = 1:3012; id2 = c((min(locaij)-3):max(locaij)); yr1 = 2100-1850+1
                    }
          
          #tranij_raw_1  = tran_ij_tas$trend
          albe_ij_tas_1 = albe_ij_tas$random
          albe_ij_raw   = albe_ij_tas$trend
          
          dare1 = getre(albe_ij_tas_1,tasij_raw, sw_ij_raw, albe_ij_raw, id1, id2, yr1)
          
          albe_ijk_tas_sd = dare1$albe_ij_tas
          tas_raw_ijk     = dare1$tasij_raw 
          sw__raw_ijk     = dare1$sw_ij_raw 
          albe_ij_raw     = dare1$albe_ij_raw 
          
          ####tas mean####
          rsdsctr = mean(sw__raw_ijk[4:6,],na.rm=T)
          tas_ctr = mean(tas_raw_ijk[4:6,1:window_yr],na.rm=T)
          
          f = 1/(4*0.95*5.67*10^-8*tas_ctr^3)
          if(ncol(albe_ij_raw)!=(2100-1850+1)){ albe_ij_raw = cbind(rep(NA,6), albe_ij_raw) }
          
          lai_tas = albe_ij_raw[4:6,]*f*(-rsdsctr)
          alb_tas_mean_c[i,j,] = rowMeans(embed(colMeans(lai_tas), window_yr), na.rm=T)
          #####lai tr wind#######################
          getalbedo = function(albe_bgc,lai__bgc,pr___bgc,wind_bgc,id1=09:3008, id2=c(1:6), yr=251){
            return(list(albe_ijk = matrix(albe_bgc[id1], ncol=yr)[id2,],
                        lai__ijk = matrix(lai__bgc[id1], ncol=yr)[id2,], 
                        pr___ijk = matrix(pr___bgc[id1], ncol=yr)[id2,],
                        wind_ijk = matrix(wind_bgc[id1], ncol=yr)[id2,])) }
          
          albedo = getalbedo(da_albedo$trend,
                             da_alblai_$trend,
                             da_albpr__$trend,
                             da_albwind$trend, id1, id2, yr1)
          
          albe_ijk = albedo$albe_ijk
          lai__ijk = albedo$lai__ijk
          pr___ijk = albedo$pr___ijk
          wind_ijk = albedo$wind_ijk
          
          if(ncol(albe_ijk) != yr){ albe_ijk = cbind(rep(NA, 6), albe_ijk)}
          if(ncol(lai__ijk) != yr){ lai__ijk = cbind(rep(NA, 6), lai__ijk)}
          if(ncol(pr___ijk) != yr){ pr___ijk = cbind(rep(NA, 6), pr___ijk)}
          if(ncol(wind_ijk) != yr){ wind_ijk = cbind(rep(NA, 6), wind_ijk)}
          
          lai__alb_tas = lai__ijk[4:6,]*f*(-rsdsctr)
          pr___alb_tas = pr___ijk[4:6,]*f*(-rsdsctr)
          wind_alb_tas = wind_ijk[4:6,]*f*(-rsdsctr)
          
          lai_alb_tas_mean_c[i,j,] = rowMeans(embed(colMeans(lai__alb_tas), window_yr), na.rm=T)
          pr__alb_tas_mean_c[i,j,] = rowMeans(embed(colMeans(pr___alb_tas), window_yr), na.rm=T)
          win_alb_tas_mean_c[i,j,] = rowMeans(embed(colMeans(wind_alb_tas), window_yr), na.rm=T)
        }
      }
    }
    
    out = list(alb_tas_mean_c     = alb_tas_mean_c    ,
               lai_alb_tas_mean_c = lai_alb_tas_mean_c,
               pr__alb_tas_mean_c = pr__alb_tas_mean_c,
               win_alb_tas_mean_c = win_alb_tas_mean_c)
    return(out)
  }
  
  t1 = proc.time()
  core=core1
  cl = makeCluster(core)
  registerDoParallel(cl)
  lmf_contribution = foreach (kk = 1:17, .packages=c('forecast','pracma','abind','pls')) %dopar% {#kk=1
    out_re = get_re_sm_mean_g(fpall, modelna[kk], scenario, window_yr=window_yr)
    out_re
  }
  stopCluster(cl)
  t2 = proc.time()
  print((t2 - t1)[3]/60)
  
  
  modelna = c("BCC-CSM2-MR","_CanESM5_","_CanESM5-1_","_CanESM5-CanOE_",
              "_CESM2_","CESM2-WACCM","CMCC-CM2-SR5","CMCC-ESM2","CNRM-CM6-1",
              "CNRM-ESM2-1","GISS-E2-1-G","GISS-E2-1-H","GISS-E2-2-G",
              "IPSL-CM6A-LR","MIROC-ES2L","MIROC-ES2H","MPI-ESM1-2-LR",
              "NorESM2-LM","NorESM2-MM","TaiESM1","UKESM1-0-L",'CAS-ESM2-0', 'EC-Earth3-Veg')
  outna = 'LAI_Pr_Wind_Albedo_induce_change_3mon_'
  setwd('/GPUFS/ygo_fwf_1/00junli/01next_step/10results_Albedo_induce_change')
  saveRDS(lmf_contribution[[01]], paste(outna, scenario, '_',modelna[01], window_yr,sep=''))
  saveRDS(lmf_contribution[[02]], paste(outna, scenario, '_',modelna[02], window_yr,sep=''))
  saveRDS(lmf_contribution[[03]], paste(outna, scenario, '_',modelna[03], window_yr,sep=''))
  saveRDS(lmf_contribution[[04]], paste(outna, scenario, '_',modelna[04], window_yr,sep=''))
  saveRDS(lmf_contribution[[05]], paste(outna, scenario, '_',modelna[05], window_yr,sep=''))
  saveRDS(lmf_contribution[[06]], paste(outna, scenario, '_',modelna[06], window_yr,sep=''))
  saveRDS(lmf_contribution[[07]], paste(outna, scenario, '_',modelna[07], window_yr,sep=''))
  saveRDS(lmf_contribution[[08]], paste(outna, scenario, '_',modelna[08], window_yr,sep=''))
  saveRDS(lmf_contribution[[09]], paste(outna, scenario, '_',modelna[09], window_yr,sep=''))
  saveRDS(lmf_contribution[[10]], paste(outna, scenario, '_',modelna[10], window_yr,sep=''))
  saveRDS(lmf_contribution[[11]], paste(outna, scenario, '_',modelna[11], window_yr,sep=''))
  saveRDS(lmf_contribution[[12]], paste(outna, scenario, '_',modelna[12], window_yr,sep=''))
  saveRDS(lmf_contribution[[13]], paste(outna, scenario, '_',modelna[13], window_yr,sep=''))
  saveRDS(lmf_contribution[[14]], paste(outna, scenario, '_',modelna[14], window_yr,sep=''))
  saveRDS(lmf_contribution[[15]], paste(outna, scenario, '_',modelna[15], window_yr,sep=''))
  saveRDS(lmf_contribution[[16]], paste(outna, scenario, '_',modelna[16], window_yr,sep=''))
  saveRDS(lmf_contribution[[17]], paste(outna, scenario, '_',modelna[17], window_yr,sep=''))
  saveRDS(lmf_contribution[[18]], paste(outna, scenario, '_',modelna[18], window_yr,sep=''))
  saveRDS(lmf_contribution[[19]], paste(outna, scenario, '_',modelna[19], window_yr,sep=''))
  saveRDS(lmf_contribution[[20]], paste(outna, scenario, '_',modelna[20], window_yr,sep=''))
  saveRDS(lmf_contribution[[21]], paste(outna, scenario, '_',modelna[21], window_yr,sep=''))
  saveRDS(lmf_contribution[[22]], paste(outna, scenario, '_',modelna[22], window_yr,sep=''))
  saveRDS(lmf_contribution[[23]], paste(outna, scenario, '_',modelna[23], window_yr,sep=''))
  
  rm(lmf_contribution)
  gc()
}
#alb_ano_plsr(fpall, modelna, '_ssp585_', 23, window_yr=30)
#alb_ano_plsr(fpall, modelna, '_ssp370_', 23, window_yr=30)
#alb_ano_plsr(fpall, modelna, '_ssp245_', 23, window_yr=30)
#alb_ano_plsr(fpall, modelna, '_ssp126_', 23, window_yr=30)

##
###################Tran alb variability and cov 3 co2 sensitive######################

tran_ano = function(core1){
  
  #modelna 
  #modelnai = modelna[01];window_yr=30
  get_re_sm_mean_g = function(modelnai, window_yr=30){
    
    fphis = '/GPUFS/ygo_fwf_1/00junli/01next_step/Ctr_bgc_rad_data'
    setwd(fphis)
    fp_file = list.files(fphis)
    for (k in 1:length(fp_file)) {#k=10
      con_tas   = grepl(c('tas'        ), fp_file[k])
      con__lai  = grepl(c('lai_L'      ), fp_file[k])
      con_et    = grepl(c('evspsbl_A'  ), fp_file[k])#
      con_tran  = grepl(c('tran_L'     ), fp_file[k])#evspsbl_A
      con_hfss  = grepl(c('hfss_A'     ), fp_file[k])
      con_pr    = grepl(c('pr_A'       ), fp_file[k])
      con_hurs  = grepl(c('hurs_A'     ), fp_file[k])
      con_mod   = grepl(c(modelnai     ), fp_file[k])
      con_ctr   = grepl(c('1pctCO2_'   ), fp_file[k])
      con_bgc   = grepl(c('1pctCO2-bgc'), fp_file[k])
      con_rad   = grepl(c('1pctCO2-rad'), fp_file[k])
      
      if(con_tas  & con_bgc & con_mod) {tas__his_bgcna = fp_file[k]}
      if(con__lai & con_bgc & con_mod) {lai__his_bgcna = fp_file[k]}
      if(con_et   & con_bgc & con_mod) {et_his_bgcna   = fp_file[k]}
      if(con_tran & con_bgc & con_mod) {tranhis_bgcna  = fp_file[k]}
      if(con_hfss & con_bgc & con_mod) {hfss_his_bgcna = fp_file[k]}
      if(con_pr   & con_bgc & con_mod) {pr_his_bgcna   = fp_file[k]}
      if(con_hurs & con_bgc & con_mod) {hurs_his_bgcna = fp_file[k]}
      
      if(con_tas  & con_rad & con_mod) {tas__his_radna = fp_file[k]}
      if(con__lai & con_rad & con_mod) {lai__his_radna = fp_file[k]}
      if(con_et   & con_rad & con_mod) {et_his_radna   = fp_file[k]}
      if(con_tran & con_rad & con_mod) {tranhis_radna  = fp_file[k]}
      if(con_hfss & con_rad & con_mod) {hfss_his_radna = fp_file[k]}
      if(con_pr   & con_rad & con_mod) {pr_his_radna   = fp_file[k]}
      if(con_hurs & con_rad & con_mod) {hurs_his_radna = fp_file[k]}
      
      if(con_tas  & con_ctr & con_mod) {tas__his_ctrna = fp_file[k]}
      if(con__lai & con_ctr & con_mod) {lai__his_ctrna = fp_file[k]}
      if(con_et   & con_ctr & con_mod) {et_his_ctrna   = fp_file[k]}
      if(con_tran & con_ctr & con_mod) {tranhis_ctrna  = fp_file[k]}
      if(con_hfss & con_ctr & con_mod) {hfss_his_ctrna = fp_file[k]}
      if(con_pr   & con_ctr & con_mod) {pr_his_ctrna   = fp_file[k]}
      if(con_hurs & con_ctr & con_mod) {hurs_his_ctrna = fp_file[k]}
    }
    tas__bgc = readRDS(tas__his_bgcna)[,,1:1680]
    lai__bgc = readRDS(lai__his_bgcna)[,,1:1680]
    et___bgc = readRDS(et_his_bgcna  )[,,1:1680]
    tran_bgc = readRDS(tranhis_bgcna )[,,1:1680]
    nrad_bgc = readRDS(hfss_his_bgcna)[,,1:1680] + 2.54*1000000*et___bgc#net radiation (w m-2)
    pr___bgc = readRDS(pr_his_bgcna  )[,,1:1680] 
    hurs_bgc = readRDS(hurs_his_bgcna)[,,1:1680]
    
    tas__rad = readRDS(tas__his_radna)[,,1:1680]
    lai__rad = readRDS(lai__his_radna)[,,1:1680]
    et___rad = readRDS(et_his_radna  )[,,1:1680]
    tran_rad = readRDS(tranhis_radna )[,,1:1680]
    nrad_rad = readRDS(hfss_his_radna)[,,1:1680] + 2.54*1000000*et___rad#net radiation (w m-2)
    pr___rad = readRDS(pr_his_radna  )[,,1:1680] 
    hurs_rad = readRDS(hurs_his_radna)[,,1:1680]
    
    tas__ctr = readRDS(tas__his_ctrna)[,,1:1680]
    lai__ctr = readRDS(lai__his_ctrna)[,,1:1680]
    et___ctr = readRDS(et_his_ctrna  )[,,1:1680]
    tran_ctr = readRDS(tranhis_ctrna )[,,1:1680]
    nrad_ctr = readRDS(hfss_his_ctrna)[,,1:1680] + 2.54*1000000*et___ctr#net radiation (w m-2)
    pr___ctr = readRDS(pr_his_ctrna  )[,,1:1680] 
    hurs_ctr = readRDS(hurs_his_ctrna)[,,1:1680]
    
    gs___ctr = tran_ctr/lai__ctr; gs___ctr[is.infinite(gs___ctr)] = NA
    gs___rad = tran_rad/lai__rad; gs___rad[is.infinite(gs___rad)] = NA
    gs___bgc = tran_bgc/lai__bgc; gs___bgc[is.infinite(gs___bgc)] = NA
    
    hurs_ctr[hurs_ctr > 100] = 100 
    hurs_rad[hurs_rad > 100] = 100 
    hurs_bgc[hurs_bgc > 100] = 100 
    vpd__ctr = (610.8*exp(17.29*(tas__ctr-273.15)/((tas__ctr-273.15)+237.3))/1000)*(1-hurs_ctr/100)# (Pa)
    vpd__rad = (610.8*exp(17.29*(tas__rad-273.15)/((tas__rad-273.15)+237.3))/1000)*(1-hurs_rad/100)# (Pa)
    vpd__bgc = (610.8*exp(17.29*(tas__bgc-273.15)/((tas__bgc-273.15)+237.3))/1000)*(1-hurs_bgc/100)# (Pa)
    rm(hurs_ctr); gc()
    rm(hurs_rad); gc()
    rm(hurs_bgc); gc()
    rm(et___ctr)
    rm(et___rad)
    rm(et___bgc)
    
    sca=3
    yr=140
    len = yr-window_yr+1
    nrow1=240
    ncol1=93
    tran_sd_bgc          = matrix(NA, nrow=240, ncol=93)
    lai__tran_sd_bgc     = matrix(NA, nrow=240, ncol=93)
    gs___tran_sd_bgc     = matrix(NA, nrow=240, ncol=93)
    vpd__tran_sd_bgc     = matrix(NA, nrow=240, ncol=93)
    pr___tran_sd_bgc     = matrix(NA, nrow=240, ncol=93)
    nrad_tran_sd_bgc     = matrix(NA, nrow=240, ncol=93)
    cov_tran_01_bgc      = matrix(NA, nrow=240, ncol=93)
    cov_lai__tran_01_bgc = matrix(NA, nrow=240, ncol=93)
    cov_gs___tran_01_bgc = matrix(NA, nrow=240, ncol=93)
    cov_vpd__tran_01_bgc = matrix(NA, nrow=240, ncol=93)
    cov_pr___tran_01_bgc = matrix(NA, nrow=240, ncol=93)
    cov_nrad_tran_01_bgc = matrix(NA, nrow=240, ncol=93)
    tran_sd_rad          = matrix(NA, nrow=240, ncol=93)
    lai__tran_sd_rad     = matrix(NA, nrow=240, ncol=93)
    gs___tran_sd_rad     = matrix(NA, nrow=240, ncol=93)
    vpd__tran_sd_rad     = matrix(NA, nrow=240, ncol=93)
    pr___tran_sd_rad     = matrix(NA, nrow=240, ncol=93)
    nrad_tran_sd_rad     = matrix(NA, nrow=240, ncol=93)
    cov_tran_01_rad      = matrix(NA, nrow=240, ncol=93)
    cov_lai__tran_01_rad = matrix(NA, nrow=240, ncol=93)
    cov_gs___tran_01_rad = matrix(NA, nrow=240, ncol=93)
    cov_vpd__tran_01_rad = matrix(NA, nrow=240, ncol=93)
    cov_pr___tran_01_rad = matrix(NA, nrow=240, ncol=93)
    cov_nrad_tran_01_rad = matrix(NA, nrow=240, ncol=93)
    tran_sd_ctr          = matrix(NA, nrow=240, ncol=93)
    lai__tran_sd_ctr     = matrix(NA, nrow=240, ncol=93)
    gs___tran_sd_ctr     = matrix(NA, nrow=240, ncol=93)
    vpd__tran_sd_ctr     = matrix(NA, nrow=240, ncol=93)
    pr___tran_sd_ctr     = matrix(NA, nrow=240, ncol=93)
    nrad_tran_sd_ctr     = matrix(NA, nrow=240, ncol=93)
    cov_tran_01_ctr      = matrix(NA, nrow=240, ncol=93)
    cov_lai__tran_01_ctr = matrix(NA, nrow=240, ncol=93)
    cov_gs___tran_01_ctr = matrix(NA, nrow=240, ncol=93)
    cov_vpd__tran_01_ctr = matrix(NA, nrow=240, ncol=93)
    cov_pr___tran_01_ctr = matrix(NA, nrow=240, ncol=93)
    cov_nrad_tran_01_ctr = matrix(NA, nrow=240, ncol=93)
    lai_sd_bgc           = matrix(NA, nrow=240, ncol=93) 
    lai_sd_rad           = matrix(NA, nrow=240, ncol=93) 
    lai_sd_ctr           = matrix(NA, nrow=240, ncol=93) 
    lai_cov_bgc          = matrix(NA, nrow=240, ncol=93) 
    lai_cov_rad          = matrix(NA, nrow=240, ncol=93) 
    lai_cov_ctr          = matrix(NA, nrow=240, ncol=93) 
    lucc = readRDS('/GPUFS/ygo_fwf_1/00junli/01next_step/lucc2016_1.5degree_vegetated_area.rds')
    for (i in 1:nrow(tas__bgc)) {
      for (j in 1:ncol(tas__bgc)) {
        if(is.na(lucc[i,j])  | is.na(tran_bgc[i,j,10]) | is.na(lai__bgc[i,j,10]) |
           is.na(gs___bgc[i,j,10]) | is.na(vpd__bgc[i,j,10]) | is.na(pr___bgc[i,j,10]) |
           is.na(nrad_bgc[i,j,10])) next#i=200;j=85;image.plot(1:240,1:93,tas__bgc[,,8]);abline(v=i,h=j)
        
        tas__ij  = rowMeans(embed(c(tas__ctr[i,j,],rep(NA, sca-1)), sca), na.rm=T)
        loca_bgc = which.max(rowMeans(matrix(tas__ij, ncol=yr), na.rm=T))
        locaij = sort(embed(c(1:12,1,2),3)[loca_bgc,])
        loca_sel = embed(1:yr, window_yr)
        
        tran_bgc_ij = try(decompose_own(ts(tran_bgc[i,j,]*86400*30, frequency=12),window_yr=30, type = "additive"), silent = T)
        lai__bgc_ij = try(decompose_own(ts(lai__bgc[i,j,]         , frequency=12),window_yr=30, type = "additive"), silent = T)
        gs___bgc_ij = try(decompose_own(ts(gs___bgc[i,j,]         , frequency=12),window_yr=30, type = "additive"), silent = T)
        vpd__bgc_ij = try(decompose_own(ts(vpd__bgc[i,j,]         , frequency=12),window_yr=30, type = "additive"), silent = T)
        pr___bgc_ij = try(decompose_own(ts(pr___bgc[i,j,]*86400*30, frequency=12),window_yr=30, type = "additive"), silent = T)
        nrad_bgc_ij = try(decompose_own(ts(nrad_bgc[i,j,]         , frequency=12),window_yr=30, type = "additive"), silent = T)
        
        tran_rad_ij = try(decompose_own(ts(tran_rad[i,j,]*86400*30, frequency=12),window_yr=30, type = "additive"), silent = T)
        lai__rad_ij = try(decompose_own(ts(lai__rad[i,j,]         , frequency=12),window_yr=30, type = "additive"), silent = T)
        gs___rad_ij = try(decompose_own(ts(gs___rad[i,j,]         , frequency=12),window_yr=30, type = "additive"), silent = T)
        vpd__rad_ij = try(decompose_own(ts(vpd__rad[i,j,]         , frequency=12),window_yr=30, type = "additive"), silent = T)
        pr___rad_ij = try(decompose_own(ts(pr___rad[i,j,]*86400*30, frequency=12),window_yr=30, type = "additive"), silent = T)
        nrad_rad_ij = try(decompose_own(ts(nrad_rad[i,j,]         , frequency=12),window_yr=30, type = "additive"), silent = T)
        
        tran_ctr_ij = try(decompose_own(ts(tran_ctr[i,j,]*86400*30, frequency=12),window_yr=30, type = "additive"), silent = T)
        lai__ctr_ij = try(decompose_own(ts(lai__ctr[i,j,]         , frequency=12),window_yr=30, type = "additive"), silent = T)
        gs___ctr_ij = try(decompose_own(ts(gs___ctr[i,j,]         , frequency=12),window_yr=30, type = "additive"), silent = T)
        vpd__ctr_ij = try(decompose_own(ts(vpd__ctr[i,j,]         , frequency=12),window_yr=30, type = "additive"), silent = T)
        pr___ctr_ij = try(decompose_own(ts(pr___ctr[i,j,]*86400*30, frequency=12),window_yr=30, type = "additive"), silent = T)
        nrad_ctr_ij = try(decompose_own(ts(nrad_ctr[i,j,]         , frequency=12),window_yr=30, type = "additive"), silent = T)
        
        if(class(tran_bgc_ij)[1] != "try-error" & class(tran_rad_ij)[1] != "try-error" & class(tran_ctr_ij)[1] != "try-error" &
           class(lai__bgc_ij)[1] != "try-error" & class(lai__rad_ij)[1] != "try-error" & class(lai__ctr_ij)[1] != "try-error" &
           class(gs___bgc_ij)[1] != "try-error" & class(gs___rad_ij)[1] != "try-error" & class(gs___ctr_ij)[1] != "try-error" &
           class(vpd__bgc_ij)[1] != "try-error" & class(vpd__rad_ij)[1] != "try-error" & class(vpd__ctr_ij)[1] != "try-error" &
           class(pr___bgc_ij)[1] != "try-error" & class(pr___rad_ij)[1] != "try-error" & class(pr___ctr_ij)[1] != "try-error" &
           class(nrad_bgc_ij)[1] != "try-error" & class(nrad_rad_ij)[1] != "try-error" & class(nrad_ctr_ij)[1] != "try-error"){
          
          getre = function(tran_bgc_ij,lai__bgc_ij,gs___bgc_ij,vpd__bgc_ij,pr___bgc_ij,nrad_bgc_ij, id1=09:3008, id2=c(1:6), yr=251){
            return(list(tran_ijk = matrix(tran_bgc_ij[id1], ncol=yr)[id2,],
                        lai__ijk = matrix(lai__bgc_ij[id1], ncol=yr)[id2,],
                        gs___ijk = matrix(gs___bgc_ij[id1], ncol=yr)[id2,],
                        vpd__ijk = matrix(vpd__bgc_ij[id1], ncol=yr)[id2,],
                        pr___ijk = matrix(pr___bgc_ij[id1], ncol=yr)[id2,],
                        nrad_ijk = matrix(nrad_bgc_ij[id1], ncol=yr)[id2,])) }
          
          if(loca_bgc==12){
            id1 = 09:1676; id2 = c(1:6); yr1 = 139 }else if(loca_bgc==11){
              id1 = 08:1675; id2 = c(1:6); yr1 = 139 }else if(loca_bgc==1){
                id1 = 10:1677; id2 = c(1:6); yr1 = 139 }else if(loca_bgc==2){
                  id1 = 11:1678; id2 = c(1:6); yr1 = 139 }else if(loca_bgc==3){
                    id1 = 12:1679; id2 = c(1:6); yr1 = 139 }else{
                      id1 = 1:1680 ; id2 = c((min(locaij)-3):max(locaij)); yr1 = 140}
          
          dare_bgc = getre(tran_bgc_ij$random,lai__bgc_ij$random,gs___bgc_ij$random,vpd__bgc_ij$random,pr___bgc_ij$random,nrad_bgc_ij$random,id1, id2, yr1)
          dare_rad = getre(tran_rad_ij$random,lai__rad_ij$random,gs___rad_ij$random,vpd__rad_ij$random,pr___rad_ij$random,nrad_rad_ij$random,id1, id2, yr1)
          dare_ctr = getre(tran_ctr_ij$random,lai__ctr_ij$random,gs___ctr_ij$random,vpd__ctr_ij$random,pr___ctr_ij$random,nrad_ctr_ij$random,id1, id2, yr1)
          
          cal_rc = function(dare1){
            tran_ijk = dare1$tran_ijk
            lai__ijk = dare1$lai__ijk
            gs___ijk = dare1$gs___ijk
            vpd__ijk = dare1$vpd__ijk
            pr___ijk = dare1$pr___ijk  
            nrad_ijk = dare1$nrad_ijk
            
            da1 = na.omit(data.frame(tran=tran_ijk[1,],lai=lai__ijk[1,],gs=gs___ijk[1,],vpd=vpd__ijk[1,],pr=pr___ijk[1,],nrad=nrad_ijk[1,]))
            da2 = na.omit(data.frame(tran=tran_ijk[2,],lai=lai__ijk[2,],gs=gs___ijk[2,],vpd=vpd__ijk[2,],pr=pr___ijk[2,],nrad=nrad_ijk[2,]))
            da3 = na.omit(data.frame(tran=tran_ijk[3,],lai=lai__ijk[3,],gs=gs___ijk[3,],vpd=vpd__ijk[3,],pr=pr___ijk[3,],nrad=nrad_ijk[3,]))
            da4 = na.omit(data.frame(tran=tran_ijk[4,],lai=lai__ijk[4,],gs=gs___ijk[4,],vpd=vpd__ijk[4,],pr=pr___ijk[4,],nrad=nrad_ijk[4,]))
            da5 = na.omit(data.frame(tran=tran_ijk[5,],lai=lai__ijk[5,],gs=gs___ijk[5,],vpd=vpd__ijk[5,],pr=pr___ijk[5,],nrad=nrad_ijk[5,]))
            da6 = na.omit(data.frame(tran=tran_ijk[6,],lai=lai__ijk[6,],gs=gs___ijk[6,],vpd=vpd__ijk[6,],pr=pr___ijk[6,],nrad=nrad_ijk[6,]))
            daall = na.omit(data.frame(tran = matrix(dare1$tran_ijk, ncol=1),
                                       lai  = matrix(dare1$lai__ijk, ncol=1),
                                       gs   = matrix(dare1$gs___ijk, ncol=1),
                                       vpd  = matrix(dare1$vpd__ijk, ncol=1),
                                       pr   = matrix(dare1$pr___ijk, ncol=1),
                                       nrad = matrix(dare1$nrad_ijk, ncol=1)))
            
            rc1   = try(plsr(tran ~ lai + gs + vpd + pr + nrad, data=da1  ), silent = T)
            rc2   = try(plsr(tran ~ lai + gs + vpd + pr + nrad, data=da2  ), silent = T)
            rc3   = try(plsr(tran ~ lai + gs + vpd + pr + nrad, data=da3  ), silent = T)
            rc4   = try(plsr(tran ~ lai + gs + vpd + pr + nrad, data=da4  ), silent = T)
            rc5   = try(plsr(tran ~ lai + gs + vpd + pr + nrad, data=da5  ), silent = T)
            rc6   = try(plsr(tran ~ lai + gs + vpd + pr + nrad, data=da6  ), silent = T)
            rcall = try(plsr(tran ~ lai + gs + vpd + pr + nrad, data=daall), silent = T)
            
            allda = list(rc1   = rc1  ,rc2   = rc2  ,rc3   = rc3  ,
                         rc4   = rc4  ,rc5   = rc5  ,rc6   = rc6  ,rcall = rcall)
            return(allda)
          }
          
          rc_bgc = cal_rc(dare_bgc)
          rc_rad = cal_rc(dare_rad)
          rc_ctr = cal_rc(dare_ctr)
          
          if(ncol(dare_bgc$tran_ijk) != yr){ tran_bgc_ijk = cbind(rep(NA, 6), dare_bgc$tran_ijk)}else{tran_bgc_ijk = dare_bgc$tran_ijk}
          if(ncol(dare_bgc$lai__ijk) != yr){ lai__bgc_ijk = cbind(rep(NA, 6), dare_bgc$lai__ijk)}else{lai__bgc_ijk = dare_bgc$lai__ijk}
          if(ncol(dare_bgc$gs___ijk) != yr){ gs___bgc_ijk = cbind(rep(NA, 6), dare_bgc$gs___ijk)}else{gs___bgc_ijk = dare_bgc$gs___ijk}
          if(ncol(dare_bgc$vpd__ijk) != yr){ vpd__bgc_ijk = cbind(rep(NA, 6), dare_bgc$vpd__ijk)}else{vpd__bgc_ijk = dare_bgc$vpd__ijk}
          if(ncol(dare_bgc$pr___ijk) != yr){ pr___bgc_ijk = cbind(rep(NA, 6), dare_bgc$pr___ijk)}else{pr___bgc_ijk = dare_bgc$pr___ijk}
          if(ncol(dare_bgc$nrad_ijk) != yr){ nrad_bgc_ijk = cbind(rep(NA, 6), dare_bgc$nrad_ijk)}else{nrad_bgc_ijk = dare_bgc$nrad_ijk}
          
          if(ncol(dare_rad$tran_ijk) != yr){ tran_rad_ijk = cbind(rep(NA, 6), dare_rad$tran_ijk)}else{tran_rad_ijk = dare_rad$tran_ijk}
          if(ncol(dare_rad$lai__ijk) != yr){ lai__rad_ijk = cbind(rep(NA, 6), dare_rad$lai__ijk)}else{lai__rad_ijk = dare_rad$lai__ijk}
          if(ncol(dare_rad$gs___ijk) != yr){ gs___rad_ijk = cbind(rep(NA, 6), dare_rad$gs___ijk)}else{gs___rad_ijk = dare_rad$gs___ijk}
          if(ncol(dare_rad$vpd__ijk) != yr){ vpd__rad_ijk = cbind(rep(NA, 6), dare_rad$vpd__ijk)}else{vpd__rad_ijk = dare_rad$vpd__ijk}
          if(ncol(dare_rad$pr___ijk) != yr){ pr___rad_ijk = cbind(rep(NA, 6), dare_rad$pr___ijk)}else{pr___rad_ijk = dare_rad$pr___ijk}
          if(ncol(dare_rad$nrad_ijk) != yr){ nrad_rad_ijk = cbind(rep(NA, 6), dare_rad$nrad_ijk)}else{nrad_rad_ijk = dare_rad$nrad_ijk}
          
          if(ncol(dare_ctr$tran_ijk) != yr){ tran_ctr_ijk = cbind(rep(NA, 6), dare_ctr$tran_ijk)}else{tran_ctr_ijk = dare_ctr$tran_ijk}
          if(ncol(dare_ctr$lai__ijk) != yr){ lai__ctr_ijk = cbind(rep(NA, 6), dare_ctr$lai__ijk)}else{lai__ctr_ijk = dare_ctr$lai__ijk}
          if(ncol(dare_ctr$gs___ijk) != yr){ gs___ctr_ijk = cbind(rep(NA, 6), dare_ctr$gs___ijk)}else{gs___ctr_ijk = dare_ctr$gs___ijk}
          if(ncol(dare_ctr$vpd__ijk) != yr){ vpd__ctr_ijk = cbind(rep(NA, 6), dare_ctr$vpd__ijk)}else{vpd__ctr_ijk = dare_ctr$vpd__ijk}
          if(ncol(dare_ctr$pr___ijk) != yr){ pr___ctr_ijk = cbind(rep(NA, 6), dare_ctr$pr___ijk)}else{pr___ctr_ijk = dare_ctr$pr___ijk}
          if(ncol(dare_ctr$nrad_ijk) != yr){ nrad_ctr_ijk = cbind(rep(NA, 6), dare_ctr$nrad_ijk)}else{nrad_ctr_ijk = dare_ctr$nrad_ijk}
          
          if(class(rc_bgc$rc1  ) != "try-error" & class(rc_rad$rc1  ) != "try-error" & class(rc_ctr$rc1  ) != "try-error" &
             class(rc_bgc$rc2  ) != "try-error" & class(rc_rad$rc2  ) != "try-error" & class(rc_ctr$rc2  ) != "try-error" &
             class(rc_bgc$rc3  ) != "try-error" & class(rc_rad$rc3  ) != "try-error" & class(rc_ctr$rc3  ) != "try-error" &
             class(rc_bgc$rc4  ) != "try-error" & class(rc_rad$rc4  ) != "try-error" & class(rc_ctr$rc4  ) != "try-error" &
             class(rc_bgc$rc5  ) != "try-error" & class(rc_rad$rc5  ) != "try-error" & class(rc_ctr$rc5  ) != "try-error" &
             class(rc_bgc$rc6  ) != "try-error" & class(rc_rad$rc6  ) != "try-error" & class(rc_ctr$rc6  ) != "try-error" &
             class(rc_bgc$rcall) != "try-error" & class(rc_rad$rcall) != "try-error" & class(rc_ctr$rcall) != "try-error"){
            
            lai_bgc_ijk_raw = lai__bgc_ijk
            lai_rad_ijk_raw = lai__rad_ijk
            lai_ctr_ijk_raw = lai__ctr_ijk
            
            lai__bgc_ijk = lai__bgc_ijk*c(coef(rc_bgc$rc1)[1],coef(rc_bgc$rc2)[1],coef(rc_bgc$rc3)[1],coef(rc_bgc$rc4)[1],coef(rc_bgc$rc5)[1],coef(rc_bgc$rc6)[1])
            gs___bgc_ijk = gs___bgc_ijk*c(coef(rc_bgc$rc1)[2],coef(rc_bgc$rc2)[2],coef(rc_bgc$rc3)[2],coef(rc_bgc$rc4)[2],coef(rc_bgc$rc5)[2],coef(rc_bgc$rc6)[2])
            vpd__bgc_ijk = vpd__bgc_ijk*c(coef(rc_bgc$rc1)[3],coef(rc_bgc$rc2)[3],coef(rc_bgc$rc3)[3],coef(rc_bgc$rc4)[3],coef(rc_bgc$rc5)[3],coef(rc_bgc$rc6)[3])
            pr___bgc_ijk = pr___bgc_ijk*c(coef(rc_bgc$rc1)[4],coef(rc_bgc$rc2)[4],coef(rc_bgc$rc3)[4],coef(rc_bgc$rc4)[4],coef(rc_bgc$rc5)[4],coef(rc_bgc$rc6)[4])
            nrad_bgc_ijk = nrad_bgc_ijk*c(coef(rc_bgc$rc1)[5],coef(rc_bgc$rc2)[5],coef(rc_bgc$rc3)[5],coef(rc_bgc$rc4)[5],coef(rc_bgc$rc5)[5],coef(rc_bgc$rc6)[5])
            
            lai__rad_ijk = lai__rad_ijk*c(coef(rc_rad$rc1)[1],coef(rc_rad$rc2)[1],coef(rc_rad$rc3)[1],coef(rc_rad$rc4)[1],coef(rc_rad$rc5)[1],coef(rc_rad$rc6)[1])
            gs___rad_ijk = gs___rad_ijk*c(coef(rc_rad$rc1)[2],coef(rc_rad$rc2)[2],coef(rc_rad$rc3)[2],coef(rc_rad$rc4)[2],coef(rc_rad$rc5)[2],coef(rc_rad$rc6)[2])
            vpd__rad_ijk = vpd__rad_ijk*c(coef(rc_rad$rc1)[3],coef(rc_rad$rc2)[3],coef(rc_rad$rc3)[3],coef(rc_rad$rc4)[3],coef(rc_rad$rc5)[3],coef(rc_rad$rc6)[3])
            pr___rad_ijk = pr___rad_ijk*c(coef(rc_rad$rc1)[4],coef(rc_rad$rc2)[4],coef(rc_rad$rc3)[4],coef(rc_rad$rc4)[4],coef(rc_rad$rc5)[4],coef(rc_rad$rc6)[4])
            nrad_rad_ijk = nrad_rad_ijk*c(coef(rc_rad$rc1)[5],coef(rc_rad$rc2)[5],coef(rc_rad$rc3)[5],coef(rc_rad$rc4)[5],coef(rc_rad$rc5)[5],coef(rc_rad$rc6)[5])
            
            lai__ctr_ijk = lai__ctr_ijk*c(coef(rc_ctr$rc1)[1],coef(rc_ctr$rc2)[1],coef(rc_ctr$rc3)[1],coef(rc_ctr$rc4)[1],coef(rc_ctr$rc5)[1],coef(rc_ctr$rc6)[1])
            gs___ctr_ijk = gs___ctr_ijk*c(coef(rc_ctr$rc1)[2],coef(rc_ctr$rc2)[2],coef(rc_ctr$rc3)[2],coef(rc_ctr$rc4)[2],coef(rc_ctr$rc5)[2],coef(rc_ctr$rc6)[2])
            vpd__ctr_ijk = vpd__ctr_ijk*c(coef(rc_ctr$rc1)[3],coef(rc_ctr$rc2)[3],coef(rc_ctr$rc3)[3],coef(rc_ctr$rc4)[3],coef(rc_ctr$rc5)[3],coef(rc_ctr$rc6)[3])
            pr___ctr_ijk = pr___ctr_ijk*c(coef(rc_ctr$rc1)[4],coef(rc_ctr$rc2)[4],coef(rc_ctr$rc3)[4],coef(rc_ctr$rc4)[4],coef(rc_ctr$rc5)[4],coef(rc_ctr$rc6)[4])
            nrad_ctr_ijk = nrad_ctr_ijk*c(coef(rc_ctr$rc1)[5],coef(rc_ctr$rc2)[5],coef(rc_ctr$rc3)[5],coef(rc_ctr$rc4)[5],coef(rc_ctr$rc5)[5],coef(rc_ctr$rc6)[5])
            
            loca_sel1 = rbind(loca_sel[1,], loca_sel[111,])
            cal_sd_cov = function(lai_bgc_ijk_raw, tran_bgc_ijk,lai__bgc_ijk,gs___bgc_ijk,vpd__bgc_ijk,pr___bgc_ijk,nrad_bgc_ijk,loca_sel1){
              
              lai_bgc_0   = rbind(colMeans(lai_bgc_ijk_raw[1:3,]),colMeans(lai_bgc_ijk_raw[2:4,]),colMeans(lai_bgc_ijk_raw[3:5,]))
              tran_bgc_0  = rbind(colMeans(tran_bgc_ijk   [1:3,]),colMeans(tran_bgc_ijk   [2:4,]),colMeans(tran_bgc_ijk   [3:5,]))
              lai__tran_0 = rbind(colMeans(lai__bgc_ijk   [1:3,]),colMeans(lai__bgc_ijk   [2:4,]),colMeans(lai__bgc_ijk   [3:5,]))
              gs___tran_0 = rbind(colMeans(gs___bgc_ijk   [1:3,]),colMeans(gs___bgc_ijk   [2:4,]),colMeans(gs___bgc_ijk   [3:5,]))
              vpd__tran_0 = rbind(colMeans(vpd__bgc_ijk   [1:3,]),colMeans(vpd__bgc_ijk   [2:4,]),colMeans(vpd__bgc_ijk   [3:5,]))
              pr___tran_0 = rbind(colMeans(pr___bgc_ijk   [1:3,]),colMeans(pr___bgc_ijk   [2:4,]),colMeans(pr___bgc_ijk   [3:5,]))
              nrad_tran_0 = rbind(colMeans(nrad_bgc_ijk   [1:3,]),colMeans(nrad_bgc_ijk   [2:4,]),colMeans(nrad_bgc_ijk   [3:5,]))
              
              lai_bgc_1   = rbind(lai_bgc_ijk_raw[4,],lai_bgc_ijk_raw[5,],lai_bgc_ijk_raw[6,])
              tran_bgc_1  = rbind(tran_bgc_ijk   [4,],tran_bgc_ijk   [5,],tran_bgc_ijk   [6,])
              lai__tran_1 = rbind(lai__bgc_ijk   [4,],lai__bgc_ijk   [5,],lai__bgc_ijk   [6,])
              gs___tran_1 = rbind(gs___bgc_ijk   [4,],gs___bgc_ijk   [5,],gs___bgc_ijk   [6,])
              vpd__tran_1 = rbind(vpd__bgc_ijk   [4,],vpd__bgc_ijk   [5,],vpd__bgc_ijk   [6,])
              pr___tran_1 = rbind(pr___bgc_ijk   [4,],pr___bgc_ijk   [5,],pr___bgc_ijk   [6,])
              nrad_tran_1 = rbind(nrad_bgc_ijk   [4,],nrad_bgc_ijk   [5,],nrad_bgc_ijk   [6,])
              
              cov_lai_bgc_01_ijk   = c()
              cov_tran_bgc_01_ijk  = c()
              cov_lai__tran_01_ijk = c()
              cov_gs___tran_01_ijk = c()
              cov_vpd__tran_01_ijk = c()
              cov_pr___tran_01_ijk = c()
              cov_nrad_tran_01_ijk = c()
              for (ijk in 1:nrow(loca_sel1)) {#ijk=1
                cov_lai_bgc_01_ijk   = c(cov_lai_bgc_01_ijk  , cov(na.omit(cbind(matrix(lai_bgc_0  [,loca_sel1[ijk,]], ncol=1), matrix(lai_bgc_1  [,loca_sel1[ijk,]], ncol=1))))[1,2])
                
                cov_tran_bgc_01_ijk  = c(cov_tran_bgc_01_ijk , cov(na.omit(cbind(matrix(tran_bgc_0 [,loca_sel1[ijk,]], ncol=1), matrix(tran_bgc_1 [,loca_sel1[ijk,]], ncol=1))))[1,2])
                cov_lai__tran_01_ijk = c(cov_lai__tran_01_ijk, cov(na.omit(cbind(matrix(lai__tran_0[,loca_sel1[ijk,]], ncol=1), matrix(lai__tran_1[,loca_sel1[ijk,]], ncol=1))))[1,2])
                cov_gs___tran_01_ijk = c(cov_gs___tran_01_ijk, cov(na.omit(cbind(matrix(gs___tran_0[,loca_sel1[ijk,]], ncol=1), matrix(gs___tran_1[,loca_sel1[ijk,]], ncol=1))))[1,2])
                cov_vpd__tran_01_ijk = c(cov_vpd__tran_01_ijk, cov(na.omit(cbind(matrix(vpd__tran_0[,loca_sel1[ijk,]], ncol=1), matrix(vpd__tran_1[,loca_sel1[ijk,]], ncol=1))))[1,2])
                cov_pr___tran_01_ijk = c(cov_pr___tran_01_ijk, cov(na.omit(cbind(matrix(pr___tran_0[,loca_sel1[ijk,]], ncol=1), matrix(pr___tran_1[,loca_sel1[ijk,]], ncol=1))))[1,2])
                cov_nrad_tran_01_ijk = c(cov_nrad_tran_01_ijk, cov(na.omit(cbind(matrix(nrad_tran_0[,loca_sel1[ijk,]], ncol=1), matrix(nrad_tran_1[,loca_sel1[ijk,]], ncol=1))))[1,2])}
              
              lai_sd_ijk_bgc       = c()
              tran_sd_ijk_bgc      = c()
              lai__tran_sd_ijk_bgc = c()
              gs___tran_sd_ijk_bgc = c()
              vpd__tran_sd_ijk_bgc = c()
              pr___tran_sd_ijk_bgc = c()
              nrad_tran_sd_ijk_bgc = c()
              for (ijk in 1:nrow(loca_sel1)) {#ij=222
                lai_sd_ijk_bgc       = c(lai_sd_ijk_bgc      , sd(lai_bgc_0  [,loca_sel1[ijk,]], na.rm=T)*sd(lai_bgc_1  [,loca_sel1[ijk,]],na.rm=T))
                
                tran_sd_ijk_bgc      = c(tran_sd_ijk_bgc     , sd(tran_bgc_0 [,loca_sel1[ijk,]], na.rm=T)*sd(tran_bgc_1 [,loca_sel1[ijk,]],na.rm=T))
                lai__tran_sd_ijk_bgc = c(lai__tran_sd_ijk_bgc, sd(lai__tran_0[,loca_sel1[ijk,]], na.rm=T)*sd(lai__tran_1[,loca_sel1[ijk,]],na.rm=T))
                gs___tran_sd_ijk_bgc = c(gs___tran_sd_ijk_bgc, sd(gs___tran_0[,loca_sel1[ijk,]], na.rm=T)*sd(gs___tran_1[,loca_sel1[ijk,]],na.rm=T))
                vpd__tran_sd_ijk_bgc = c(vpd__tran_sd_ijk_bgc, sd(vpd__tran_0[,loca_sel1[ijk,]], na.rm=T)*sd(vpd__tran_1[,loca_sel1[ijk,]],na.rm=T))
                pr___tran_sd_ijk_bgc = c(pr___tran_sd_ijk_bgc, sd(pr___tran_0[,loca_sel1[ijk,]], na.rm=T)*sd(pr___tran_1[,loca_sel1[ijk,]],na.rm=T))
                nrad_tran_sd_ijk_bgc = c(nrad_tran_sd_ijk_bgc, sd(nrad_tran_0[,loca_sel1[ijk,]], na.rm=T)*sd(nrad_tran_1[,loca_sel1[ijk,]],na.rm=T))}
              
              return(list(tran_sd_ijk          = tran_sd_ijk_bgc     [2] - tran_sd_ijk_bgc     [1],
                          lai__tran_sd_ijk     = lai__tran_sd_ijk_bgc[2] - lai__tran_sd_ijk_bgc[1],
                          gs___tran_sd_ijk     = gs___tran_sd_ijk_bgc[2] - gs___tran_sd_ijk_bgc[1],
                          vpd__tran_sd_ijk     = vpd__tran_sd_ijk_bgc[2] - vpd__tran_sd_ijk_bgc[1],
                          pr___tran_sd_ijk     = pr___tran_sd_ijk_bgc[2] - pr___tran_sd_ijk_bgc[1],
                          nrad_tran_sd_ijk     = nrad_tran_sd_ijk_bgc[2] - nrad_tran_sd_ijk_bgc[1],
                          cov_tran_01_ijk      = cov_tran_bgc_01_ijk [2] - cov_tran_bgc_01_ijk [1],
                          cov_lai__tran_01_ijk = cov_lai__tran_01_ijk[2] - cov_lai__tran_01_ijk[1],
                          cov_gs___tran_01_ijk = cov_gs___tran_01_ijk[2] - cov_gs___tran_01_ijk[1],
                          cov_vpd__tran_01_ijk = cov_vpd__tran_01_ijk[2] - cov_vpd__tran_01_ijk[1],
                          cov_pr___tran_01_ijk = cov_pr___tran_01_ijk[2] - cov_pr___tran_01_ijk[1],
                          cov_nrad_tran_01_ijk = cov_nrad_tran_01_ijk[2] - cov_nrad_tran_01_ijk[1],
                          lai_sd_ijk_bgc       = lai_sd_ijk_bgc      [2] - lai_sd_ijk_bgc      [1],
                          cov_lai_bgc_01_ijk   = cov_lai_bgc_01_ijk  [2] - cov_lai_bgc_01_ijk  [1]))
            }
            
            sd_cov_bgc = cal_sd_cov(lai_bgc_ijk_raw, tran_bgc_ijk,lai__bgc_ijk,gs___bgc_ijk,vpd__bgc_ijk,pr___bgc_ijk,nrad_bgc_ijk,loca_sel1)
            sd_cov_rad = cal_sd_cov(lai_rad_ijk_raw, tran_rad_ijk,lai__rad_ijk,gs___rad_ijk,vpd__rad_ijk,pr___rad_ijk,nrad_rad_ijk,loca_sel1)
            sd_cov_ctr = cal_sd_cov(lai_ctr_ijk_raw, tran_ctr_ijk,lai__ctr_ijk,gs___ctr_ijk,vpd__ctr_ijk,pr___ctr_ijk,nrad_ctr_ijk,loca_sel1)
            
            #(sd_cov_bgc$tran_sd_ijk[2] - sd_cov_bgc$tran_sd_ijk[1])
            #(sd_cov_bgc$tran_sd_ijk[2:111] - sd_cov_bgc$tran_sd_ijk[1]     )[110]
            tran_sd_bgc         [i,j] = sd_cov_bgc$tran_sd_ijk
            lai__tran_sd_bgc    [i,j] = sd_cov_bgc$lai__tran_sd_ijk 
            gs___tran_sd_bgc    [i,j] = sd_cov_bgc$gs___tran_sd_ijk 
            vpd__tran_sd_bgc    [i,j] = sd_cov_bgc$vpd__tran_sd_ijk 
            pr___tran_sd_bgc    [i,j] = sd_cov_bgc$pr___tran_sd_ijk 
            nrad_tran_sd_bgc    [i,j] = sd_cov_bgc$nrad_tran_sd_ijk 
            cov_tran_01_bgc     [i,j] = sd_cov_bgc$cov_tran_01_ijk   
            cov_lai__tran_01_bgc[i,j] = sd_cov_bgc$cov_lai__tran_01_ijk    
            cov_gs___tran_01_bgc[i,j] = sd_cov_bgc$cov_gs___tran_01_ijk    
            cov_vpd__tran_01_bgc[i,j] = sd_cov_bgc$cov_vpd__tran_01_ijk    
            cov_pr___tran_01_bgc[i,j] = sd_cov_bgc$cov_pr___tran_01_ijk    
            cov_nrad_tran_01_bgc[i,j] = sd_cov_bgc$cov_nrad_tran_01_ijk 
            tran_sd_rad         [i,j] = sd_cov_rad$tran_sd_ijk
            lai__tran_sd_rad    [i,j] = sd_cov_rad$lai__tran_sd_ijk 
            gs___tran_sd_rad    [i,j] = sd_cov_rad$gs___tran_sd_ijk 
            vpd__tran_sd_rad    [i,j] = sd_cov_rad$vpd__tran_sd_ijk 
            pr___tran_sd_rad    [i,j] = sd_cov_rad$pr___tran_sd_ijk 
            nrad_tran_sd_rad    [i,j] = sd_cov_rad$nrad_tran_sd_ijk 
            cov_tran_01_rad     [i,j] = sd_cov_rad$cov_tran_01_ijk   
            cov_lai__tran_01_rad[i,j] = sd_cov_rad$cov_lai__tran_01_ijk    
            cov_gs___tran_01_rad[i,j] = sd_cov_rad$cov_gs___tran_01_ijk    
            cov_vpd__tran_01_rad[i,j] = sd_cov_rad$cov_vpd__tran_01_ijk    
            cov_pr___tran_01_rad[i,j] = sd_cov_rad$cov_pr___tran_01_ijk    
            cov_nrad_tran_01_rad[i,j] = sd_cov_rad$cov_nrad_tran_01_ijk 
            tran_sd_ctr         [i,j] = sd_cov_ctr$tran_sd_ijk
            lai__tran_sd_ctr    [i,j] = sd_cov_ctr$lai__tran_sd_ijk 
            gs___tran_sd_ctr    [i,j] = sd_cov_ctr$gs___tran_sd_ijk 
            vpd__tran_sd_ctr    [i,j] = sd_cov_ctr$vpd__tran_sd_ijk 
            pr___tran_sd_ctr    [i,j] = sd_cov_ctr$pr___tran_sd_ijk 
            nrad_tran_sd_ctr    [i,j] = sd_cov_ctr$nrad_tran_sd_ijk 
            cov_tran_01_ctr     [i,j] = sd_cov_ctr$cov_tran_01_ijk   
            cov_lai__tran_01_ctr[i,j] = sd_cov_ctr$cov_lai__tran_01_ijk    
            cov_gs___tran_01_ctr[i,j] = sd_cov_ctr$cov_gs___tran_01_ijk    
            cov_vpd__tran_01_ctr[i,j] = sd_cov_ctr$cov_vpd__tran_01_ijk    
            cov_pr___tran_01_ctr[i,j] = sd_cov_ctr$cov_pr___tran_01_ijk    
            cov_nrad_tran_01_ctr[i,j] = sd_cov_ctr$cov_nrad_tran_01_ijk 
            lai_sd_bgc          [i,j] = sd_cov_bgc$lai_sd_ijk
            lai_sd_rad          [i,j] = sd_cov_rad$lai_sd_ijk
            lai_sd_ctr          [i,j] = sd_cov_ctr$lai_sd_ijk
            lai_cov_bgc         [i,j] = sd_cov_bgc$cov_lai_bgc_01_ijk   
            lai_cov_rad         [i,j] = sd_cov_rad$cov_lai_bgc_01_ijk   
            lai_cov_ctr         [i,j] = sd_cov_ctr$cov_lai_bgc_01_ijk   
          }
        }
      }
    }
    
    out = list(tran_sd_bgc          = tran_sd_bgc         ,
               lai__tran_sd_bgc     = lai__tran_sd_bgc    ,
               gs___tran_sd_bgc     = gs___tran_sd_bgc    ,
               vpd__tran_sd_bgc     = vpd__tran_sd_bgc    ,
               pr___tran_sd_bgc     = pr___tran_sd_bgc    ,
               nrad_tran_sd_bgc     = nrad_tran_sd_bgc    ,
               cov_tran_01_bgc      = cov_tran_01_bgc     ,
               cov_lai__tran_01_bgc = cov_lai__tran_01_bgc,
               cov_gs___tran_01_bgc = cov_gs___tran_01_bgc,
               cov_vpd__tran_01_bgc = cov_vpd__tran_01_bgc,
               cov_pr___tran_01_bgc = cov_pr___tran_01_bgc,
               cov_nrad_tran_01_bgc = cov_nrad_tran_01_bgc,
               tran_sd_rad          = tran_sd_rad         ,
               lai__tran_sd_rad     = lai__tran_sd_rad    ,
               gs___tran_sd_rad     = gs___tran_sd_rad    ,
               vpd__tran_sd_rad     = vpd__tran_sd_rad    ,
               pr___tran_sd_rad     = pr___tran_sd_rad    ,
               nrad_tran_sd_rad     = nrad_tran_sd_rad    ,
               cov_tran_01_rad      = cov_tran_01_rad     ,
               cov_lai__tran_01_rad = cov_lai__tran_01_rad,
               cov_gs___tran_01_rad = cov_gs___tran_01_rad,
               cov_vpd__tran_01_rad = cov_vpd__tran_01_rad,
               cov_pr___tran_01_rad = cov_pr___tran_01_rad,
               cov_nrad_tran_01_rad = cov_nrad_tran_01_rad,
               tran_sd_ctr          = tran_sd_ctr         ,
               lai__tran_sd_ctr     = lai__tran_sd_ctr    ,
               gs___tran_sd_ctr     = gs___tran_sd_ctr    ,
               vpd__tran_sd_ctr     = vpd__tran_sd_ctr    ,
               pr___tran_sd_ctr     = pr___tran_sd_ctr    ,
               nrad_tran_sd_ctr     = nrad_tran_sd_ctr    ,
               cov_tran_01_ctr      = cov_tran_01_ctr     ,
               cov_lai__tran_01_ctr = cov_lai__tran_01_ctr,
               cov_gs___tran_01_ctr = cov_gs___tran_01_ctr,
               cov_vpd__tran_01_ctr = cov_vpd__tran_01_ctr,
               cov_pr___tran_01_ctr = cov_pr___tran_01_ctr,
               cov_nrad_tran_01_ctr = cov_nrad_tran_01_ctr,
               lai_sd_bgc  = lai_sd_bgc ,
               lai_sd_rad  = lai_sd_rad ,
               lai_sd_ctr  = lai_sd_ctr ,
               lai_cov_bgc = lai_cov_bgc,
               lai_cov_rad = lai_cov_rad,
               lai_cov_ctr = lai_cov_ctr)
    return(out)
  }
  
  t1 = proc.time()
  core=core1
  cl = makeCluster(core)
  registerDoParallel(cl)
  lmf_contribution = foreach (kk = 1:10, .packages=c('forecast','pracma','abind','pls')) %dopar% {#kk=1
    out_re = get_re_sm_mean_g(modelna[kk], window_yr=30)
    out_re
  }
  stopCluster(cl)
  t2 = proc.time()
  print((t2 - t1)[3]/60)
  
  outna = 'Tran_by_lai_gs_vpd_pr_nrad_ctr_bgc_rad_'
  setwd('/GPUFS/ygo_fwf_1/00junli/01next_step/10results/tran_ano_explain')
  saveRDS(lmf_contribution[[01]], paste(outna, modelna[01],sep=''))
  saveRDS(lmf_contribution[[02]], paste(outna, modelna[02],sep=''))
  saveRDS(lmf_contribution[[03]], paste(outna, modelna[03],sep=''))
  saveRDS(lmf_contribution[[04]], paste(outna, modelna[04],sep=''))
  saveRDS(lmf_contribution[[05]], paste(outna, modelna[05],sep=''))
  saveRDS(lmf_contribution[[06]], paste(outna, modelna[06],sep=''))
  saveRDS(lmf_contribution[[07]], paste(outna, modelna[07],sep=''))
  saveRDS(lmf_contribution[[08]], paste(outna, modelna[08],sep=''))
  saveRDS(lmf_contribution[[09]], paste(outna, modelna[09],sep=''))
  saveRDS(lmf_contribution[[10]], paste(outna, modelna[10],sep=''))

  rm(lmf_contribution)
  gc()
}
#tran_ano(10)

albedo_ano = function(core1){
  
  #modelna 
  get_re_sm_mean_g = function(modelnai, window_yr=30){
    
    fphis = '/GPUFS/ygo_fwf_1/00junli/01next_step/Ctr_bgc_rad_data'
    setwd(fphis)
    fp_file = list.files(fphis)
    for (k in 1:length(fp_file)) {#k=10
      con_tas   = grepl(c('tas'        ), fp_file[k])
      con__lai  = grepl(c('lai_L'      ), fp_file[k])
      con_et    = grepl(c('evspsbl_A'  ), fp_file[k])#
      con_tran  = grepl(c('tran_L'     ), fp_file[k])#evspsbl_A
      con_hfss  = grepl(c('hfss_A'     ), fp_file[k])
      con_pr    = grepl(c('pr_A'       ), fp_file[k])
      con_hurs  = grepl(c('hurs_A'     ), fp_file[k])
      con_rsus  = grepl(c('rsus_A'     ), fp_file[k])
      con_rsds  = grepl(c('rsds_A'     ), fp_file[k])
      con_wind  = grepl(c('clt_A'      ), fp_file[k])
      con_mod   = grepl(c(modelnai     ), fp_file[k])
      con_ctr   = grepl(c('1pctCO2_'   ), fp_file[k])
      con_bgc   = grepl(c('1pctCO2-bgc'), fp_file[k])
      con_rad   = grepl(c('1pctCO2-rad'), fp_file[k])
      
      if(con_tas  & con_bgc & con_mod) {tas__his_bgcna = fp_file[k]}
      if(con__lai & con_bgc & con_mod) {lai__his_bgcna = fp_file[k]}
      if(con_et   & con_bgc & con_mod) {et_his_bgcna   = fp_file[k]}
      if(con_tran & con_bgc & con_mod) {tranhis_bgcna  = fp_file[k]}
      if(con_hfss & con_bgc & con_mod) {hfss_his_bgcna = fp_file[k]}
      if(con_pr   & con_bgc & con_mod) {pr_his_bgcna   = fp_file[k]}
      if(con_hurs & con_bgc & con_mod) {hurs_his_bgcna = fp_file[k]}
      
      if(con_tas  & con_rad & con_mod) {tas__his_radna = fp_file[k]}
      if(con__lai & con_rad & con_mod) {lai__his_radna = fp_file[k]}
      if(con_et   & con_rad & con_mod) {et_his_radna   = fp_file[k]}
      if(con_tran & con_rad & con_mod) {tranhis_radna  = fp_file[k]}
      if(con_hfss & con_rad & con_mod) {hfss_his_radna = fp_file[k]}
      if(con_pr   & con_rad & con_mod) {pr_his_radna   = fp_file[k]}
      if(con_hurs & con_rad & con_mod) {hurs_his_radna = fp_file[k]}
      
      if(con_tas  & con_ctr & con_mod) {tas__his_ctrna = fp_file[k]}
      if(con__lai & con_ctr & con_mod) {lai__his_ctrna = fp_file[k]}
      if(con_et   & con_ctr & con_mod) {et_his_ctrna   = fp_file[k]}
      if(con_tran & con_ctr & con_mod) {tranhis_ctrna  = fp_file[k]}
      if(con_hfss & con_ctr & con_mod) {hfss_his_ctrna = fp_file[k]}
      if(con_pr   & con_ctr & con_mod) {pr_his_ctrna   = fp_file[k]}
      if(con_hurs & con_ctr & con_mod) {hurs_his_ctrna = fp_file[k]}
      
      if(con_rsus & con_bgc & con_mod) {rsus_his_bgcna = fp_file[k]}
      if(con_rsus & con_rad & con_mod) {rsus_his_radna = fp_file[k]}
      if(con_rsus & con_ctr & con_mod) {rsus_his_ctrna = fp_file[k]}
      
      if(con_rsds & con_bgc & con_mod) {rsds_his_bgcna = fp_file[k]}
      if(con_rsds & con_rad & con_mod) {rsds_his_radna = fp_file[k]}
      if(con_rsds & con_ctr & con_mod) {rsds_his_ctrna = fp_file[k]}
      
      if(con_wind & con_bgc & con_mod) {wind_his_bgcna = fp_file[k]}
      if(con_wind & con_rad & con_mod) {wind_his_radna = fp_file[k]}
      if(con_wind & con_ctr & con_mod) {wind_his_ctrna = fp_file[k]}
    }
    tas__bgc = readRDS(tas__his_bgcna)[,,1:1680]
    lai__bgc = readRDS(lai__his_bgcna)[,,1:1680]
    pr___bgc = readRDS(pr_his_bgcna  )[,,1:1680] 
    
    tas__rad = readRDS(tas__his_radna)[,,1:1680]
    lai__rad = readRDS(lai__his_radna)[,,1:1680]
    pr___rad = readRDS(pr_his_radna  )[,,1:1680] 
    
    tas__ctr = readRDS(tas__his_ctrna)[,,1:1680]
    lai__ctr = readRDS(lai__his_ctrna)[,,1:1680]
    pr___ctr = readRDS(pr_his_ctrna  )[,,1:1680] 
    
    rsus_bgc = readRDS(rsus_his_bgcna)[,,1:1680] 
    rsus_rad = readRDS(rsus_his_radna)[,,1:1680] 
    rsus_ctr = readRDS(rsus_his_ctrna)[,,1:1680] 
    rsds_bgc = readRDS(rsds_his_bgcna)[,,1:1680] 
    rsds_rad = readRDS(rsds_his_radna)[,,1:1680] 
    rsds_ctr = readRDS(rsds_his_ctrna)[,,1:1680] 
    wind_bgc = readRDS(wind_his_bgcna)[,,1:1680] 
    wind_rad = readRDS(wind_his_radna)[,,1:1680] 
    wind_ctr = readRDS(wind_his_ctrna)[,,1:1680] 
    
    albe_bgc = rsus_bgc/rsds_bgc
    albe_rad = rsus_rad/rsds_rad
    albe_ctr = rsus_ctr/rsds_ctr
    albe_bgc[is.infinite(albe_bgc)] = NA
    albe_rad[is.infinite(albe_rad)] = NA
    albe_ctr[is.infinite(albe_ctr)] = NA
    rm(rsus_bgc)
    rm(rsus_rad)
    rm(rsus_ctr)
    rm(rsds_bgc)
    rm(rsds_rad)
    rm(rsds_ctr); gc()
    
    sca=3
    yr=140
    len = yr-window_yr+1
    nrow1=240
    ncol1=93
    albe_bgc_c      = matrix(NA, nrow=240, ncol=93)
    lai__alb_bgc_c  = matrix(NA, nrow=240, ncol=93)
    pr___alb_bgc_c  = matrix(NA, nrow=240, ncol=93)
    wind_alb_bgc_c  = matrix(NA, nrow=240, ncol=93)
    albe_rad_c      = matrix(NA, nrow=240, ncol=93)
    lai__alb_rad_c  = matrix(NA, nrow=240, ncol=93)
    pr___alb_rad_c  = matrix(NA, nrow=240, ncol=93)
    wind_alb_rad_c  = matrix(NA, nrow=240, ncol=93)
    albe_ctr_c      = matrix(NA, nrow=240, ncol=93)
    lai__alb_ctr_c  = matrix(NA, nrow=240, ncol=93)
    pr___alb_ctr_c  = matrix(NA, nrow=240, ncol=93)
    wind_alb_ctr_c  = matrix(NA, nrow=240, ncol=93)
    lucc = readRDS('/GPUFS/ygo_fwf_1/00junli/01next_step/lucc2016_1.5degree_vegetated_area.rds')
    for (i in 1:nrow(tas__bgc)) {
      for (j in 1:ncol(tas__bgc)) {
        if(is.na(lucc[i,j])  |  is.na(lai__bgc[i,j,10]) | is.na(pr___bgc[i,j,10])) next#i=50;j=80;image.plot(1:240,1:93,tas__bgc[,,8]);abline(v=i,h=j)
        
        tas__ij  = rowMeans(embed(c(tas__ctr[i,j,],rep(NA, sca-1)), sca), na.rm=T)
        loca_bgc = which.max(rowMeans(matrix(tas__ij, ncol=yr), na.rm=T))
        locaij = sort(embed(c(1:12,1,2),3)[loca_bgc,])
        loca_sel = embed(1:yr, window_yr)
        
        calbgcradctr = function(albe_bgc,lai__bgc,pr___bgc,wind_bgc, i, j){
          yr=140
          daalbe = matrix(albe_bgc[i,j,], ncol=yr)
          dalai_ = matrix(lai__bgc[i,j,], ncol=yr)
          dapr__ = matrix(pr___bgc[i,j,], ncol=yr)
          dawind = matrix(wind_bgc[i,j,], ncol=yr)
          
          cal_att = function(daalbe,dalai_,dapr__,dawind){#k=7
            
            alblai_ = matrix(NA, nrow=12, ncol=140)
            albpr__ = matrix(NA, nrow=12, ncol=140)
            albwind = matrix(NA, nrow=12, ncol=140)
            albfit  = matrix(NA, nrow=12, ncol=140)
            for (k in 1:12) {
              rc=try(plsr(daalbe[k,]~
                          dalai_[k,]+
                          dapr__[k,]+
                          dawind[k,]), silent = T)
              if(class(rc)!= "try-error"){
                if(length(rc$fitted.values) == 140){
                  #print(k)
                  alblai_[k,] = dalai_[k,]*coef(rc)[1]
                  albpr__[k,] = dapr__[k,]*coef(rc)[2]
                  albwind[k,] = dawind[k,]*coef(rc)[3]
                  albfit [k,] = rc$fitted.values
                }
              }
            }
            return(list(alblai_ = alblai_,
                        albpr__ = albpr__,
                        albwind = albwind,
                        albfit  = albfit))
          }
          re = cal_att(daalbe,dalai_,dapr__,dawind)
          
          da_albedo  = try(decompose_own(ts(albe_bgc            [i,j,], frequency = 12), window_yr=30), silent = T)
          da_alblai_ = try(decompose_own(ts(matrix(re$alblai_, ncol=1), frequency = 12), window_yr=30), silent = T)
          da_albpr__ = try(decompose_own(ts(matrix(re$albpr__, ncol=1), frequency = 12), window_yr=30), silent = T)
          da_albwind = try(decompose_own(ts(matrix(re$albwind, ncol=1), frequency = 12), window_yr=30), silent = T)
          
          if(class(da_albedo )[1] != "try-error" &
             class(da_alblai_)[1] != "try-error" &
             class(da_albpr__)[1] != "try-error" &
             class(da_albwind)[1] != "try-error") {
            
            getalbedo = function(albe_bgc,lai__bgc,pr___bgc,wind_bgc,id1=09:3008, id2=c(1:6), yr=251){
              return(list(albe_ijk = matrix(albe_bgc[id1], ncol=yr)[id2,],
                          lai__ijk = matrix(lai__bgc[id1], ncol=yr)[id2,], 
                          pr___ijk = matrix(pr___bgc[id1], ncol=yr)[id2,],
                          wind_ijk = matrix(wind_bgc[id1], ncol=yr)[id2,])) }
            if(loca_bgc==12){
              id1 = 09:3008; id2 = c(1:6); yr1 = 2100-1850 }else if(loca_bgc==11){
                id1 = 08:3007; id2 = c(1:6); yr1 = 2100-1850 }else if(loca_bgc==1){
                  id1 = 10:3009; id2 = c(1:6); yr1 = 2100-1850 }else if(loca_bgc==2){
                    id1 = 11:3010; id2 = c(1:6); yr1 = 2100-1850 }else if(loca_bgc==3){
                      id1 = 12:3011; id2 = c(1:6); yr1 = 2100-1850 }else{
                        id1 = 1:3012; id2 = c((min(locaij)-3):max(locaij)); yr1 = 2100-1850+1
                      }
            
            albedo = getalbedo(da_albedo $trend,
                               da_alblai_$trend,
                               da_albpr__$trend,
                               da_albwind$trend, id1, id2, yr1)
            
            albe_ijk = albedo$albe_ijk
            lai__ijk = albedo$lai__ijk
            pr___ijk = albedo$pr___ijk
            wind_ijk = albedo$wind_ijk
            
            if(ncol(albe_ijk) != yr){ albe_ijk = cbind(rep(NA, 6), albe_ijk)}
            if(ncol(lai__ijk) != yr){ lai__ijk = cbind(rep(NA, 6), lai__ijk)}
            if(ncol(pr___ijk) != yr){ pr___ijk = cbind(rep(NA, 6), pr___ijk)}
            if(ncol(wind_ijk) != yr){ wind_ijk = cbind(rep(NA, 6), wind_ijk)}
            
            albe_1 = rbind(albe_ijk[4,], albe_ijk[5,], albe_ijk[6,])
            lai__1 = rbind(lai__ijk[4,], lai__ijk[5,], lai__ijk[6,])
            pr___1 = rbind(pr___ijk[4,], pr___ijk[5,], pr___ijk[6,])
            wind_1 = rbind(wind_ijk[4,], wind_ijk[5,], wind_ijk[6,])
            
            id1=1;id2=111
            albe_c = mean(albe_1[,loca_sel[id2,]], na.rm=T) - mean(albe_1[,loca_sel[id1,]], na.rm=T)
            lai__c = mean(lai__1[,loca_sel[id2,]], na.rm=T) - mean(lai__1[,loca_sel[id1,]], na.rm=T)
            pr___c = mean(pr___1[,loca_sel[id2,]], na.rm=T) - mean(pr___1[,loca_sel[id1,]], na.rm=T)
            wind_c = mean(wind_1[,loca_sel[id2,]], na.rm=T) - mean(wind_1[,loca_sel[id1,]], na.rm=T)
            
            return(c(albe_c,lai__c,pr___c,wind_c))
          }else{
            return(c(NA,NA,NA,NA))
          }
          
        }
        
        re_bgc = calbgcradctr(albe_bgc,lai__bgc,pr___bgc,wind_bgc, i, j)
        re_rad = calbgcradctr(albe_rad,lai__rad,pr___rad,wind_rad, i, j)
        re_ctr = calbgcradctr(albe_ctr,lai__ctr,pr___ctr,wind_ctr, i, j)
        
        albe_bgc_c    [i,j] = re_bgc[1]  
        lai__alb_bgc_c[i,j] = re_bgc[2]      
        pr___alb_bgc_c[i,j] = re_bgc[3]      
        wind_alb_bgc_c[i,j] = re_bgc[4]      
        albe_rad_c    [i,j] = re_rad[1]  
        lai__alb_rad_c[i,j] = re_rad[2]      
        pr___alb_rad_c[i,j] = re_rad[3]      
        wind_alb_rad_c[i,j] = re_rad[4]      
        albe_ctr_c    [i,j] = re_ctr[1] 
        lai__alb_ctr_c[i,j] = re_ctr[2]     
        pr___alb_ctr_c[i,j] = re_ctr[3]     
        wind_alb_ctr_c[i,j] = re_ctr[4]     
        
      }
    }
    
    out = list(albe_bgc_c     = albe_bgc_c    ,
               lai__alb_bgc_c = lai__alb_bgc_c,
               pr___alb_bgc_c = pr___alb_bgc_c,
               wind_alb_bgc_c = wind_alb_bgc_c,
               albe_rad_c     = albe_rad_c    ,
               lai__alb_rad_c = lai__alb_rad_c,
               pr___alb_rad_c = pr___alb_rad_c,
               wind_alb_rad_c = wind_alb_rad_c,
               albe_ctr_c     = albe_ctr_c    ,
               lai__alb_ctr_c = lai__alb_ctr_c,
               pr___alb_ctr_c = pr___alb_ctr_c,
               wind_alb_ctr_c = wind_alb_ctr_c)
    return(out)
  }
  
  t1 = proc.time()
  core=core1
  cl = makeCluster(core)
  registerDoParallel(cl)
  lmf_contribution = foreach (kk = 1:10, .packages=c('forecast','pracma','abind','pls')) %dopar% {#kk=1
    out_re = get_re_sm_mean_g(modelna[kk], window_yr=30)
    out_re
  }
  stopCluster(cl)
  t2 = proc.time()
  print((t2 - t1)[3]/60)
  
  #modelna = 
  setwd('/GPUFS/ygo_fwf_1/00junli/01next_step/10results_Albedo_induce_change/alb_tr_attribution_ctr_bgc_rad')
  saveRDS(lmf_contribution[[01]], paste(outna, modelna[01],sep=''))
  saveRDS(lmf_contribution[[02]], paste(outna, modelna[02],sep=''))
  saveRDS(lmf_contribution[[03]], paste(outna, modelna[03],sep=''))
  saveRDS(lmf_contribution[[04]], paste(outna, modelna[04],sep=''))
  saveRDS(lmf_contribution[[05]], paste(outna, modelna[05],sep=''))
  saveRDS(lmf_contribution[[06]], paste(outna, modelna[06],sep=''))
  saveRDS(lmf_contribution[[07]], paste(outna, modelna[07],sep=''))
  saveRDS(lmf_contribution[[08]], paste(outna, modelna[08],sep=''))
  saveRDS(lmf_contribution[[09]], paste(outna, modelna[09],sep=''))
  saveRDS(lmf_contribution[[10]], paste(outna, modelna[10],sep=''))

  rm(lmf_contribution)
  gc()
}
#albedo_ano(12)