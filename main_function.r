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
        
        ta__thre_detrend = cdfnor(quanor(threshold, para=c(0, sd(tas_1 , na.rm=T))), para=c(0, sd(tas_2 , na.rm=T)))
        sm__thre_detrend = cdfnor(quanor(threshold, para=c(0, sd(mrso_1, na.rm=T))), para=c(0, sd(mrso_2, na.rm=T)))
        
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

#for (kk in 1:17) {
#  tryCatch({get_lmf_all(fpall, '_ssp585_', modelna[kk], 0.90, 30);print(kk)}, error=function(e){print(paste('Error:', kk))})}
#for (kk in 1:17) {
#  tryCatch({get_lmf_all(fpall, '_ssp126_', modelna[kk], 0.90, 30);print(kk)}, error=function(e){print(paste('Error:', kk))})}
#for (kk in 1:17) {
#  tryCatch({get_lmf_all(fpall, '_ssp245_', modelna[kk], 0.90, 30);print(kk)}, error=function(e){print(paste('Error:', kk))})}
#for (kk in 1:17) {
#  tryCatch({get_lmf_all(fpall, '_ssp370_', modelna[kk], 0.90, 30);print(kk)}, error=function(e){print(paste('Error:', kk))})}
#
#
#for (kk in 1:17) {
#  tryCatch({get_lmf_all(fpall, '_ssp585_', modelna[kk], 0.90, 20);print(kk)}, error=function(e){print(paste('Error:', kk))})}
#for (kk in 1:17) {
#  tryCatch({get_lmf_all(fpall, '_ssp585_', modelna[kk], 0.90, 40);print(kk)}, error=function(e){print(paste('Error:', kk))})}
#for (kk in 1:17) {
#  tryCatch({get_lmf_all(fpall, '_ssp585_', modelna[kk], 0.95, 30);print(kk)}, error=function(e){print(paste('Error:', kk))})}
#for (kk in 1:17) {
#  tryCatch({get_lmf_all(fpall, '_ssp585_', modelna[kk], 0.80, 30);print(kk)}, error=function(e){print(paste('Error:', kk))})}
#for (kk in 1:17) {
#  tryCatch({get_lmf_all(fpall, '_ssp585_', modelna[kk], 0.85, 30);print(kk)}, error=function(e){print(paste('Error:', kk))})}

###################tran-induced change##########################################################


#cal_tran_induced_change(fpall, modelna, '_ssp585_', 17, window_yr=30)
#cal_tran_induced_change(fpall, modelna, '_ssp370_', 17, window_yr=30)
#cal_tran_induced_change(fpall, modelna, '_ssp245_', 17, window_yr=30)
#cal_tran_induced_change(fpall, modelna, '_ssp126_', 17, window_yr=30)
#cal_tran_induced_change(fpall, modelna, '_ssp585_', 17, window_yr=20)
#cal_tran_induced_change(fpall, modelna, '_ssp585_', 17, window_yr=40)


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

##
###################Tran alb variability and cov ############################


#Tran_ano_plsr('_ssp585_', window_yr=30, 23)
#Tran_ano_plsr('_ssp370_', window_yr=30, 23)
#Tran_ano_plsr('_ssp245_', window_yr=30, 23)
#Tran_ano_plsr('_ssp126_', window_yr=30, 23)


##
###################Tran alb variability and cov 3 co2 sensitive######################

#tran_ano(10)

