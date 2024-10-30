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

#Please note that:
#The overarching function is composed of one or two constituent functions; the intended outcome can be computed directly from the input variables.
###################lmf calculation###############################

#This function is designed for the calculation of the LMF and its decomposition across all datasets utilized in CMIP6, employing parallel computation.
#The parameter "fpall" specifies the file path containing the CMIP6 simulation data. The dataset comprises three-dimensional representations for historical (1850-2014) and future (2015-2100) periods for each model.
#The parameter "scenario" refers to the CMIP6 scenarios, such as "_ssp585_", "_ssp126_", and "_ssp245_", which are employed to read the relevant data.
#The parameter "modelnai" denotes the model name, which is also used to load the CMIP6 data.
#The "threshold" parameter defines the criteria for identifying compound dry-hot events, with a typical value of 0.9.
#The "window_yr" parameter specifies the length of the moving window used to calculate compound events, commonly set to 30 years.
#The "get_lmf_all" function encompasses a nested function, i.e., get_lmf_ij_sd, which facilitates LMF calculation and decomposition at each grid point. In this function, i and j represent grid locations, while the other parameters correspond to those in the get_lmf_all function.

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
    foreach (k2 = 1:dim(tas__bgc)[2],.combine='rbind') %dopar% {  # 
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

  lmf_bgc_timing = array(NA, c(dim(tas__bgc)[1], dim(tas__bgc)[2], da_len))
  for (k3 in 1:length(lmf_bgc1)) {#k3=1
    lmf_bgc_timing[k3,,] = lmf_bgc1[[k3]]
  }
  
  out_na = 'lmf-sm-ta-mean-sm-ta-var-sm-ta-dep-total_'
  setwd('/GPUFS/ygo_fwf_1/00junli/01next_step/7lmf_all_and_ensemble')
  saveRDS(lmf_bgc_timing, file = paste(out_na, scenario, modelnai, '_', threshold, '_', window_yr, sep=''))
} 

#This is an example.
#for (kk in 1:17) {
#  tryCatch({get_lmf_all(fpall, '_ssp585_', modelna[kk], 0.90, 30);print(kk)}, error=function(e){print(paste('Error:', kk))})}

#
###################tran-induced change##########################################################

#The function "decompose_own" is employed to detrend and deseasonalize the data, while "x" represents the time series at each grid derived from CMIP6 data. Additionally, "window_yr" refers to the moving window utilized in the analysis.
decompose_own = function(x, window_yr){
  #We eliminated the mean seasonality from the first 30-year period, along with the long-term trends.
  type = "additive"
  filter = NULL
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


#The functions "cal_tran_induced_change" and "cal_albedo_induced_change" are utilized to compute the Tran-induced and albedo-induced changes, respectively, across all datasets employed, leveraging parallel computation.
#The variable "fpall" denotes the file path containing the CMIP6 simulation data, which is structured as three-dimensional arrays representing both historical (1850-2014) and future (2015-2100) periods for each model.
#The variable "scenario" refers to the CMIP6 scenarios, such as '_ssp585_'. 
#The variable "modelna" identifies the specific model name required for loading the CMIP6 data.
#The variable "core1" is designated for defining the computational core for parallel processing. 
#The variable "window_yr" represents a moving window utilized for calculating compound events.
#The function "get_re_sm_mean_g", which is integrated within "cal_tran_induced_change", shares the same parameters and is employed to derive tran-induced changes for each CMIP6 model based on water and energy balance equations. 
#Similarly, the function "get_re_sm_alb_g", contained within "cal_albedo_induced_change", utilizes the same parameters to obtain albedo-induced changes for each CMIP6 model, also grounded in the water and energy balance equations.
#For a comprehensive understanding of the water and energy balance equations, please refer to our published paper.
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
    
    out = list(da_tran_sm           = da_tran_sm        ,#
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
  
  outna = 'Tran_induced_changes_'
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
cal_albedo_induced_change = function(fpall, modelna, scenario, core1, window_yr){
  
  get_re_sm_alb_g = function(fpall, modelnai, scenario, window_yr=window_yr){
    
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
    out_re = get_re_sm_alb_g(fpall, modelna[kk], scenario, window_yr=window_yr)
    out_re
  }
  stopCluster(cl)
  t2 = proc.time()
  print((t2 - t1)[3]/60)
  
  outna = 'Albedo_induce_change_'
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
                 
#There are examples.                
#cal_tran_induced_change(fpall, modelna, '_ssp585_', 17, window_yr=30)
##cal_albedo_induced_change(fpall, modelna, '_ssp585_', 17, window_yr=30)             
##
###################lmf Prediction####################################################
#This function is designed for predicting LMF related to vegetation greening across all modelled datasets. 
#The variable "rena" represents the data output from the get_lmf_all function, which includes each LMF and the associated changes. 
#The variable "lai_gs_indc" refers to the data on tran-induced or albedo-induced changes derived from the aforementioned functions. 
#The "threshold" parameter defines the criteria for identifying compound dry-hot events, such as a value of 0.9. 
#The "window_yr" parameter represents a moving window utilized for calculating these compound events. 
#The variable "out_na" specifies the output file name.
#Additionally, two auxiliary functions are incorporated: get_da and cal_re. 
#The get_da function is employed for LMF prediction at each grid cell, utilizing inputs from the cal_re function. Few fitting equations yield intercepts marginally exceeding or falling zero, and the Tran-induced changes are comparatively modest, we opt to deduct the prediction model intercept to avoid the potential underestimation or overestimation in predictions.              

lmf_pre1 = function(rena, lai_gs_indc, modelna, threshold=0.9, window_yr=30, out_na){
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

#This is an example.
#lmf_pre1(rena585_0.95_30, tran_indc585_30yr, modelna, threshold=0.95, window_yr=30, 'Pre_sm_ta_sd_together_ssp585_')
                 

