read_year_data=function(yi) {
  library(imputeTS)
  library(data.table)
  cat("Year:",yi)
  omni_features=c('B','Bz','Sigma_M','Sigma_V','Vs','Pd','Pt','Fp','Ef','Rap','Pbeta','Amn','Mmn')
  interpret_option="linear"
  
  kp_file  = paste0("Data/Origin/Kp/Kp_ap_",yi,".txt")
  omni_file= paste0("Data/Origin/omni-high-resolution-1min/",yi,"-0101-1231-omni_min_Data.txt")
  
  Kp_ap_i=read.table(kp_file, quote="\"")
  omni_raw_i= read.table(omni_file, quote="\"", comment.char="")
  
  # assign names (do it once via list instead of repeating)
  names(Kp_ap_i)=c('Year','Month','Day','Hour.head','Hour.mid','Days','Days.mid','Kp','ap','Definitive')
  
  if(yi=="2026"){
    names(omni_raw_i)=c('Year','Day','Hour','Minute','B','Bx','By','Bz','Sigma_M','Sigma_V','Vs','Pd','Pt',
                        'Fp','Ef','Pbeta','Amn','Mmn','SYMd','SYMh')
  }else{
    names(omni_raw_i)=c('Year','Day','Hour','Minute','B','Bx','By','Bz','Sigma_M','Sigma_V','Vs','Pd','Pt',
                        'Rap','Fp','Ef','Pbeta','Amn','Mmn','SYMd','SYMh')
  }
  omni_raw_i=omni_raw_i[,1:(ncol(omni_raw_i)-2)] # drop last 2 cols for SYMd, SYMh
  
  # interpolate omni features
  omni_i=omni_raw_i
  for(col_id in omni_features){
    if(col_id %in%  names(omni_raw_i) ){
      omni_raw_i[which(omni_raw_i[col_id]==max(omni_raw_i[col_id],na.rm=TRUE)),col_id]=NA
      omni_i[which(omni_i[col_id]==max(omni_i[col_id],na.rm=TRUE)),col_id]=NA
      omni_i[col_id]=na_interpolation(omni_i[col_id],option=interpret_option)
    }else{
      omni_i[col_id]=NA
      omni_raw_i[col_id]=NA
    }
    
  }
  
  list(Kp_ap=Kp_ap_i,omni=omni_i, omni_raw=omni_raw_i,days=nrow(Kp_ap_i)/8)
}


