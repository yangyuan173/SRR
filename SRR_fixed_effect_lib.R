#SRR fixed effects model
library(Matrix)
SRR_fixed <- function(data0, Y_name, z_names, Fac_name, test = "Score", refer_fac = NULL,
                      criterion = 1e-8, max.iter = 500, bound = 8, size_cut = 10){
  Yid <- match(Y_name, colnames(data0))
  if(is.na(Yid)){
    stop(paste("Can not find ", Y_name, " in the data!", sep=""))
  }
  
  zid <- match(z_names, colnames(data0))
  if(length(which(is.na(zid)))>0){
    stop(paste("Error: Can not find ", z_names[which(is.na(zid))], " in the data!", sep=""))
  }
  
  Facid <- match(Fac_name, colnames(data0))
  if(is.na(Facid)){
    stop(paste("Error: Can not find ", Fac_name, " in the data!", sep=""))
  }
  if(sum(is.na(data0[,Facid]))>0){
    stop("Error: You have missingness in your Facility variable")
  }
  Zmis_index=which(apply(data0[, zid],2,function(x) sum(is.na(x)))>0)
  if(sum(is.na(data0[, zid]))>0){
    print(paste0("Warning: You have missingness in your explanatory covariate: ",z_names[Zmis_index]))
  }
  var0_index<-which(apply(data0[, zid],2,var)==0)
  if(length(var0_index)>0){
    print(paste0("Warning: You have no variation in your variable: ",z_names[var0_index]))
  }

  Zmis_row_index<-which(apply(as.matrix(data0[, zid][,Zmis_index]),1,function(x) sum(is.na(x))>0))
  if(rankMatrix(as.matrix(data0[, zid][-Zmis_row_index,]))<length(z_names)){
    print(paste0("Warning: rank of the variables is less than the predictor numbers"))
  }
  
  z <- as.matrix(data0[,zid])
  tmp <- apply(z, 1, sum)
  NAFacid <- which(is.na(tmp))
  
  ######### check facility/hospital with 0/all events
  #check facility with 0 events
  if(length(NAFacid)!=0){
    Y <- as.vector(data0[-NAFacid,Yid])
    z <- as.matrix(data0[-NAFacid,zid])
    Fac <- as.vector(data0[-NAFacid,Facid])
    m2 <- table(Fac)
  }else{
    Y <- as.vector(data0[,Yid])
    z <- as.matrix(data0[,zid])
    Fac <- as.vector(data0[,Facid])
    m2 <- table(Fac)
  }
  
  
  Y_sort = Y[order(Fac)]
  z_sort = z[order(Fac),]
  Fac_sort = Fac[order(Fac)]
  
  ### 
  
  size_long = rep(m2,m2)
  
  # all 0 or all 1 or  facility size smaller than 10
  
  Y_sum = sapply(split(Y_sort,Fac_sort),sum)
  Y_sum_long = rep(Y_sum,m2)
  include = !(Y_sum_long == 0 | Y_sum_long == size_long | size_long < size_cut)
  
  exclude_fac = unique(Fac_sort[!include])
  if(length(exclude_fac)>0){
    warning(paste("Facility",exclude_fac,"has 0 events or all events or small size (size <", size_cut,").\n"))
  }
  
  Y_cut = Y_sort[include]
  z_cut = z_sort[include,]
  Fac_cut = Fac_sort[include]
  m2_cut = table(Fac_cut)
  
  
  Est = SRR_profile(Y = Y_cut, facility = Fac_cut, z = z_cut, criterion = criterion, max.iter = max.iter, bound = bound)
  #list(beta.est = beta.est, gamma.est = gamma.est, llk_all = llk_all, rep.n = rep.n-1)
  beta.est = Est$beta.est
  
  gamma.est = Est$gamma.est
  iterN = Est$rep.n
  
  gamma_cut_long = rep(gamma.est,m2_cut)
  
  all_Fac_names = names(m2)
  all_Fac_effect =  (Y_sum_long/size_long >= 0.5)*bound - (Y_sum_long/size_long <  0.5)*bound
  all_Fac_effect[include] = gamma_cut_long
  all_gamma = sapply(split(all_Fac_effect,Fac_sort), mean) 
  
  if(is.null(refer_fac)){
    proposed_gamma = median(gamma.est)
  }else{
    proposed_gamma = gamma.est[names(m2_cut) == refer_fac]
  }
  test_rslt = P_value_score_test(proposed_gamma, Y = Y_cut, z = z_cut, Fac_ind = Fac_cut, beta = beta.est)
  
  Observed = as.numeric( sapply( split(Y_cut, factor(Fac_cut)),sum)) 
  Expected =  as.numeric( sapply( split(plogis(proposed_gamma + z_cut%*%beta.est), factor(Fac_cut)),sum))#Amy changed Median gamma to reference gamma
  
  srr = Observed/Expected
  names(srr) = sort(unique(Fac_cut))#split() with factor() will sort it-self
  smry = list(srr=srr, test_rslt = test_rslt, beta=beta.est, gamma = all_gamma, iterN=iterN, Observed=Observed, Expected=Expected, Fnames=all_Fac_names)
  return(smry)
}




# to use this function, data should be sorted by facilities.
SRR_profile <- function(Y = Y, facility = facility, z = z, criterion = 1e-8, max.iter = 200, bound = 8){
  
  N = length(Y)
  nF = length(unique(facility))
  size = sapply(split(Y,facility), length)
  gamma.est = rep(0, nF)  #dim(gamma.est)=1*F
  beta.length = ncol(z)
  beta.est = rep(0,beta.length)
  # Method 1 
  rep.n = 0
  llk_all = NULL
  repeat{ #repeat2
    rep.n = rep.n+1
    if(rep.n > max.iter){
      break
    }
    print(rep.n)
    gamma.est.long = rep(gamma.est, size)  
    p = plogis( gamma.est.long+z%*%beta.est )  #dim=n*1
    llk = sum(log(p))
    llk_all = c(llk_all,llk)
    q = 1-p
    p.times.q = p*q
    gamma.update = sapply(split(Y-p,facility), sum)/sapply(split(p.times.q, facility), sum)
    gamma.est = gamma.est + gamma.update  #dim(gamma.est)=1*F 
    gamma.est.long = rep(gamma.est, times = size)  #dim(gamma.est.long)=n*1
    #####  beta
    p = plogis( gamma.est.long+z%*%beta.est )
    q = 1-p
    p.times.q = p*q
    beta.U = t(Y-p)%*%z
    beta.I = t(z)%*%(z*c(p.times.q))
    beta.update = solve(beta.I) %*% t(beta.U)
    beta.est = beta.est +  beta.update#dim=3*1
    gamma.est[gamma.est > bound] = bound
    gamma.est[gamma.est < -bound] = bound
    beta.est[beta.est > bound] = bound
    beta.est[beta.est < -bound] = bound
    
    if( max(abs(c(beta.update))) < criterion ){
      break
    } ## end if
    
  } ## end repeat2
  
  return(list(beta.est = beta.est, gamma.est = gamma.est, llk_all = llk_all, rep.n = rep.n-1))
}


P_value_score_test <- function(proposed_gamma, Y, z, Fac_ind, beta){
  p_ik = 1-1/(exp(z%*%beta+proposed_gamma)+1)
  p_i = sapply(split(p_ik,Fac_ind),sum)
  y_i = sapply(split(Y,Fac_ind),sum)
  sigma_i = sqrt(sapply(split(p_ik*(1-p_ik),Fac_ind),sum))
  z_stat = (y_i-p_i)/sigma_i
  pv1 = pnorm(z_stat,lower.tail = FALSE)
  pv2 = pnorm(z_stat,lower.tail = TRUE)
  Fac_flag =  rep(0, length(unique(Fac_ind)))
  Fac_flag[pv1<0.025] = 1
  Fac_flag[pv2<0.025] = -1
  pvalue = (pv1<pv2)*2*pv1 + (pv1>=pv2)*2*pv2  # 2*min
  rslt = cbind(pvalue,Fac_flag,z_stat)
  colnames(rslt) = c("Pvalue","flag","z_stat")
  rownames(rslt) = unique(Fac_ind)
  return(rslt)
}

