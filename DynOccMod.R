######################################################################
###
### Bayesian hierachical dynamic occupancy model for raccoon rabies
###
### March 22, 2019
### Amy J. Davis
###
######################################################################

RDOcc.gcovs.negts.mcmc<-function(y,n,Xgam,Xlocs,gridnames,catnames,xgamnames,gam.tune,GWtQu=FALSE,Gtimesince=TRUE,
                                 a.psi=1,b.psi=1,a.e=1,b.e=1,n.mcmc){
  
  ### y = array of number of rabies positive animals by site, season, and method
  ### n = array of number of animals samples by site, season, and method
  ### Xlocs = data fram of Grid ID values and associated latitude and longitude centroids
  ### catnames = vector of names for the surveillance categories
  ### xgamnames = vector names of covariates on colonization parameter (gamma)
  ### xgam = matrix of covariates for colonization, including at minimum an intercept, columns should match xgamnames, dimensions should be the same as y but one fewer columns reflecting transitions
  ### gam.tune = tuning parameter for gamma, can be a single value or a vector of the same length as xgamnames
  ### GWtQU = TRUE/FALSE value as to whether the neighbor effect should be considered (Queen's neighbor)
  ### Gtimesince = TRUE/FALSE value as to whether the time since a cell was occupied should be included as a covariate
  ### a.psi = prior value for alpha for psi1 (initial occupancy)
  ### b.psi = prior value for beta for psi1 (initial occupancy)
  ### a.e = prior value for alpha for epsilon (extinction probability)
  ### b.e = prior value for beta for epsilon (extinction probability) 
  ### n.mcmc = number of MCMC iterations
  
  
  ### Libraries
  library(igraph)
  library(akima)
  library(boot)
  library(MASS)
  library(mvtnorm)
  library(fields)
  library(abind)
  library(geosphere)
  library(Hmisc)
  
  gettimesince<-function(z){
    zlen=apply(z,1,rle)
    zts=t(sapply(zlen,FUN=function(y)unlist(lapply(y$lengths,FUN=function(x)seq(1,x,1)))))
    zts[which(z==1,arr.ind = TRUE)]=0
    zts
  }
  
  
  ###
  ###  Data manipulation
  ###
  mths=dim(y)[2]
  sites=dim(y)[1]
  ncats=dim(y)[3]
  Xgam2=Xgam
  tpos=apply(y,c(1,2),sum)
  tneg=apply(n-y,c(1,2),sum)
  xgamnames2=xgamnames
  catrep=rep(1:ncats,each=sites)
  siterep=rep(1:sites,ncats)
  sitemrep=rep(1:sites,n.mcmc)
  nw=which(n>0,arr.ind = TRUE)
  
  ### Get n and y into a long format to reduce memory needs
  nl=apply(n,2L,c)
  yl=apply(y,2L,c)
  
  xgamvals=dim(Xgam2)[3]
  Xlocs=as.matrix(Xlocs[,c("Longitude","Latitude")])
  XlocsUsedB=Xlocs[gridnames,]
  
  #####################
  ###
  ### Starting values
  ###
  #####################
  y0=as.matrix(aggregate(yl,by=list(cat=siterep),FUN=sum)[,-1])
  y0=ifelse(y0>0,1,0)
  z=ifelse(y0>0,1,rbinom(sites*mths,1,0.1))
  
  ### Neighbor effects
  sitedist=distm(Xlocs[,c(1,2)])
  ndist=sort(unique(c(round(sitedist/1000,0))))[2]*2000
  siteneighb=ifelse(sitedist>ndist,0,1)
  diag(siteneighb)=1
  siteneighb=siteneighb[gridnames,gridnames]
  
  
  ###
  ### Neighbor effects
  ###
  gneighb=apply(z[,-mths],2,function(w)colSums(w*siteneighb))/rowSums(siteneighb)
  zneb=ifelse(gneighb>0,1,0)
  
  ###
  ### Get the time lag x matrix
  ###
  xgamts=array(gettimesince(zneb),c(sites,mths-1,1))
  xgamts2=xgamts/(mths-1)
  
  ###
  ### Location colonization neighbor effects
  ###
  Xgam=Xgam2
  if(GWtQu==TRUE){
    Xgam=abind(Xgam2,gneighb)
    xgamnames=c(xgamnames,"Neighbors")
  }
  if(Gtimesince==TRUE){
    Xgam=abind(Xgam,xgamts2)
    xgamnames=c(xgamnames,"TimeSince")
  }
  xgamvals=dim(Xgam)[3]
  gammabetas=rnorm(xgamvals,0,1)
  gam=inv.logit(apply(Xgam*array(rep(gammabetas,each=(sites*(mths-1))),c(sites,(mths-1),xgamvals)),c(1,2),sum))  
  Xgammean=apply(Xgam,3,mean)
  
  
  ### Epsilon start
  eps=rbeta(1,a.e,b.e)
  
  ## Start for p 
  if(dim(n)[3]>1){
    pint=which(array(apply(n,3,function(x)x*z),c(sites,mths,ncats))!=0,arr.ind = TRUE)
    p=tapply((y/array(apply(n,3,function(x)x*z),c(sites,mths,ncats)))[pint],pint[,3],mean)
  }else{
    p=as.matrix(sum(y)/sum(n))
  }
  
  ## Start for Psi
  psi=matrix(rbeta(sites*mths,z+1,1-z+1),sites,mths)
  
  ###
  n.burn=round(n.mcmc/2)           ### Create a burn in portion 
  
  ### 
  ### Save samples
  ###
  psi.save=matrix(0,sites*n.mcmc,mths)
  z.save=matrix(0,sites*n.mcmc,mths)
  eps.save=rep(0,n.mcmc)
  gam.save=matrix(0,sites*n.mcmc,mths-1)
  betagam.save=matrix(0,xgamvals,n.mcmc)
  p.save=matrix(0,ncats,n.mcmc)
  elim.save=matrix(0,sites*n.mcmc,mths)
  
  waicsPDmatrix=matrix(0,dim(nw)[1],n.mcmc)
  
  ll1=rep(1,n.mcmc)
  accept.gam <- matrix(0,xgamvals,n.mcmc)
  gam.tunesave=matrix(gam.tune,xgamvals,n.mcmc)
  
  
  ###
  ### MCMC loop
  ###
  for(k in 1:n.mcmc){
    if(k%%100==0)cat(k," ");flush.console()
    
    #########################################################################
    ###
    ### Sample from psi 1 
    ###
    #########################################################################
    
    psi[,1]=rbeta(1,sum(z[,1])+a.psi,sum(1-z[,1])+b.psi)
    
    
    #########################################################################
    ###
    ### Sample from epsilon, calculate summary stats
    ###
    #########################################################################
    
    pind=which(z[,-mths]==1,arr.ind = TRUE)
    p2ind=pind+matrix(c(0,1),dim(pind)[1],dim(pind)[2],byrow=TRUE)
    
    table(z[p2ind])
    length(z[p2ind]==1)
    
    eps=rbeta(1,sum(1-z[p2ind])+a.e,sum(z[p2ind])+b.e)
    eps
    
    #########################################################################
    ###
    ### Sample from gamma, calculate summary stats
    ###
    #########################################################################
    ### Uncomment this chunk to autotune to a given acceptance rate
    # if(k%%100==0){
    #   if(dim(accept.gam)[1]>1){
    #     a.r=rowMeans(accept.gam[,(k-100):(k)])
    #     gam.tune=a.r/0.7*gam.tune
    #   }else{
    #     a.r=rowMeans(accept.gam[,(k-100):(k)])
    #     gam.tune=a.r/0.7*gam.tune
    #   }
    # }
    ###
    ### Neighbor effects
    ###
    gneighb=apply(z[,-mths],2,function(w)colSums(w*siteneighb))/rowSums(siteneighb)
    zneb=ifelse(gneighb>0,1,0)
    
    ###
    ### Get the time lag x matrix
    ###
    xgamts=array(gettimesince(zneb),c(sites,mths-1,1))
    xgamts2=(xgamts/39)
    
    Xgam=Xgam2
    ###
    ### Location colonization neighbor effects
    ###
    if(GWtQu==TRUE){
      Xgam=abind(Xgam2,gneighb)
    }
    if(Gtimesince==TRUE){
      Xgam=abind(Xgam,xgamts2)
    }
    
    gind=which((1-z[,-mths])==1,arr.ind = TRUE)
    g2ind=gind+matrix(c(0,1),dim(gind)[1],dim(gind)[2],byrow=TRUE)
    gsp=which(z[g2ind]==1,arr.ind = TRUE)
    a.g=rep(1,dim(gind)[1])
    a.g[gsp]=2
    b.g=rep(1,dim(gind)[1])
    b.g[-gsp]=2 
    
    ### 
    ### Sample from gamma0 using Metropolis-Hastings
    ###    
    gambetastar=rnorm(xgamvals,gammabetas,gam.tune)
    gammastar=inv.logit(apply(Xgam*array(rep(gambetastar,each=(sites*(mths-1))),c(sites,(mths-1),xgamvals)),c(1,2),sum,na.rm=TRUE))
    gam=inv.logit(apply(Xgam*array(rep(gammabetas,each=(sites*(mths-1))),c(sites,(mths-1),xgamvals)),c(1,2),sum,na.rm=TRUE))
    
    gammastar[gammastar>0.999]=0.999
    gammastar[gammastar<0.001]=0.001
    gam[gam>0.999]=0.999
    gam[gam<0.001]=0.001
    
    gamMHratio=exp(sum(dbinom(a.g-1,1,gammastar[gind],log=TRUE))+dnorm(gambetastar,0,1,log=TRUE)-
                     sum(dbinom(a.g-1,1,gam[gind],log=TRUE))-dnorm(gammabetas,0,1,log=TRUE))
    
    tmp.keep=gamMHratio>runif(xgamvals)
    gammabetas[tmp.keep]=gambetastar[tmp.keep]
    gammabetas
    accept.gam[,k]=tmp.keep*1
    gam.tunesave[,k]=gam.tune
    
    gam=inv.logit(apply(Xgam*array(rep(gammabetas,each=(sites*(mths-1))),c(sites,(mths-1),xgamvals)),c(1,2),sum))
    
    
    #########################################################################
    ###
    ### Sample from psi* and z
    ###
    #########################################################################
    ### Detection by site and time
    parray=array(rep(p,each=sites*mths),c(sites,mths,ncats))
    pe=(1-parray)^n
    ptot=1-apply(pe,c(1,2),prod)
    epsm=matrix(eps,sites,mths-1)
    for(t in 1:mths){
      ###
      ### Sample from z (only when y=0)
      ###
      indy0=which(y0[,(t)]==0)
      
      ### For t=1
      if(t==1){
        psicond=z[,2]*(1-epsm[,1])/((1-epsm[,1])+gam[,1])+(1-z[,2])*(1-(1-epsm[,1]))/((1-(1-epsm[,1]))+(1-gam[,1]))
        psicond[is.na(psicond)]=0
        z[indy0,1]=rbinom(length(indy0),1,(1-ptot[indy0,1])*psicond[indy0]/((1-ptot[indy0,1])*psicond[indy0]+(1-psicond[indy0])))
        
      } else 
        if(t<mths){
          psicond=z[,(t-1)]*z[,(t+1)]*(1-epsm[,t-1])*(1-epsm[,t])/((1-epsm[,t-1])*(1-epsm[,t])+(1-(1-epsm[,t-1]))*gam[,t])+
            z[,(t-1)]*(1-z[,(t+1)])*(1-epsm[,t-1])*(1-(1-epsm[,t]))/((1-epsm[,t-1])*(1-(1-epsm[,t]))+(1-(1-epsm[,t-1]))*(1-gam[,t]))+
            (1-z[,(t-1)])*(z[,(t+1)])*gam[,t-1]*(1-epsm[,t])/(gam[,t-1]*(1-epsm[,t])+(1-gam[,t-1])*gam[,t])+
            (1-z[,(t-1)])*(1-z[,(t+1)])*gam[,t-1]*(1-(1-epsm[,t]))/(gam[,t-1]*(1-(1-epsm[,t]))+(1-gam[,t-1])*(1-gam[,t]))
          psicond[is.na(psicond)]=0
          z[indy0,t]=rbinom(length(indy0),1,(1-ptot[indy0,t])*psicond[indy0]/((1-ptot[indy0,t])*psicond[indy0]+(1-psicond[indy0])))
        }else {
          psicond=ifelse(z[,mths]==1,(1-epsm[,t-1])/((1-epsm[,t-1])+(1-(1-epsm[,t-1]))),gam[,t-1]/(gam[,t-1]+(1-gam[,t-1])))
          psicond[is.na(psicond)]=0
          z[indy0,mths]=rbinom(length(indy0),1,(1-ptot[indy0,t])*psicond[indy0]/((1-ptot[indy0,t])*psicond[indy0]+(1-psicond[indy0])))
        }
    } 
    
    ###
    ### Calculate psi for each run
    ### 
    psi[,-1]=(1-eps)*z[,-mths]+gam*(1-z[,-mths])
    
    
    ###
    ###  Sample from p (prevelance parameter for disease if present), when z = 1
    ###
    z1ind=which(z==1,arr.ind = TRUE)
    p=rbeta(ncats,colSums(apply(y,3,function(w)w[z1ind]))+1,colSums(apply(n,3,function(w)w[z1ind])-apply(y,3,function(w)w[z1ind]))+1)
    
    
    ###
    ### Save samples
    ###
    betagam.save[,k]=gammabetas
    
    eps.save[k]=eps
    gam.save[(sites*(k-1)+1):(sites*(k-1)+sites),]=gam
    
    p.save[,k]=p
    psi.save[(sites*(k-1)+1):(sites*(k-1)+sites),]=psi
    z.save[(sites*(k-1)+1):(sites*(k-1)+sites),]=z
    pstar=1-apply(n,2,function(x) apply((1-matrix(p,sites,ncats,byrow=TRUE))^x,1,prod))
    elim.save[(sites*(k-1)+1):(sites*(k-1)+sites),]=(1-psi)/((1-psi)+psi*(1-pstar))
    
    waicsPDmatrix[,k]=dbinom(y[cbind(nw[,1],nw[,2],nw[,3])],n[cbind(nw[,1],nw[,2],nw[,3])],
                             p[nw[,3]]*z[cbind(nw[,1],nw[,2])],log=TRUE)
  }
  
  
  cat("\n")
  
  ####
  ####  Write Output 
  ####
  
  list(z=z.save,psi=psi.save,p=p.save,eps=eps.save,gammas=gam.save,
       betagamma=betagam.save)
  
}










