# This R code contains two functions:

# 1) "Nuclear" function that fits the integrative latent class model using training data
# 2) "Nuclear_predict" function for predicting the obtruction status status of new kidney

# written by Changgee Chang and Jeong Hoon Jang
# ver 20230710

#### Install Dependent Packages #####
if ( !require(msm) )
{
  install.packages("msm")
  library(msm)
}

if ( !require(tmvtnorm) )
{
  install.packages("tmvtnorm")
  library(tmvtnorm)
}

if ( !require(abind) )
{
  install.packages("abind")
  library(abind)
}

if ( !require(mvtnorm) )
{
  install.packages("mvtnorm")
  library(mvtnorm)
}


##################################################################################################
#### Function "Nuclear" that trains the model
# Arguments:
# y1, y2: J1 (and J2) by n observation matrix of pre-furosemide (and post-furosemide)
# sbj: n by 1 vector of subject index from 1 to M.
# X1: r1 by M covariates matrix for subject specified features
# X2: r2 by n covariates matrix for kidney specific features
# X3: r3 by n covariates matrix for latent factors
# W: L by n matrix; Expert Ratings; Values 0, 1, or 2
# t1, t2: J1 (and J2) time points at which the data are collected for pre- (and post-) furosemide
# p, nu, psi: spline basis function parameters; p is the number of basis functions, and see the technical documents for nu and psi
# q: c(q1,q2)
# ups : tuning parameter upsilon
# a: tuning parameters a = c(a1,a2). If NULL, adaptive with initial c(2,2)
# a_zeta, b_zeta: 
# a_iota, b_iota: 
# a_xi, b_xi: 
# tau: expert variability parameter
# kappa:
# theta: sparsity tuning parameter


Nuclear <- function(nsample,sbj,y1,y2,X1,X2,X3,W,t1,t2,p,q,nu,psi,ups,a,a_sigma,b_sigma,a_zeta,b_zeta,a_iota,b_iota,a_xi,b_xi,tau,kappa,theta)
{
  # number of subjects
  M = ncol(X1)
  
  # number of kidneys
  n = ncol(X2)
  n2M = matrix(0,n,M)
  for ( i in 1:n )
    n2M[i,sbj[i]] = 1
  nkid = colSums(n2M)
  
  # number of clinical features for subjects
  r1 = nrow(X1)
  
  # number of clinical features for kidneys
  r2 = nrow(X2)
  
  # number of features for latent factors
  r3 = nrow(X3)
  XX3 = X3 %*% t(X3)
  
  # Concatenated X
  X = rbind(X1[,sbj],X2)
  XX = X %*% t(X)
  r = r1+r2
  
  # number of experts
  L = nrow(W)
  
  # number of time points
  J1 = length(t1)
  J2 = length(t2)
  
  # matrices for kernel
  A1 = matrix(1,J1,p)
  A2 = matrix(1,J2,p)
  for ( k in 2:p )
  {
    A1[,k] = exp(-nu*(t1-(k-1)*psi)^2)
    A2[,k] = exp(-nu*(t2-(k-1)*psi)^2)
  }
  
  AA = array(0,c(p,p,2))
  AA[,,1] = t(A1)%*%A1
  AA[,,2] = t(A2)%*%A2
  
  Ay = array(0,c(p,n,2))
  Ay[,,1] = t(A1)%*%y1
  Ay[,,2] = t(A2)%*%y2
  
  
  K = diag(c(1,rep(2,p-2),1),p)
  K[0:(p-2)*p+2:p] = -1
  K[1:(p-1)*(p+1)] = -1
  
  # FF is the matrix F
  tmp = rep(1:p,p)-rep(1:p,each=p)
  FF = matrix(sqrt(pi/nu/2)*exp(-nu*psi^2/2*tmp^2),p,p)
  tF = t(FF)
  
  # adaptive a
  if ( is.null(a) )
  {
    ada = TRUE
    a = c(2,2)
  }
  else
    ada = FALSE
  
  # initialize C
  sumW = apply(W,2,sum)
  C = as.integer(sumW>L)
  U = (sumW==0 | sumW==2*L) # unanimous cases. prevents flip of two classes C = 0,1 during burn-in
  pC = rep(0,n)
  
  # initialize O
  O = rep(0,M)
  for ( m in 1:M )
    O[m] = sum(C[sbj==m])
  
  # initialize H and N
  H1 = matrix(0,M,3)
  H1[O==0,1] = 1
  H1[O==1,2] = 1
  H1[O==2,3] = 1
  N1 = apply(H1,2,sum)
  
  H2 = matrix(c(1-C,C),n)
  Cwh = H2==1
  N2 = apply(H2,2,sum)
  
  # initialize iota
  iota = array(0,c(p,M,2))
  
  # initialize sigma2
  beta = array(0,c(p,n,2))
  for ( s in 1:2 )
    beta[,,s] = chol2inv(chol(AA[,,s]+diag(0.01,p)))%*%Ay[,,s]
  eps1 = y1 - A1%*%beta[,,1]
  eps2 = y2 - A2%*%beta[,,2]
  sigma2 = (1+sum(eps1^2)+sum(eps2^2))/(1+n*(J1+J2))
  
  # initialize Lambda and gamma
  lam = array(0,c(p,max(q),2,2))
  gamma = array(0,c(max(q),r3,2))
  iXX30 = chol2inv(chol(X3%*%(t(X3)*H2[,1])+diag(0.01,r3)))
  iXX31 = chol2inv(chol(X3%*%(t(X3)*H2[,2])+diag(0.01,r3)))
  for ( s in 1:2 )
  {
    tmp0 = beta[,,s]%*%(t(X3)*H2[,1])%*%iXX30
    tmp1 = beta[,,s]%*%(t(X3)*H2[,2])%*%iXX31
    tmp = svd(rbind(tmp0,tmp1),q[s],q[s])
    lam[,1:q[s],1,s] = tmp$u[1:p,] %*% diag(tmp$d/sqrt(r3),q[s])
    lam[,1:q[s],2,s] = tmp$u[1:p+p,] %*% diag(tmp$d/sqrt(r3),q[s])
    if ( q[s] <= r3 )
      gamma[1:q[s],,s] = t(tmp$v)*sqrt(r3)
    else
    {
      gamma[1:r3,,s] = t(tmp$v)*sqrt(r3)
      gamma[(r3+1):q[s],,s] = rnorm((q[s]-r3)*r3)
    }    
  }    
  
  # initialize eta
  eta = array(0,c(max(q),n,2))
  for ( s in 1:2 )
    eta[1:q[s],,s] = gamma[1:q[s],,s]%*%X3
  
  # initialize beta
  for ( s in 1:2 )
  {
    V = chol2inv(chol(AA[,,s]/sigma2+diag(a_sigma/b_sigma,p)))
    for ( i in 1:n )
      beta[,i,s] =  V %*% (Ay[,i,s]/sigma2+a_sigma*lam[,1:q[s],C[i]+1,s]%*%eta[1:q[s],i,s]/b_sigma)
  }
  
  # initialize mu  
  mu = array(0,c(max(q),2,2))
  mu[1,,] = a[1]
  pi = array(a[1],c(max(q),2,2))
  for ( s in 1:2 )
    if ( q[s]>1 )
    {
      mu[2:q[s],,s] = a[2]
      for ( h in 2:q[s] )
        pi[h,,s] = pi[h-1,,s]*mu[h,,s]
    }
  
  # initialize phi    
  phi = (ups+1) / (ups+lam^2*rep(pi,each=p))
  
  # initialize omega
  omega = 2/(gamma^2+1)
  
  # initialize sig2_zeta
  sig2_zeta = 0
  for ( s in 1:2 )
    for ( c in 0:1 )
    {
      zeta = beta[,,s] - lam[,1:q[s],c+1,s] %*% eta[1:q[s],,s]
      sig2_zeta = sig2_zeta + sum(t(zeta)^2*H2[,c+1])
    }
  sig2_zeta = (b_zeta+sig2_zeta) / (a_zeta+2*n*p)
  
  # initialize sig2_iota
  sig2_iota = b_iota/a_iota
  
  # initialize Z
  Z = W-1/2
  
  # initialize b
  b = matrix(0,L,4)
  b[,1] = -Inf
  b[,3] = 1 # initial value
  b[,4] = Inf
  
  mu_b = 1
  sig2_b = 1/tau
  
  # initialize rho
  iXX0 = chol2inv(chol(X%*%(t(X)*H2[,1])+diag(0.01,r)))
  iXX1 = chol2inv(chol(X%*%(t(X)*H2[,2])+diag(0.01,r)))
  rho = array(0,c(L,r,2))
  rho[,,1] = (Z %*% (t(X)*H2[,1])) %*% iXX0
  rho[,,2] = (Z %*% (t(X)*H2[,2])) %*% iXX1
  rhoX = array(0,c(L,n,2))
  rhoX[,,1] = rho[,,1] %*% X
  rhoX[,,2] = rho[,,2] %*% X
  
  mu_rho = apply(rho,c(2,3),mean)
  sig2_rho = matrix(1/tau,r,2)
  
  # initialize xi
  xi = array(0,c(p,L,2,2))
  iSig_xi = matrix(0,2*p,2*p)
  mu_xi = matrix(0,2*p,L)
  sig2_xi = b_xi / a_xi
  
  Fbeta1 = FF%*%beta[,,1]
  Fbeta2 = FF%*%beta[,,2]
  for ( c in 0:1 )
  {
    iSig_xi[1:p,1:p] = Fbeta1%*%(t(Fbeta1)*H2[,c+1])
    iSig_xi[1:p,1:p+p] = Fbeta1%*%(t(Fbeta2)*H2[,c+1])
    iSig_xi[1:p+p,1:p] = Fbeta2%*%(t(Fbeta1)*H2[,c+1])
    iSig_xi[1:p+p,1:p+p] = Fbeta2%*%(t(Fbeta2)*H2[,c+1])
    iSig_xi[1:p,1:p] = iSig_xi[1:p,1:p] + K/sig2_xi
    iSig_xi[1:p+p,1:p+p] = iSig_xi[1:p+p,1:p+p] + K/sig2_xi
    mu_xi[1:p,] = Fbeta1%*%(t(Z-rhoX[,,c+1])*H2[,c+1])
    mu_xi[1:p+p,] = Fbeta2%*%(t(Z-rhoX[,,c+1])*H2[,c+1])
    ciSig_xi = chol(iSig_xi)
    tmp = backsolve(ciSig_xi,forwardsolve(t(ciSig_xi),mu_xi))
    xi[,,c+1,1] = tmp[1:p,]
    xi[,,c+1,2] = tmp[1:p+p,]
  }
  
  # initialize sig2_xi
  sig2_xi = 0
  for ( s in 1:2 )
    for ( c in 0:1 )
      sig2_xi = sig2_xi + sum(diag(t(xi[,,c+1,s])%*%K%*%xi[,,c+1,s]))
  sig2_xi = (b_xi+sig2_xi) / (a_xi+4*(p-1)*L)
  
  # initialize delta and S
  delta1 = X1%*%H1/rep(N1,each=r1)
  delta2 = X2%*%H2/rep(N2,each=r2)
  S1 = (diag(kappa,r1) + (X1-delta1[,O+1])%*%t(X1-delta1[,O+1])) / (kappa+M)
  S2 = (diag(kappa,r2) + (X2-delta2[,C+1])%*%t(X2-delta2[,C+1])) / (kappa+n)
  iS1 = chol2inv(chol(S1))
  iS2 = chol2inv(chol(S2))
  
  # repositories
  rep_beta = array(0,c(p,n,2,nsample))
  rep_lam = array(0,c(p,max(q),2,2,nsample))
  rep_eta = array(0,c(max(q),n,2,nsample))
  rep_iota = array(0,c(p,M,2,nsample))
  rep_gamma = array(0,c(max(q),r3,2,nsample))
  rep_C = matrix(1,n,nsample)
  rep_pC = matrix(1,n,nsample)
  rep_phi = array(1,c(p,max(q),2,2,nsample))
  rep_pi = array(1,c(max(q),2,2,nsample))
  rep_omega = array(1,c(max(q),r3,2,nsample))
  rep_sig2_zeta = rep(0,nsample)
  rep_sig2_iota = rep(0,nsample)
  rep_sig2_xi = rep(0,nsample)
  rep_sigma2 = rep(0,nsample)
  rep_mu = array(1,c(max(q),2,2,nsample))
  
  rep_Z = array(0,c(L,n,nsample))
  rep_b = matrix(0,L,nsample)
  rep_mu_b = rep(0,nsample)
  rep_sig2_b = rep(0,nsample)
  rep_rho = array(0,c(L,r,2,nsample))
  rep_mu_rho = array(0,c(r,2,nsample))
  rep_sig2_rho = array(0,c(r,2,nsample))
  rep_xi = array(0,c(p,L,2,2,nsample))
  rep_delta1 = array(0,c(r1,3,nsample))
  rep_delta2 = array(0,c(r2,2,nsample))
  rep_S1 = array(0,c(r1,r1,nsample))
  rep_S2 = array(0,c(r2,r2,nsample))
  rep_iS1 = rep_S1
  rep_iS2 = rep_S2
  
  if ( ada )
    rep_a = matrix(0,2,nsample)
  else
    rep_a = a
  
  
  for ( it in 1:nsample )
  {
    # beta
    for ( c in 0:1 )
    {
      Fxi = rbind(tF%*%xi[,,c+1,1],tF%*%xi[,,c+1,2])
      iV = Fxi%*%t(Fxi)
      iV[1:p,1:p] = iV[1:p,1:p] + AA[,,1]/sigma2 + diag(1/sig2_zeta,p)
      iV[1:p+p,1:p+p] = iV[1:p+p,1:p+p] + AA[,,2]/sigma2 + diag(1/sig2_zeta,p)
      ciV = chol(iV)
      iVM1 = Ay[,Cwh[,c+1],1]/sigma2 + (lam[,,c+1,1]%*%eta[,Cwh[,c+1],1] + iota[,sbj[Cwh[,c+1]],1])/sig2_zeta
      iVM2 = Ay[,Cwh[,c+1],2]/sigma2 + (lam[,,c+1,2]%*%eta[,Cwh[,c+1],2] + iota[,sbj[Cwh[,c+1]],2])/sig2_zeta
      iVM = rbind(iVM1,iVM2)
      iVM = iVM + Fxi%*%(Z-rhoX[,,c+1])[,Cwh[,c+1]]
      MM = backsolve(ciV,forwardsolve(t(ciV),iVM))
      tmp = MM + backsolve(ciV,matrix(rnorm(2*p*N2[c+1]),2*p))
      beta[,Cwh[,c+1],1] = tmp[1:p,]
      beta[,Cwh[,c+1],2] = tmp[p+1:p,]
    }
    Fbeta1 = FF%*%beta[,,1]
    Fbeta2 = FF%*%beta[,,2]
    
    
    # eta
    for ( s in 1:2 )
      for ( c in 0:1 )
      {
        iV = t(lam[,1:q[s],c+1,s]) %*% (lam[,1:q[s],c+1,s]/sig2_zeta) + diag(1,q[s])
        ciV = chol(iV)
        iVM = t(lam[,1:q[s],c+1,s]) %*% (beta[,Cwh[,c+1],s] - iota[,sbj[Cwh[,c+1]],s])/sig2_zeta + gamma[1:q[s],,s]%*%X3[,Cwh[,c+1]]
        MM = backsolve(ciV,forwardsolve(t(ciV),iVM))
        eta[1:q[s],Cwh[,c+1],s] = MM + backsolve(ciV,matrix(rnorm(q[s]*N2[c+1]),q[s]))
      }
    
    # iota
    for ( s in 1:2 )
    {
      iVM = beta[,,s] - lam[,1:q[s],1,s]%*%eta[1:q[s],,s]
      iVM[,Cwh[,2]] = beta[,Cwh[,2],s] - lam[,1:q[s],2,s]%*%eta[1:q[s],Cwh[,2],s]
      iVM = iVM%*%n2M / sig2_zeta
      MM = iVM / rep(nkid/sig2_zeta+1/sig2_iota,each=p)
      iota[,,s] = MM + rnorm(p*M) / rep(sqrt(nkid/sig2_zeta+1/sig2_iota),each=p)
    }
    
    # lambda
    for ( s in 1:2 )
      for ( c in 0:1 )
        for ( k in 1:p )
        {
          iV = eta[1:q[s],,s]%*%(t(eta[1:q[s],,s])*H2[,c+1])/sig2_zeta + diag(phi[k,1:q[s],c+1,s]*pi[1:q[s],c+1,s],q[s])
          iVM = eta[1:q[s],,s]%*%((beta[k,,s]-iota[k,sbj,s])*H2[,c+1])/sig2_zeta
          ciV = chol(iV)
          MM = backsolve(ciV,forwardsolve(t(ciV),iVM))
          lam[k,1:q[s],c+1,s] = MM + backsolve(ciV,rnorm(q[s]))
        }
    
    # gamma
    for ( s in 1:2 )
      for ( h in 1:q[s] )
      {
        iV = XX3 + diag(omega[h,,s],r3)
        iVM = X3%*%eta[h,,s]
        ciV = chol(iV)
        MM = backsolve(ciV,forwardsolve(t(ciV),iVM))
        gamma[h,,s] = MM + backsolve(ciV,rnorm(r3))
      }
    
    # phi
    for ( s in 1:2 )
      phi[,1:q[s],,s] = rgamma(p*q[s]*2,(ups+1)/2,(ups+lam[,1:q[s],,s]^2*rep(pi[1:q[s],,s],each=p))/2)
    
    # mu
    for ( s in 1:2 )
      for ( c in 0:1 )
      {
        if ( q[s] > 1 )
          philam2 = apply(phi[,1:q[s],c+1,s]*lam[,1:q[s],c+1,s]^2,2,sum)
        else
          philam2 = sum(phi[,1,c+1,s]*lam[,1,c+1,s]^2)
        
        for ( d in 1:q[s] )
        {
          pibar = pi[d:q[s],c+1,s]/mu[d,c+1,s]
          if ( ada == FALSE )
          {
            if ( d==1 )
              mu[d,c+1,s] = rgamma(1,a[1]+p*q[s]/2,1+sum(pibar*philam2)/2)
            else
              mu[d,c+1,s] = rgamma(1,a[2]+p*(q[s]-d+1)/2,1+sum(pibar*philam2[d:q[s]])/2)
          }
          else
          {
            if ( d==1 )
              mu[d,c+1,s] = rgig(1,4*a[1],4/a[1]+sum(pibar*philam2),p*q[s]/2-1/2)
            else
              mu[d,c+1,s] = rgig(1,4*a[2],4/a[2]+sum(pibar*philam2[d:q[s]]),p*(q[s]-d+1)/2-1/2)
          }
          pi[d:q[s],c+1,s] = pibar*mu[d,c+1,s]
        }
      }
    
    
    # a
    if ( ada )
    {
      mu1 = mu[1,,]
      a[1] = rgig(1,8+4*sum(mu1),2+4*sum(1/mu1),3/2)
      mu2 = mu[-1,,]
      a[2] = rgig(1,8+4*sum(mu2),2+4*sum(1/mu2[mu2!=0]),sum(q)-5/2)
    }
    
    
    # omega
    for ( s in 1:2 )
      omega[1:q[s],,s] = rexp(q[s]*r3,(gamma[1:q[s],,s]^2+1)/2)
    
    # sig2_zeta
    sig2_zeta = 0
    for ( s in 1:2 )
      for ( c in 0:1 )
      {
        zeta = beta[,,s] - lam[,1:q[s],c+1,s] %*% eta[1:q[s],,s] - iota[,sbj,s]
        sig2_zeta = sig2_zeta + sum(t(zeta)^2*H2[,c+1])
      }
    sig2_zeta = 1/rgamma(1,(a_zeta+2*n*p)/2,(b_zeta+sig2_zeta)/2)
    
    # sig2_iota
    sig2_iota = 1/rgamma(1,(a_iota+2*M*p)/2,(b_iota+sum(iota^2))/2)
    
    # sigma2
    eps1 = y1 - A1%*%beta[,,1]
    eps2 = y2 - A2%*%beta[,,2]
    sigma2 = 1/rgamma(1,(a_sigma+n*(J1+J2))/2,(b_sigma+sum(eps1^2)+sum(eps2^2))/2)
    
    # Z
    MZ = rhoX[,,1] + t(xi[,,1,1])%*%Fbeta1 + t(xi[,,1,2])%*%Fbeta2
    MZ[,Cwh[,2]] = rhoX[,Cwh[,2],2] + t(xi[,,2,1])%*%Fbeta1[,Cwh[,2]] + t(xi[,,2,2])%*%Fbeta2[,Cwh[,2]]
    for ( l in 1:L )
      Z[l,] = rtnorm(n,MZ[l,],1,b[l,W[l,]+1],b[l,W[l,]+2])
    
    # xi
    for ( c in 0:1 )
    {
      iSig_xi[1:p,1:p] = Fbeta1%*%(t(Fbeta1)*H2[,c+1])
      iSig_xi[1:p,1:p+p] = Fbeta1%*%(t(Fbeta2)*H2[,c+1])
      iSig_xi[1:p+p,1:p] = Fbeta2%*%(t(Fbeta1)*H2[,c+1])
      iSig_xi[1:p+p,1:p+p] = Fbeta2%*%(t(Fbeta2)*H2[,c+1])
      iSig_xi[1:p,1:p] = iSig_xi[1:p,1:p] + K/sig2_xi
      iSig_xi[1:p+p,1:p+p] = iSig_xi[1:p+p,1:p+p] + K/sig2_xi
      mu_xi[1:p,] = Fbeta1%*%(t(Z-rhoX[,,c+1])*H2[,c+1])
      mu_xi[1:p+p,] = Fbeta2%*%(t(Z-rhoX[,,c+1])*H2[,c+1])
      ciSig_xi = chol(iSig_xi)
      MM = backsolve(ciSig_xi,forwardsolve(t(ciSig_xi),mu_xi))
      tmp = MM + backsolve(ciSig_xi,rnorm(2*p))
      xi[,,c+1,1] = tmp[1:p,]
      xi[,,c+1,2] = tmp[1:p+p,]
    }
    
    # sig2_xi
    tmp = 0
    for ( s in 1:2 )
      for ( c in 0:1 )
        for ( l in 1:L )
          tmp = tmp + crossprod(xi[,l,c+1,s],K%*%xi[,l,c+1,s])   
    sig2_xi = 1/rgamma(1,(a_xi+4*(p-1)*L)/2,(b_xi+tmp)/2)
    
    # b
    for ( l in 1:L )
    {
      lb = max(0,Z[l,W[l,]==1])
      ub = min(Inf,Z[l,W[l,]==2])
      b[l,3] = rtnorm(1,mu_b,sqrt(sig2_b),lb,ub)
    }
    
    # mu_b
    mu_b = rtnorm(1,sum(b[,3])/(sig2_b+L),sqrt(sig2_b/(sig2_b+L)),0)
    while ( TRUE )
    {
      if ( runif(1) < (2*pnorm(mu_b/sqrt(sig2_b)))^(-L) )
        break
      mu_b = rtnorm(1,sum(b[,3])/(sig2_b+L),sqrt(sig2_b/(sig2_b+L)),0)
    }
    
    # sig2_b
    sig2_b = 1/rgamma(1,(tau+L)/2,(1+sum((b[,3]-mu_b)^2))/2)
    while ( TRUE )
    {
      if ( runif(1) < (2*pnorm(mu_b/sqrt(sig2_b)))^(-L) )
        break
      sig2_b = 1/rgamma(1,(tau+L)/2,(1+sum((b[,3]-mu_b)^2))/2)
    }
    
    # rho
    for ( c in 0:1 )
    {
      iV = X%*%(t(X)*H2[,c+1]) + diag(1/sig2_rho[,c+1],r)
      V = chol2inv(chol(iV))
      MiV =  t(X%*%(t(Z - t(xi[,,c+1,1])%*%Fbeta1 - t(xi[,,c+1,2])%*%Fbeta2)*H2[,c+1]) + mu_rho[,c+1]/sig2_rho[,c+1])
      MM = MiV %*% V
      cV = kronecker(t(chol(V)),diag(1,L))
      rho[,,c+1] = MM + as.vector(cV %*% rnorm(L*r))
      rhoX[,,c+1] = rho[,,c+1] %*% X
    }
    
    # mu_rho
    for ( c in 0:1 )
    {
      MM = apply(rho[,,c+1],2,sum)/(sig2_rho[,c+1]+L)
      V = sig2_rho[,c+1]/(sig2_rho[,c+1]+L)
      mu_rho[,c+1] = rnorm(r,MM,sqrt(V))
    }
    
    # sig2_rho
    for ( c in 0:1 )
      sig2_rho[,c+1] = 1/rgamma(r,(tau+L)/2,(1+apply((t(rho[,,c+1])-mu_rho[,c+1])^2,1,sum))/2)
    
    # delta
    for ( o in 0:2 )
    {
      iV1 = diag(1,r1) + N1[o+1]*iS1
      ciV1 = chol(iV1)
      iVM1 = iS1 %*% (X1%*%H1[,o+1])
      M1 = backsolve(ciV1,forwardsolve(t(ciV1),iVM1))
      delta1[,o+1] = M1 + backsolve(ciV1,rnorm(r1))
    }
    for ( c in 0:1 )
    {
      iV2 = diag(1,r2) + N2[c+1]*iS2
      ciV2 = chol(iV2)
      iVM2 = iS2 %*% (X2%*%H2[,c+1])
      M2 = backsolve(ciV2,forwardsolve(t(ciV2),iVM2))
      delta2[,c+1] = M2 + backsolve(ciV2,rnorm(r2))
    }
    
    # S
    V1 = chol2inv(chol(diag(kappa,r1) + (X1-delta1[,O+1]) %*% t(X1-delta1[,O+1])))
    V2 = chol2inv(chol(diag(kappa,r2) + (X2-delta2[,C+1]) %*% t(X2-delta2[,C+1])))
    iS1 = rWishart(1,kappa+M,V1)[,,1]
    iS2 = rWishart(1,kappa+n,V2)[,,1]
    S1 = chol2inv(chol(iS1))
    S2 = chol2inv(chol(iS2))
    
    # C
    term1 = matrix(rep(c(0,-theta),each=n),n)
    term2 = matrix(0,n,2)
    term3 = matrix(0,n,2)
    term4 = matrix(0,n,2)
    term6 = matrix(0,n,2)
    for ( c in 0:1 )
    {
      term2[,c+1] = -colSums((lam[,,c+1,1]%*%eta[,,1] + iota[,sbj,1])^2)/sig2_zeta/2
      term2[,c+1] = term2[,c+1] - colSums((lam[,,c+1,2]%*%eta[,,2] + iota[,sbj,2])^2)/sig2_zeta/2
      
      term3[,c+1] = -colSums((Z-rhoX[,,c+1])^2)/2
      
      Fxi = rbind(tF%*%xi[,,c+1,1],tF%*%xi[,,c+1,2])
      iV = Fxi%*%t(Fxi)
      iV[1:p,1:p] = iV[1:p,1:p] + AA[,,1]/sigma2 + diag(1/sig2_zeta,p)
      iV[1:p+p,1:p+p] = iV[1:p+p,1:p+p] + AA[,,2]/sigma2 + diag(1/sig2_zeta,p)
      ciV = chol(iV)
      
      iVM = Fxi%*%(Z-rhoX[,,c+1])
      iVM[1:p,] = iVM[1:p,] + Ay[,,1]/sigma2 + (lam[,,c+1,1]%*%eta[,,1] + iota[,sbj,1])/sig2_zeta
      iVM[1:p+p,] = iVM[1:p+p,] + Ay[,,2]/sigma2 + (lam[,,c+1,2]%*%eta[,,2] + iota[,sbj,2])/sig2_zeta
      
      tmp = forwardsolve(t(ciV),iVM)
      term4[,c+1] = colSums(tmp^2)/2 - sum(log(diag(ciV)))
      
      term6[,c+1] = -colSums((X2-delta2[,c+1])*(iS2%*%(X2-delta2[,c+1])))/2
    }
    allterm = term1 + term2 + term3 + term4 + term6
    
    for ( i in 1:n )
    {
      if ( it<500 & U[i] )
        next
      llk = rep(0,2)
      Omi = O[sbj[i]] - C[i]
      for ( c in 0:1 )
        llk[c+1] = allterm[i,c+1] - sum((X1[,sbj[i]]-delta1[,Omi+c+1])*(iS1%*%(X1[,sbj[i]]-delta1[,Omi+c+1])))/2
      pC[i] = 1/(1+exp(llk[1]-llk[2]))
      C[i] = rbinom(1,1,pC[i])
      O[sbj[i]] = Omi + C[i]
    }
    H1 = matrix(0,M,3)
    H1[O==0,1] = 1
    H1[O==1,2] = 1
    H1[O==2,3] = 1
    N1 = apply(H1,2,sum)
    
    H2 = matrix(c(1-C,C),n)
    Cwh = H2==1
    N2 = apply(H2,2,sum)
    
    # store
    rep_beta[,,,it] = beta
    rep_lam[,,,,it] = lam
    rep_eta[,,,it] = eta
    rep_gamma[,,,it] = gamma
    rep_C[,it] = C
    rep_pC[,it] = pC
    rep_phi[,,,,it] = phi
    rep_pi[,,,it] = pi
    rep_omega[,,,it] = omega
    rep_sig2_zeta[it] = sig2_zeta
    rep_sig2_iota[it] = sig2_iota
    rep_sig2_xi[it] = sig2_xi
    rep_sigma2[it] = sigma2
    rep_mu[,,,it] = mu
    if ( ada )
      rep_a[,it] = a
    
    rep_Z[,,it] = Z
    rep_b[,it] = b[,3]
    rep_mu_b[it] = mu_b
    rep_sig2_b[it] = sig2_b
    rep_rho[,,,it] = rho
    rep_mu_rho[,,it] = mu_rho
    rep_sig2_rho[,,it] = sig2_rho
    rep_xi[,,,,it] = xi
    
    rep_delta1[,,it] = delta1
    rep_delta2[,,it] = delta2
    rep_S1[,,it] = S1
    rep_iS1[,,it] = iS1
    rep_S2[,,it] = S2
    rep_iS2[,,it] = iS2
    
    print(it)
  }
  
  ret = list(nsample=nsample,y1=y1,y2=y2,J=c(J1,J2),n=n,r1=r1,r2=r2,r3=r3,L=L,X1=X1,X2=X2,X3=X3,W=W,t1=t1,t2=t2,
             p=p,q=q,nu=nu,psi=psi,ups=ups,tau=tau,kappa=kappa,theta=theta,
             A1=A1,A2=A2,AA=AA,
             beta=rep_beta,lam=rep_lam,eta=rep_eta,gamma=rep_gamma,C=rep_C,pC=rep_pC,phi=rep_phi,pi=rep_pi,omega=rep_omega,
             sig2_zeta=rep_sig2_zeta,sig2_iota=rep_sig2_iota,sig2_xi=rep_sig2_xi,sigma2=rep_sigma2,mu=rep_mu,a=rep_a,
             Z=rep_Z,b=rep_b,mu_b=rep_mu_b,sig2_b=rep_sig2_b,rho=rep_rho,mu_rho=rep_mu_rho,sig2_rho=rep_sig2_rho,xi=rep_xi,
             delta1=rep_delta1,delta2=rep_delta2,S1=rep_S1,iS1=rep_iS1,S2=rep_S2,iS2=rep_iS2)
}









#############################################################################
# Function "Nuclear_predict" for predicting the obstruction status of a new kidney
# Arguments:
# bi: burn-ins; the first bi samples are ignored
# y1: J1 by 2 (left and right) new observation vector of pre-furosemide curve
# y2: J2 by 2 (left and right) new observation vector of post-furosemide curve
# x1: r1 by 1 covariates vector
# x2: r2 by 2 (left and right) covariates vector
# x3: r3 by 2 (left and right) covariates vector
#

Nuclear_predict <- function(fit,bi,y1,y2,x1,x2,x3)
{
  p = fit$p
  q = fit$q
  N = fit$nsample
  Ay = abind(t(fit$A1)%*%y1,t(fit$A2)%*%y2,along=3)
  
  llk = array(0,c(N-bi,2,2))
  llk[,,2] = llk[,,2] - fit$theta
  llk[,2,] = llk[,2,] - fit$theta
  
  for ( it in (bi+1):N )
  {
    llk[it-bi,1,1] = llk[it-bi,1,1] - sum((x1-fit$delta1[,1,it])*(fit$iS1[,,it]%*%(x1-fit$delta1[,1,it])))/2
    llk[it-bi,1,2] = llk[it-bi,1,2] - sum((x1-fit$delta1[,2,it])*(fit$iS1[,,it]%*%(x1-fit$delta1[,2,it])))/2
    llk[it-bi,2,1] = llk[it-bi,2,1] - sum((x1-fit$delta1[,2,it])*(fit$iS1[,,it]%*%(x1-fit$delta1[,2,it])))/2
    llk[it-bi,2,2] = llk[it-bi,2,2] - sum((x1-fit$delta1[,3,it])*(fit$iS1[,,it]%*%(x1-fit$delta1[,3,it])))/2
    
    for ( c in 0:1 )
    {
      llk[it-bi,c+1,] = llk[it-bi,c+1,] - sum((x2[,1]-fit$delta2[,c+1,it])*(fit$iS2[,,it]%*%(x2[,1]-fit$delta2[,c+1,it])))/2
      llk[it-bi,,c+1] = llk[it-bi,,c+1] - sum((x2[,2]-fit$delta2[,c+1,it])*(fit$iS2[,,it]%*%(x2[,2]-fit$delta2[,c+1,it])))/2
    }
    
    
    
    for ( s in 1:2 )
    {
      m = rep(0,3*p+2*q[s])
      m[1:p] = Ay[,1,s]/fit$sigma2[it]
      m[p+1:q[s]] = fit$gamma[1:q[s],,s,it]%*%x3[,1]
      m[p+q[s]+1:p] = Ay[,2,s]/fit$sigma2[it]
      m[2*p+q[s]+1:q[s]] = fit$gamma[1:q[s],,s,it]%*%x3[,2]
      
      for ( c1 in 0:1 )
        for ( c2 in 0:1 )
        {
          V = matrix(0,3*p+2*q[s],3*p+2*q[s])
          V[1:p,1:p] = fit$AA[,,s]/fit$sigma2[it] + diag(1/fit$sig2_zeta[it],p)
          V[1:p,p+1:q[s]] = -fit$lam[,1:q[s],c1+1,s,it] / fit$sig2_zeta[it]
          V[p+1:q[s],1:p] = t(V[1:p,p+1:q[s]])
          V[p+1:q[s],p+1:q[s]] = -t(fit$lam[,1:q[s],c1+1,s,it]) %*% V[1:p,p+1:q[s]] + diag(1,q[s])
          
          V[p+q[s]+1:p,p+q[s]+1:p] = fit$AA[,,s]/fit$sigma2[it] + diag(1/fit$sig2_zeta[it],p)
          V[p+q[s]+1:p,2*p+q[s]+1:q[s]] = -fit$lam[,1:q[s],c2+1,s,it] / fit$sig2_zeta[it]
          V[2*p+q[s]+1:q[s],p+q[s]+1:p] = t(V[p+q[s]+1:p,2*p+q[s]+1:q[s]])
          V[2*p+q[s]+1:q[s],2*p+q[s]+1:q[s]] = -t(fit$lam[,1:q[s],c2+1,s,it]) %*% V[p+q[s]+1:p,2*p+q[s]+1:q[s]] + diag(1,q[s])
          
          V[2*p+2*q[s]+1:p,1:p] = -diag(1/fit$sig2_zeta[it],p)
          V[2*p+2*q[s]+1:p,p+1:q[s]] = -V[1:p,p+1:q[s]]
          V[2*p+2*q[s]+1:p,p+q[s]+1:p] = V[2*p+2*q[s]+1:p,1:p]
          V[2*p+2*q[s]+1:p,2*p+q[s]+1:q[s]] = -V[p+q[s]+1:p,2*p+q[s]+1:q[s]]
          V[2*p+2*q[s]+1:p,2*p+2*q[s]+1:p] = diag(2/fit$sig2_zeta[it]+1/fit$sig2_iota[it],p)
          V[1:(2*p+2*q[s]),2*p+2*q[s]+1:p] = t(V[2*p+2*q[s]+1:p,1:(2*p+2*q[s])])
          
          cV = chol(V)
          llk[it-bi,c1+1,c2+1] = llk[it-bi,c1+1,c2+1] - sum(log(diag(cV))) + sum(forwardsolve(t(cV),m)^2)/2
        }
    }
  }
  
  llk = llk - apply(llk,1,max)
  exp(llk)/apply(exp(llk),1,sum)
}

