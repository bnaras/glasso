"
c
c           Lasso regularized covariance matrix estimate
c
c                         version (9/13/11)
c
c call glasso(n,ss,rho,ia,is,itr,ipen,thr,maxit,ww,wwi,niter,del,jerr)
c
c Input:
c    n = dimension of matrix
c    ss(n,n) = data covariance matrix
c    rho(n,n) = regularization strength parameters for each element
c              (must be symmetric: rho(i,j)=rho(j,i))
c    ia = approximation flag
c       ai =  0 => exact solution
c       ia != 0 => Meinhausen-Buhlmann approximation
c    is = initialization flag
c       is  = 0 => cold start: initialize using ss
c       is != 0 => warm start: initialize with previous solution
c                  stored in ww and wwi (see below)
c    itr = trace flag
c       itr != 0 => trace information printed
c       itr =  0 => trace information not printed
c    ipen = diagonal penalty flag
c       ipen != 0 => diagonal is penalized
c       ipen =  0 => diagonal is not penalized
c    thr = convergence threshold: iterations stop when average absolute
c          parameter change is less than thr * ave(abs(offdiag(ss)))
c          (suggested value 1.0e-4)
c    maxit = maximum number of iterations (no effect for ia ! = 0)
c
c Output:
c    ww(n,n) = solution covariance matrix estimate (ia = 0)
c               (not used for ia != 0)
c    wwi(n,n) = solution inverse covariance matrix estimate (ia = 0)
c             = off-diagonal lasso coefficients (ia != 0)
c    niter = number of iterations
c    del = average absolute parameter change at termination
c             (not used for ia != 0)
c    jerr = memory allocation error flag
c      jerr = 0 => no error
c      jerr != 0 => memory allocation error - no output returned
c
c
"
subroutine
   glasso(nn,sss,rrho,ia,is,itr,ipen,thr,maxit,www,wwwi,nniter,ddel,jerr);
implicit double precision(a-h,o-z);
double precision sss(nn,nn),rrho(nn,nn),www(nn,nn),wwwi(nn,nn);
%fortran
      double precision, dimension (:), allocatable :: ss,rho,ww,wwi
      integer, dimension (:), allocatable :: ir,ie
      integer, dimension (:,:), allocatable :: ic
%mortran
if ia.ne.0 <
  call lasinv1(nn,sss,rrho,ia,is,itr,ipen,thr,maxit,www,wwwi,nniter,ddel,jerr);
  return;
>
%fortran
      allocate(ic(1:2,1:nn),stat=jerr)
%mortran
if (jerr.ne.0) return;
allocate(ir(1:nn),stat=jerr);
if (jerr.ne.0) return;
allocate(ie(1:nn),stat=jerr);
if(jerr.ne.0) return;
call connect(nn,sss,rrho,nc,ic,ir,ie);
nnq=0; <kc=1,nc; nnq=max(ic(2,kc)-ic(1,kc)+1,nnq);> nnq=nnq**2;
allocate(ss(1:nnq),stat=jerr);
if(jerr.ne.0) return;
allocate(rho(1:nnq),stat=jerr);
if(jerr.ne.0) return;
allocate(ww(1:nnq),stat=jerr);
if(jerr.ne.0) return;
allocate(wwi(1:nnq),stat=jerr);
if(jerr.ne.0) return;
nniter=0; ddel=0.0; l=0;
<kc=1,nc; n=ic(2,kc)-ic(1,kc)+1;
   if n.le.1 < k=ir(ic(1,kc));
      www(:,k)=0.0; www(k,:)=0.0; wwwi(:,k)=0.0; wwwi(k,:)=0.0; next;
   >
   kb=ic(1,kc); ke=ic(2,kc); l=0;
   <k=kb,ke; ik=ir(k); <j=kb,ke; ij=ir(j); l=l+1;
      ss(l)=sss(ij,ik); rho(l)=rrho(ij,ik);
      ww(l)=www(ij,ik); wwi(l)=wwwi(ij,ik);
   >>
   call lasinv1(n,ss,rho,ia,is,itr,ipen,thr,maxit,ww,wwi,niter,del,jerr);
   if(jerr.ne.0) return; nniter=nniter+niter; ddel=ddel+del;
   <j=kb,ke; k=ir(j);
      www(:,k)=0.0; www(k,:)=0.0; wwwi(:,k)=0.0; wwwi(k,:)=0.0;
   >
   l=0; <k=kb,ke; ik=ir(k); <j=kb,ke; l=l+1; wwwi(ir(j),ik)=wwi(l);>>
   if ia.eq.0 < l=0;
      <k=kb,ke; ik=ir(k); <j=kb,ke; l=l+1; www(ir(j),ik)=ww(l);>>
   >
>
ddel=ddel/nc;
if ia.eq.0 <   
   <j=1,nn; if(www(j,j).ne.0.0) next;
      if ipen.eq.0 < www(j,j)=sss(j,j);> else < www(j,j)=sss(j,j)+rrho(j,j);>
      wwwi(j,j)=1.0/www(j,j);
   >
>
return;
end;
subroutine connect(n,ss,rho,nc,ic,ir,ie);
implicit double precision(a-h,o-z);
double precision ss(n,n),rho(n,n); integer ic(2,n),ir(n),ie(n);
ie=0; nc=0; is=1;
<k=1,n; if(ie(k).gt.0) next;
   ir(is)=k; nc=nc+1; ie(k)=nc; ic(1,nc)=is; is=is+1;
   call row(nc,1,ir((is-1):n),n,ss,rho,ie,na,ir(is:n));
   if na.eq.0 < ic(2,nc)=is-1; next;>
   loop < nas=na; iss=is; il=iss+nas-1;
      if(il.ge.n) exit; is=is+na;
      call row(nc,nas,ir(iss:n),n,ss,rho,ie,na,ir(is:n));
   > until na.eq.0;
   ic(2,nc)=il;
>
return;
end;
subroutine row(nc,nr,jr,n,ss,rho,ie,na,kr);
implicit double precision(a-h,o-z);
double precision ss(n,n),rho(n,n); integer jr(nr),ie(n),kr(*"na");
na=0;
<l=1,nr; k=jr(l);
   <j=1,n; if(ie(j).gt.0) next; if(j.eq.k) next;
      if(abs(ss(j,k)).le.rho(j,k)) next;
      na=na+1; kr(na)=j; ie(j)=nc;
   >
>
return;
end;  
subroutine lasinv1(n,ss,rho,ia,is,itr,ipen,thr,maxit,ww,wwi,niter,del,jerr);
implicit double precision(a-h,o-z);
parameter(eps=1.0e-7);
double precision ss(n,n),rho(n,n),ww(n,n),wwi(n,n);
%fortran
      double precision, dimension (:,:), allocatable :: vv,xs
      double precision, dimension (:), allocatable :: s,x,z,ws,ro,so
      integer, dimension (:), allocatable :: mm
      nm1=n-1
      allocate(vv(1:nm1,1:nm1),stat=jerr)
      if(jerr.ne.0) return
      if(ia.eq.0) allocate(xs(1:nm1,1:n),stat=jerr)
      if(jerr.ne.0) return
%mortran
allocate(s(1:nm1),stat=jerr);
if(jerr.ne.0) return;
allocate(so(1:nm1),stat=jerr);
if(jerr.ne.0) return;
allocate(x(1:nm1),stat=jerr);
if(jerr.ne.0) return;
allocate(z(1:nm1),stat=jerr);
if(jerr.ne.0) return;
allocate(mm(1:nm1),stat=jerr);
if(jerr.ne.0) return;
allocate(ro(1:nm1),stat=jerr);
if ia.eq.0 < allocate(ws(1:n),stat=jerr); >
if(jerr.ne.0) return;
shr=0.0; <j=1,n; <k=1,n; if(j.eq.k) next; shr=shr+abs(ss(j,k));>>
if shr.eq.0.0 < ww=0.0; wwi=0.0;
   <j=1,n; if ipen.eq.0 < ww(j,j)=ss(j,j);> else < ww(j,j)=ss(j,j)+rho(j,j);>
      wwi(j,j)=1.0/max(ww(j,j),eps);
   >
   return;
>
shr=thr*shr/nm1;
if ia.ne.0 < if(is.eq.0) wwi=0.0;
   <m=1,n; call setup(m,n,ss,rho,ss,vv,s,ro); l=0;
      <j=1,n; if(j.eq.m) next; l=l+1; x(l)=wwi(j,m);>
      call lasso(ro,nm1,vv,s,shr/n,x,z,mm);
      l=0; <j=1,n; if(j.eq.m) next; l=l+1; wwi(j,m)=x(l);>
   >
   niter=1; return;
>     
if is.eq.0 < ww=ss; xs=0.0;>
else <
   <j=1,n; xjj=-wwi(j,j); l=0;
      <k=1,n; if(k.eq.j) next; l=l+1; xs(l,j)=wwi(k,j)/xjj;>
   >
>
<j=1,n; if ipen.eq.0 < ww(j,j)=ss(j,j);> else < ww(j,j)=ss(j,j)+rho(j,j);>>
niter=0;
loop < dlx=0.0;
   <m=1,n;
      if itr.ne.0 < call intpr('m',1,m,1); >
      x=xs(:,m); ws=ww(:,m);
      call setup(m,n,ss,rho,ww,vv,s,ro); so=s;
      call lasso(ro,nm1,vv,s,shr/sum(abs(vv)),x,z,mm);
      l=0;
      <j=1,n; if(j.eq.m) next;
         l=l+1; /ww(j,m),ww(m,j)/=so(l)-s(l);
      >
      dlx=max(dlx,sum(abs(ww(:,m)-ws)));
      xs(:,m)=x;
   >
   niter=niter+1; if(niter.ge.maxit) exit;
> until dlx.lt.shr;
del=dlx/nm1; call inv(n,ww,xs,wwi);
return;
end;
subroutine setup(m,n,ss,rho,ww,vv,s,r);
implicit double precision(a-h,o-z);
double precision ss(n,n),rho(n,n),ww(n,n),vv(n-1,n-1),s(n-1),r(n-1);
l=0;
<j=1,n; if(j.eq.m) next; l=l+1; r(l)=rho(j,m); s(l)=ss(j,m); i=0;
   <k=1,n; if(k.eq.m) next; i=i+1; vv(i,l)=ww(k,j);>
>
return;
end;  
subroutine lasso(rho,n,vv,s,thr,x,z,mm);
implicit double precision(a-h,o-z);
double precision rho(n),vv(n,n),s(n),x(n),z(n); integer mm(n);
call fatmul(2,n,vv,x,s,z,mm);
loop < dlx=0.0;
   <j=1,n; xj=x(j); x(j)=0.0; t=s(j)+vv(j,j)*xj;
      if (abs(t)-rho(j).gt.0.0) x(j)=sign(abs(t)-rho(j),t)/vv(j,j);
      if(x(j).eq.xj) next; del=x(j)-xj;
      dlx=max(dlx,abs(del)); s=s-del*vv(:,j);
   >
> until dlx.lt.thr;
return;
end;
subroutine fatmul(it,n,vv,x,s,z,m);
implicit double precision(a-h,o-z);
parameter(fac=0.2);
double precision vv(n,n),x(n),s(n),z(n); integer m(n);
l=0;
<j=1,n; if(x(j).eq.0.0) next; l=l+1; m(l)=j; z(l)=x(j);>
if l.gt.int(fac*n) <
   if it.eq.1 < s=matmul(vv,x);> else < s=s-matmul(x,vv);>
>
elseif it.eq.1 < <j=1,n; s(j)=dot_product(vv(j,m(1:l)),z(1:l));>>
else < <j=1,n; s(j)=s(j)-dot_product(vv(m(1:l),j),z(1:l));>>
return;
end;
subroutine inv(n,ww,xs,wwi);
implicit double precision(a-h,o-z);
double precision ww(n,n),xs(n-1,n),wwi(n,n);
nm1=n-1; xs=-xs;
wwi(1,1)=1.0/(ww(1,1)+dot_product(xs(:,1),ww(2:n,1)));
wwi(2:n,1)=wwi(1,1)*xs(:,1);
wwi(n,n)=1.0/(ww(n,n)+dot_product(xs(:,n),ww(1:nm1,n)));
wwi(1:nm1,n)=wwi(n,n)*xs(:,n);
<j=2,nm1; jm1=j-1; jp1=j+1;
   wwi(j,j)=1.0/(ww(j,j)+dot_product(xs(1:jm1,j),ww(1:jm1,j))
      +dot_product(xs(j:nm1,j),ww(jp1:n,j)));
   wwi(1:jm1,j)=wwi(j,j)*xs(1:jm1,j);
   wwi(jp1:n,j)=wwi(j,j)*xs(j:nm1,j);
>
return;
end;

subroutine glassopath(beta,what,jerrs,rholist,nrho,n,ss,rho,ia,itr,
ipen,thr,maxit,ww,wwi,niter,del,jerr);
implicit double precision(a-h,o-z);
integer nrho,n,jerrs(nrho);
double precision rholist(nrho),beta(n,n,nrho),what(n,n,nrho);
double precision ss(n,n),rho(n,n),ww(n,n),wwi(n,n);

is=0;
<
  j=1,n;
  <
    k=1,n;
    rho(j,k) = rholist(nrho);
  >
>
call glasso(n,ss,rho,ia,is,itr,ipen,thr,maxit,ww,wwi,niter,del,jerr);
jerrs(1)=jerr;
<
  j=1,n;
  <
    k=1,n;
    beta(j,k,nrho)=wwi(j,k);
    what(j,k,nrho)=ww(j,k);
  >
>
is=1;
<
  i=nrho,1,-1;
  <
    j=1,n;
    <
      k=1,n;
      rho(j,k)=rholist(i);
    >
  >
  if(itr.gt.0) call dblepr('rho=', -1, rholist(i), 1);
  itr2=itr;
  if(itr2.gt.0) itr2=itr-1;
  call glasso(n,ss,rho,ia,is,itr2,ipen,thr,maxit,ww,wwi,niter,del,jerr);
  jerrs(i)=jerr;
  <
    j=1,n;
    <
      k=1,n;
      beta(j,k,i)=wwi(j,k);
      what(j,k,i)=ww(j,k);
    >
  >
>
return;
end;
%%
   
