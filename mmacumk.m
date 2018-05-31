function [betahat,w,yhat,ehat,r2,cn,wstar] = mmacumk(y,x,k)

%  This Matlab function computes the Mallows Model Average (MMA) 
%  least-squares estimates, as described in the paper "The Asymptotic 
%  Risk of the Least Squares Averaging estimator".
%
%  written by
%  Chu-An Liu and Bruce E. Hansen
%  Department of Economics
%  University of Wisconsin
%
%  Format: 
%
%  [betahat,w,yhat,ehat,r2,cn,wstar]=mmacumk(y,x,k)
%  
%  Inputs:
%  y           nx1      dependent variable
%  x           nxp      regressor matrix
%  k           mx1      number of regressors for submodels 
%
%  Outputs:
%  betahat     px1      parameter estimate
%  w           mx1      weight vector
%  yhat        nx1      fitted values   
%  ehat        nx1      fitted residuals   
%  r2          1x1      R-squared
%  cn          1x1      value of Mallows criterion
%  wstar       (m-1)x1  cummulative weight vector
%
%  Note:
%  The regressors columns should be ordered, with the intercept first 
%  and then in order of relevance. 
% 

[n,p]=size(x);
m=length(k);
xx=x'*x;
sxy=x'*y;
bbeta=zeros(p,m);
kc=cumsum(k);
for j=1:m;
    bbeta(1:kc(j),j)=(xx(1:kc(j),1:kc(j))\sxy(1:kc(j)))';
end
ee=y*ones(1,m)-x*bbeta;
ehat=y-x*bbeta(:,m);
sighat=(ehat'*ehat)/(n-p);
ff=zeros(m-1,1);
for j=1:m-1;
    ff(j)=(ee(:,j)'*ee(:,j)-ee(:,j+1)'*ee(:,j+1))/sighat;
end
kk=k(2:m,1);
wbar=kk./ff;
indx=wbar(1:end-1)>=wbar(2:end);
indx1=[0;indx];
while sum(indx)>0;
    kbar=kk;
    fbar=ff;
    for j=1:m-1;
        if indx1(j)==1;
            kbar(j)=kbar(j-1)+kbar(j);
            fbar(j)=fbar(j-1)+fbar(j);
            kbar(j-1)=0;
            fbar(j-1)=0;
        end
    end
    wbar=kbar./fbar;   
    wbar=wbar(isnan(wbar)==0);
    indx=wbar(1:end-1)>=wbar(2:end);
    indx2=[0;indx];
    indx1(indx1==0)=indx2;
end
wbar(wbar>1)=1;
wstar=zeros(m-1,1);
wstar(indx1==0)=wbar;
w=(indx1==0)+0;
ww=[wbar(1);wbar(2:end)-wbar(1:end-1)];
w(w>0)=ww;
w=[w;1-wbar(end)];
betahat=bbeta*w;
ybar=mean(y);
yhat=x*betahat;
ehat=y-yhat;
r2=sum((yhat-ybar).^2)/sum((y-ybar).^2);
a1=ee'*ee;
a2=kc*sighat;
cn=(w'*a1*w+2*a2'*w)/n;
end

