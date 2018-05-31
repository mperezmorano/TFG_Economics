%MMA, JMA and GC
%Written by Qingfeng Liu 
%Based on the program of JMA by Chu-An Liu and Bruce E. Hansen
% y           nx1   dependent variable
% x           nxp   regressor matrix
% s           mxp   selection matrix. optional. m: number of models
%                    Example:
%                       Suppose there are 3 candidate models.
%                       Model 1: y=beta1*x1+beta2*x2+e
%                       Model 2: y=beta1*x1+beta3*x3+e
%                       Model 3: y=beta1*x1+beta2*x2+beta4*x4+e
%                       Then s=[1,1,0,0;
%                               1,0,1,0;
%                               1,1,0,1]
function [bhat_hrcp,bhat_jma,bhat_mma,muhat_hrcp,muhat_jma,muhat_mma,whrcp,wjma,wmma]=ma(y,x,s)

%%The following statement is the original note by Chu-An Liu and Bruce E. Hansen

%This Matlab function computes the Mallows Model Average (MMA) and
%  the Jackknife Model Average (JMA) least-squares estimates.
%
%  written by
%  Chu-An Liu and Bruce E. Hansen
%  Department of Economics
%  University of Wisconsin
%
%  Format:
%
%  [betahat,w,muhat,ehat,R2,Cn]=gma(y,x,method,subset,s)
%  or  [betahat,w,muhat,ehat,R2,Cn]=gma(y,x,method,subset)
%
%  Inputs:
%  y           nx1   dependent variable
%  x           nxp   regressor matrix
%  method      1x1   set to 1 for Mallows model average estimates
%                    set to 2 for Jackknife model average estimates
%                    set to 3 for Generalized Cp model average estimates
%  subset      1x1   set to 1 for pure nested subsets
%                    set to 2 for all combinations of subsets
%                    set to 3 for using the selection matrix
%  s           mxp   selection matrix. optional. m: number of models
%                    Example:
%                       Suppose there are 3 candidate models.
%                       Model 1: y=beta1*x1+beta2*x2+e
%                       Model 2: y=beta1*x1+beta3*x3+e
%                       Model 3: y=beta1*x1+beta2*x2+beta4*x4+e
%                       Then s=[1,1,0,0;
%                               1,0,1,0;
%                               1,1,0,1]
%
%  Outputs:
%  betahat     px1   parameter estimate
%  w           mx1   weight vector
%  muhat        nx1   fitted values
%  ehat        nx1   fitted residuals
%  r2          1x1   R-squared
%  cn          1x1   Value of Mallows criterion or Cross-Validation criterion
%
%  Note:
%  For pure nested subsets, the regressors columns should be ordered, with the
%  intercept first and then in order of relevance.
%  For all combinations of subsets, p is less than about 20.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
[n,p]=size(x);
m=length(s(:,1));
bbeta=zeros(p,m);
eejma=zeros(n,m);
for j=1:m;
    ss=ones(n,1)*s(j,:)>0;
    xs=x(ss);
    xs=reshape(xs,n,length(xs)/n);
    betas=(xs'*xs)\(xs'*y);
    sj=s(j,:)>0;
    bbeta(sj,j)=betas;
    ei=y-xs*betas;
    hi=diag(xs/(xs'*xs)*xs');
    eejma(:,j)=ei.*(1./(1-hi));
end
ee=y*ones(1,m)-x*bbeta;
ehat=y-x*bbeta(:,m);
sighat=(ehat'*ehat)/(n-p);
a1=ee'*ee;
a1jma=eejma'*eejma;
a2mma=sum(s,2)*sighat;
a2jma=zeros(m,1);

%% from here, HRCp, by LIU
Omega_hat=diag(ee(:,end).^2);
a2hrcp=zeros(m,1);
for i=1:m;
    ss=ones(n,1)*s(i,:)>0;
    xs=x(ss);
    xs=reshape(xs,n,length(xs)/n);
    a2hrcp(i,1)=length(y)/(length(y)-p)*trace(Omega_hat*(xs/(xs'*xs)*xs'));
end;

w0=ones(m,1)/m;
options=optimset('LargeScale','off','Display','off', 'Algorithm','interior-point-convex');
w=quadprog(a1,a2hrcp,zeros(1,m),0,ones(1,m),1,zeros(m,1),ones(m,1),w0,options);
w=w.*(w>0);
whrcp=w/sum(w);
bhat_hrcp=bbeta*whrcp;
muhat_hrcp=x*bhat_hrcp;

% to here, by LIU

%%
w=quadprog(a1,a2mma,zeros(1,m),0,ones(1,m),1,zeros(m,1),ones(m,1),w0,options);
w=w.*(w>0);
wmma=w/sum(w);
bhat_mma=bbeta*wmma;
muhat_mma=x*bhat_mma;

w=quadprog(a1jma,a2jma,zeros(1,m),0,ones(1,m),1,zeros(m,1),ones(m,1),w0,options);
w=w.*(w>0);
wjma=w/sum(w);
bhat_jma=bbeta*wjma;
muhat_jma=x*bhat_jma;
