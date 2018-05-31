%Written by Qingfeng Liu based on Chu-An Liu and Bruce E. Hansen
%For \bar{C1}, with original y and x.
function [mu_hat_glsma,bhat_glsma,wglsma]=ma_C1_bar(y,x,s,sigma2)

[n,p]=size(x);
m=length(s(:,1));
bbeta=zeros(p,m);
a2glsma=zeros(m,1);
for j=1:m;
    ss=ones(n,1)*s(j,:)>0;
    xs=x(ss);
    xs=reshape(xs,n,length(xs)/n);
    xstar=diag(diag(sigma2).^(-1/2))*xs;
    ystar=diag(diag(sigma2).^(-1/2))*y;
    betas=(xstar'*xstar)\xstar'*ystar;
    a2glsma(j,1)=length(y)/(length(y)-p)*trace(xs/(xstar'*xstar)*xs');
    %a2glsma(j,1)=trace(xs/(xstar'*xstar)*xs');
    sj=s(j,:)>0;
    bbeta(sj,j)=betas;
end
ee=y*ones(1,m)-x*bbeta;
a1=ee'*ee;

w0=ones(m,1)/m;
options=optimset('LargeScale','off','Display','off', 'Algorithm','active-set');
w=quadprog(a1,a2glsma,zeros(1,m),0,ones(1,m),1,zeros(m,1),ones(m,1),w0,options);
w=w.*(w>0);
wglsma=w/sum(w);
bhat_glsma=bbeta*wglsma;
mu_hat_glsma=x*bhat_glsma;






