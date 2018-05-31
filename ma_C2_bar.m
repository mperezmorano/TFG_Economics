%For \bar{C2}, with weighted y and x.
function [muhat_glsma,bhat_glsma,wglsma]=ma_C2_bar(y,x,s)
[n,p]=size(x);
m=length(s(:,1));
bbeta=zeros(p,m);
for j=1:m;
    ss=ones(n,1)*s(j,:)>0;
    xs=x(ss);
    xs=reshape(xs,n,length(xs)/n);
    betas=(xs'*xs)\(xs'*y);
    sj=s(j,:)>0;
    bbeta(sj,j)=betas;
end
ee=y*ones(1,m)-x*bbeta;
a1=ee'*ee;
a2glsma=sum(s,2);

w0=ones(m,1)/m;
options=optimset('LargeScale','off','Display','off', 'Algorithm','active-set');
w=quadprog(a1,a2glsma,zeros(1,m),0,ones(1,m),1,zeros(m,1),ones(m,1),w0,options);
w=w.*(w>0);
wglsma=w/sum(w);
bhat_glsma=bbeta*wglsma;
muhat_glsma=x*bhat_glsma;








