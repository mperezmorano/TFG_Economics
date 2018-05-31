function s=candidate(x,subset)
p=size(x,2);
if subset==1;
    s=tril(ones(p+1,p+1),-1);
    s=s(:,1:end-1);
    
    %%by Liu
    s(1,:)=[];
    s(:,1)=1;
    %%
    
elseif subset==2;
    p=p-1;
    s=zeros(2^p,p);
    s0=[1,zeros(1,p-1)];
    s1=zeros(1,p);
    for i=2:2^p;
        s1=s0+s1;
        for j=1:p;
            if s1(j)==2;
                s1(j+1)=s1(j+1)+1;
                s1(j)=0;
            end
        end
        s(i,:)=s1;
    end
    
    %%by Liu
    s=[ones(size(s,1),1) s];
    %%
end
end