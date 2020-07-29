function [din,dout]=ba_model(n, m, a, gamma, rho)
rng('shuffle');
T=(rand(m,m)<rho).*(ones(m,m)-eye(m));
targets=[1:m];
dout=sum(T');
din=sum(T);
source=m+1;
D = cumsum(gamma);
f = waitbar(0,'Please wait...');
while source<=n
    cc=0;
    for i=targets
        dout(i)=dout(i)+1;
        cc=cc+1;
    end
    dout=[dout,0];
    din=[din, cc];
    tr=dout+a; str=sum(tr); tr=tr./str;
    bins=[0 cumsum(tr)];
    mt=source;
    while mt>source-1
        r=rand();
        mt=find(r<D,1,'first')-1;
    end
    j=1; targets=zeros(1,mt); alt=0;
    while j<=mt
       r=rand();
       cb=discretize(r,bins);
       if max(targets==cb)==0
           targets(j)=cb;
           j=j+1;
       end
       alt=alt+1;
       if alt>source && source>10000
           figure(4); plot([0 1],[0 0],'--r');
       end
    end
    source=source+1;
    if mod(source,500)==0
        waitbar(source/n,f,'Please wait...');
    end
end
close(f)