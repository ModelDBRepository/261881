function [sim_indeg,sim_outdeg,gamma_true]=conv_model(data_indeg,data_outdeg,Ek,l,rho,m0,pup,pdw,pc)

N=length(data_indeg);
binvi=[0:5:max(data_indeg)];
binvo=[0:5:max(data_outdeg)];

figure(10); 
subplot(2,2,1); hold on;
title(strcat('in-degree for N=',num2str(N)))
hdi = histogram(data_indeg,binvi,'Normalization','pdf');
xlabel('in-degree'); ylabel('probability');

subplot(2,2,2); hold on;
title(strcat('out-degree for N=',num2str(N)))
hdo = histogram(data_outdeg,binvo,'Normalization','pdf');
xlabel('out-degree'); ylabel('probability');


h=hdi;

ed=h.BinEdges;
data_vedi=[0]; data_vhi=[0]; c=2;
for i=1:h.NumBins
    if h.Values(i)>0
        data_vedi(c)=(ed(i)+ed(i+1))*0.5;
        data_vhi(c)=h.Values(i);
        c=c+1;
    end
end

h=hdo;

ed=h.BinEdges;
data_vedo=[0]; data_vho=[0]; c=2;
for i=1:h.NumBins
    if h.Values(i)>0
        data_vedo(c)=(ed(i)+ed(i+1))*0.5;
        data_vho(c)=h.Values(i);
        c=c+1;
    end
end

n=[1:floor(data_vedo(end))]; ni=[1:floor(data_vedi(end))];
binsi=[-0.5:1:max(data_indeg)];
binso=[-0.5:1:max(data_outdeg)];

m=mean(data_indeg);
N1=N/2; N2=N/2; M=floor(N1/l)
cutoff=0.001;
nf=floor(data_vedi(end));
hi2=interp1(data_vedi,data_vhi,[0:nf]);
T=hi2.*[0:nf];

Es=sum(T); E=cumsum(T);
s=sum(hi2); S=cumsum(hi2);
V=((Es-E)-[1:nf+1].*(s-S))./(s-S);
figure(10); hold on;
plot(V); plot([1 nf+1],[m-Ek m-Ek]);
del=find(V<(m-Ek),1,'first')
c=(m-Ek-m0/N1*(m0-1)*rho)*N1/(N1-m0)
a=pc*c;

f=beta([0 n]+c,3)./(beta(c,2));
ho1=interp1(data_vedo,data_vho,n); go1=impl_inf_deconv(f,[0 ho1]);
ho2=interp1(data_vedo,data_vho,n,'spline'); go2=impl_inf_deconv(f,[0 ho2]);
p=(Ek/N1-pdw)/(pup-pdw)

hi1=zeros(1,floor(data_vedi(end)));
hi1(1:data_vedi(end))=interp1(data_vedi,data_vhi,ni);



alpha_true=zeros(1,ni(end));
alpha_true([1:ni(end)-del+1])=hi1([del:ni(end)]);
alpha_true=alpha_true/sum(alpha_true);
sa=sum(alpha_true.*[0:ni(end)-1]);
la=length(alpha_true);

iota=binopdf([0:la-1],m0-1,rho);
gamma_true=N1/(N1-m0)*alpha_true-m0/(N1-m0)*iota;
gamma_true=gamma_true/sum(gamma_true);
ga=sum(gamma_true.*[0:la-1])


rng('shuffle');
nsim=10;


[F1_in,F1_out]=ba_model(N1, m0, a, gamma_true, rho);
[F2_in,F2_out]=ba_model(N2, m0, a, gamma_true, rho);
upper_graph_in=(rand(M,M)<p);
upper_graph_out=(rand(M,M)<p);
C1_in=zeros(1,N1); C2_in=zeros(1,N2);
C1_out=zeros(1,N1); C2_out=zeros(1,N2);
part1=randperm(N1); part2=randperm(N2);
for li=1:M
    for lj=1:M
        if upper_graph_in(li,lj)==1
            Ma=(rand(l,l)<pup);
            C1_in(part2((lj-1)*l+1:lj*l))=C1_in(part2((lj-1)*l+1:lj*l))+sum(Ma);
            C1_out(part1((li-1)*l+1:li*l))=C1_out(part1((li-1)*l+1:li*l))+sum(Ma');
        else
            Ma=(rand(l,l)<pdw);
            C1_in(part2((lj-1)*l+1:lj*l))=C1_in(part2((lj-1)*l+1:lj*l))+sum(Ma);
            C1_out(part1((li-1)*l+1:li*l))=C1_out(part1((li-1)*l+1:li*l))+sum(Ma');
        end
    end
end

for li=1:M
    for lj=1:M
        if upper_graph_out(li,lj)==1
            Ma=(rand(l,l)<pup);
            C2_in(part1((li-1)*l+1:li*l))=C2_in(part1((li-1)*l+1:li*l))+sum(Ma);
            C2_out(part2((lj-1)*l+1:lj*l))=C2_out(part2((lj-1)*l+1:lj*l))+sum(Ma');
        else
            Ma=(rand(l,l)<pdw);
            C2_in(part1((li-1)*l+1:li*l))=C2_in(part1((li-1)*l+1:li*l))+sum(Ma);
            C2_out(part2((lj-1)*l+1:lj*l))=C2_out(part2((lj-1)*l+1:lj*l))+sum(Ma');
        end
    end
end

sim_indeg=[F1_in+C2_in,F2_in+C1_in];
sim_outdeg=[F1_out+C1_out,F2_out+C2_out];

