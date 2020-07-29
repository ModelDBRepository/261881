function plot_figure(dinEE,doutEE,dinII,doutII,sim_dinEE,sim_doutEE,sim_dinII,sim_doutII,dinEE4,doutEE4,dinII4,doutII4,dinEE3,doutEE3,dinII3,doutII3)

nE=length(dinEE); nI=length(dinII);
N=nE+nI; fsz=20; pprt=10; ntype='pdf';

figure(1); hold on; 
subplot(2,2,1); hold on;
title(strcat('NC in-degree EE'))
hdEEi = histogram(dinEE,pprt,'Normalization',ntype);
xlabel('in-degree'); ylabel('probability density');
legend({'data','sim'},'Location','northeast');

subplot(2,2,2); hold on;
title(strcat('NC out-degree EE'))
hdEEo = histogram(doutEE,pprt,'Normalization',ntype);
xlabel('out-degree'); ylabel('probability density');
legend({'data','sim'},'Location','northeast');

subplot(2,2,3); hold on;
title(strcat('NC in-degree II'))
hdIIi = histogram(dinII,pprt,'Normalization',ntype);
xlabel('in-degree'); ylabel('probability density');
legend({'data','sim'},'Location','northeast');

subplot(2,2,4); hold on;
title(strcat('NC out-degree II'))
hdIIo = histogram(doutII,pprt,'Normalization',ntype);
xlabel('out-degree'); ylabel('probability density');
legend({'data','sim'},'Location','northeast');

figure(1); hold on;
subplot(2,2,1); hold on;
hsEEi = histogram(sim_dinEE,pprt,'Normalization',ntype);
legend({'data','sim'},'Location','northeast');

subplot(2,2,2); hold on;
hsEEo = histogram(sim_doutEE,pprt,'Normalization',ntype);
legend({'data','sim'},'Location','northeast');

subplot(2,2,3); hold on;
hsIIi = histogram(sim_dinII,pprt,'Normalization',ntype);
legend({'data','sim'},'Location','northeast');

subplot(2,2,4); hold on;
hsIIo = histogram(sim_doutII,pprt,'Normalization',ntype);
legend({'data','sim'},'Location','northeast');

Nf=5180;
ks=[0:1:Nf-1]; tc=mean(dinII)
a=100*tc;
f=beta(ks+a,2+a/tc)./beta(a,1+a/tc);
sum(f)
ct=cumsum(f); data_vhi=sum(f)-ct(1:end-1);
lfei=log10(ks(1:end-1)); lfhi=log10(data_vhi);
lfeo=lfei; lfho=lfhi;

pbi=mean(dinII)/Nf;
pbo=mean(doutII)/Nf;

exi=binopdf(ks,Nf,pbi);
ct=cumsum(exi); data_vhi=sum(exi)-ct(1:end-1);
lxei=log10(ks(1:end-1)); lxhi=log10(data_vhi);

exo=binopdf(ks,Nf,pbo);
ct=cumsum(exo); data_vho=sum(exo)-ct(1:end-1);
lxeo=log10(ks(1:end-1)); lxho=log10(data_vho);

h=hdEEi;
clear data_vedi;
clear data_vhi;

ed=h.BinEdges;
step=ed(2)-ed(1);
data_vedi=[]; data_vhi=[]; c=1;
for i=1:h.NumBins
    if h.Values(i)>0
        data_vedi(c)=(ed(i)+ed(i+1))*0.5;
        data_vhi(c)=h.Values(i);
        c=c+1;
    end
end
W=[[0 data_vedi]', [0 data_vhi]'];
ct=cumsum(data_vhi); data_vhi=sum(data_vhi)-ct(1:end-1);
ldei=log10(data_vedi(1:end-1)); ldhi=log10(data_vhi*step);

figure(9); hold on;
subplot(1,2,1); hold on;
plot([0 ldei],[0 ldhi],'linewidth',3,'color',[0 0 1]);
set(gca,'linewidth',2); set(gca,'FontSize',fsz);

figure(12); hold on;
subplot(1,2,1); hold on;
plot([0 ldei],[0 ldhi],'linewidth',3,'color',[0 0 1]);
plot(lfei,lfhi,'linewidth',3,'color',[0 0 0]);
plot(lxei,lxhi,'linewidth',3,'color',[1 0 0]);
set(gca,'linewidth',2); set(gca,'FontSize',fsz);

h=hsEEi;
clear sim_vedi;
clear sim_vhi;

ed=h.BinEdges;
step=ed(2)-ed(1);
sim_vedi=[]; sim_vhi=[]; c=1;
for i=1:h.NumBins
    if h.Values(i)>0
        sim_vedi(c)=(ed(i)+ed(i+1))*0.5;
        sim_vhi(c)=h.Values(i);
        c=c+1;
    end
end
W=[sim_vedi', sim_vhi'];
ct=cumsum(sim_vhi); sim_vhi=sum(sim_vhi)-ct(1:end-1);
lsei=log10(sim_vedi(1:end-1)); lshi=log10(sim_vhi*step);


figure(9); hold on;
subplot(1,2,1); hold on;
plot([0,lsei],[0,lshi],'--','linewidth',3,'color',[0 0 1]);
set(gca,'linewidth',2); set(gca,'FontSize',fsz);

h=hdEEo; 
clear data_vedo;
clear data_vho;

ed=h.BinEdges;
step=ed(2)-ed(1);
data_vedi=[]; data_vhi=[]; c=1;
for i=1:h.NumBins
    if h.Values(i)>0
        data_vedo(c)=(ed(i)+ed(i+1))*0.5;
        data_vho(c)=h.Values(i);
        c=c+1;
    end
end
W=[[0 data_vedo]', [0 data_vho]'];
ct=cumsum(data_vho); data_vho=sum(data_vho)-ct(1:end-1);
ldeo=log10(data_vedo(1:end-1)); ldho=log10(data_vho*step);

figure(9); hold on;
subplot(1,2,2); hold on;
plot([0,ldeo],[0,ldho],'linewidth',3,'color',[0 0 1]);
xlabel('log_{10}(NC outdegree)'); ylabel('log_{10}(NC survival function)');
set(gca,'linewidth',2); set(gca,'FontSize',fsz);

figure(12); hold on;
subplot(1,2,2); hold on;
plot([0 ldeo],[0 ldho],'linewidth',3,'color',[0 0 1]);
xlabel('log_{10}(NC outdegree)'); ylabel('log_{10}(NC survival function)');
set(gca,'linewidth',2); set(gca,'FontSize',fsz);


h=hsEEo;
clear sim_vedo;
clear sim_vho;

ed=h.BinEdges;
step=ed(2)-ed(1);
sim_vedo=[]; sim_vho=[]; c=1;
for i=1:h.NumBins
    if h.Values(i)>0
        sim_vedo(c)=(ed(i)+ed(i+1))*0.5;
        sim_vho(c)=h.Values(i);
        c=c+1;
    end
end
W=[sim_vedo', sim_vho'];
ct=cumsum(sim_vho); sim_vho=sum(sim_vho)-ct(1:end-1);
lseo=log10(sim_vedo(1:end-1)); lsho=log10(sim_vho*step);

figure(9); hold on;
subplot(1,2,2); hold on;
plot([0,lseo],[0,lsho],'--','linewidth',3,'color',[0 0 1]);

inh_col=[1 0.5 0.25];

h=hdIIi;
clear data_vedi;
clear data_vhi;

ed=h.BinEdges;
step=ed(2)-ed(1);
data_vedi=[]; data_vhi=[]; c=1;
for i=1:h.NumBins
    if h.Values(i)>0
        data_vedi(c)=(ed(i)+ed(i+1))*0.5;
        data_vhi(c)=h.Values(i);
        c=c+1;
    end
end
W=[[0 data_vedi]', [0 data_vhi]'];
ct=cumsum(data_vhi); data_vhi=sum(data_vhi)-ct(1:end-1);
ldei=log10(data_vedi(1:end-1)); ldhi=log10(data_vhi*step);

figure(9); hold on;
subplot(1,2,1); hold on;
plot([0,ldei],[0,ldhi],'linewidth',3,'color',inh_col,'DisplayName','Inhibitory');
xlabel('log_{10}(NC indegree)'); ylabel('log_{10}(NC survival function)');

figure(12); hold on;
subplot(1,2,1); hold on;
plot([0 ldei],[0 ldhi],'linewidth',3,'color',inh_col,'DisplayName','Inhibitory');
xlabel('log_{10}(NC indegree)'); ylabel('log_{10}(NC survival function)');

h=hsIIi;
clear sim_vedi;
clear sim_vhi;

ed=h.BinEdges;
step=ed(2)-ed(1);
sim_vedi=[]; sim_vhi=[]; c=1;
for i=1:h.NumBins
    if h.Values(i)>0
        sim_vedi(c)=(ed(i)+ed(i+1))*0.5;
        sim_vhi(c)=h.Values(i);
        c=c+1;
    end
end
W=[sim_vedi', sim_vhi'];
ct=cumsum(sim_vhi); sim_vhi=sum(sim_vhi)-ct(1:end-1);
lsei=log10(sim_vedi(1:end-1)); lshi=log10(sim_vhi*step);

figure(9); hold on;
subplot(1,2,1); hold on;
plot([0,lsei],[0,lshi],'--','linewidth',3,'color',inh_col,'DisplayName','Inhibitory');

h=hdIIo; 
clear data_vedo;
clear data_vho;

ed=h.BinEdges;
step=ed(2)-ed(1);
data_vedi=[]; data_vhi=[]; c=1;
for i=1:h.NumBins
    if h.Values(i)>0
        data_vedo(c)=(ed(i)+ed(i+1))*0.5;
        data_vho(c)=h.Values(i);
        c=c+1;
    end
end
W=[[0 data_vedo]', [0 data_vho]'];
ct=cumsum(data_vho); data_vho=sum(data_vho)-ct(1:end-1);
ldeo=log10(data_vedo(1:end-1)); ldho=log10(data_vho*step);

figure(9); hold on;
subplot(1,2,2); hold on;
plot([0,ldeo],[0,ldho],'linewidth',3,'color',inh_col,'DisplayName','Inhibitory');

figure(12); hold on;
subplot(1,2,2); hold on;
plot([0 ldeo],[0 ldho],'linewidth',3,'color',inh_col,'DisplayName','Inhibitory');


h=hsIIo;
clear sim_vedo;
clear sim_vho;

ed=h.BinEdges;
step=ed(2)-ed(1);
sim_vedo=[]; sim_vho=[]; c=1;
for i=1:h.NumBins
    if h.Values(i)>0
        sim_vedo(c)=(ed(i)+ed(i+1))*0.5;
        sim_vho(c)=h.Values(i);
        c=c+1;
    end
end
W=[sim_vedo', sim_vho'];
ct=cumsum(sim_vho); sim_vho=sum(sim_vho)-ct(1:end-1);
lseo=log10(sim_vedo(1:end-1)); lsho=log10(sim_vho*step);

figure(9); hold on;
subplot(1,2,2); hold on;
plot([0,lseo],[0,lsho],'--','linewidth',3,'color',inh_col,'DisplayName','Inhibitory');

figure(10); hold on; 
subplot(2,2,1); hold on;
title(strcat('NC indegree EE from 6'))
hdEEi3 = histogram(dinEE3,pprt,'Normalization','pdf');
xlabel('indegree'); ylabel('probability density');
legend({'data','sim'},'Location','northeast');

subplot(2,2,2); hold on;
title(strcat('NC outdegree EE from 6'))
hdEEo3 = histogram(doutEE3,pprt,'Normalization','pdf');
xlabel('outdegree'); ylabel('probability density');
legend({'data','sim'},'Location','northeast');

subplot(2,2,3); hold on;
title(strcat('NC indegree II from 6'))
hdIIi3 = histogram(dinII3,pprt,'Normalization','pdf');
xlabel('indegree'); ylabel('probability density');
legend({'data','sim'},'Location','northeast');

subplot(2,2,4); hold on;
title(strcat('NC outdegree II from 6'))
hdIIo3 = histogram(doutII3,pprt,'Normalization','pdf');
xlabel('out-degree'); ylabel('probability density');
legend({'data','sim'},'Location','northeast');


figure(10); hold on; 
subplot(2,2,1); hold on;
title(strcat('NC indegree EE from 4'))
hdEEi4 = histogram(dinEE4,pprt,'Normalization','pdf');
xlabel('indegree'); ylabel('probability density');
legend({'data','sim'},'Location','northeast');

subplot(2,2,2); hold on;
title(strcat('NC outdegree EE from 4'))
hdEEo4 = histogram(doutEE4,pprt,'Normalization','pdf');
xlabel('outdegree'); ylabel('probability density');
legend({'data','sim'},'Location','northeast');

subplot(2,2,3); hold on;
title(strcat('NC indegree II from 4'))
hdIIi4 = histogram(dinII4,pprt,'Normalization','pdf');
xlabel('indegree'); ylabel('probability density');
legend({'data','sim'},'Location','northeast');

subplot(2,2,4); hold on;
title(strcat('NC outdegree II from 4'))
hdIIo4 = histogram(doutII4,pprt,'Normalization','pdf');
xlabel('out-degree'); ylabel('probability density');
legend({'data','sim'},'Location','northeast');

h=hdIIi3;
clear data_vedi;
clear data_vhi;

ed=h.BinEdges;
step=ed(2)-ed(1);
data_vedi=[]; data_vhi=[]; c=1;
for i=1:h.NumBins
    if h.Values(i)>0
        data_vedi(c)=(ed(i)+ed(i+1))*0.5;
        data_vhi(c)=h.Values(i);
        c=c+1;
    end
end
W=[[0 data_vedi]', [0 data_vhi]'];
ct=cumsum(data_vhi); data_vhi=sum(data_vhi)-ct(1:end-1);
ldei=log10(data_vedi(1:end-1)); ldhi=log10(data_vhi*step);

figure(12); hold on;
subplot(1,2,1); hold on;
plot([0 ldei],[0 ldhi],'linewidth',3,'color',inh_col,'DisplayName','Inhibitory');

h=hdIIi4;
clear data_vedi;
clear data_vhi;

ed=h.BinEdges;
step=ed(2)-ed(1);
data_vedi=[]; data_vhi=[]; c=1;
for i=1:h.NumBins
    if h.Values(i)>0
        data_vedi(c)=(ed(i)+ed(i+1))*0.5;
        data_vhi(c)=h.Values(i);
        c=c+1;
    end
end
W=[[0 data_vedi]', [0 data_vhi]'];
ct=cumsum(data_vhi); data_vhi=sum(data_vhi)-ct(1:end-1);
ldei=log10(data_vedi(1:end-1)); ldhi=log10(data_vhi*step);

figure(12); hold on;
subplot(1,2,1); hold on;
plot([0 ldei],[0 ldhi],'linewidth',3,'color',inh_col,'DisplayName','Inhibitory');

h=hdIIo3; 
clear data_vedo;
clear data_vho;

ed=h.BinEdges;
step=ed(2)-ed(1);
data_vedi=[]; data_vhi=[]; c=1;
for i=1:h.NumBins
    if h.Values(i)>0
        data_vedo(c)=(ed(i)+ed(i+1))*0.5;
        data_vho(c)=h.Values(i);
        c=c+1;
    end
end
W=[[0 data_vedo]', [0 data_vho]'];
ct=cumsum(data_vho); data_vho=sum(data_vho)-ct(1:end-1);
ldeo=log10(data_vedo(1:end-1)); ldho=log10(data_vho*step);

figure(12); hold on;
subplot(1,2,2); hold on;
plot([0 ldeo],[0 ldho],'linewidth',3,'color',inh_col,'DisplayName','Inhibitory');

h=hdIIo4; 
clear data_vedo;
clear data_vho;

ed=h.BinEdges;
step=ed(2)-ed(1);
data_vedi=[]; data_vhi=[]; c=1;
for i=1:h.NumBins
    if h.Values(i)>0
        data_vedo(c)=(ed(i)+ed(i+1))*0.5;
        data_vho(c)=h.Values(i);
        c=c+1;
    end
end
W=[[0 data_vedo]', [0 data_vho]'];
ct=cumsum(data_vho); data_vho=sum(data_vho)-ct(1:end-1);
ldeo=log10(data_vedo(1:end-1)); ldho=log10(data_vho*step);

figure(12); hold on;
subplot(1,2,2); hold on;
plot([0 ldeo],[0 ldho],'linewidth',3,'color',inh_col,'DisplayName','Inhibitory');

h=hdEEi4;
clear data_vedi;
clear data_vhi;

ed=h.BinEdges;
step=ed(2)-ed(1);
data_vedi=[]; data_vhi=[]; c=1;
for i=1:h.NumBins
    if h.Values(i)>0
        data_vedi(c)=(ed(i)+ed(i+1))*0.5;
        data_vhi(c)=h.Values(i);
        c=c+1;
    end
end
W=[[0 data_vedi]', [0 data_vhi]'];
ct=cumsum(data_vhi); data_vhi=sum(data_vhi)-ct(1:end-1);
ldei=log10(data_vedi(1:end-1)); ldhi=log10(data_vhi*step);

figure(12); hold on;
subplot(1,2,1); hold on;
plot([0 ldei],[0 ldhi],'linewidth',3,'color',[0 0 1]);
set(gca,'linewidth',2); set(gca,'FontSize',fsz);

h=hdEEi3;
clear data_vedi;
clear data_vhi;

ed=h.BinEdges;
step=ed(2)-ed(1);
data_vedi=[]; data_vhi=[]; c=1;
for i=1:h.NumBins
    if h.Values(i)>0
        data_vedi(c)=(ed(i)+ed(i+1))*0.5;
        data_vhi(c)=h.Values(i);
        c=c+1;
    end
end
W=[[0 data_vedi]', [0 data_vhi]'];
ct=cumsum(data_vhi); data_vhi=sum(data_vhi)-ct(1:end-1);
ldei=log10(data_vedi(1:end-1)); ldhi=log10(data_vhi*step);

figure(12); hold on;
subplot(1,2,1); hold on;
plot([0 ldei],[0 ldhi],'linewidth',3,'color',[0 0 1]);
set(gca,'linewidth',2); set(gca,'FontSize',fsz);


h=hdEEo3; 
clear data_vedo;
clear data_vho;

ed=h.BinEdges;
step=ed(2)-ed(1);
data_vedi=[]; data_vhi=[]; c=1;
for i=1:h.NumBins
    if h.Values(i)>0
        data_vedo(c)=(ed(i)+ed(i+1))*0.5;
        data_vho(c)=h.Values(i);
        c=c+1;
    end
end
W=[[0 data_vedo]', [0 data_vho]'];
ct=cumsum(data_vho); data_vho=sum(data_vho)-ct(1:end-1);
ldeo=log10(data_vedo(1:end-1)); ldho=log10(data_vho*step);

figure(12); hold on;
subplot(1,2,2); hold on;
plot([0 ldeo],[0 ldho],'linewidth',3,'color',[0 0 1]);
set(gca,'linewidth',2); set(gca,'FontSize',fsz);

h=hdEEo4; 
clear data_vedo;
clear data_vho;

ed=h.BinEdges;
step=ed(2)-ed(1);
data_vedi=[]; data_vhi=[]; c=1;
for i=1:h.NumBins
    if h.Values(i)>0
        data_vedo(c)=(ed(i)+ed(i+1))*0.5;
        data_vho(c)=h.Values(i);
        c=c+1;
    end
end
W=[[0 data_vedo]', [0 data_vho]'];
ct=cumsum(data_vho); data_vho=sum(data_vho)-ct(1:end-1);
ldeo=log10(data_vedo(1:end-1)); ldho=log10(data_vho*step);

figure(12); hold on;
subplot(1,2,2); hold on;
plot([0 ldeo],[0 ldho],'linewidth',3,'color',[0 0 1]);
set(gca,'linewidth',2); set(gca,'FontSize',fsz);

figure(9); hold on;
subplot(1,2,1); axis([0 3 -4 0]);
subplot(1,2,2); axis([0 3 -4 0]);

figure(12); hold on;
subplot(1,2,1); axis([0 3 -4 0]);
subplot(1,2,2); axis([0 3 -4 0]);

figure(1); close(gcf);
figure(10); close(gcf);




