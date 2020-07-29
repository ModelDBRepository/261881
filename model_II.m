function [sim_dinII,sim_doutII]=model_II(dinEE,doutEE,dinII,doutII)

nE=length(dinEE); nI=length(dinII);

m0=100; rho=0.1; Ek=10; l=15;
pup=0.9; pdw=0.0005;

[sim_dinII,sim_doutII]=conv_model(dinII(1:end-1),doutII(1:end-1),Ek,l,rho,m0,pup,pdw,2.75);
