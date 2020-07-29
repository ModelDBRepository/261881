function [sim_dinEE,sim_doutEE]=model_EE(dinEE,doutEE,dinII,doutII)

nE=length(dinEE); nI=length(dinII);

m0=100; rho=0.1; Ek=115; l=100;
pup=0.9; pdw=0.0005;

[sim_dinEE,sim_doutEE]=conv_model(dinEE(1:end-1),doutEE(1:end-1),Ek,l,rho,m0,pup,pdw,300);

