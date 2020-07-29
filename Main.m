
clear all 
close all

load('data.mat')

[sim_dinEE,sim_doutEE]=model_EE(dinEE,doutEE,dinII,doutII);
[sim_dinII,sim_doutII]=model_II(dinEE,doutEE,dinII,doutII);
close all
plot_figure(dinEE,doutEE,dinII,doutII,sim_dinEE,sim_doutEE,sim_dinII,sim_doutII,dinEE4,doutEE4,dinII4,doutII4,dinEE3,doutEE3,dinII3,doutII3);


