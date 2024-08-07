clear;clc
addpath 'tool\' 'connect\' 'func_analysis\' 'main_fun\' 'fun\'
%% Parameters 
Naver=1; step=0.01; rng(1436,'threefry'); % Random number control network structure
tic
NR_node=100;  NS_node=100;
%% Fixed adjacency matrix parameters
KEtoE=20; KEtoI=20; KItoE=20; KItoI=20;
[matrixR] = func_WS_network(NR_node,KEtoE,KEtoI,KItoE,KItoI);
KEtoE1 = 5;
[matrixS] = func_WS_network(NS_node,KEtoE1,KEtoI,KItoE,KItoI);
[matrixRtoS] = func_RtoS(NR_node,NS_node,0.05);
[matrixStoR] = func_StoR(NS_node,NR_node,0.2);

gNMDA_R = 0.13;   gAMPA_R = 0.11;   gGABA_R = 0.47; gStoR = 0.22; gRtoS = 0.22; gPoisson=0.01;

%%
[time,outputS,outputR,Nt,break_time,sum_Energy_ions,time_fir_R,n_fir_R,time_fir_S,n_fir_S] = ...
    trial_R(gPoisson,NS_node,NR_node,step,gRtoS,gAMPA_R,gNMDA_R,gGABA_R,gStoR,...
    matrixS,matrixR,matrixRtoS,matrixStoR);
%%
[time_firR,n_firR,output_spike_R] = single_text(time,outputR,Nt,NS_node);
[time_firS,n_firS,output_spike_S] = single_text(time,outputS,Nt,NR_node);

figure (20)
subplot(2,1,1),plot(time_firS,n_firS,'ro'); axis([0 6500,-inf inf]);  % %%%%%%%%%%%%%%%%%%%%%%
subplot(2,1,2),plot(time_firR,n_firR,'ro'); axis([0 6500,-inf inf]);  % %%%%%%%%%%%%%%%%%%%%%%
title ('Spiking evolution diagram')

%% firing rate
delta_t=round(10/step); t_min=0; t_max=round(6500/step);

[tS,num_spikeS] = fun_Nspike(output_spike_S,delta_t,t_min,t_max);
[tR,num_spikeR] = fun_Nspike(output_spike_R,delta_t,t_min,t_max);
mean_firingrate_S = num_spikeS/(delta_t*step*0.001*(NS_node));
mean_firingrate_R = num_spikeR/(delta_t*step*0.001*(NR_node));

figure (22)
subplot(4,1,1),plot(tS*step,mean_firingrate_S); axis([-inf inf,-inf inf]);  %
subplot(4,1,2),plot(tR*step,mean_firingrate_R); axis([-inf inf,-inf inf]);  %
