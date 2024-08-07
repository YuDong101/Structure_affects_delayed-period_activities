clear;clc 
addpath 'tool\' 'connect\' 'func_analysis\' 'main_fun\'

%% Parameters
Nj = 48; Np=48; step=0.01; rng(1);
NR_node=100;  NS_node=100;
%%

gNMDA_R = 0.12;   gAMPA_R = 0.11;   gGABA_R = 0.5; gStoR = 0.22; gRtoS = 0.22; gPoisson=0.01;
gmax=0.13; gmin=0.11; gstep=(gmax-gmin)/(Np-1);
gNMDA_R=gmin:gstep:gmax;
% gmax=0.2; gmin=0.05; gstep=(gmax-gmin)/(Np-1);
% gAMPA_R=gmin:gstep:gmax;
gmax=0.6; gmin=0.4; gstep=(gmax-gmin)/(Np-1);
gGABA_R = gmin:gstep:gmax;
%% Fixed network structure
KEtoE=20; KEtoI=20; KItoE=20; KItoI=20;
[matrixR] = func_WS_network(NR_node,KEtoE,KEtoI,KItoE,KItoI);
KEtoE1 = 5;
[matrixS] = func_WS_network(NS_node,KEtoE1,KEtoI,KItoE,KItoI);
[matrixRtoS] = func_RtoS(NR_node,NS_node,0.05);
[matrixStoR] = func_StoR(NS_node,NR_node,0.2);
%% 
tic
for kkk = 1:1
    for ppp = 1:Np
        parfor jjj = 1:Nj
            rng(1);
            [time,outputS(ppp,jjj),outputR(ppp,jjj),Nt(ppp,jjj),break_time(ppp,jjj),sum_Energy_ions(ppp,jjj),time_fir_R(ppp,jjj),n_fir_R(ppp,jjj),time_fir_S(ppp,jjj),n_fir_S(ppp,jjj)] = ...
            trial_R(gPoisson,NS_node,NR_node,step,gRtoS,gAMPA_R,gNMDA_R(jjj),gGABA_R(ppp),gStoR,matrixS,matrixR,matrixRtoS,matrixStoR);
            
            T_persistent(ppp,jjj) = break_time(ppp,jjj);
            mean_Energy(ppp,jjj) = sum_Energy_ions(ppp,jjj);
            fprintf('ppp %d/%d jjj %d/%d Duration= %f\n',ppp,Np,jjj,Nj,T_persistent(ppp,jjj));
        end
        toc
        save Duration-NMDA=0.11-0.13_AMPA=0.11_GABA=0.4-0.6.mat
    end
end

%%
heatmap(break_time)