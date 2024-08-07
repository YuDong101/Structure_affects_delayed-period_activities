clear;clc 
addpath 'tool\' 'connect\' 'func_analysis\' 'main_fun\'

%% Parameters
Naver = 1000; Np=1; step=0.01; rng(1);
NR_node=100;  NS_node=100;
%%
gNMDA_R = 0.13;   gAMPA_R = 0.11;   gGABA_R = 0.47; gStoR = 0.22; gRtoS = 0.22; gPoisson=0.01;

%% 
tic

for kkk = 1:1
    for ppp = 1:Np
        parfor jjj = 1:Naver
            rng(jjj+1000);
            KEtoE=20; KEtoI=20; KItoE=20; KItoI=20;
            [matrixR] = func_WS_network(NR_node,KEtoE,KEtoI,KItoE,KItoI);
            KEtoE1 = 5;
            [matrixS] = func_WS_network(NS_node,KEtoE1,KEtoI,KItoE,KItoI);
            [matrixRtoS] = func_RtoS(NR_node,NS_node,0.05);
            [matrixStoR] = func_StoR(NS_node,NR_node,0.2);
%             rng(jjj);
            [time,outputS(jjj),outputR(jjj),Nt(jjj),break_time(jjj),sum_Energy_ions(jjj),time_fir_R(jjj),n_fir_R(jjj),time_fir_S(jjj),n_fir_S(jjj)] = ...
            trial_R(gPoisson,NS_node,NR_node,step,gRtoS,gAMPA_R,gNMDA_R,gGABA_R,gStoR,matrixS,matrixR,matrixRtoS,matrixStoR);
            
            T_persistent(ppp,jjj) = break_time(jjj);
            mean_Energy(ppp,jjj) = sum_Energy_ions(jjj);
            matrix_rec_R(kkk,ppp,jjj,:,:)=matrixR;
            matrix_rec_S(kkk,ppp,jjj,:,:)=matrixS;
            matrix_rec_RtoS(kkk,ppp,jjj,:,:)=matrixRtoS;
            matrix_rec_StoR(kkk,ppp,jjj,:,:)=matrixStoR;
            
            fprintf('ppp %d/%d jjj %d/%d Duration= %f\n',ppp,Np,jjj,Naver,T_persistent(ppp,jjj));
        end

        mean_T_persistent=sum(T_persistent)/Naver;
        toc
    end
end
%%
[T_p,P_Tp] = fun_histogram_Tdur(T_persistent);
figure (8)
subplot (1,1,1), bar(log10(T_p),P_Tp/Naver); axis([log10(100) log10(5100),-inf inf]);
save A-NMDA=0.13_AMPA=0.11_GABA=0.47_gStoR=0.22_gPoisson=0.01-1000runs-ran=1000.mat