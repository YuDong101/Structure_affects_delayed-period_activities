clear;clc 
addpath 'tool\' 'connect\' 'func_analysis\' 'main_fun\'

%% Parameters
Naver = 50; Np=241; step=0.01; rng(1);
NR_node=100;  NS_node=100;
%%
gNMDA_R = 0.12;   gAMPA_R = 0.11;   gGABA_R = 0.5; gStoR = 0.22; gRtoS = 0.22; gPoisson=0.01;
gmax=0.18; gmin=0.1; gstep=(gmax-gmin)/(Np-1);
gNMDA_R=gmin:gstep:gmax;
%% 
tic
for kkk = 1:1
    for ppp = 122:Np
        v_rand=rand(Naver);
        KEtoE=20;
        KEtoI=20; KItoE=20; KItoI=20;

        TTT=T_persistent(ppp-1,:);

        parfor jjj = 1:Naver
            rng(jjj+1000);
            [matrixR] = func_WS_network(NR_node,KEtoE,KEtoI,KItoE,KItoI);
            KEtoE1 = 5;
            [matrixS] = func_WS_network(NS_node,KEtoE1,KEtoI,KItoE,KItoI);            
            [matrixRtoS] = func_RtoS(NR_node,NS_node,0.05);
            [matrixStoR] = func_StoR(NS_node,NR_node,0.2);
            
            if (TTT(jjj)<4000)
                [time,outputS(jjj),outputR(jjj),Nt(jjj),break_time(jjj),sum_Energy_ions(jjj),time_fir_R(jjj),n_fir_R(jjj),time_fir_S(jjj),n_fir_S(jjj)] = ...
                trial_R(gPoisson,NS_node,NR_node,step,gRtoS,gAMPA_R,gNMDA_R(ppp),gGABA_R,gStoR,matrixS,matrixR,matrixRtoS,matrixStoR);
                
                T_persistent(ppp,jjj) = break_time(jjj);
                mean_Energy(ppp,jjj) = sum_Energy_ions(jjj);
                matrix_rec_R(kkk,ppp,jjj,:,:)=matrixR;
                matrix_rec_S(kkk,ppp,jjj,:,:)=matrixS;
                matrix_rec_RtoS(kkk,ppp,jjj,:,:)=matrixRtoS;
                matrix_rec_StoR(kkk,ppp,jjj,:,:)=matrixStoR;                
            else
                T_persistent(ppp,jjj) = 4000;            
            end
            fprintf('ppp %d/%d jjj %d/%d 持续时间= %f\n',ppp,Np,jjj,Naver,T_persistent(ppp,jjj));
        end

        mean_T_persistent=sum(T_persistent)/Naver;
        
        toc
    end
end
save Duration-NMDA=0.1-0.14_AMPA=0.11_GABA=0.47.mat
%%
% gNMDA_R=0.001:gstep:0.001+gstep*40; nnn=length(gNMDA_R);

% T_persistent(end:80,:)=4000; 
% T_persistent(:,[31,44])=[]; T_persistent(:,[3,10,13,23])=[];

meanT=mean(T_persistent,2); errorT=std(T_persistent');
%%
figure (7)
subplot (3,1,1), errorbar (gNMDA_R,meanT,errorT'); axis([-inf inf,-inf inf]);  % semilogy
subplot (3,1,2), plot (gNMDA_R,T_persistent(:,38)); axis([-inf inf,-inf inf]);  % semilogy
subplot (3,1,3), plot (gNMDA_R,errorT); axis([-inf inf,-inf inf]);  % semilogy
%%
select_arr=[1,5,6,7,8,12,15,22,37,40];
figure (9)
subplot (2,1,1), plot (gNMDA_R,T_persistent(:,select_arr)); axis([-inf inf,-inf inf]);  % semilogy
subplot (2,1,2), plot (gNMDA_R,T_persistent(:,11)); axis([-inf inf,-inf inf]);  % semilogy
%%
figure (8)
subplot (1,1,1), plot (gNMDA_R,T_persistent,'ro'); axis([-inf inf,-inf inf]);  % semilogy
