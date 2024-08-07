function [time,outputS,outputR,Nt,break_time,sum_Energy_ions,time_fir_R,n_fir_R,time_fir_S,n_fir_S] = trial_R(gPoisson,NS_node,NR_node,step,gRtoS,gAMPA_RR,gNMDA_RR,gGABA_RR,gStoR...
    ,matrixS,matrixR,matrixRtoS,matrixStoR)
%% Fixed adjacency matrix parameters
% clear; clc
% Naver=1; step=0.01; rng(1);
% NR_node=100;  NS_node=100;
% KEtoE=20; KEtoI=20; KItoE=20; KItoI=20;
% [matrixR] = func_WS_network(NR_node,KEtoE,KEtoI,KItoE,KItoI);
% KEtoE1 = 5;
% [matrixS] = func_WS_network(NS_node,KEtoE1,KEtoI,KItoE,KItoI);
% Per=0.2;
% [matrixRtoS] = func_RtoS(NR_node,NS_node,0.2);
% [matrixStoR] = func_StoR(NS_node,NR_node,0.05);
% gNMDA_RR = 0.14; gAMPA_RR = 0.12; gGABA_RR = 0.36; gStoR = 0.2; gRtoS = 0.2; gPoisson = 0.01;
addpath 'fun\' 'tool\' 'connect\' 'func_analysis\' 'main_fun\'
%% Par
T_tot=6500; N_print=1;
Nt=round(T_tot/step);

tao_AMPA=2; tao_GABA=10;  tao_NMDA_rise=2;  tao_NMDA_decay=100; Alpha=0.5;
gAMPA_SS = gAMPA_RR; gNMDA_SS = gNMDA_RR; gGABA_SS = gGABA_RR;

I0_R(1:NR_node)=0.0;  I0_R(round(0.8*NR_node+1):NR_node)=2;
I0_S(1:NS_node)=2; I0_S((0.8*NS_node+1):NS_node)=0;
r_poissonS=2300; r_poissonR=2300; % unite Hz
P_poisson_S=1-exp(-r_poissonS*step/1000);
P_poisson_R=1-exp(-r_poissonR*step/1000);
poissonSpikeS(1:NS_node)=0.0;
poissonSpikeR(1:NR_node)=0.0;

matrixSfromE=matrixS; matrixSfromE(round(0.8*NS_node+1):NS_node,:)=0;
matrixSfromI=matrixS; matrixSfromI(1:round(0.8*NS_node),:)=0;
matrixRtoS(round(0.8*NS_node+1):NS_node,:)=0;

matrixRfromE=matrixR; matrixRfromE(round(0.8*NR_node+1):NR_node,:)=0;
matrixRfromI=matrixR; matrixRfromI(1:round(0.8*NR_node),:)=0;
matrixStoR(round(0.8*NR_node+1):NR_node,:)=0;  matrixStoR(:,round(0.8*NR_node+1):NR_node)=0;

KSE=1; KSI=1; KRE=1; KRI=1;

gAMPA_S=gAMPA_SS.*KSE;gNMDA_S=gNMDA_SS.*KSE;gGABA_S=gGABA_SS.*KSI;
gAMPA_R=gAMPA_RR.*KRE;gNMDA_R=gNMDA_RR.*KRE;gGABA_R=gGABA_RR.*KRI;

Energy_ions(1:NR_node)=0; sum_Energy_ions=0;
time_firR=0; n_firR=0; time_firS=0; n_firS=0;
%% Var
flag_firing=0.0; break_time=0.0;

vR(1:NR_node)=-65; mR(1:NR_node)=0.0024; hR(1:NR_node)=0.9963; nR(1:NR_node)=0.0254;
vS(1:NS_node)=-65; mS(1:NS_node)=0.0024; hS(1:NS_node)=0.9963; nS(1:NS_node)=0.0254;

spike_R(1:NR_node)=0; spike_S(1:NS_node)=0; spike_Rr(1:NR_node)=0;
spike_RtoS(1:NR_node)=0; spike_StoR(1:NS_node)=0;

rsynAMPA_R(1:NR_node)=0.0; rsynGABA_R(1:NR_node)=0.0; x_NMDA_R(1:NR_node) = 0.0; rsynNMDA_R(1:NR_node) = 0.0;
rsynAMPA_S(1:NS_node)=0.0; rsynGABA_S(1:NS_node)=0.0; x_NMDA_S(1:NS_node) = 0.0; rsynNMDA_S(1:NS_node) = 0.0;
rsynAMPA_RtoS(1:NS_node)=0.0; rsynAMPA_StoR(1:NR_node)=0.0;

outputS(1:NR_node,1:Nt/N_print)=0.0; outputR(1:NR_node,1:Nt/N_print)=0.0; jjjj=1;
output_spike_R(1:NR_node,1:Nt)=0.0;
time=step:step:T_tot;

CmR(1:round(0.8*NR_node))=0.5;
CmR(round(0.8*NR_node+1):NR_node)=0.25;
CmS(1:round(0.8*NS_node))=0.5;
CmS(round(0.8*NS_node+1):NS_node)=0.25;

ss=1; sss=1; is_spikeR(1:NR_node)=0;is_spikeS(1:NS_node)=0;
%% time loop
for iii  = 1:Nt
    spike_Rr(vR>0.0)=1;    spike_Rr(vR<=0.0)=0;
    
    Rand_poissonS=rand (1,NS_node); Rand_poissonR=rand (1,NR_node);
    poissonSpikeS(P_poisson_S>Rand_poissonS)=1; poissonSpikeS(P_poisson_S<=Rand_poissonS)=0;
    poissonSpikeR(P_poisson_R>Rand_poissonR)=1; poissonSpikeR(P_poisson_R<=Rand_poissonR)=0;
    %% Sender

    delta_spike_S_E=matrixSfromE.*spike_S';  delta_spike_S_I=matrixSfromI.*spike_S';
    delta_spike_RtoS=matrixRtoS.*spike_RtoS';

    drsynAMPA_S = rsynAMPA_S + step*(-rsynAMPA_S/tao_AMPA+delta_spike_S_E);
    drsynAMPA_RtoS = rsynAMPA_RtoS + step*(-rsynAMPA_RtoS/tao_AMPA+delta_spike_RtoS+gPoisson.*poissonSpikeS); 
    drsynGABA_S = rsynGABA_S + step*(-rsynGABA_S/tao_GABA+delta_spike_S_I);
    dx_NMDA_S   = x_NMDA_S + step*(-x_NMDA_S/tao_NMDA_rise+delta_spike_S_E);
    drsynNMDA_S = rsynNMDA_S + step*(-rsynNMDA_S/tao_NMDA_decay+Alpha*x_NMDA_S.*(1-rsynNMDA_S));
    
    sum_rsynAMPA_S=sum(drsynAMPA_S,1);
    sum_rsynAMPA_RtoS=sum(drsynAMPA_RtoS,1);
    sum_rsynNMDA_S=sum(drsynNMDA_S,1);
    sum_rsynGABA_S=sum(drsynGABA_S,1);

    Irec_AMPA_S = gAMPA_S.*sum_rsynAMPA_S.*(0-vS); Irec_NMDA_S = gNMDA_S.*sum_rsynNMDA_S.*(0-vS)./(1+exp((-0.062*vS)/3.57));
    Irec_GABA_S = gGABA_S.*sum_rsynGABA_S.*(-70-vS); IextRtoS = gRtoS.*sum_rsynAMPA_RtoS.*(0-vS);
    Isny_S = Irec_AMPA_S + Irec_NMDA_S + Irec_GABA_S + IextRtoS;
    IextS=0;
    if iii>=round(0.5*1000/step) && iii<= round(1.0*1000/step)
        IextS = I0_S;
    end
    if iii>=round((T_tot-1500)/step) && iii<= round((T_tot-1000)/step)
        IextS = I0_S;
    end

    dvS = vS + step*fun_v(vS,mS,hS,nS,Isny_S,IextS,CmS);
    dmS = mS + step*fun_m(vS,mS);
    dhS = hS + step*fun_h(vS,hS);
    dnS = nS + step*fun_n(vS,nS);
    %% Receiver    
    delta_spike_R_E=matrixRfromE.*spike_R';    delta_spike_R_I=matrixRfromI.*spike_R'; 
    delta_spike_StoR=matrixStoR.*spike_StoR'; 

    drsynAMPA_R = rsynAMPA_R + step*(-rsynAMPA_R/tao_AMPA+delta_spike_R_E);
    drsynAMPA_StoR = rsynAMPA_StoR + step*(-rsynAMPA_StoR/tao_AMPA+delta_spike_StoR+gRtoS/gStoR*gPoisson.*poissonSpikeR); 
    drsynGABA_R = rsynGABA_R + step*(-rsynGABA_R/tao_GABA+delta_spike_R_I);
    dx_NMDA_R   = x_NMDA_R + step*(-x_NMDA_R/tao_NMDA_rise+delta_spike_R_E);
    drsynNMDA_R = rsynNMDA_R + step*(-rsynNMDA_R/tao_NMDA_decay+Alpha*x_NMDA_R.*(1-rsynNMDA_R));
    
    sum_rsynAMPA_R=sum(drsynAMPA_R,1);
    sum_rsynAMPA_StoR=sum(drsynAMPA_StoR,1);
    sum_rsynNMDA_R=sum(drsynNMDA_R,1);
    sum_rsynGABA_R=sum(drsynGABA_R,1);

    Irec_AMPA_R = gAMPA_R.*sum_rsynAMPA_R.*(0-vR); Irec_NMDA_R = gNMDA_R.*sum_rsynNMDA_R.*(0-vR)./(1+exp((-0.062*vR)/3.57));
    Irec_GABA_R = gGABA_R.*sum_rsynGABA_R.*(-70-vR); IextStoR = gStoR.*sum_rsynAMPA_StoR.*(0-vR);
    Isny_R = Irec_AMPA_R + Irec_NMDA_R + Irec_GABA_R + IextStoR;
    IextR = 0;
    if time(iii)>(T_tot-1000) && time(iii)<(T_tot-700)
        IextR=I0_R;
%         IextR=2;
    end

    [funv,Energy_ion] = fun_v(vR,mR,hR,nR,Isny_R,IextR,CmR);
    dvR = vR + step*funv;
    dmR = mR + step*fun_m(vR,mR);
    dhR = hR + step*fun_h(vR,hR);
    dnR = nR + step*fun_n(vR,nR);
    
    if time(iii)>1000, Energy_ions=Energy_ions+Energy_ion.*step; end

    if mod(iii,N_print)==0
        outputS(:,jjjj)=vS; outputR(:,jjjj)=vR;
        jjjj=jjjj+1;
    end

    output_spike_R(:,iii)=spike_Rr;

%     LFP_S(iii)=sum(vS(1:round(0.8*NS_node)))/(0.8*NS_node); LFP_R(iii)=sum(vR(1:round(0.8*NR_node)))/(0.8*NR_node);
%     LFP_S(iii)=sum(Isny_S(1:round(0.8*NS_node)))/(0.8*NS_node); LFP_R(iii)=sum(Isny_R(1:round(0.8*NR_node)))/(0.8*NR_node);
    %% return value
    rsynAMPA_S    = drsynAMPA_S;
    rsynGABA_S    = drsynGABA_S;
    rsynAMPA_R    = drsynAMPA_R;
    rsynGABA_R    = drsynGABA_R;
    rsynAMPA_RtoS = drsynAMPA_RtoS;
    rsynAMPA_StoR = drsynAMPA_StoR;

    x_NMDA_S      = dx_NMDA_S;
    rsynNMDA_S    = drsynNMDA_S;
    x_NMDA_R      = dx_NMDA_R;
    rsynNMDA_R    = drsynNMDA_R;

    is_spikeR(dvR>0.0 & vR<=0.0)=1;
    index_spikeR=find(is_spikeR==1);
    if isempty(index_spikeR)==0
        for jjj=1:length(index_spikeR)
            time_firR(ss)=time(iii);
            n_firR(ss)=index_spikeR(jjj);
            ss=ss+1;
        end
    end
    is_spikeR(1:NR_node)=0;

    is_spikeS(dvS>0.0 & vS<=0.0)=1;
    index_spikeS=find(is_spikeS==1);
    if isempty(index_spikeS)==0
        for jjj=1:length(index_spikeS)
            time_firS(sss)=time(iii);
            n_firS(sss)=index_spikeS(jjj);
            sss=sss+1;
        end
    end
    is_spikeS(1:NS_node)=0;
    
    spike_R(1:NR_node)=0; spike_S(1:NS_node)=0;
    spike_RtoS(1:NR_node)=0; spike_StoR(1:NS_node)=0;

    spike_R(dvR>0.0 & vR<=0.0)=1; spike_S(dvS>0.0 & vS<=0.0)=1;
    spike_RtoS(dvR(1:round(0.8*NR_node))>0.0 & vR(1:round(0.8*NR_node))<=0.0)=1;
    spike_StoR(dvS(1:round(0.8*NS_node))>0.0 & vS(1:round(0.8*NS_node))<=0.0)=1;

    vS = dvS; mS = dmS; hS = dhS; nS = dnS;
    vR = dvR; mR = dmR; hR = dhR; nR = dnR;
%% Determine whether there is spontaneous discharge
    if iii==round(500/step) && ...
            sum(output_spike_R(1:round(0.8*NR_node),iii-round(100/step):iii),'all')~=0 %判断有没有激起自发震荡
        fprintf('Not Silent\n'); flag_Spontaneous = 1; break_time=5000; 
%         clear output_spike_R time
%         time=0;
        return;
    end
%% Calculated duration
    if iii>1000/step      
        flag_firing=flag_firing+sum(output_spike_R(1:round(0.8*NR_node),iii),'all');
        if iii>1100/step
            flag_firing=flag_firing-sum(output_spike_R(1:round(0.8*NR_node),iii-round(100/step)),'all');
            if(flag_firing==0), break_time=time(iii)-1100; break; end
        end
        if iii==Nt, break_time=time(Nt)-1000; flag_persitent = 1; end        
    end
end

time_fir_R = {time_firR}; n_fir_R = {n_firR};
time_fir_S = {time_firS}; n_fir_S = {n_firS};

sum_Energy_ions = sum(Energy_ions,"all")/(break_time*NR_node);

clear output_spike_R

%%
% figure (1)
% subplot(4,1,1),plot(time(1:N_print:end),outputS); axis([0 inf,-inf inf]);  %
% subplot(4,1,2),plot(time(1:N_print:end),outputR); axis([0 inf,-inf inf]);  %
% subplot(4,1,3),plot(time(1:N_print:end),outputS(1,:),time(1:N_print:end),outputS(NS_node,:)); axis([0 inf,-inf inf]);  %
% subplot(4,1,4),plot(time(1:N_print:end),outputR(1,:),time(1:N_print:end),outputR(NS_node,:)); axis([0 inf,-inf inf]);  %