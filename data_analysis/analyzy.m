
clc;clear; load A-NMDA=0.13_AMPA=0.11_GABA=0.47_gStoR=0.22_gPoisson=0.01-1000runs-ran=1000
addpath 'tool\' 'connect\' 'func_analysis\' 'main_fun\' 'BCT' 'BCT\data_and_demos' ResultAnal\

kkk=1; sss=1; tic
for ppp = 1:Np
    for sss = 1:Naver
        [matrixx] = matrix_trans(matrix_rec_R,sss,ppp);
        matrix=matrixx;
        matrixx_E=matrix(1:round(0.8*NR_node),1:round(0.8*NR_node));
        matrix_E=matrixx_E;
        [Dds_out(kkk,:),Dds_out_to_E(kkk,:),Dds_out_to_I(kkk,:),Dds_in(kkk,:),Dds_in_from_E(kkk,:),Dds_in_from_I(kkk,:),Dds_avg_out(kkk,:),Dds_avg_in(kkk,:)]...
            = func_Degree_Distribution(matrix); % 出入度的计算 不在乎来源
        [Dds_EfromE(kkk,:),Dds_EfromI(kkk,:),Dds_IfromE(kkk,:),Dds_IfromI(kkk,:)]=func_Degree_Distribution_E_from_I(matrix); % 来源为 E
        [Dds_fromE(kkk,:),Dds_fromI(kkk,:)]=func_Degree_Distribution_EI(matrix);
        [Cc,Cc_avg(kkk,:)]          = func_Cluster_Coeff(matrix);
        [Cc_E,Cc_avg_E(kkk,:)]      = func_Cluster_Coeff(matrix_E);
        [Lens_E,Lens_avg_E(kkk,:)]  = func_Path_Length(matrix_E);
        [Lens,Lens_avg(kkk,:)]      = func_Path_Length(matrix);
        [sum_Num_Cycle(kkk,:),sum_Num_Cycle_E(kkk,:)] = Calculate_Basis_Cycle(matrix);
        [sum_Num_All_Cycle(kkk,:),sum_Num_All_Cycle_E(kkk,:)] = Calculate_All_Cycle(matrix);
        [f_motifs(kkk,:),F_motifs(kkk,:,:)]=motif3struct_bin(matrix);
        [f_motifs_E(kkk,:),F_motifs_E(kkk,:,:)]=motif3struct_bin(matrix_E);
        T_persistent_c(kkk) = T_persistent(ppp,sss);
        mean_Energy_c(kkk) = mean_Energy(ppp,sss);
        kkk=kkk+1;
%         kkk
    end
end
[T_p,P_Tp] = fun_histogram_Tdur(T_persistent_c);
figure (8)
subplot (1,1,1), bar(log10(T_p),P_Tp/Naver); axis([log10(100) log10(5100),-inf inf]);  % semilogy
title 'Duration distribution'

%% motif
arry_persist=find(T_persistent_c<=4000);
P_motifs=f_motifs_E./f_motifs;
 P_motifs(isnan(P_motifs))=0;%% The number of motifs 13 is relatively small, and there may be a situation where f_motifs is 0. Exclude it
d_motifs=f_motifs-f_motifs_E; d_Fmotifs=F_motifs(:,:,1:80)-F_motifs_E;

flag_X=P_motifs;
figure (7)
subplot (4,4,1), plot (flag_X(arry_persist,1),T_persistent_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-1  % semilogy
subplot (4,4,2), plot (flag_X(arry_persist,2),T_persistent_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-2  % semilogy
subplot (4,4,3), plot (flag_X(arry_persist,3),T_persistent_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-3  % semilogy
subplot (4,4,4), plot (flag_X(arry_persist,4),T_persistent_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-4  % semilogy
subplot (4,4,5), plot (flag_X(arry_persist,5),T_persistent_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-5  % semilogy
subplot (4,4,6), plot (flag_X(arry_persist,6),T_persistent_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-6  % semilogy
subplot (4,4,7), plot (flag_X(arry_persist,7),T_persistent_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-7  % semilogy
subplot (4,4,8), plot (flag_X(arry_persist,8),T_persistent_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-8  % semilogy
subplot (4,4,9), plot (flag_X(arry_persist,9),T_persistent_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-9  % semilogy
subplot (4,4,10), plot (flag_X(arry_persist,10),T_persistent_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-10  % semilogy
subplot (4,4,11), plot (flag_X(arry_persist,11),T_persistent_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-11  % semilogy
subplot (4,4,12), plot (flag_X(arry_persist,12),T_persistent_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-12  % semilogy
subplot (4,4,14), plot (flag_X(arry_persist,13),T_persistent_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-13  % semilogy
subplot (4,4,16), plot (mean_Energy_c,T_persistent_c,'ro'); axis([-inf inf,-inf inf]); title T-E  % semilogy
figure (8)
subplot (4,4,1), plot (flag_X(arry_persist,1),mean_Energy_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-1  % semilogy
subplot (4,4,2), plot (flag_X(arry_persist,2),mean_Energy_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-2  % semilogy
subplot (4,4,3), plot (flag_X(arry_persist,3),mean_Energy_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-3  % semilogy
subplot (4,4,4), plot (flag_X(arry_persist,4),mean_Energy_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-4  % semilogy
subplot (4,4,5), plot (flag_X(arry_persist,5),mean_Energy_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-5  % semilogy
subplot (4,4,6), plot (flag_X(arry_persist,6),mean_Energy_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-6  % semilogy
subplot (4,4,7), plot (flag_X(arry_persist,7),mean_Energy_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-7  % semilogy
subplot (4,4,8), plot (flag_X(arry_persist,8),mean_Energy_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-8  % semilogy
subplot (4,4,9), plot (flag_X(arry_persist,9),mean_Energy_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-9  % semilogy
subplot (4,4,10), plot (flag_X(arry_persist,10),mean_Energy_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-10  % semilogy
subplot (4,4,11), plot (flag_X(arry_persist,11),mean_Energy_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-11  % semilogy
subplot (4,4,12), plot (flag_X(arry_persist,12),mean_Energy_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-12  % semilogy
subplot (4,4,14), plot (flag_X(arry_persist,13),mean_Energy_c(arry_persist),'ro'); axis([-inf inf,-inf inf]); title class-13  % semilogy

aaaa=corr(T_persistent_c',mean_Energy_c',"type","Spearman");
corr_Emotifs=corr(f_motifs_E(arry_persist,:),T_persistent_c(arry_persist)',"type","Spearman");
corr_Emotifs_Energy=corr(f_motifs_E(arry_persist,:),mean_Energy_c(arry_persist)',"type","Spearman");

% [~, p_value] = ttest(T_persistent_c' - corr_Emotifs'.* f_motifs_E);
aaa=isnan(P_motifs(:,13))==0; %% The number of motifs 13 is relatively small, and there may be a situation where f_motifs is 0. Exclude it
corr_Pmotifs=corr(P_motifs(aaa,:),T_persistent_c(aaa)',"type","Spearman");
corr_Pmotifs_Energy=corr(P_motifs(aaa,:),mean_Energy_c(aaa)',"type","Spearman");

corr_motifs=corr(f_motifs(arry_persist,:),T_persistent_c(arry_persist)',"type","Spearman");
corr_motifs_Energy=corr(f_motifs(arry_persist,:),mean_Energy_c(arry_persist)',"type","Spearman");

corr_dmotifs=corr(d_motifs(arry_persist,:),T_persistent_c(arry_persist)',"type","Spearman");
corr_dmotifs_Energy=corr(d_motifs(arry_persist,:),mean_Energy_c(arry_persist)',"type","Spearman");

%%
ppp=[corr_dmotifs,corr_motifs,corr_Emotifs,corr_Pmotifs];
pppp=[corr_dmotifs_Energy,corr_motifs_Energy,corr_Emotifs_Energy,corr_Pmotifs_Energy];
figure (13)
subplot(211), plot (1:4,ppp); axis([-inf inf,-inf inf]);
subplot(212), plot (1:4,pppp); axis([-inf inf,-inf inf]);

[hhh,pmotif]=ttest2(corr_Emotifs,corr_Pmotifs)
%% EI-balance
aaEE=sum(Dds_EfromE(:,:),2);
aaEI=sum(Dds_EfromI(:,:),2);
aaIE=sum(Dds_IfromE(:,:),2);
aaII=sum(Dds_IfromI(:,:),2);

P_Dds_IE=(aaII+aaEE)./(aaEI+aaIE); % (aa4+aa1)./(aa2+aa3)

figure (12)
subplot (2,2,1), plot (aaEE,T_persistent_c,'ro'); axis([-inf inf,-inf inf]); title 'E to E'  % semilogy
subplot (2,2,2), plot (aaII,T_persistent_c,'ro'); axis([-inf inf,-inf inf]); title 'I to I' % semilogy
subplot (2,2,3), plot (aaEI,T_persistent_c,'ro'); axis([-inf inf,-inf inf]); title 'I to E'  % semilogy
subplot (2,2,4), plot (aaIE,T_persistent_c,'ro'); axis([-inf inf,-inf inf]); title 'E to I' % semilogy

clear rr pp
[rr(1),pp(1)]=corr(aaEE,T_persistent_c',"type","Spearman");
[rr(2),pp(2)]=corr(aaII,T_persistent_c',"type","Spearman");
[rr(3),pp(3)]=corr(aaEI,T_persistent_c',"type","Spearman");
[rr(4),pp(4)]=corr(aaIE,T_persistent_c',"type","Spearman");
[rr(5),pp(5)]=corr(P_Dds_IE,T_persistent_c',"type","Spearman");
rr,pp

figure (13)
subplot (1,1,1), semilogx (P_Dds_IE,T_persistent_c,'ro'); axis([-inf inf,-inf inf]);
%% Rich_Hub
NE=round(0.8*NR_node);
P_Dds_EI_single=Dds_fromE./(Dds_fromE+Dds_fromI);

[ai,aj] = size(P_Dds_EI_single); aa(1:ai,1:aj) = 0; aa(P_Dds_EI_single(:,1:NE)>0.8)=1;  bb(1:ai,1:aj) = 0; bb(P_Dds_EI_single(:,NE+1:100)>0.8)=1;
P_Dds_E_large=sum(aa,2)./NE; P_Dds_I_large=sum(bb,2)./20;
Raito_richHub=P_Dds_E_large./P_Dds_I_large;

figure (6)
subplot (311), plot (P_Dds_E_large,T_persistent_c,'ro'); axis([-inf inf,-inf inf]);  % semilogy
subplot (312), plot (P_Dds_I_large,T_persistent_c,'ro'); axis([-inf inf,-inf inf]);  % semilogy
subplot (313), plot (Raito_richHub,T_persistent_c,'ro'); axis([-inf inf,-inf inf]);  % semilogy

clear rr pp
[rr(1),pp(1)]=corr(P_Dds_E_large,T_persistent_c',"type","Spearman");
[rr(2),pp(2)]=corr(P_Dds_I_large,T_persistent_c',"type","Spearman");
[rr(3),pp(3)]=corr(Raito_richHub,T_persistent_c',"type","Spearman");
rr,pp

%% Small_World
P_lens=Lens_avg_E./Lens_avg; P_Cc=Cc_avg_E./Cc_avg;
Small_World=P_lens./P_Cc;

figure (13)
subplot (2,3,1), plot (Cc_avg,T_persistent_c,'ro'); axis([-inf inf,-inf inf]);
subplot (2,3,2), plot (Cc_avg_E,T_persistent_c,'ro'); axis([-inf  inf,-inf inf]);
subplot (2,3,3), plot (P_Cc,T_persistent_c,'ro'); axis([-inf inf,-inf inf]);
subplot (2,3,4), plot (Lens_avg,T_persistent_c,'ro'); axis([-inf inf,-inf inf]);
subplot (2,3,5), plot (Lens_avg_E,T_persistent_c,'ro'); axis([-inf inf,-inf inf]);
subplot (2,3,6), plot (P_lens,T_persistent_c,'ro'); axis([-inf inf,-inf inf]);

clear rr pp
[rr(1),pp(1)]=corr(Cc_avg,T_persistent_c',"type","Spearman");
[rr(2),pp(2)]=corr(Cc_avg_E,T_persistent_c',"type","Spearman");
[rr(3),pp(3)]=corr(P_Cc,T_persistent_c',"type","Spearman");
[rr(4),pp(4)]=corr(Lens_avg,T_persistent_c',"type","Spearman");
[rr(5),pp(5)]=corr(Lens_avg_E,T_persistent_c',"type","Spearman");
[rr(6),pp(6)]=corr(P_lens,T_persistent_c',"type","Spearman");
rr,pp

%%
Cc_corr_All=corr(Cc_avg,T_persistent_c',"type","Spearman"); Cc_corr_E=corr(Cc_avg_E,T_persistent_c',"type","Spearman"); Cc_corr_P=corr(P_Cc,T_persistent_c',"type","Spearman");
lens_corr_All=corr(Lens_avg,T_persistent_c',"type","Spearman"); lens_corr_E=corr(Lens_avg_E,T_persistent_c',"type","Spearman"); lens_corr_P=corr(P_lens,T_persistent_c',"type","Spearman");
corr_Cc=[Cc_corr_All,Cc_corr_E,Cc_corr_P];
corr_lens=[lens_corr_All,lens_corr_E,lens_corr_P];
figure (14)
subplot(211), plot (1:3,corr_Cc); axis([-inf inf,-inf inf]);
subplot(212), plot (1:3,corr_lens); axis([-inf inf,-inf inf]);

%% Cycle

P_Num_All_Cycle=sum_Num_All_Cycle_E(:,3)./(sum_Num_All_Cycle(:,3)+sum_Num_All_Cycle_E(:,3));
nor_P_Num_All_Cycle=(P_Num_All_Cycle-min(P_Num_All_Cycle))./(max(P_Num_All_Cycle)-min(P_Num_All_Cycle));
ratio_Cycle=sum_Num_All_Cycle_E(:,3)./(sum_Num_All_Cycle(:,3));

figure (14)
subplot (2,1,1), plot (P_Num_All_Cycle,T_persistent_c,'ro'); axis([-inf inf,-inf inf]);  % semilogy
subplot (2,1,2), plot (ratio_Cycle,T_persistent_c,'ro'); axis([-inf inf,-inf inf]);  % semilogy

figure (15)
subplot (4,1,2), plot (sum_Num_All_Cycle_E(:,3),T_persistent_c,'ro'); axis([-inf inf,-inf inf]);  % semilogy
subplot (4,1,3), plot (sum_Num_All_Cycle(:,3),T_persistent_c,'ro'); axis([-inf inf,-inf inf]);  % semilogy
subplot (4,1,1), plot (sum_Num_All_Cycle(:,3)+sum_Num_All_Cycle_E(:,3),T_persistent_c,'ro'); axis([-inf inf,-inf inf]);  % semilogy
subplot (4,1,4), plot (P_Num_All_Cycle,T_persistent_c,'ro'); axis([-inf inf,-inf inf]);  % semilogy


Num_All_Cycle_E=sum_Num_All_Cycle_E(:,3);
Num_All_Cycle_d=sum_Num_All_Cycle(:,3);
Num_All_Cycle=sum_Num_All_Cycle(:,3)+sum_Num_All_Cycle_E(:,3);


clear rr pp
[rr(1),pp(1)]=corr(sum_Num_All_Cycle(:,3)+sum_Num_All_Cycle_E(:,3),T_persistent_c',"type","Spearman");
[rr(2),pp(2)]=corr(sum_Num_All_Cycle_E(:,3),T_persistent_c',"type","Spearman");
[rr(3),pp(3)]=corr(sum_Num_All_Cycle(:,3),T_persistent_c',"type","Spearman");
[rr(4),pp(4)]=corr(P_Num_All_Cycle,T_persistent_c',"type","Spearman");
rr,pp

%%  single neuron level *******************************************************
clear corrEI_r corrEI_p corrCycle_r corrCycle_p
TatolDuration=[];
TatolDds_fromE=[]; TatolDds_fromI=[]; TatolEI_ratio=[];
lens_ap=length(find(T_persistent_c==4000));  lens_at=length(find(T_persistent_c<4000));

for iii=1:lens_ap
    index_singel=iii+lens_at;
% for iii=1:1000
%     index_singel=iii;
    analyzy_singel
    corrEI_r(iii,:)=EI_r;corrEI_p(iii,:)=EI_p;
    corrCycle_r(iii,:)=Cycle_r;corrCycle_p(iii,:)=Cycle_p;
end
%%
mean_corrEI=mean(corrEI_r,1); std_corrEI=std(corrEI_r,0,1);
mean_corrEI=mean(corrEI_r,1); std_corrEI=std(corrEI_r,0,1);

figure(3);
subplot(311),bar(mean_corrEI);
hold on;
errorbar(1:length(mean_corrEI), mean_corrEI, std_corrEI, 'k', 'linestyle', 'none');

mean_corrCycle=mean(corrCycle_r,1);
std_corrCycle=std(corrCycle_r,0,1);
subplot(312),bar(mean_corrCycle);
hold on;
errorbar(1:length(mean_corrCycle), mean_corrCycle, std_corrCycle, 'k', 'linestyle', 'none');

figure (7)
subplot(221); plot(TatolDds_fromE,TatolDuration,'ro');
subplot(222); plot(TatolDds_fromI,TatolDuration,'ro');
subplot(223); plot(TatolEI_ratio,TatolDuration,'ro');
subplot(4,2,6); histogram(TatolDuration, 'Normalization', 'probability', 'NumBins', 200);
subplot(4,2,8); histogram(TatolEI_ratio, 'Normalization', 'probability', 'NumBins', 200);
%%

index_singel=950;
analyzy_singel
figure (2100)
subplot(221); plot(Dds_fromE(NN_FLAG,N_OB),fir_rate(N_OB),'ro');
subplot(222); plot(Dds_fromI(NN_FLAG,N_OB),fir_rate(N_OB),'ro');
subplot(223); plot(EI_ratio(N_OB),fir_rate(N_OB),'ro');
EI_r ,EI_p
index_singel=904;
analyzy_singel
figure (2200)
subplot(221); plot(Cycle_i(N_OB),fir_rate(N_OB),'ro');
subplot(222); plot(CycleE_i(N_OB),fir_rate(N_OB),'ro');
subplot(223); plot(Cycle_d_i(N_OB),fir_rate(N_OB),'ro');
Cycle_r ,Cycle_p
