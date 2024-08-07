close all
clc
ap=find(T_persistent_c==4000);  at=find(T_persistent_c<4000);
%% small-world Cc
[x_Cc,ht_Cc,hp_Cc,fRxt_Cc,fRyt_Cc,fRxp_Cc,fRyp_Cc,kl_d_Cc]=gaussfit_and_comparison(Cc_avg,T_persistent_c,20);
[x_Cc_E,ht_Cc_E,hp_Cc_E,fRxt_Cc_E,fRyt_Cc_E,fRxp_Cc_E,fRyp_Cc_E,kl_d_Cc_E]=gaussfit_and_comparison(Cc_avg_E,T_persistent_c,20);
[x_Cc_P,ht_Cc_P,hp_Cc_P,fRxt_Cc_P,fRyt_Cc_P,fRxp_Cc_P,fRyp_Cc_P,kl_d_Cc_P]=gaussfit_and_comparison(P_Cc,T_persistent_c,20);

figure (1);
subplot (2,2,1), bar(x_Cc,[ht_Cc' hp_Cc']); hold on; title('Cc_{All}')
plot(fRxt_Cc,fRyt_Cc, 'b'); plot(fRxp_Cc,fRyp_Cc, 'r'); hold off;

subplot (2,2,2), bar(x_Cc_E,[ht_Cc_E' hp_Cc_E']); hold on; title('Cc_{E}')
plot(fRxt_Cc_E,fRyt_Cc_E, 'b'); plot(fRxp_Cc_E,fRyp_Cc_E, 'r'); hold off;

subplot (2,2,3), bar(x_Cc_P,[ht_Cc_P' hp_Cc_P']); hold on; title('Cc_{P}')
plot(fRxt_Cc_P,fRyt_Cc_P, 'b'); plot(fRxp_Cc_P,fRyp_Cc_P, 'r'); hold off;

subplot (2,2,4), plot(1:3,[kl_d_Cc,kl_d_Cc_E,kl_d_Cc_P], 'b');

%% small-world lens
[x_lens,ht_lens,hp_lens,fRxt_lens,fRyt_lens,fRxp_lens,fRyp_lens,kl_d_lens]=gaussfit_and_comparison(Lens_avg,T_persistent_c,20);
[x_lens_E,ht_lens_E,hp_lens_E,fRxt_lens_E,fRyt_lens_E,fRxp_lens_E,fRyp_lens_E,kl_d_lens_E]=gaussfit_and_comparison(Lens_avg_E,T_persistent_c,20);
[x_lens_P,ht_lens_P,hp_lens_P,fRxt_lens_P,fRyt_lens_P,fRxp_lens_P,fRyp_lens_P,kl_d_lens_P]=gaussfit_and_comparison(P_lens,T_persistent_c,20);
figure (2);
subplot (2,2,1), bar(x_lens,[ht_lens' hp_lens']); hold on; title('lens_{All}')
plot(fRxt_lens,fRyt_lens, 'b'); plot(fRxp_lens,fRyp_lens, 'r'); hold off;

subplot (2,2,2), bar(x_lens_E,[ht_lens_E' hp_lens_E']); hold on; title('lens_{E}')
plot(fRxt_lens_E,fRyt_lens_E, 'b'); plot(fRxp_lens_E,fRyp_lens_E, 'r'); hold off;

subplot (2,2,3), bar(x_lens_P,[ht_lens_P' hp_lens_P']); hold on; title('lens_{P}')
plot(fRxt_lens_P,fRyt_lens_P, 'b'); plot(fRxp_lens_P,fRyp_lens_P, 'r'); hold off;

subplot (2,2,4), plot(1:3,[kl_d_lens,kl_d_lens_E,kl_d_lens_P], 'b');

%% E/I balance
[x_aaEE,ht_aaEE,hp_aaEE,fRxt_aaEE,fRyt_aaEE,fRxp_aaEE,fRyp_aaEE,kl_d_aaEE]=gaussfit_and_comparison(aaEE,T_persistent_c,20);
[x_aaEI,ht_aaEI,hp_aaEI,fRxt_aaEI,fRyt_aaEI,fRxp_aaEI,fRyp_aaEI,kl_d_aaEI]=gaussfit_and_comparison(aaEI,T_persistent_c,20);
[x_aaIE,ht_aaIE,hp_aaIE,fRxt_aaIE,fRyt_aaIE,fRxp_aaIE,fRyp_aaIE,kl_d_aaIE]=gaussfit_and_comparison(aaIE,T_persistent_c,20);
[x_aaII,ht_aaII,hp_aaII,fRxt_aaII,fRyt_aaII,fRxp_aaII,fRyp_aaII,kl_d_aaII]=gaussfit_and_comparison(aaII,T_persistent_c,20);

[x_IE,ht_IE,hp_IE,fRxt_IE,fRyt_IE,fRxp_IE,fRyp_IE,kl_d_IE]=gaussfit_and_comparison(P_Dds_IE,T_persistent_c,20);

figure (3);
subplot (3,2,1), bar(x_aaEE,[ht_aaEE' hp_aaEE']); hold on; title('E/I Balance')
plot(fRxt_aaEE,fRyt_aaEE, 'b'); plot(fRxp_aaEE,fRyp_aaEE, 'r'); hold off;
subplot (3,2,2), bar(x_aaII,[ht_aaII' hp_aaII']); hold on; title('E/I Balance')
plot(fRxt_aaII,fRyt_aaII, 'b'); plot(fRxp_aaII,fRyp_aaII, 'r'); hold off;
subplot (3,2,3), bar(x_aaEI,[ht_aaEI' hp_aaEI']); hold on; title('E/I Balance')
plot(fRxt_aaEI,fRyt_aaEI, 'b'); plot(fRxp_aaEI,fRyp_aaEI, 'r'); hold off;
subplot (3,2,4), bar(x_aaIE,[ht_aaIE' hp_aaIE']); hold on; title('E/I Balance')
plot(fRxt_aaIE,fRyt_aaIE, 'b'); plot(fRxp_aaIE,fRyp_aaIE, 'r'); hold off;
subplot (3,2,5), bar(x_IE,[ht_IE' hp_IE']); hold on; title('E/I Balance')
plot(fRxt_IE,fRyt_IE, 'b'); plot(fRxp_IE,fRyp_IE, 'r'); hold off;

subplot (3,2,6), plot(1:5,[kl_d_aaEE,kl_d_aaII,kl_d_aaEI,kl_d_aaIE,kl_d_IE], 'b');

%% Rich_Hub
[x_LargeE,ht_LargeE,hp_LargeE,fRxt_LargeE,fRyt_LargeE,fRxp_LargeE,fRyp_LargeE,kl_d_LargeE]=gaussfit_and_comparison(P_Dds_E_large,T_persistent_c,12);
[x_LargeI,ht_LargeI,hp_LargeI,fRxt_LargeI,fRyt_LargeI,fRxp_LargeI,fRyp_LargeI,kl_d_LargeI]=gaussfit_and_comparison(P_Dds_I_large,T_persistent_c,12);
[x_RH,ht_RH,hp_RH,fRxt_RH,fRyt_RH,fRxp_RH,fRyp_RH,kl_d_RH]=gaussfit_and_comparison(Raito_richHub,T_persistent_c,20);

figure (4);
subplot (2,2,1), bar(x_LargeE,[ht_LargeE' hp_LargeE']); hold on; title('Large E')
plot(fRxt_LargeE,fRyt_LargeE, 'b'); plot(fRxp_LargeE,fRyp_LargeE, 'r'); hold off;
subplot (2,2,2), bar(x_LargeI,[ht_LargeI' hp_LargeI']); hold on; title('Large E')
plot(fRxt_LargeI,fRyt_LargeI, 'b'); plot(fRxp_LargeI,fRyp_LargeI, 'r'); hold off;
subplot (2,2,3), bar(x_RH,[ht_RH' hp_RH']); hold on; title('Large E')
plot(fRxt_RH,fRyt_RH, 'b'); plot(fRxp_RH,fRyp_RH, 'r'); hold off;

%% CYCLE
[x_cycle_All,ht_cycle_All,hp_cycle_All,fRxt_cycle_All,fRyt_cycle_All,fRxp_cycle_All,fRyp_cycle_All,kl_d_cycle_All]=gaussfit_and_comparison(Num_All_Cycle,T_persistent_c,20);
[x_cycle_d,ht_cycle_d,hp_cycle_d,fRxt_cycle_d,fRyt_cycle_d,fRxp_cycle_d,fRyp_cycle_d,kl_d_cycle_d]=gaussfit_and_comparison(Num_All_Cycle_d,T_persistent_c,20);
[x_cycle_E,ht_cycle_E,hp_cycle_E,fRxt_cycle_E,fRyt_cycle_E,fRxp_cycle_E,fRyp_cycle_E,kl_d_cycle_E]=gaussfit_and_comparison(Num_All_Cycle_E,T_persistent_c,20);
[x_cycle_P,ht_cycle_P,hp_cycle_P,fRxt_cycle_P,fRyt_cycle_P,fRxp_cycle_P,fRyp_cycle_P,kl_d_cycle_P]=gaussfit_and_comparison(P_Num_All_Cycle,T_persistent_c,20);

figure (5);
subplot (2,2,1), bar(x_cycle_All,[ht_cycle_All' hp_cycle_All']); hold on; title('cycle_{All}')
plot(fRxt_cycle_All,fRyt_cycle_All, 'b'); plot(fRxp_cycle_All,fRyp_cycle_All, 'r'); hold off;
subplot (2,2,2), bar(x_cycle_d,[ht_cycle_d' hp_cycle_d']); hold on; title('cycle_{d}')
plot(fRxt_cycle_d,fRyt_cycle_d, 'b'); plot(fRxp_cycle_d,fRyp_cycle_d, 'r'); hold off;
subplot (2,2,3), bar(x_cycle_E,[ht_cycle_E' hp_cycle_E']); hold on; title('cycle_{E}')
plot(fRxt_cycle_E,fRyt_cycle_E, 'b'); plot(fRxp_cycle_E,fRyp_cycle_E, 'r'); hold off;
subplot (2,2,4), bar(x_cycle_P,[ht_cycle_P' hp_cycle_P']); hold on; title('cycle P')
plot(fRxt_cycle_P,fRyt_cycle_P, 'b'); plot(fRxp_cycle_P,fRyp_cycle_P, 'r'); hold off;

figure (6);
plot(1:4,[kl_d_cycle_All,kl_d_cycle_d,kl_d_cycle_E,kl_d_cycle_P]);

%%
function [x_data,ht_data,hp_data,fRxt,fRyt,fRxp,fRyp,kl_divergence]...
    =gaussfit_and_comparison(topology_data,T_persistent_c,N_bin)
    ap=find(T_persistent_c==4000);  at=find(T_persistent_c<4000);
    [x_data,ht_data,hp_data]=htp(topology_data ,at,ap,N_bin); 
    ht_data=ht_data/length(at); hp_data=hp_data/length(ap);
    [fRxt,fRyt,fitResult] = gaussfitt(x_data,ht_data);
    [fRxp,fRyp,fitResult] = gaussfitt(x_data,hp_data);
    
    kl_divergence = sum(fRyt .* log(fRyt ./ fRyp), 'omitnan');
end