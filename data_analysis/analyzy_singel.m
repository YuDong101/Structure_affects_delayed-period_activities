clear tabulate_nfir Duration_i

[~, index_mlg] = sort(T_persistent_c);
% 
NN_FLAG=index_mlg(index_singel); % 906 index_singel
matrixR(:,:) = matrix_rec_R(1,1,NN_FLAG,:,:); matrixS(:,:) = matrix_rec_S(1,1,NN_FLAG,:,:); matrixRtoS(:,:) = matrix_rec_RtoS(1,1,NN_FLAG,:,:); matrixStoR(:,:) = matrix_rec_StoR(1,1,NN_FLAG,:,:);
% figure (300)
% G=digraph(matrixR);
% plot(G)
time_fir_c=cell2mat(time_fir_R(NN_FLAG));
n_fir_c=cell2mat(n_fir_R(NN_FLAG));

for ii=1:100
    indices = find(n_fir_c == ii);
    LastIndex = max(indices);
    if isempty(LastIndex), Duration_i(ii)=0;
    else, Duration_i(ii)=time_fir_c(LastIndex)-1000;
    end
end
Duration_i(Duration_i<0)=0;

TatolDuration=[TatolDuration Duration_i];
% figure (4000)
% bar(Duration_i);

%%  network reordering
tabulate_nfir=tabulate(n_fir_c(find(time_fir_c>=1000,1):end));
tabulate_nfir(end+1:100,1:3)=0;

fir_rate=tabulate_nfir(:,2)'./(Duration_i*0.001);
fir_rate(isnan(fir_rate))=0; fir_rate(tabulate_nfir(:,2)<2)=0;
[~,arr_recong]=sort(fir_rate);
arr_recong_E=arr_recong(arr_recong<=80); %%%%

matrixR_reco_E=matrixR(arr_recong_E,arr_recong_E);

%% EI ratio
N_OB=1:80;

[EI_r(1), EI_p(1)]=corr(fir_rate(N_OB)',Dds_fromE(NN_FLAG,N_OB)', 'Type', 'Pearson');
[EI_r(2), EI_p(2)]=corr(fir_rate(N_OB)',Dds_fromI(NN_FLAG,N_OB)', 'Type', 'Pearson');
EI_ratio=Dds_fromE(NN_FLAG,:)./(Dds_fromI(NN_FLAG,:)+Dds_fromE(NN_FLAG,:));
[EI_r(3), EI_p(3)]=corr(fir_rate(N_OB)',EI_ratio(N_OB)', 'Type', 'Pearson');
% EI_r ,EI_p


%% cycle

Cycle_i(:)=F_motifs(NN_FLAG,7,:)+F_motifs(NN_FLAG,10,:)+F_motifs(NN_FLAG,12,:)+F_motifs(NN_FLAG,13,:);
CycleE_i(:)=F_motifs_E(NN_FLAG,7,:)+F_motifs_E(NN_FLAG,10,:)+F_motifs_E(NN_FLAG,12,:)+F_motifs_E(NN_FLAG,13,:);  
Cycle_d_i(:)=d_Fmotifs(NN_FLAG,7,:)+d_Fmotifs(NN_FLAG,10,:)+d_Fmotifs(NN_FLAG,12,:)+d_Fmotifs(NN_FLAG,13,:);
[Cycle_r(1), Cycle_p(1)]=corr(fir_rate(1:80)',Cycle_i(1:80)', 'Type', 'Pearson');
[Cycle_r(2), Cycle_p(2)]=corr(fir_rate(1:80)',CycleE_i', 'Type', 'Pearson');
[Cycle_r(3), Cycle_p(3)]=corr(fir_rate(1:80)',Cycle_d_i', 'Type', 'Pearson');

%% motifs
clear motifs_r motifs_p
motifs_i(:,:)=F_motifs(NN_FLAG,:,:); motifsE_i(:,:)=F_motifs_E(NN_FLAG,:,:); motifsd_i(:,:)=d_Fmotifs(NN_FLAG,:,:);

[motifs_r(1,:), motifs_p(1,:)]=corr(fir_rate(1:80)',motifs_i(:,1:80)', 'Type', 'Pearson');
[motifs_r(2,:), motifs_p(2,:)]=corr(fir_rate(1:80)',motifsE_i(:,1:80)', 'Type', 'Pearson');
[motifs_r(3,:), motifs_p(3,:)]=corr(fir_rate(1:80)',motifsd_i(:,1:80)', 'Type', 'Pearson');

EI_ratio=Dds_fromE(NN_FLAG,:)./(Dds_fromI(NN_FLAG,:)+Dds_fromE(NN_FLAG,:));
TatolDds_fromE=[TatolDds_fromE Dds_fromE(NN_FLAG,:)];
TatolDds_fromI=[TatolDds_fromI Dds_fromI(NN_FLAG,:)]; 
TatolEI_ratio=[TatolEI_ratio EI_ratio];
