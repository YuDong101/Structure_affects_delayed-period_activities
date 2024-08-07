
function [time_fir,n_fir] = single_text(time,outputS,outputR,Nt,NR_node)

ss=1;
for iii  = 1:Nt
    for nnn = 1:NR_node
        if  iii>1
%             if (outputS(nnn,iii)>0.0 && outputS(nnn,iii-1)<=0.0)
%                 output_spike_S(nnn,iii) = 1;
% %                 fprintf(fp1,'%f %f\n',time(iii),nnn);
%             end
            if (outputR(nnn,iii)>0.0 && outputR(nnn,iii-1)<=0.0)
                output_spike_R(nnn,iii) = 1;
                time_fir(ss)=time(iii);
                n_fir(ss)=nnn;
                ss=ss+1;
%                 fprintf(fp3,'%f %f\n',time(iii),nnn);
            end
        end
    end
end

% figure (11)
% G = digraph(aaa_flag_matrix_R);
% plot(G)
% [x_E_flag,N_loop_E_flag,x_I_flag,N_loop_I_flag,Cc_avg_flag,Lens_avg_flag...
%     ] = matrix_analysis(aaa_flag_matrix,K_R);
% figure (12)
% subplot (3,1,1), plot (x_flag,N_loop_flag);
% subplot (3,1,2), plot (x_E_flag,N_loop_E_flag);
% subplot (3,1,3), plot (x_I_flag,N_loop_I_flag);

%  disp(['平均路径长度为：',num2str(Lens_avg)]);
% matrix1(:,:)=matrix_rec_R(kkk,ppp,1,:,:);
% matrix2(:,:)=matrix_rec_R(kkk,ppp,2,:,:);
% figure (11)
% G = digraph(matrix1);
% plot(G)
% figure (12)
% G = digraph(matrix2);
% plot(G)