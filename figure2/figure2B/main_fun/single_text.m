
function [time_fir,n_fir,output_spike_R] = single_text(time,outputR,Nt,NR_node)

output_spike_R(NR_node,Nt) = 0; % time_fir=0; n_fir=0;
ss=1;
for iii  = 1:Nt
    for nnn = 1:NR_node
        if  iii>1
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
