function [out_x,ht,hp]=htp(flage_var,at,ap,N_bin)

x_min=min(flage_var);x_max=max(flage_var);dx=(x_max-x_min)/N_bin;
[~,ht] = fun_histogram_one_par(flage_var(at),x_min,x_max,dx);
[out_x,hp] = fun_histogram_one_par(flage_var(ap),x_min,x_max,dx);
[out_x,hp] = fun_histogram_one_par(flage_var(ap),x_min,x_max,dx);