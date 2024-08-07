function [out_x,output_histogram] = fun_histogram_one_par(input_x,x_min,x_max,dx)
% clear out_x output_histogram num_histogram

% input_x = TatolDuration;
% x_min=min(input_x);x_max=max(input_x);dx=(x_max-x_min)/20;
x=x_min-dx/2:dx:x_max+dx/2;
Nx_bin = round((x_max-x_min)/dx);

num_histogram(1:Nx_bin-1)=0;

for ii=1:Nx_bin+1
    ccc = find(input_x>x(ii)&input_x<=x(ii+1));
    num_histogram(ii)=length(ccc);
    clear ccc
end

for ss = 1:Nx_bin+1
    out_x(ss) = (x(ss)+x(ss+1))/2;
end

output_histogram=num_histogram; %sum(num_histogram,'all')
% figure (25)
% bar (out_x,output_histogram);
