
function [x,output] = fun_Nspike(input,dx,x_min,x_max)

% input = sum(matrixS,1); x_min=0; x_max=300; dx=1; %sum(matrixS,1)
% load histogram_par.mat;

[m,N_trail] = size (input);
x_range=x_max-x_min;
N_bin = round (x_range/dx);

num_histogram(1:N_bin)=0;

x=x_min:dx:x_max;

for ii = 1:N_bin  % 分成 N_range 个区间
    for jj = 1:dx
        for iii = 1:m
            if input(iii,(ii-1)*dx+jj)==1
                num_histogram(ii)=num_histogram(ii)+1; end
        end
    end
end

x=x(2:end)-dx/2;
output=num_histogram; %sum(num_histogram,'all')
