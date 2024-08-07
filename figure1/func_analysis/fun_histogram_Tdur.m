
function [x_out,output] = fun_histogram_Tdur(T_persistent)
% clear; clc
% load 100次.mat

input = T_persistent; x_min=10; N_bin = 50;

num_histogram(1:N_bin-1)=0;

x(1:N_bin)=0;
for ii=1:N_bin
    x(ii)=x_min*1.15^ii;
end

N_tot = size (input,2);

for ss = 1:N_tot
    for jj = 1:N_bin-1
        if input(ss)>=x(jj) && input(ss)<x(jj+1)
            num_histogram(jj) = num_histogram(jj)+1;
        end
    end
end
for ss = 1:N_bin-1
    x_out(ss) = (x(ss)+x(ss+1))/2;
end
output=num_histogram; %sum(num_histogram,'all')

% save histogram_par.mat;

% subplot (1,1,1), bar (log10(x),output);

% [x2,y2,aaa,bbb] = Gaussian_fit(num_histogram/sum(matrixS,'all'),x,dx);
% 
% figure (2)
% plot(x2,y2,'ro');
% hold on
% plot(aaa,bbb,'b');

