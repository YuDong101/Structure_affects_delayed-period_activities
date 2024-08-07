
function [x,output] = fun_Nspike(input,dx,x_min,x_max)
N_ave=4; % N_ave 倍的 dx 为积分窗口

[m,N_trail] = size (input(1:80,:));
x_range=x_max-x_min-N_ave*dx;
N_bin = round (x_range/dx);

num_histogram(1:N_bin)=0;

x=x_min:dx:x_max-N_ave*dx;

for ii = 1:N_bin  % 分成 N_range 个区间
    for jj = 1:dx*N_ave
        for iii = 1:m
            if input(iii,(ii-1)*dx+jj)==1
                num_histogram(ii)=num_histogram(ii)+1; end
        end
    end
end

x=x(2:end)+N_ave*dx/2;
output=num_histogram/N_ave; %sum(num_histogram,'all')

% save histogram_par.mat;

% subplot (1,1,1), plot (x,output);

% [x2,y2,aaa,bbb] = Gaussian_fit(num_histogram/sum(matrixS,'all'),x,dx);
% 
% figure (2)
% plot(x2,y2,'ro');
% hold on
% plot(aaa,bbb,'b');

