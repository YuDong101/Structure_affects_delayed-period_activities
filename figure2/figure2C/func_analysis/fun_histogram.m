
function [x,output] = fun_histogram(N_loop)

input_s = N_loop; x_min=-0.5; x_max=10.5; dx=1; %sum(matrixS,1)
% load histogram_par.mat;

[m,N] = size (input_s);
ss=1;
for iii = 1:m
    for jjj = 1:N
        input(ss) = input_s(iii,jjj);
        ss=ss+1;
    end
end

x_range=x_max-x_min;
N_bin = round (x_range/dx);

num_histogram(1:N_bin)=0;

x=x_min:dx:x_max;
N_tot = size (input,2);

for ss = 1:N_tot  % 分成 N_range 个区间
    for jj = 1:N_bin
        if input(ss)==jj
            num_histogram(jj) = num_histogram(jj)+1;
        end
    end
end

x=x(2:end)-dx/2;
output=num_histogram; %sum(num_histogram,'all')

% save histogram_par.mat;

% subplot (1,1,1), plot (x,output);

% [x2,y2,aaa,bbb] = Gaussian_fit(num_histogram/sum(matrixS,'all'),x,dx);
% 
% figure (2)
% plot(x2,y2,'ro');
% hold on
% plot(aaa,bbb,'b');

