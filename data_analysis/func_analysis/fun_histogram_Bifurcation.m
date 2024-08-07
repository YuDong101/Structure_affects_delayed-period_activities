function [out_x,out_y,output_histogram] = fun_histogram_Bifurcation(input_x,input_y,x_min,x_max,dx)

% input_x = P_Dds_EfromE; input_y = ssT_persistent;
% x_min=0; x_max=1; dx=0.1;
x=x_min:dx:x_max;
Nx_bin = round((x_max-x_min)/dx);

% y_min=10; Ny_bin = 50;
% y(1:Ny_bin)=0;
% for ii=1:Ny_bin
%     y(ii)=y_min*1.13^ii;
% end

y_min=-250; y_max=4500; dy=(y_max-y_min)/20;
y=y_min:dy:y_max;
Ny_bin = round((y_max-y_min)/dy);

num_histogram(1:Nx_bin-1,1:Ny_bin-1)=0;

for ii=1:Nx_bin-1
    input_yy = [];
    ccc = find(input_x>x(ii)&input_x<x(ii+1));
    for iii = 1:length(ccc)        
        input_yy(iii) = input_y(ccc(iii));
    end
    
    if isempty(input_yy)        
    else
        N_tot = length (input_yy);
        
        for ss = 1:N_tot
            for jj = 1:Ny_bin-1
                if input_yy(ss)>=y(jj) && input_yy(ss)<y(jj+1)
                    num_histogram(ii,jj) = num_histogram(ii,jj)+1;
                end
            end
        end
    end
end

for ss = 1:Nx_bin-1
    out_x(ss) = (x(ss)+x(ss+1))/2;
end

for ss = 1:Ny_bin-1
    out_y(ss) = (y(ss)+y(ss+1))/2;
end

% 列归一
v=sum(num_histogram,2);
D=diag(v);
norm_num_histogram=D^-1*num_histogram;

% 全归一
% D=sum(num_histogram,'all');
% norm_num_histogram=D^-1*num_histogram;

output_histogram=norm_num_histogram; %sum(num_histogram,'all')

% fp6 = fopen('三维 0.028 补充.dat','w');
% for ppp = 1:length(out_x)
%     for jjj= 1:length(out_y)
%         fprintf(fp6,'%f %f %f\n',out_x(ppp),out_y(jjj),output_histogram(ppp,jjj));
%     end
% end
% fclose (fp6);
% subplot(2,2,2),contour(out_x,out_y,output_histogram',8);

