function [Lens,Lens_avg,flag_connect]=func_Path_Length(matrix)
 flag_connect=1;
 Num = size(matrix,2);
 Lens   = matrix;
 Lens(find(Lens==0))=inf;    %将邻接矩阵变为邻接距离矩阵，两点无边相连时赋值为inf，自身到自身的距离为0.
 for i=1:Num           
     Lens(i,i)=0;       
 end   
 for k=1:Num            %Floyd算法求解任意两点的最短距离
     for i=1:Num
         for j=1:Num
             if Lens(i,j)>Lens(i,k)+Lens(k,j)
                Lens(i,j)=Lens(i,k)+Lens(k,j);
             end
         end
     end
 end
 Lens_avg=sum(sum(Lens))/(Num*(Num-1));  %平均路径长度
 if Lens_avg==inf
     disp('该网络图不是连通图');
     flag_connect = 0;
 end
         
 
 
 