function [Lens,Lens_avg]=func_Path_Length(matrix)
 
 Num = size(matrix,2);
 Lens   = matrix;
 Lens(find(Lens==0))=inf;    %���ڽӾ����Ϊ�ڽӾ�����������ޱ�����ʱ��ֵΪinf����������ľ���Ϊ0.
 for i=1:Num           
     Lens(i,i)=0;       
 end   
 for k=1:Num            %Floyd�㷨��������������̾���
     for i=1:Num
         for j=1:Num
             if Lens(i,j)>Lens(i,k)+Lens(k,j)
                Lens(i,j)=Lens(i,k)+Lens(k,j);
             end
         end
     end
 end
 Lens_avg=sum(sum(Lens))/(Num*(Num-1));  %ƽ��·������
 if Lens_avg==inf
     disp('������ͼ������ͨͼ');
 end
         
 
 
 