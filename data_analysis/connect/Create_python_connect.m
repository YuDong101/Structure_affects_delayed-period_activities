% function [matrix] = Create_python_connect(matrix)

G = digraph(matrix);

CCC=G.Edges(:,1);

aaa=allcycles(G);

% save Endnodes.txt -ascii CCC
% for kkk = 1:100
%     End_nodes(kkk) = CCC(kkk,1);
% end

% fp1=fopen('Endnodes.txt','w');
% fprintf(fp1,'%d\n',CCC);