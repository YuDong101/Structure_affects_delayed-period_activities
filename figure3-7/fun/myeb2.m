function h=myeb2(x,ymean,ysigma,color)


revx=fliplr(x);
yminus=ymean-ysigma;
yplus=fliplr(ymean+ysigma);
fillerror=fill([x revx],[yminus yplus],0.25*color+[0.75 0.75 0.75]);
hold on;set(fillerror,'EdgeColor','None');
h=plot(x,ymean,'Color',color,'LineWidth',3);
set(gca, 'Layer','top');z=0.01;
hold off;
h=fillerror;



