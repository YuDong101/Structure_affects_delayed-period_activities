
function [dv,Energy_ions] = fun_v(v,m,h,n,Isny,Iext,Cm)
gNa=12.5; gK=4.74; gL=0.025;
ENa=40; EK=-80; EL=-65;
iNa=gNa.*m.^3.*h.*(v-ENa);
iK=gK.*n.^4.*(v-EK);
iL=gL.*(v-EL);
dv=(-iNa-iK-iL+Isny+Iext)./Cm;

Energy_ions=iNa.*(v-ENa)+iK.*(v-EK)+iL.*(v-EL);



% 
%         AAA(iS,iii)=Cm;
%         bbb(iS,iii)=iS;
% 
% %     dv = v + step*fun_v(v,m,h,n,Isny,Iext,Cm);
% %     dm = m + step*fun_m(v,m);
% %     dh = h + step*fun_h(v,h);
% %     dn = n + step*fun_n(v,n);
% save fun_data.mat