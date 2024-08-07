
function dn = fun_h(v,h)
alpga_h=0.07*exp(-(v+30)/20);
beta_h=1./(1+exp(-v/10));
dn=alpga_h.*(1-h)-beta_h.*h;