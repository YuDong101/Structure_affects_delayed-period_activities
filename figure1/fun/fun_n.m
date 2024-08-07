
function dn = fun_n(v,n)
alpga_n=0.01*(v+20)./(1-exp(-(v+20)/10));
beta_n=0.125*exp(-(v+30)/80);
dn=alpga_n.*(1-n)-beta_n.*n;