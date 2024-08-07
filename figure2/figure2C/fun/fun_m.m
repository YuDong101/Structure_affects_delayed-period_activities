
function dm = fun_m(v,m)
alpga_m=0.1*(v+16)./(1-exp(-(v+16)/10));
beta_m=4*exp(-(v+41)/18);
dm=alpga_m.*(1-m)-beta_m.*m;

