
function [pxx,fx]=fun_FFT_periodogram(par,re,minfreq)
  
  
  
  %power for the selected layer:
  restate=re(1,round((par.transient+par.dt)/par.dt):end);
  [pxx1,fx1]=periodogram(restate',[],[],1./(par.dt));
  %mean(restate)
  
  bin=5;%compress the data (pick up one in every 'bin' points):
  remaining=mod(length(pxx1),bin);
  pxx0=pxx1(1:end-remaining,:);fx0=fx1(1:end-remaining,:);
  pxx=binTraces(pxx0,bin)';fx=binTraces(fx0,bin)';
  

  save pard.mat


  %%