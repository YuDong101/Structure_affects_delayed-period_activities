
function [f,psdx] = fun_FFT(input,dt)

XX_FFT=input';

[N,m]=size(XX_FFT);

fs=1/dt;
n=0:N-1;   
y=fft(XX_FFT,N);    %对信号进行快速Fourier变换
% mag=abs(y').^2/N^2;     %求得Fourier变换后的振幅
% f=n*fs/N;

y = y(1:round(N/2),:);

psdx = (1/(fs*N)).*abs(y).^2;

psdx(2:end-1) = 2*psdx(2:end-1);

psdx=psdx';

f = 0:fs/length(XX_FFT):fs/2;


save fun_FFT

