%% data sets
function [x,y,aaa,bbb] = Gaussian_fit(num_histogram,px,dx)

y=num_histogram;
x=px(2:end)-dx/2;

y1=fliplr(y)/fliplr(max(y));
x1=fliplr(x)/fliplr(max(x));
ymax=max(y);

%% Calculate coefficient
fun=fittype('A*exp(-(x-mu)^2/(2*sigma^2))');

[cf,gof]=fit(x1(:),y1(:),fun,'Start',[])

%% Interpolate the data

k=1:length(x1);
ki=linspace(1,length(x1),2000);
xi=interp1(k,x1,ki,'linear');
yi=interp1(k,y1,ki,'linear');
Yi=cf.A*exp(-(xi-cf.mu).^2/(2*cf.sigma^2))*ymax;

%% plot data
% plot(x,y,'ro');
% hold on
% plot(xi*max(x),Yi,'b');
aaa=xi*max(x); bbb= Yi;