clc
clear all
close all
t=0:0.005:10;
%original signal
y=(exp(-0.2.*t).*cos((6.*pi.*t)+(0.5.*sin(6.*pi.*t)))+0.5.*sin(20.*pi.*t)); 
%plotting original signal
figure(1)
plot(t,y);


%genetic algorithm for parameter estimation in eemd
yn=transpose(y);                     
aim=7;                     % numbers of IMF
NR = 100;                      % value of ensemble
Nstd = 5;                   % param to white noise
[imf]=eemd(yn,aim,NR,Nstd);


a=zeros(2001,7);
a=imf(1:2001,3:7);
a=transpose(a);
b=transpose(imf(1:2001,1:2));
imf1=zeros(1,2001);
imf2=zeros(1,2001);
imf3=zeros(1,2001);
imf4=zeros(1,2001);
imf5=zeros(1,2001);
imf6=zeros(1,2001);
imf7=zeros(1,2001);
imf1=b(1,1:2001);
imf2=b(2,1:2001);
imf3=a(1,1:2001);
imf4=a(2,1:2001);
imf5=a(3,1:2001);
imf6=a(4,1:2001);
imf7=a(5,1:2001);
imf_final=imf3+imf4+imf5+imf6+imf7;

%plotting results for eemd
figure(2)
plot(t,y,'ob',t,imf_final,'+r');

%wavelet thresholding 
x1=imf1;
x2=imf2;
x=x1+x2;
XDEN=wdenoise(x,3);

%reconstruction using ssa

original_signal=imf_final+XDEN;

figure(3)

plot(t,y,'ob',t,original_signal,'+r');
[final]=ssa(original_signal,1);
final_result=zeros(1,2001);
final_result=transpose(final);
snr_ssa=snr(y,y-final_result)
mse_ssa=immse(final_result,y)
snr_eemd_wavelet=snr(y,y-original_signal)
snr_eemd=snr(y,y-imf_final)
mse_eemd_wavelet = immse(original_signal,y)
mse_eemd=immse(imf_final,y)

