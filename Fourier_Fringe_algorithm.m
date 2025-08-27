clc
close all
clear all 

lamda=632.8e-9*1000;%wavelength (mm)
f=500;%focal length

dx=5e-3;%mm
N=256;
x=[-N/2:N/2-1]*dx;
y=x;
[X,Y]=meshgrid(x,y);

a=1/(2*dx);
df=2*a/N;
fx=[-a:2*a/N:a-2*a/N];
fy=[-a:2*a/N:a-2*a/N];
[Fx,Fy]=meshgrid(fx,fy);

temp_square=(Fx).^2+(Fy).^2; %X^2+Y^2 of object
temp=1-(Fx/f).^2-(Fy/f).^2;
temp(find(temp<0))=0;
kz=2*pi/lamda*sqrt(temp);

alpha=1.1;%Sagnac scaling                               
alpha_=(alpha-1/alpha); %=(alpha-1/alpha)^-1

%% Generate complex object with amplitude and phase
load('object/2GT_phase.mat')
load('object/2GT_abs.mat')

abs_fieldU=(abs_fieldU+1)/2*255;%convert to 8-bit value 
ang_fieldU=(ang_fieldU+1)/2*pi;%convert to 0-pi value 
field=abs_fieldU.*exp(1i*ang_fieldU);
field(1:100,1:256)=0; % eliminate the conjugate object.

figure
mesh(angle(field))
%% at the hologram plane
% complex amplitude of the off-axis object
Fresnel_term=exp(1i*pi/lamda/f*temp_square);
UOB_z=Fresnel_term.*(fftshift(fft2(fftshift((field)))))*dx*dx;

% Reference wave
Ur=1*exp(1i*pi/lamda/f*temp_square);

Holo=Ur+UOB_z;%interference
Intensity_holo=(abs(Holo).^2);%Recorded hologram
figure
imagesc(fx,fy,Intensity_holo)
axis off
title('Hologram')

%% The coherence function 
Coherence_f=fftshift(ifft2(fftshift(Intensity_holo)))/dx/dx;
%make small numbers to be zeros
Coherence_f(abs(Coherence_f)<1e-11)=0;
Coherence_f_temp=Coherence_f;

Coherence_f=fftshift(ifft2(fftshift(Intensity_holo)))/dx/dx;

figure
imagesc((abs(Coherence_f)))
%% image captured by CCD with the introduction of carrier frequency
off_ang=2.15;
R=sum(sum(Intensity_holo))*(df*f*lamda)*(df*f*lamda); %first term of eq (7)
Re=Coherence_f; %Coherence function in term of delta x, delta y, delta z from eq(6)
Re_=Re;%introduce carrier frequency 
I_true=2*R+2*real(Re_.*exp(1i*2*pi*1/lamda*sind(off_ang)*(X+Y))); %captured image 
num=500000;
I=add_noise_CCD(I_true,num,5,14);
% imagesc(log(I))%using log() to rescale the image for better visualization
im(I)
%%
% Fourier Fringe analysis 
I_FFT=fftshift(fft2(fftshift(I)))*dx*dx;  %FFT of the CCD image

%simulate Notch filter
xsize=100; %size of the filter is 256 by 256 pixels
ysize=100;
filter=zeros(N);
% the position of the carrier frequency in pixel
pixel=round((max(-fx)-(1/lamda*sind(off_ang)))/df)+1;
tic
filtered_sig=I_FFT(pixel-xsize/2+1:pixel+xsize/2,pixel-ysize/2+1:pixel+ysize/2);%simulate Notch filter

%one padding the filtered signal becasue the amplitude of the spectrum of
%the coherence function is always larger than 1 becasue of the interfernce
%of the reference wave and the object wave 
In_recon=ones(N); 
center=floor(N/2)+1;
In_recon(center-xsize/2:center+xsize/2-1,center-ysize/2:center+ysize/2-1)=filtered_sig-1; % one padding

%make small numbers to be zeros
recon_real_img=fftshift(ifft2(fftshift((In_recon))))/dx/dx;
recon_real_img(find(abs(imag(recon_real_img))<10e-0))=real(recon_real_img(find(abs(imag(recon_real_img))<10e-0))); %zero imagine part
recon_real_img(find(abs(recon_real_img)<10e-0))=0;%zero real part 
toc
im(abs(recon_real_img))
title('amplitude of the reconstructed coherence function')
axis off 
im(angle(recon_real_img)) 
colormap('gray')
title('phase of the reconstructed coherence function')
axis off
figure
imagesc(angle(Coherence_f_temp))
colormap('gray')
title('phase of the coherence function')
axis off