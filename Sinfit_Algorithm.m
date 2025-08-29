%referenece:
%(1) Coherence holography by achromatic 3-D field correlation of generic thermal light with an imaging Sagnac shearing interferometer 
%       Dinesh N. Naik,Takahiro Ezawa,Rakesh Kumar Singh,Yoko Miyamoto,and Mitsuo Takeda
% (2)Properties of the IEEE-STD-1057 four parameter sine wave fit algorithm
%                       Peter H¨andel, Senior Member, IEEE
clc
close all
clear all 


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

temp_square=(Fx*lamda).^2-(Fy*lamda).^2; %X^2+Y^2 of object
temp=1-temp_square;
temp(find(temp<0))=0;
kz=2*pi/lamda.*sqrt(temp);

alpha=1.1;%Sagnac scaling                               
alpha_=(alpha-1/alpha); %=(alpha-1/alpha)^-1
n=4; %number of CCD images to reconstruct a subobject
w=2*pi/n; %angular velocity of phase shift for each CCD image 
% deltaz_a=0;


%% Collect the CCD images to reconstruct the sub-object
clc
in=1;%index of the subobject to reconstruct 
clear I
z=deltaz_a(in);
% z=z_rec;
for m=1:n
    G_m=G*exp(1i*w*(m));%equa (2) adding phase shift to object spectrum
    
    H=abs(G)+1/2*G.*exp(1i*w*(m))+1/2*conj(G).*exp(-1i*w*(m));%equation (3) phase-shifted holograms
    background=sum(sum((H)))*(df)*(df); %Equation (6) 
    Cohere=fftshift(ifft2(fftshift(H.*exp(1i*kz*z))))/dx/dx;
    I_true=(2*background+2*real(Cohere));
    input_irrad_photons=I_true;
    num=500000;
    I{m}=add_noise_CCD(I_true,num,5,14);%captured image 
    temp=I{m}-I_true;
    SNR_=10*log10(mean(I_true(:))/var(temp(:)))
end
%%
%% display one of the CCD images and holograms
close all
figure %Hologram figure 10
imagesc(H)
colormap('hot')

axis off
 figure
imagesc(I{1}) %first CCD image
colormap('gray')
axis off
%% (2)Properties of the IEEE-STD-1057 four parameter sine wave fit algorithm
Num=length(I);
t=linspace(1,Num,Num);
w=2*pi/Num;
tic
for i=1:Num
yy(i,:,:)=I{i}; % equation (3)
D(i,:)=[cos(w*t(i)) sin(w*t(i)) 1];  % equation (5)
end
L=inv(D'*D)*D'; %(D'D)^-1D'
X=pagetimes(L,yy);% equation (6) MATLAB function will be pagemtimes
for i=1:length(X) 
    AA(:,i)=X(1,:,i)';
    BB(:,i)=X(2,:,i)';
    CC(:,i)=X(3,:,i)';
end
A_new=AA;
A_new(find(abs(AA)<10))=0;% small values are assigned to 0

B_new=BB;
B_new(find(abs(BB)<10))=0;% small values are assigned to 0

amp=sqrt(A_new.^2+B_new.^2); 
phase=atan2(B_new,A_new);

toc
figure %figure 12, 15, 18
imagesc(x/alpha_,y/alpha,amp) 
colormap(gray(256))
title('amplitude of the reconstructed subobject')
axis off
figure %figure 13, 16, 19
imagesc(x/alpha_,y/alpha,phase)

colormap(gray(256))
axis off
title('Phase of the reconstructed object')