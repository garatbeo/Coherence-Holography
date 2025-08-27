function [noisy_image]=add_noise_CCD(Input_img,fw,bit_level,rand_noise)
%% Define FullWell and sensor Resolution
% fw = 2000000;
% bit_level = 14;
%% Load input image and allocate output image
% im = double((imread('capitol.jpg')));
rng(1)
photo_electrons_number= graylevel2photoelectrons(Input_img,fw,bit_level );
new_photo_electrons_number=poissrnd(photo_electrons_number);
noisy_image = photoelectrons2graylevel(new_photo_electrons_number,fw,bit_level);
noisy_image=noisy_image+rand_noise.*randn(size(noisy_image));
% temp=noisy_image-I_true;
% snr=10*log10(mean(I_true(:))/var(temp(:)))
% figure
% subplot(1,2,1)
% imagesc(noisy_image)
% colormap('gray')
% subplot(1,2,2)
% imagesc(I_true)
% colormap('gray')
% %%
% figure; imshow(uint16(noisy_image)); colormap gray; axis equal
% title('Left - clean image ; Right - Noisy Image')