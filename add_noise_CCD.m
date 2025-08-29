function [noisy_image]=add_noise_CCD(Input_img,fw,bit_level,rand_noise)
%% Define FullWell and sensor Resolution
%Input_img: input non-noise image
% fw: FullWell;
% bit_level: sensor Resolution;
%rand_noise: std of experiment noise
% reference: Quang Trieu and George Nehmetallah,"Deep Learning Based Coherence Holography (DCH) Reconstruction of 3D Objects "
%           at section 5. Noise Simulation mimicking experimental data
rng(1)
photo_electrons_number= graylevel2photoelectrons(Input_img,fw,bit_level );
new_photo_electrons_number=poissrnd(photo_electrons_number);
noisy_image = photoelectrons2graylevel(new_photo_electrons_number,fw,bit_level);
noisy_image=noisy_image+rand_noise.*randn(size(noisy_image));
