%RUN_DEMO run through this project. 
%input image is set in `generate_init_data` file
clear all;close all
%function [ vsnr vpsnr ] = run_demo(sigma, kernelSize, lenth )
sigma = 50, length = 40.1;%%43.6
addpath('./l1_ls_matlab');
addpath('./cho_code');
verbose = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% if (abs(length - 43.6))<= 2e-6 %==0
%     kernelSize = 41;
% else
%     kernelSize = 45;
% end
%%
%%
if mod((ceil(length)),2) == 0
    kernelSize = ceil(length) + 1;
else
    kernelSize = ceil(length);
end
%%
opts.kernel_size = kernelSize;  saturation = 0;
lambda_pixel = 4e-3; lambda_grad = 4e-3; opts.gamma_correct = 2.2;
lambda_tv = 0.002; lamda_l0 = 2e-4; weight_ring = 1;
opts.lambda = 5;
opts.beta = 1.0;
opts.niter = 20;
opts.k_thresh = 9;%值越大，收缩模糊核能力越弱，只在最后一步收缩！
opts.xk_iter = 1;
opts.M_high = 0.1;
opts.M_low = 0.05;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  generate blur image
len = length;
theta = 45;
y = im2double(imread('./image_map/map.png')); %skd_map 768x1024x3 double
k = fspecial('motion', len, 180); % used for the example in report

if mod((ceil(len)),2) == 0
    len = ceil(len) + 1;
else
    len = ceil(len);
end

Ktmp = zeros(len);
Ktmp(int8(len)/2,:) = k;
k = Ktmp;
%H = fspecial('motion',LEN,THETA)为运动模糊算子，有两个参数，表示摄像物体逆时针方向以theta角度运动了len个像素，len的默认值为9，theta的默认值为0；
imageBlur = imfilter(y, k);%skd_map 768x1024x3 double
imageBlur = imageBlur(51:end-50,251:end-50,:);%skd_map_clip 668x724x3 double
y = y(51:end-50,251:end-50,:);% do need ,:

[psnr_blur, snr_blur] = psnr(imageBlur,y)
%ssim_noise=ssim(imageBlur,y)
%skd_map_clip psnr_blur =28.8827nsnr_blur =19.9948
% figure,imshow(imageBlur);
% title('imageBlur');
%% generate noise image
addpath('./BM3D');
randn('seed', 0);
%sigma = 75;
z = y + (sigma/255)*randn(size(y));
%z = im2double(imread('noise.jpg'));
% tic;
[NA, y_est] = CBM3D(1, z, sigma);
% elapsed_time = toc %
 [psnr_noise, snr_noise] = psnr(z,y)
[psnr_denoise, snr_denoise] = psnr(y_est,y)
%ssim_noise=ssim(y_est,y)
%ssim_denoise=ssim(y_est,y)

% show the noisy image 'z' and the denoised 'y_est'
%figure; imshow(y);  
% figure; 
% subplot(121);imshow(z);  
% subplot(122);imshow(y_est);
% title('noisy image and denoised image');
%% 
imageNoise = 'noise_enhenceRegi';
director = ['./image_map/' 'map_' num2str(len) '_' num2str(sigma) '_kernelsize' num2str(kernelSize)];%such as './image_map/map_5_5'
%s = num2str(A)converts a numeric array into a character array that represents the numbers.
mkdir(director);
B = imageBlur;
Nd = y_est; % denoised
% %%% 加入一个双边滤波器，对噪声图像再滤波
% % B = B(200:600,200:700,:);
% % Nd = Nd(200:600,200:700,:);
% w     = 5;       % bilateral filter half-width
% sigma = [3 0.1]; % bilateral filter standard deviations
% Nd_before = Nd;
% Nd = bfilter2( Nd_before, w, sigma );
% figure;subplot(1,2,2);imshow(Nd_before);title('仅仅去噪的图像')
% subplot(1,2,1);imshow(Nd);title('加双边滤波器的图像')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matlab builtin regularized filter
% [Ireg ,K1] = deblur( Nd, B, true, 'reg', verbose,opts );
%[K] = ms_blind_deconv(Nd,B, opts)

% %% get 41 pixel's kernel
% yy = im2double(imread('./image_map/map.png')); %skd_map 768x1024x3 double
% kk = fspecial('motion', len-3, 180); % used for the example in report
% %H = fspecial('motion',LEN,THETA)为运动模糊算子，有两个参数，表示摄像物体逆时针方向以theta角度运动了len个像素，len的默认值为9，theta的默认值为0；
% imageBlurr = imfilter(yy, kk);%skd_map 768x1024x3 double
% imageBlurr = imageBlurr(51:end-50,251:end-50,:);%skd_map_clip 668x724x3 double
% yy = yy(51:end-50,251:end-50,:);% do need ,:
% randn('seed', 0);
% %sigma = 75;
% z = y + (sigma/255)*randn(size(y));
% %z = im2double(imread('noise.jpg'));
% [NAA, y_estt] = CBM3D(1, z, sigma);
% BB = imageBlurr;
% Ndd = y_estt;
% [Ireg111 ,K] = deblur( Ndd, BB, true, 'reg', verbose,opts );
%% 得出与blur不同的另外一个blur 不需要estimate kernel
temp = fspecial('motion', kernelSize, 180);
K = zeros(kernelSize);
K(int8(kernelSize)/2,:) = temp;


%%
write_kernel(K, [director '/map' ]);
%imwrite(Ireg, 'images/deblurred_reg.jpg');

% call deconv directly since we already have K

% matlab builtin Richardson-Lucy
Ilucy = deconv(Nd, B, K, 'lucy', verbose);
imwrite(Ilucy, [director '/'  'map_deblur_RL.png']);

% residual RL algorithm
Ir = deconv(Nd, B, K, 'resRL', verbose);
imwrite(Ir, [director '/'  'map_deblur_residual_RL.png']);

% gain-controlled RL
Ig = deconv(Nd, B, K, 'gcRL', verbose);
imwrite(Ig, [director '/'  'map_deblur_gain_controlled_RL.png']);

% gain-controlled Rl with detail
% no need to call deconv again, because we already have Ir and Ig
% Idetailed = deconv(double(Nd), double(B), double(K), 'detailedRL', verbose);
disp('Calculating Ibar with joint/cross bilateral filtering...');
Ibar = zeros(size(Ir), 'double');
[~, ~, d] = size(Ir);
for i = 1:d
    % joint/cross bilateral filtering
    Ibar(:,:,i) = jbfilter2(Ir(:,:,i), Ig(:,:,i), 2, [1.6, 0.08]);
end
Id = Ir - Ibar; % detail layer Id
imwrite(Id+0.8, [director '/'  'map_detail_layer.png']);
Idetailed = Ig + Id;    % final result
imwrite(Idetailed,[director '/'  'map_deblurred_detail_RL.png']);
%% calculate SNR
[psnr_res_RL, snr_res_RL] = psnr(Ir,y)
[psnr_gain_RL, snr_gain_RL] = psnr(Ig,y)
[psnr_detail_RL, snr_detail_RL] = psnr(Idetailed,y)

imwrite(z, [director '/'  'map_noise.png']);
imwrite(y_est, [director '/'  'map_denoise.png']);
imwrite(imageBlur, [director '/'  'map_blur.png']);

vsnr(6)=0;
vsnr(1)=snr_noise;
vsnr(2)=snr_blur;
vsnr(3)=snr_denoise;
vsnr(4)=snr_res_RL;
vsnr(5)=snr_gain_RL;
vsnr(6)=snr_detail_RL;

vpsnr(6)=0;
vpsnr(1)=psnr_noise;
vpsnr(2)=psnr_blur;
vpsnr(3)=psnr_denoise;
vpsnr(4)=psnr_res_RL;
vpsnr(5)=psnr_gain_RL;
vpsnr(6)=psnr_detail_RL;

disp('ALL CLEAR.')