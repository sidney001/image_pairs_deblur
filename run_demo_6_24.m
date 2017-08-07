%RUN_DEMO run through this project. 
%   input image is set in `generate_init_data` file
clear;close all;
addpath('./l1_ls_matlab');
addpath('./cho_code');
verbose = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.kernel_size = 39;  saturation = 0;
lambda_pixel = 4e-3; lambda_grad = 4e-3; opts.gamma_correct = 2.2;
lambda_tv = 0.002; lambda_l0 = 2e-4; weight_ring = 1;
opts.lambda = 5;
opts.beta = 1.0;
opts.niter = 20;
opts.k_thresh = 9;%值越大，收缩模糊核能力越弱，只在最后一步收缩！
opts.xk_iter = 1;
opts.M_high = 0.1;
opts.M_low = 0.05;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('loading...');
% B = im2double(imread('images/blurred.jpg'));
% Nd = im2double(imread('images/denoised.jpg'));
imageBlur = 'IMG_3515_Regi_SURF';
sub_director = 'campus/3515/'; % 'map/' or ''
imageNoise = 'IMG_3514_Regi_SURF';
director = ['./images/' sub_director imageBlur];
mkdir(director);
B = im2double(imread(['images/' sub_director imageBlur '.png']));
Nd = im2double(imread(['images/' sub_director imageNoise '.png']));
%%% 加入一个双边滤波器，对噪声图像再滤波
% B = B(200:600,200:700,:);
% Nd = Nd(200:600,200:700,:);
w     = 2;       % bilateral filter half-width     w     = 3;  
sigma = [1.5 0.05]; % bilateral filter standard deviations    sigma = [2 0.08];
Nd_before = Nd;
Nd = bfilter2( Nd_before, w, sigma );
figure;subplot(1,2,2);imshow(Nd_before);title('仅仅去噪的图像')
subplot(1,2,1);imshow(Nd);title('加双边滤波器的图像')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matlab builtin regularized filter
[Ireg,K] = deblur( Nd, B, true, 'reg', verbose,opts );
%[K] = ms_blind_deconv(Nd,B, opts)


write_kernel(K, [director '/' imageBlur]);
%imwrite(Ireg, 'images/deblurred_reg.jpg');

% call deconv directly since we already have K

% matlab builtin Richardson-Lucy
Ilucy = deconv(Nd, B, K, 'lucy', verbose);
imwrite(Ilucy, [director '/' imageBlur '_deblur_lucy.png']);

% residual RL algorithm
Ir = deconv(Nd, B, K, 'resRL', verbose);
imwrite(Ir, [director '/' imageBlur '_deblur_residual_RL.png']);

% gain-controlled RL
Ig = deconv(Nd, B, K, 'gcRL', verbose);
imwrite(Ig, [director '/' imageBlur '_deblur_gain_controlled_RL.png']);

% gain-controlled Rl with detail
% no need to call deconv again, because we already have Ir and Ig
% Idetailed = deconv(double(Nd), double(B), double(K), 'detailedRL', verbose);
disp('Calculating Ibar with joint/cross bilateral filtering...');
Ibar = zeros(size(Ir), 'double');
[~, ~, d] = size(Ir);
for i = 1:d
    % joint/cross bilateral filtering
    Ibar(:,:,i) = jbfilter2(Nd(:,:,i), Ig(:,:,i), 2, [1.6, 0.08]);
end
Id = Nd - Ibar; % detail layer Id
imwrite(Id+0.2, [director '/' imageBlur '_detail_layer.png']);
Idetailed = Ig + Id*1.5;    % final result
imwrite(Idetailed,[director '/' imageBlur '_deblurred_detail_RL.png']);


disp('ALL CLEAR.')