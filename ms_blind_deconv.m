function [kernel] = ms_blind_deconv(Nd,B, opts)

% 
% Do multi-scale blind deconvolution given input file name and options
% structure opts. Returns a double deblurred image along with estimated
% kernel. Following the kernel estimation, a non-blind deconvolution is run.
%
% Copyright (2011): Dilip Krishnan, Rob Fergus, New York University.
%

  %y = im2double(imread(fn));
  y = B;
  noise = Nd;

% prescale the image if it's too big; kernel size is defined for the SCALED
% image(如果图像太大就按照比例缩小)
for k = 1:size(y, 3)
  y1(:, :, k) = imresize(y(:, :, k), opts.prescale, 'bilinear');
  noise1(:, :, k) = imresize(noise(:, :, k), opts.prescale, 'bilinear');
end;
y = y1;
noise = noise1;
% save off for non-blind deconvolution
yorig = y;
noise_orig = noise;
% gamma correct（伽马校正）
y = y.^opts.gamma_correct;
noise = noise.^opts.gamma_correct;

% use a window to estimate kernel（用选择的窗口去估计图像，否则就是选择整个图像去估计，并转成灰度图）
if (~isempty(opts.kernel_est_win))
  w = opts.kernel_est_win;
  if (size(y, 3) == 3)
    y = rgb2gray(y([w(1):w(3)], [w(2):w(4)], :));
    noise = rgb2gray(noise([w(1):w(3)], [w(2):w(4)], :));
  end;
else
  if (size(y, 3) == 3)
    y = rgb2gray(y);
    noise = rgb2gray(noise);
  end;
end;

b = zeros(opts.kernel_size);%b=模糊核大小
bhs = floor(size(b, 1)/2);%bhs=模糊核一半大小

% set kernel size for coarsest level - must be odd（设置最粗层的核大小）
minsize = max(3, 2*floor(((opts.kernel_size - 1)/16)) + 1);
fprintf('Kernel size at coarsest level is %d\n', minsize);

% derivative filters
dx = [-1 1; 0 0];
dy = [-1 0; 1 0];

% l2 norm of gradient images(为啥等于6？？？？？？？？？？？？？？？？？？？）
l2norm = 6;
  
resize_step = sqrt(2);
% determine number of scales（设置模糊核的尺度，从小到大，依次递增）
num_scales = 1;
tmp = minsize;
while(tmp < opts.kernel_size)
  ksize(num_scales) = tmp;
  num_scales = num_scales + 1;
  tmp = ceil(tmp * resize_step);
  if (mod(tmp, 2) == 0) 
    tmp = tmp + 1;
  end;
end;
ksize(num_scales) = opts.kernel_size;
%例如 ksize = [3,5,9,13,19,27,31]

% blind deconvolution - multiscale processing（多尺度的盲估计过程）
for s = 1:num_scales
  if (s == 1)
    % at coarsest level, initialize kernel（最粗层的核初始化）
    ks{s} = init_kernel(ksize(1));
    k1 = ksize(1);
    k2 = k1; % always square kernel assumed
  else
    % upsample kernel from previous level to next finer level
    k1 = ksize(s);
    k2 = k1; % always square kernel assumed
    
    % resize kernel from previous level
    tmp = ks{s-1};
    tmp(tmp<0) = 0;
    tmp = tmp/sum(tmp(:));
    %将上一层的k（用tmp）表示，线性插值到[k1*k2]大小，比如上一层k大小为3*3
    %k1*k2=5*5,那么这一层就把3*3的核线性插值到5*5大小
    ks{s} = imresize(tmp, [k1 k2], 'bilinear');
    % bilinear interpolantion not guaranteed to sum to 1 - so renormalize
    % 线性插值不能保证k之和为1，所以要归一化
    ks{s}(ks{s} < 0) = 0;
    sumk = sum(ks{s}(:));
    ks{s} = ks{s}./sumk;
  end;
  
  % image size at this level（在模糊核的这层尺度上 的 原图像大小）
  r = floor(size(y, 1) * k1 / size(b, 1));%图像的行
  c = floor(size(y, 2) * k2 / size(b, 2));%图像的列
  
  if (s == num_scales)
    r = size(y, 1);
    c = size(y, 2);
  end;
  
  fprintf('Processing scale %d/%d; kernel size %dx%d; image size %dx%d\n', ...
            s, num_scales, k1, k2, r, c);
  
    
  % resize y according to the ratio of filter sizes（根据上一步重新调整图像大小）
  %把y调整到r*c大小，利用线性插值
  ys = imresize(y, [r c], 'bilinear');
  yx = conv2(ys, dx, 'valid'); %图像的横向梯度
  yy = conv2(ys, dy, 'valid'); %图像的竖直梯度
  
  c = min(size(yx, 2), size(yy, 2));%设置本层梯度图像的列
  r = min(size(yx, 1), size(yy, 1));%设置本层梯度图像的行
  
  g = [yx yy];%横向和竖直梯度的集合
  
  % normalize to have l2 norm of a certain size
  %目的应该是，把yx和yy转化为具有指定l2模大小的矩阵，比如之前设置的l2norm=6
  %下面的目的就是把yx和yy归一化为模大小为6的矩阵，为啥归一化为模为6而不是1？？？
  tmp1 = g(:, 1:c);%图像的横向梯度yx
  %在第一层，norm(tmp1(:))竟然等于6，l2norm/norm(tmp1(:))=1，为啥这么巧？？？
  tmp1 = tmp1*l2norm/norm(tmp1(:));%？？？？？？？？？？？？？？？？？？？
  g(:, 1:c) = tmp1;
  tmp1 = g(:, c+1:end);
  tmp1 = tmp1*l2norm/norm(tmp1(:));
  g(:, c+1:end) = tmp1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    % image size at this level（在模糊核的这层尺度上 的 原图像大小）
  r = floor(size(y, 1) * k1 / size(b, 1));%图像的行
  c = floor(size(y, 2) * k2 / size(b, 2));%图像的列
  
  if (s == num_scales)
    r = size(y, 1);
    c = size(y, 2);
  end;
% resize y according to the ratio of filter sizes（根据上一步重新调整图像大小）
  %把y调整到r*c大小，利用线性插值
  noise_s = imresize(noise, [r c], 'bilinear');
  noise_x = conv2(noise_s, dx, 'valid'); %图像的横向梯度
  noise_y = conv2(noise_s, dy, 'valid'); %图像的竖直梯度
  
  c = min(size(yx, 2), size(yy, 2));%设置本层梯度图像的列
  r = min(size(yx, 1), size(yy, 1));%设置本层梯度图像的行
 
  noise_grd = [noise_x noise_y];%横向和竖直梯度的集合

  % normalize to have l2 norm of a certain size
  %目的应该是，把yx和yy转化为具有指定l2模大小的矩阵，比如之前设置的l2norm=6
  %下面的目的就是把yx和yy归一化为模大小为6的矩阵，为啥归一化为模为6而不是1？？？
  tmp1 = noise_grd(:, 1:c);%图像的横向梯度yx
  %在第一层，norm(tmp1(:))竟然等于6，l2norm/norm(tmp1(:))=1，为啥这么巧？？？
  tmp1 = tmp1*l2norm/norm(tmp1(:));%？？？？？？？？？？？？？？？？？？？
  noise_grd(:, 1:c) = tmp1;
  tmp1 = noise_grd(:, c+1:end);
  tmp1 = tmp1*l2norm/norm(tmp1(:));
  noise_grd(:, c+1:end) = tmp1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ls{s} = noise_grd;
  
  % call kernel estimation for this scale（在本尺度上进行核估计）
  opts.lambda{s} = opts.min_lambda;
  
  [ls{s} ks{s} error_flag] = ss_blind_deconv(g, ls{s}, ks{s}, opts.lambda{s}, ...
                                               opts.delta, opts.x_in_iter, opts.x_out_iter, ...
                                               opts.xk_iter, opts.k_reg_wt);
 
  if (error_flag < 0)
    ks{s}(:) = 0;
    ks{s}(ceil(size(ks{s}, 1)/2), ceil(size(ks{s}, 2)/2)) = 1;
    fprintf('Bad error - just set output to delta kernel and return\n');
  end;
  
   % center the kernel（将模糊核移到中心位置，来缓解边界问题）
  c1 = (size(ls{s}, 2)) / 2;
  tmp1 = ls{s}(:, 1:c1); 
  tmp2 = ls{s}(:, c1 + 1 : end);
  [tmp1_shifted tmp2_shifted ks{s}] = center_kernel_separate(tmp1, tmp2, ks{s});
  ls{s} = [tmp1_shifted tmp2_shifted];
  
  % set elements below threshold to 0
  if (s == num_scales)
    kernel = ks{s};
    kernel(kernel(:) < opts.k_thresh * max(kernel(:))) = 0;
    kernel = kernel / sum(kernel(:));
  end;
end;

figure; imagesc(kernel); colormap gray; title('Kernel');
%
function [k] = init_kernel(minsize)
  k = zeros(minsize, minsize);
  k((minsize - 1)/2, (minsize - 1)/2:(minsize - 1)/2+1) = 1/2;
