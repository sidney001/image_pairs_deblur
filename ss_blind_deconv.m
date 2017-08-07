function [x, k, error_flag] = ss_blind_deconv(y, x, k, lambda, delta, ...
                                              x_in_iter, x_out_iter, ...
                                              xk_iter, k_reg_wt)

% --就是优化论文中的方程2
% --y表示原来模糊图像的横和竖方向的梯度集合，y的内容一直不变，只是经过不同大小的重采样，变成不同的大小而已
% --x表示上一层得到的估计清晰图像梯度，然后经过上采样到本层图像的大小，并且经过
%   模2归一化为6，经过本函数的估计过程，输出的x为本层更新后的清晰估计图像x，并
%   作为下一次迭代的基图像，然后再上采样，归6化，依次循环。。。
% --k也是类似的，k表示上一层估计得到的k，然后经过上采样到本层图像的大小，并且模归一化为1
%   经过本函数的估计，输出的k为本层更新后的模糊核，并作为下一次迭代的基核，然后上采样，归一化，依次循环。。。

% Do single-scale blind deconvolution using the input initializations
% x and k. The cost function being minimized is: min_{x,k}
% \lambda/2 |y - x \oplus k|^2 + |x|_1/|x|_2 + k_reg_wt*|k|_1
%

k_init = k;
khs = floor(size(k, 1) / 2);

[m, n] = size(y);
[k1, k2] = size(k);
m2 = n / 2;

% arrays to hold costs
lcost = [];
pcost = [];
totiter = 1;%迭代次数
error_flag = 0;

% Split y into 2 parts: x and y gradients; handle independently
% throughout（把y分成两部分，分别处理）
y1{1} = y(:, 1 : m2);%横向梯度
y1{2} = y(:, m2 + 1 : end);%竖直梯度
y2{1} = y1{1}(khs + 1 : end - khs, khs + 1 : end - khs);%大概是为了和k卷积后保持图像大小不变，所以要提前缩小一点点
y2{2} = y1{2}(khs + 1 : end - khs, khs + 1 : end - khs);

x1{1} = x(:, 1 : m2);
x1{2} = x(:, m2 + 1 : end);

lambda_orig = lambda;
delta_orig = delta;
%分别更新x和k！！！直到126行
for iter = 1:xk_iter

  %更新k！！！！
 
  % set up options for the kernel estimation（核估计）
  opts.use_fft = 1;
  opts.lambda = k_reg_wt;
  opts.pcg_tol=1e-4;
  opts.pcg_its = 1;
    
  k_prev = k;
  k = pcg_kernel_irls_conv(k_prev, x1, y2, opts); % using conv2's
  k(k < 0) = 0;
  sumk = sum(k(:));
  k = k ./ sumk;
end;

% combine back into output
x(:, 1 : m2) = x1{1};
x(:, m2 + 1 : end) = x1{2};
