function [x, k, error_flag] = ss_blind_deconv(y, x, k, lambda, delta, ...
                                              x_in_iter, x_out_iter, ...
                                              xk_iter, k_reg_wt)

% --�����Ż������еķ���2
% --y��ʾԭ��ģ��ͼ��ĺ����������ݶȼ��ϣ�y������һֱ���䣬ֻ�Ǿ�����ͬ��С���ز�������ɲ�ͬ�Ĵ�С����
% --x��ʾ��һ��õ��Ĺ�������ͼ���ݶȣ�Ȼ�󾭹��ϲ���������ͼ��Ĵ�С�����Ҿ���
%   ģ2��һ��Ϊ6�������������Ĺ��ƹ��̣������xΪ������º����������ͼ��x����
%   ��Ϊ��һ�ε����Ļ�ͼ��Ȼ�����ϲ�������6��������ѭ��������
% --kҲ�����Ƶģ�k��ʾ��һ����Ƶõ���k��Ȼ�󾭹��ϲ���������ͼ��Ĵ�С������ģ��һ��Ϊ1
%   �����������Ĺ��ƣ������kΪ������º��ģ���ˣ�����Ϊ��һ�ε����Ļ��ˣ�Ȼ���ϲ�������һ��������ѭ��������

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
totiter = 1;%��������
error_flag = 0;

% Split y into 2 parts: x and y gradients; handle independently
% throughout����y�ֳ������֣��ֱ���
y1{1} = y(:, 1 : m2);%�����ݶ�
y1{2} = y(:, m2 + 1 : end);%��ֱ�ݶ�
y2{1} = y1{1}(khs + 1 : end - khs, khs + 1 : end - khs);%�����Ϊ�˺�k����󱣳�ͼ���С���䣬����Ҫ��ǰ��Сһ���
y2{2} = y1{2}(khs + 1 : end - khs, khs + 1 : end - khs);

x1{1} = x(:, 1 : m2);
x1{2} = x(:, m2 + 1 : end);

lambda_orig = lambda;
delta_orig = delta;
%�ֱ����x��k������ֱ��126��
for iter = 1:xk_iter

  %����k��������
 
  % set up options for the kernel estimation���˹��ƣ�
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
