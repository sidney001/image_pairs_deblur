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
% image(���ͼ��̫��Ͱ��ձ�����С)
for k = 1:size(y, 3)
  y1(:, :, k) = imresize(y(:, :, k), opts.prescale, 'bilinear');
  noise1(:, :, k) = imresize(noise(:, :, k), opts.prescale, 'bilinear');
end;
y = y1;
noise = noise1;
% save off for non-blind deconvolution
yorig = y;
noise_orig = noise;
% gamma correct��٤��У����
y = y.^opts.gamma_correct;
noise = noise.^opts.gamma_correct;

% use a window to estimate kernel����ѡ��Ĵ���ȥ����ͼ�񣬷������ѡ������ͼ��ȥ���ƣ���ת�ɻҶ�ͼ��
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

b = zeros(opts.kernel_size);%b=ģ���˴�С
bhs = floor(size(b, 1)/2);%bhs=ģ����һ���С

% set kernel size for coarsest level - must be odd��������ֲ�ĺ˴�С��
minsize = max(3, 2*floor(((opts.kernel_size - 1)/16)) + 1);
fprintf('Kernel size at coarsest level is %d\n', minsize);

% derivative filters
dx = [-1 1; 0 0];
dy = [-1 0; 1 0];

% l2 norm of gradient images(Ϊɶ����6����������������������������������������
l2norm = 6;
  
resize_step = sqrt(2);
% determine number of scales������ģ���˵ĳ߶ȣ���С�������ε�����
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
%���� ksize = [3,5,9,13,19,27,31]

% blind deconvolution - multiscale processing����߶ȵ�ä���ƹ��̣�
for s = 1:num_scales
  if (s == 1)
    % at coarsest level, initialize kernel����ֲ�ĺ˳�ʼ����
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
    %����һ���k����tmp����ʾ�����Բ�ֵ��[k1*k2]��С��������һ��k��СΪ3*3
    %k1*k2=5*5,��ô��һ��Ͱ�3*3�ĺ����Բ�ֵ��5*5��С
    ks{s} = imresize(tmp, [k1 k2], 'bilinear');
    % bilinear interpolantion not guaranteed to sum to 1 - so renormalize
    % ���Բ�ֵ���ܱ�֤k֮��Ϊ1������Ҫ��һ��
    ks{s}(ks{s} < 0) = 0;
    sumk = sum(ks{s}(:));
    ks{s} = ks{s}./sumk;
  end;
  
  % image size at this level����ģ���˵����߶��� �� ԭͼ���С��
  r = floor(size(y, 1) * k1 / size(b, 1));%ͼ�����
  c = floor(size(y, 2) * k2 / size(b, 2));%ͼ�����
  
  if (s == num_scales)
    r = size(y, 1);
    c = size(y, 2);
  end;
  
  fprintf('Processing scale %d/%d; kernel size %dx%d; image size %dx%d\n', ...
            s, num_scales, k1, k2, r, c);
  
    
  % resize y according to the ratio of filter sizes��������һ�����µ���ͼ���С��
  %��y������r*c��С���������Բ�ֵ
  ys = imresize(y, [r c], 'bilinear');
  yx = conv2(ys, dx, 'valid'); %ͼ��ĺ����ݶ�
  yy = conv2(ys, dy, 'valid'); %ͼ�����ֱ�ݶ�
  
  c = min(size(yx, 2), size(yy, 2));%���ñ����ݶ�ͼ�����
  r = min(size(yx, 1), size(yy, 1));%���ñ����ݶ�ͼ�����
  
  g = [yx yy];%�������ֱ�ݶȵļ���
  
  % normalize to have l2 norm of a certain size
  %Ŀ��Ӧ���ǣ���yx��yyת��Ϊ����ָ��l2ģ��С�ľ��󣬱���֮ǰ���õ�l2norm=6
  %�����Ŀ�ľ��ǰ�yx��yy��һ��Ϊģ��СΪ6�ľ���Ϊɶ��һ��ΪģΪ6������1������
  tmp1 = g(:, 1:c);%ͼ��ĺ����ݶ�yx
  %�ڵ�һ�㣬norm(tmp1(:))��Ȼ����6��l2norm/norm(tmp1(:))=1��Ϊɶ��ô�ɣ�����
  tmp1 = tmp1*l2norm/norm(tmp1(:));%��������������������������������������
  g(:, 1:c) = tmp1;
  tmp1 = g(:, c+1:end);
  tmp1 = tmp1*l2norm/norm(tmp1(:));
  g(:, c+1:end) = tmp1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    % image size at this level����ģ���˵����߶��� �� ԭͼ���С��
  r = floor(size(y, 1) * k1 / size(b, 1));%ͼ�����
  c = floor(size(y, 2) * k2 / size(b, 2));%ͼ�����
  
  if (s == num_scales)
    r = size(y, 1);
    c = size(y, 2);
  end;
% resize y according to the ratio of filter sizes��������һ�����µ���ͼ���С��
  %��y������r*c��С���������Բ�ֵ
  noise_s = imresize(noise, [r c], 'bilinear');
  noise_x = conv2(noise_s, dx, 'valid'); %ͼ��ĺ����ݶ�
  noise_y = conv2(noise_s, dy, 'valid'); %ͼ�����ֱ�ݶ�
  
  c = min(size(yx, 2), size(yy, 2));%���ñ����ݶ�ͼ�����
  r = min(size(yx, 1), size(yy, 1));%���ñ����ݶ�ͼ�����
 
  noise_grd = [noise_x noise_y];%�������ֱ�ݶȵļ���

  % normalize to have l2 norm of a certain size
  %Ŀ��Ӧ���ǣ���yx��yyת��Ϊ����ָ��l2ģ��С�ľ��󣬱���֮ǰ���õ�l2norm=6
  %�����Ŀ�ľ��ǰ�yx��yy��һ��Ϊģ��СΪ6�ľ���Ϊɶ��һ��ΪģΪ6������1������
  tmp1 = noise_grd(:, 1:c);%ͼ��ĺ����ݶ�yx
  %�ڵ�һ�㣬norm(tmp1(:))��Ȼ����6��l2norm/norm(tmp1(:))=1��Ϊɶ��ô�ɣ�����
  tmp1 = tmp1*l2norm/norm(tmp1(:));%��������������������������������������
  noise_grd(:, 1:c) = tmp1;
  tmp1 = noise_grd(:, c+1:end);
  tmp1 = tmp1*l2norm/norm(tmp1(:));
  noise_grd(:, c+1:end) = tmp1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ls{s} = noise_grd;
  
  % call kernel estimation for this scale���ڱ��߶��Ͻ��к˹��ƣ�
  opts.lambda{s} = opts.min_lambda;
  
  [ls{s} ks{s} error_flag] = ss_blind_deconv(g, ls{s}, ks{s}, opts.lambda{s}, ...
                                               opts.delta, opts.x_in_iter, opts.x_out_iter, ...
                                               opts.xk_iter, opts.k_reg_wt);
 
  if (error_flag < 0)
    ks{s}(:) = 0;
    ks{s}(ceil(size(ks{s}, 1)/2), ceil(size(ks{s}, 2)/2)) = 1;
    fprintf('Bad error - just set output to delta kernel and return\n');
  end;
  
   % center the kernel����ģ�����Ƶ�����λ�ã�������߽����⣩
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
