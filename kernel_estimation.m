function [ks] = kernel_estimation( I, B, opts)
%KERNEL_ESTIMATION Get an estimated kernel based on benoised and blurred image
%   Solve a l1-regularized least squares problem to get the kernel
%
%   I       - is actually Nd as approxiamation //denoised
%   B       - blurred image
%   len     - kernel size; default: 64
%   lambda  _ estimation parameter mentioned in the paper: default: 25
%   method  - 
%   verbose - true to enable debug output
%
%   including a toolbox from Standford to the l1-ls problem (in 
%       l1_ls_matlab directory, https://github.com/cvxgrp/l1_ls)
%
disp('Starting kernal estimation...');
b = zeros(opts.kernel_size);
%set kernel size for coarsest level - must be odd
%minsize = max(3, 2*floor(((opts.kernel_size - 1)/16)) + 1);
%fprintf('Kernel size at coarsest level is %d\n', maxitr);
%%
ret = sqrt(0.75);    %0.8660
%% 计算需要的层数
maxitr=max(floor(log(5/min(opts.kernel_size))/log(ret)),0);
num_scales = maxitr + 1;
fprintf('Maximum iteration level is %d\n', num_scales);
%%
retv=ret.^[0:maxitr];                %对y的缩放比例
k1list=ceil(opts.kernel_size*retv);
k1list=k1list+(mod(k1list,2)==0);    %kernel的长
k2list=ceil(opts.kernel_size*retv);  %kernel的宽
k2list=k2list+(mod(k2list,2)==0);

% blind deconvolution - multiscale processing
for s = num_scales:-1:1
  if (s == num_scales)
      %%
      % at coarsest level, initialize kernel
      ks = init_kernel(k1list(s));
      k1 = k1list(s);
      k2 = k1; %always square kernel assumed 因为是方形的kernel
  else
    % upsample kernel from previous level to next finer level
    k1 = k1list(s);
    k2 = k1; % always square kernel assumed
    
    % resize kernel from previous level
    ks = resizeKer(ks,1/ret,k1list(s),k2list(s));
    % 对模糊核进行收缩，减小噪声！
    ks(ks < max(ks(:))*opts.M_low) = 0;
    ks = ks / sum(ks(:));   
  end;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%
  cret=retv(s);
  ys=downSmpImC(I,cret);    %根据缩放比例 对y进行下采样
  yB=downSmpImC(B,cret);    %根据缩放比例 对y进行下采样

  fprintf('Processing scale %d/%d; kernel size %dx%d; image size %dx%d\n', ...
            s, num_scales, k1, k2, size(ys,1), size(ys,2));
  %-----------------------------------------------------------%
  %% kernel estimation
%     K = mat2vec(ks);
%     A = mat2kmat(ys, k1);
%     At = A';
%     b = mat2vec(yB);
%     AtA = At * A;
%     Atb = At * b;
%     lambda2 = opts.lambda^2;
%     E = eye(k1^2);
%     
%     for i = 1:opts.niter
%         K = K + opts.beta * (Atb - (AtA + lambda2*E)*K);
%         K(K<0) = 0;
%         K = K / sum(K(:));
%     end
%     ks = reshape(K, [k1 k2]);
% 
%     % normalize estimated kernel
%     ks = ks / sum(ks(:));
%%

  [ks] = blind_deconv_main(ys,yB, ks,opts);
% 对模糊核进行收缩，减小噪声！
   ks(ks < max(ks(:))*opts.M_high) = 0;
   ks = ks / sum(ks(:));   

%% center the kernel
   ks = adjust_psf_center(ks);
   ks(ks(:)<0) = 0;
   sumk = sum(ks(:));
   ks = ks./sumk;
  %% set elements below threshold to 0
  if (s == 1)
    kernel = ks;
    if opts.k_thresh>0
        kernel(kernel(:) < max(kernel(:))/opts.k_thresh) = 0;
    else
        kernel(kernel(:) < 0) = 0;
    end
    ks = kernel / sum(kernel(:));
  end;
end;
figure;imshow(ks, [ ]);
disp('Kernel estimation complete...');
%% end kernel estimation
end


%% Sub-function
function [k] = init_kernel(minsize)
  k = zeros(minsize, minsize);
  k((minsize - 1)/2, (minsize - 1)/2:(minsize - 1)/2+1) = 1/2;
end

%%
function sI=downSmpImC(I,ret)
%% refer to Levin's code
if (ret==1)
    sI=I;
    return
end
%%%%%%%%%%%%%%%%%%%

sig=1/pi*ret;

g0=[-50:50]*2*pi;
sf=exp(-0.5*g0.^2*sig^2);
sf=sf/sum(sf);
csf=cumsum(sf);
csf=min(csf,csf(end:-1:1));
ii=find(csf>0.05);

sf=sf(ii);
sum(sf);

I=conv2(sf,sf',I,'valid');

[gx,gy]=meshgrid([1:1/ret:size(I,2)],[1:1/ret:size(I,1)]);

sI=interp2(I,gx,gy,'bilinear');
end
%%
function k=resizeKer(k,ret,k1,k2)
%%
% levin's code
k=imresize(k,ret);
k=max(k,0);
k=fixsize(k,k1,k2);
if max(k(:))>0
    k=k/sum(k(:));
end
end
%% 
function nf=fixsize(f,nk1,nk2)
[k1,k2]=size(f);

while((k1~=nk1)|(k2~=nk2))
    
    if (k1>nk1)
        s=sum(f,2);
        if (s(1)<s(end))
            f=f(2:end,:);
        else
            f=f(1:end-1,:);
        end
    end
    
    if (k1<nk1)
        s=sum(f,2);
        if (s(1)<s(end))
            tf=zeros(k1+1,size(f,2));
            tf(1:k1,:)=f;
            f=tf;
        else
            tf=zeros(k1+1,size(f,2));
            tf(2:k1+1,:)=f;
            f=tf;
        end
    end
    
    if (k2>nk2)
        s=sum(f,1);
        if (s(1)<s(end))
            f=f(:,2:end);
        else
            f=f(:,1:end-1);
        end
    end
    
    if (k2<nk2)
        s=sum(f,1);
        if (s(1)<s(end))
            tf=zeros(size(f,1),k2+1);
            tf(:,1:k2)=f;
            f=tf;
        else
            tf=zeros(size(f,1),k2+1);
            tf(:,2:k2+1)=f;
            f=tf;
        end
    end
    
    
    
    [k1,k2]=size(f);
    
end

nf=f;
end