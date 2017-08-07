function [ I, K ] = deblur( Nd, B, unikernel, deconvmode, verbose,opts )
%DEBLUR Deblur image B
%
%   Nd          - Denoised image used by Yuan et al.
%   B           - Blurred imaged to deblur
%   unikernel   - true to estimate a uniform kernel for all channels
%   deconvmode  - deconvolution algorithm. it can be 'reg', 'lucy',
%                 'resRL', 'gcRL', 'detailedRL'
%   verbose     - display extra info if true
%

%nargin为“number of input arguments”的缩写。 在matlab中定义一个函数时， 在函数体内部，
%nargin是用来判断输入变量个数的函数。
if nargin < 5
    verbose = true;
end
if nargin < 4
    deconvmode = 'gcRL';
end
if nargin < 3
    unikernel = true;
end

% isa(obj,ClassName);Determine if input is object of specified class
if isa(B, 'uint8')
    I = zeros(size(B), 'uint8');
elseif isa(B, 'uint16')
    I = zeros(size(B), 'uint16');
elseif isa(B, 'double')
    I = zeros(size(B), 'double');
end

[~, ~, d] = size(B); %skd 668x724x3 double d=3

if unikernel
    if d == 3 %convert RGB to gray 
        gNd = rgb2gray(Nd);
        gB = rgb2gray(B);
    else
        gNd = Nd;
        gB = B;
    end
    
    % kernel estimation
    K = kernel_estimation(gNd, gB,opts);
    
    % deconvolution
    I = deconv(Nd, B,K, deconvmode, verbose);
else
    % calculate kernel for each channel
    K = zeros([21 21 d], 'double');
    
    % loop for every channel, same as unikernel if d == 1
    for i = 1:d
        K(:,:,i) = kernel_estimation(Nd(:,:,i), B(:,:,i), 21, 5, 'l1ls', verbose);
        I(:,:,i) = deconv(double(Nd), double(B(:,:,i)), double(K(:,:,i)), deconvmode, verbose);
    end
end

end
