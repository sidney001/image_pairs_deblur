clear,close all;
name1 = '8363'
f = imread([name1 '.jpg']);
df = imresize(f, 0.20, 'bilinear');
figure,imshow( df );
imwrite(df,[name1 '.png']);

% %% denoise
% addpath('./denoise', './denoise/l1_ls_matlab');
% 
% g =ice('image',df);
% 
% df_denoise = denoise(im2double(g));
% imwrite(df_denoise,'noise_enhence_denoise.png');