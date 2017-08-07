function [ Igain ] = compute_Igain_map ( Nd, alpha )
%COMPUTE_IGAIN_MAP Compute gain map used in gain-controlled RL
%   
%   Nd      - The imaged as the source of gain map calculation
%   alpha   - The param determines the influence of gain map

lmax = 20;
gsize = [5 5];
sig = 0.5;


% g = fspecial('gaussian', [3,3], 1);
% Nd = imfilter(rgb2gray(Nd), g, 'conv');
[ h, w, d ] = size(Nd);

% Prepare Gaussian Pyramid
if d == 3
    Ndq = Gscale(rgb2gray(Nd), lmax, gsize, sig);
else
    Ndq = Gscale(Nd, lmax, gsize, sig);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sigama = [0.001, 0.005,0.01,0.05]
% for i = 1:size(Ndq,2)
%     figure;imshow(Ndq(i).img);title(['第' char(i) '层高斯金字塔']);
% end
% Preprocess the gradient of every point in the image
for k = 1:size(Ndq,2)
    [ img(k).Gx, img(k).Gy ] = imgradientxy(Ndq(k).img, 'CentralDifference');
    img(k).Gx = img(k).Gx/2^(k-1);
    img(k).Gy = img(k).Gy/2^(k-1);
    img_resize(k).Gx = imresize(img(k).Gx, [h w],'bilinear');
    img_resize(k).Gy = imresize(img(k).Gy, [h w],'bilinear');
end
% for i = 1:size(Ndq,2)
%     figure;imshow(img(i).Gx, [ ]);title(['第' char(i) '层高斯金字塔']);
% end
% for i = 1:size(Ndq,2)
%     figure;imshow(img_resize(i).Gx,[ ]);title(['第' char(i) '层高斯金字塔']);
% end
Igain = zeros([ h w ], 'double');

% Calculate Igain(i,j) for every single point in the image
for i = 1:h
    for j = 1:w
        Sum = 0;
        for k = 1:size(Ndq,2)
            %Sum = Sum + norm(get_gradient(img, , j, l));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             xx = max(fix(i/(2^l)),1);
%             yy = max(fix(j/(2^l)),1);
            Sum = Sum + sqrt((img_resize(k).Gx(i,j)) ^2 +(img_resize(k).Gy(i,j)) ^2);
        end
        
        Igain(i, j) = (1-alpha) + alpha * Sum;
    end
end

% normalize Igain
Igain = Igain / max(Igain(:));

end

