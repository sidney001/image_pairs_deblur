function [ Igain ] = compute_Igain_map ( Nd, alpha )
%COMPUTE_IGAIN_MAP Compute gain map used in gain-controlled RL
%   
%   Nd      - The imaged as the source of gain map calculation
%   alpha   - The param determines the influence of gain map

lmax = 1;
gsize = [5 5];
sig = 0.5;


% g = fspecial('gaussian', [3,3], 1);
% Nd = imfilter(rgb2gray(Nd), g, 'conv');
[ h, w, d ] = size(Nd);

% Prepare Gaussian Pyramid
if d == 3
    Ndl = Gscale(rgb2gray(Nd), lmax, gsize, sig);
else
    Ndl = Gscale(Nd, lmax, gsize, sig);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sigama = [0.001, 0.005,0.01,0.05]

% Preprocess the gradient of every point in the image
for l = 1:lmax
    [ img(l).Gx, img(l).Gy ] = imgradientxy(Ndl(l).img, 'CentralDifference');
%     g = fspecial('gaussian', [5,5], 1);
%     img(l).Gx = imfilter(img(l).Gx, g, 'conv');
%     img(l).Gy = imfilter(img(l).Gy, g, 'conv');

end

Igain = zeros([ h w ], 'double');

% Calculate Igain(i,j) for every single point in the image
for i = 1:h
    for j = 1:w
        Sum = 0;
        for l = 1:lmax
            %Sum = Sum + norm(get_gradient(img, , j, l));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             xx = max(fix(i/(2^l)),1);
%             yy = max(fix(j/(2^l)),1);
            Sum = Sum + sqrt((img(l).Gx(i,j)) ^2 +(img(l).Gy(i,j)) ^2);
        end
        
        Igain(i, j) = 0.8 + alpha * Sum*80;
    end
end

% normalize Igain
Igain = Igain / max(Igain(:));

end

function [ gx, gy ] = get_gradient ( img, x, y, l )
% calculate the corresponding point for (x,y) in level l of Gaussian
% Pyramid

gx = img(l).Gx(max(fix(x/(2^l)),1));
gy = img(l).Gy(max(fix(y/(2^l)),1));

end

