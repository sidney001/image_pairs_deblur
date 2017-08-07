function resim = gain_residual_Lucy(Nd,delta_B, k, iterations,alpha)
%Function to restore the image using Lucy-Richardson
%Inputs: ifbl, LEN, THETA, iterations.
%Returns: resim
%
%ifbl:  It is the input image.
%THETA: It is the blur angle. The angle at which the image is blurred.
%LEN:   It is the blur length. The length is the number
%       of pixels by which the image is blurred.
%iterations: It is the number of iterations.
%handle:It is the handle to the waitbar(progress bar).
%resim: It is the restored image.
%
%Example:
%       resim = Lucy(image, LEN, THETA, iterations);
%       This call takes image, blur length, blur angle & no. of iterations 
%       as input and returns the restored image

%Preprocessing
%Performing Median Filter before restoring the blurred image
% delta_B = medfilt2(abs(delta_B));
offset =1;

Igain = compute_Igain_map(Nd, alpha);
figure;imshow(Igain,[ ]);title('gain map');

%Initialising the initial estimate to the blurred image
delta_I = delta_B;

%Create PSF of degradation
PSF = k;

%Convert psf to otf of desired size
%OTF is Optical Transfer Function
OTF = psf2otf(PSF,size(delta_B));

i = 1;
while i<=iterations
    %Converting the estimate to frequency domain
    delta_I = delta_I+offset;
    delta_I = max(delta_I, 0);
    fest = fft2(delta_I);
    
    %Multiplying OTF with the estimate in frequency domain
    fblest = OTF.*fest;
    
    %Converting the blurred image estimate to spatial domain
    ifblest = ifft2(fblest);
    
    %Calculating ratio of blurred image and estimate of the deblurred image
    iratio = (delta_B+offset)./ifblest;
    
    %Converting the ratio to frequency domain
    firatio = fft2(iratio);
    
    %Calculating the correction vector
    corrvec = conj(OTF) .* firatio;
    
    %Converting the correction vector to spatial domain
    icorrvec = ifft2(corrvec);
    
    %Multiplying correction vector & estimate of deblurred image to find next estimate
    aftercorr = icorrvec.*(delta_I+offset);
    delta_I = aftercorr - offset;
    for j=1:size(delta_I,3)
       % if i ~= iterations
            delta_I(:,:,j) = Igain .* delta_I(:,:,j);
     %   end
    end
    delta_I = min(delta_I, 1);
    %Setting the waitbar indicating how much is completed
   % waitbar(i/iterations, handle);
    i = i+1;
end

%Restored image
resim = abs(delta_I);