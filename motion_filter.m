% function  [ht] = motion_filter( length , ang )
%
% Returns a filter "ht", with a length "length" and an angle "ang",
% calling funciton 搃mrotate?with bicubic interpolation.
%
% Mariana S. C. Almeida
% Instituto de Telecomunica珲es, Lisbon, Portugal 
% marianascalmeida@gmail.com
 
function  [ht] = motion_filter(length, ang)
 
%保证模糊核为奇数    
if rem(ceil(length),2) 
    length2 = ceil(length); 
else
    length2 = ceil(length)+1; 
end
 
ht = zeros(length2);
 
if length <= 1
    linha=1;
    ht(ceil(length2/2),:)=linha;
else
    naux = ceil(length/2)-1;
    linha1 = repmat(1,1,naux);
    linha2 = (length-1-(2*naux))/2;
    if linha2>0
        linha = [linha2 linha1 1 linha1 linha2];
    else
        linha = [linha1 1 linha1];
        linha(1)=linha(1)+linha2;
        linha(end)=linha(end)+linha2;
    end
    ht(ceil(length2/2),:)=linha;
end
 
ht = imrotate(ht, ang,'bicubic','crop');
ht(find(ht<0))=0;
ht=ht/sum(sum(ht));
%画出来
hf = ht/max(ht(:));
figure,imshow(hf);