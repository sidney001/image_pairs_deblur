clear;
vsnr=cell(1,10);
vpsnr=cell(1,10);
kernel = [11,21,31,41,51,61,71,81,91,101];
sigma=100;
for i=1:10
    close all;
    [ vsnr{i} , vpsnr{i} ] = run_demo(sigma, kernel(i)+4, kernel(i) );
end


vsnr9=cell(1,10);
vpsnr9=cell(1,10);
sigma=90;
for i=1:10
    close all;
    [ vsnr9{i} , vpsnr9{i} ] = run_demo(sigma, kernel(i)+4, kernel(i) );
end


vsnr6=cell(1,10);
vpsnr6=cell(1,10);
sigma=60;
for i=1:10
    close all;
    [ vsnr6{i} , vpsnr6{i} ] = run_demo(sigma, kernel(i)+4, kernel(i) );
end