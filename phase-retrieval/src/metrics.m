function [p,s] = metrics(data,Xrec,kdataset) %1 for cave dataset
[M,N,L] = size(Xrec);
p = zeros(L,1);
s = zeros(L,1);
for i=1:L
    if(kdataset==1 || kdataset==2)
        ref = data(:,:,i);
        A = Xrec(:,:,i);
    end
    p(i) = psnr(A,ref);
    s(i) = ssim(A,ref);
end
p = round(mean(p),2); % psnr
s= round(mean(s),2); % ssim
%disp("PSNR "+ num2str(p)+ " SSIM "+num2str(s));
%disp("---------------------------------------------------------------------------------------------------------")
end