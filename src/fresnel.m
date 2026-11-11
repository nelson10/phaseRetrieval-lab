%## Propagation Function
%# obj: Object with internal values.
%# opc: Forward or backward propagation.
%# dim: Problem dimension.

function out = fresnel(obj,wvl,delta1,z)
z_ = z;
pixelSize = delta1;
dim = size(obj,1);
NN=dim*dim;
N = sqrt(NN);
y1=0-(dim/2)+1; x1 = y1;
y2=dim-(dim/2); x2 = y2;
XX = x1:x2; YY = y1:y2;
[Fx,Fy]=meshgrid(XX,YY);
%#Campo lejano
% k=(obj.pixelSize^2)/(obj.wavelength*z_);
% propagation=exp(-1i*2*pi*z_/obj.wavelength)*exp(-1i*pi*k*(Fy.^2+Fx.^2));
% g=1;
%#Campo cercano
k=(wvl*z_)/((pixelSize^2)*NN);
g=exp(1i*pi*z_/wvl);
propagation=exp(1i*pi*k*(Fy.^2+Fx.^2));
R=g.*ifft2(fftshift(propagation.*fftshift(fft2(obj)/N))/N);
%R=g.*ifft2(propagation.*fft2(obj.complexObject)/N)/N;
out = R;
end