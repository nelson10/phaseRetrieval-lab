function [Uout] = ang_spec_multi_prop(Uin, wvl, delta1, deltan, z)
% split-step method for wave propagation from J. D. Schmidt (2010)
%
% Uin - input field
% wvl - wavelength
% delta1 - initial step size
% deltan - final step size
% z - propagation distance
% atmos - atmospheric model (optional)
%
% Uout - output field

%% Initial parameters
N = length(Uin);
%N = 8192;
%Uin=imresize(Uin,8192/512);
nx = linspace(-N / 2, N / 2 - 1, N);
[nx, ny] = meshgrid(nx, nx);
k = 2 * pi / wvl;
nsq = nx.^2 + ny.^2;
w = 0.47 * N;
%sg = exp(-nsq.^8 / w.^16);
sg =1;
clear nsq
n = length(z);
delta_z = z(2:end) - z(1:end - 1);


%% Propagation according to Chp. 8


alpha = z / z(end);
delta = (1 - alpha) * delta1 + alpha * deltan;
m = delta(2:end) ./ delta(1:end - 1);
x1 = nx * delta(1);
y1 = ny * delta(1);
r1sq = x1.^2 + y1.^2;

Q1 = exp(-1j*k/2*(1 - m(1))/delta_z(1).*r1sq);
Uin = Uin .* Q1;

for idx = 1:n - 1

    deltaf = 1 / (N * delta(idx));
    fX = nx * deltaf;
    fY = ny * deltaf;
    fsq = fX.^2 + fY.^2;

    Dz = delta_z(idx);

    Q2 = exp(1j * pi^2 * 2 * Dz ./ m(idx) / k .* fsq);

    Uin = sg.* ift2(Q2 .* ft2(Uin/m(idx), delta(idx)), deltaf);
    if idx == 39
        figure; imagesc(abs(Uin).^2)
    end
end

xn = nx * delta(end);
yn = ny * delta(end);

rnsq = xn.^2 + yn.^2;

Q3 = exp(-1j*k/2*(m(end-1) - 1)/(m(end-1)*Dz) .* rnsq);

Uout = Q3.*Uin;

end

function G = ft2(g, dx)

% Schmidt implementation of DFT2

G = fftshift(fft2(fftshift(g))) * dx^2;

end

function g = ift2(G, df)

% Schmidt implementation of DIFT2

N = length(G);
g = ifftshift(ifft2(ifftshift(G))) * (N * df)^2;

end
