clear;
clc;
close all;

%% Path
addpath(genpath('src'));
addpath(genpath('mask'));
kdataset = 1; % Kodak dataset

%load('kodakInput.mat')
%% Parameters
nm = 16; % number of masks
codes  = 1; % 1 Random, 2 Sphere Packing, 3 DCT, 4 Hadamard
complementary = 1; % 0 No complementary, 1 Complementary
ndataset = 24; % number of images Kodak
typePropagation  = 1; % 0  one step, 1 two steps, 2 Fresnel
iter = 10; % Number of iterations
wvl = 643e-6; % wavelength
deltan = 0.02; % DMD pixel size
delta1 = 0.02; % camera pixel size
rng(4);
result = zeros(ndataset,4);
CA = ["Random","Sphere Packing","DCT","Hadamard"];
propagator=["One-step","Multi-step","Fresnel"];
comple = ["No complementary","Complementary"];

%% Image size
if(codes ==6)
    M = 512;
    N = 512;
elseif(codes == 0 | codes ==1 |codes ==2 | codes==3 | codes==4)
    M = 384;
    N = 384;
end
dataset = zeros(M,N,ndataset);

%% Initializations
Set2med = zeros(M,N,nm);
fieldvec = ones(M,N);

%% Load Masks
[masks] = selectMasks(codes,M,N,nm,complementary);

[conditionNumber] = ComputeConditionNumber(masks);
%% Load Intensity and Phase
load('kodakInput.mat')

for k = 1:ndataset
    data = mat2gray(double(rgb2gray(inputKodak{k})));
    [d1,d2] = size(data);
    r1 = (d1/2)+(M/2);
    r2 = (d1/2)-(M/2)+1;
    c1 = (d2/2)+(M/2);
    c2 = (d2/2)-(M/2)+1;
    dataset(:,:,k) = data(r2:r1,c2:c1);
end

%% Propagation Distance
if(typePropagation==0 | typePropagation==2)
    z = 155; % propagation distance
else
    z = linspace(0,160,5); % propagation distance
end
 fig = figure();
for k = 1:ndataset
    intensity = dataset(:,:,25-k);
    phase = dataset(:,:,k);
    propagateObj = zeros(M,N,nm);

    for i=1:nm
        Obj1 = (masks(:,:,i).*intensity).*exp(1i.*phase);
        [OP] = propagation(typePropagation,Obj1,wvl, delta1, deltan, z);
        Obj(:,:,i) = abs(OP).^2;
    end

    tic
    for j = 1:iter
        for i = 1:nm
            [fieldvec,Obf2F,Phase] = PhaseLiftSimu(fieldvec,masks(:,:,i),typePropagation,Obj(:,:,i),wvl, delta1, deltan, z);
        end

        % intensityHat = uint8(mat2gray(abs(Obf2F)).*255);
        %phaseHat = uint8(mat2gray(angle(Obf2F)).*255);
        intensityHat = mat2gray(abs(Obf2F));
        phaseHat = mat2gray(angle(Obf2F));
        %uint8(intensity.*255)
        %uint8(phase.*255)
        [p1,s1] = metrics(intensity,intensityHat,kdataset);
        [p2,s2] = metrics(phase,phaseHat,kdataset);
        disp("iteration "+num2str(j) +" out of " + iter);

        %subplot(1,2,1),imagesc(intensityHat),title("PSNR= "+num2str(p1)+" SSIM= "+num2str(s1)),colormap('gray'),pbaspect([1 1 1]),axis off
        %subplot(1,2,2),imagesc(rescale(phaseHat,-pi,pi)),title("PSNR= "+num2str(p2)+" SSIM= "+num2str(s2)),colormap('gray'),pbaspect([1 1 1]),axis off;
        %montage({intensityHat,phaseHat});
        %pause(0.1)
        tlo = tiledlayout(fig,1,2,'TileSpacing','None');
        ax = nexttile(tlo);
        imshow(intensityHat,'Parent',ax),title("Intensity PSNR= "+num2str(p1)+" SSIM= "+num2str(s1)),colormap('gray'),pbaspect([1 1 1]),axis off
        ax = nexttile(tlo);
        imshow(phaseHat,'Parent',ax),title("Phase PSNR= "+num2str(p2)+" SSIM= "+num2str(s2)),colormap('gray'),pbaspect([1 1 1]),axis off;
        pause(0.01);
    end
    result(k,:) = [p1 s1 p2 s2];
    disp("--------------------------------------------------------------");
    disp("dataset "+num2str(k) +" out of " + ndataset);
    disp("--------------------------------------------------------------");
    disp("PSNR "+ num2str(p1)+ " SSIM "+num2str(s1)+ " PSNR "+ num2str(p2)+ " SSIM "+num2str(s2));
end
disp("Mask = "+CA(codes)+ ", "+comple(complementary+1)+", Propagator = "+propagator(typePropagation));
disp("--------------------------------------------------------------");
disp("dataset "+num2str(k) +" out of " + ndataset);
disp("--------------------------------------------------------------");
mean(result,1)
toc
% fig = figure();
% tlo = tiledlayout(fig,1,2,'TileSpacing','None');
% ax = nexttile(tlo);
% imshow(intensityHat,'Parent',ax),title("Intensity PSNR= "+num2str(p1)+" SSIM= "+num2str(s1)),colormap('gray'),pbaspect([1 1 1]),axis off
% ax = nexttile(tlo);
% imagesc(phaseHat,'Parent',ax),title("Phase PSNR= "+num2str(p2)+" SSIM= "+num2str(s2)),colormap('gray'),pbaspect([1 1 1]),axis off;