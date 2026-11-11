function [masks] = selectMasks(codes,M,N,nm,complementary)
masks = zeros(M,N,nm);
if(codes == 1)
    %% masks
    load('mascaras.mat')
    for i = 1:nm
        mask(:,:,i) = Set1{i};%rand(M,N)>0.25; %Set1{i};
    end
    %%Random
    for i=1:nm
        masks(1:M,1:N,i) = rand(M,N,1)> 0.5;
    end
    %% Random Complementary by pairs
    if(complementary == 1)
        [masks] = complementary_masks(masks);
    end

elseif(codes == 2)
    addpath(genpath("../2024-OOPhare-main/mask/Polar sphere packing/"+num2str(nm)));
    %% masks
    c = 1;
    for i=1:nm
        d = num2str(i);
        masks(:,:,c) = imread("image_"+d+".png");
        c = c + 1;
    end
    if(nm==8)
        load('filter384_8.mat');
    end
    %% Sphere packing Complementary by pairs
    if(complementary == 1)
        [masks] = complementary_masks(masks);
    end
    temp = sum(masks,3);
    if(max(temp(:))==1)
        disp('full complementary')
    else
        disp('no full complementary')
    end
elseif(codes == 3)
    %nm2 = 16;
    c = 1;
    addpath(genpath("../2024-OOPhare-main/mask/DCT/"+num2str(nm)));
    for i=1:nm
        d = num2str(i);
        masks(:,:,c) = imread("image_"+d+".png");
        %masks(:,:,c) = masks(:,:,c);
        c = c + 1;
    end
    %% DCT Complementary by pairs
    if(complementary == 1)
        [masks] = complementary_masks(masks);
    end
elseif(codes == 4)
    addpath(genpath("../2024-OOPhare-main/mask/Hadamard/"));
    load('Hadamard384_16_2.mat')   % 59.8158    0.9962   57.1708    0.9483
    %load('Hadamard384_16_4.mat') %  59.2008
    %load('Hadamard384_16_5.mat') % 59.9313    0.9962   57.2458    0.9558
    %load('Hadamard384_16_6.mat') % 59.4971    0.9962   56.6825    0.9475
end
end
