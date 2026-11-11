function [T] = complementary_masks(masks)
[M,N,nm] = size(masks);
T = zeros(M,N,nm);
c =0;
for i=1:2:nm
    c =c + 1;
    T(1:M,1:N,i) = masks(:,:,c);
    T(1:M,1:N,i+1) = ~masks(1:M,1:N,i);
end
end