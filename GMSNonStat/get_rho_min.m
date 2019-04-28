function [ rho_min ] = get_rho_min(K)

num_blocks = size(K,3);
rho = zeros(size(K,1),size(K,2));
for i=1:1:num_blocks
    M = K(:,:,i);
    d = (diag(M).^2);
    denom =repmat(d,1,size(M,2));
    rho = rho + (M.^2)./denom;
end

rho = rho/num_blocks;
rho_min = min(nonzeros(rho));

end

