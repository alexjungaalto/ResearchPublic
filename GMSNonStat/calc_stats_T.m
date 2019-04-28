function [Z_T] = calc_stats_T(samples, node_i, setT, block_length)
num_blocks = size(samples,1)/block_length;
Z_T = 0;
for idx_blk=1:num_blocks
    block_range = (idx_blk - 1)*block_length + (1:(idx_blk*block_length));
    X_i_b = samples(block_range, node_i);
    X_T_b = samples(block_range, setT);
    proj_T_b = (X_T_b / (X_T_b' * X_T_b)) * X_T_b';
    Z_T = Z_T + power(norm(X_i_b - proj_T_b * X_i_b),2);
end

Z_T = Z_T/(block_length*num_blocks);

end