function [ samples ] = sample_generator( cov_matrices, block_length )
%=============================================
%SAMPLE_GENERATOR 
% Generate a sequence of samples from covariance matrices
% input arguments: 
%       covariance_matrices - Sequence of covariance matrices.
% output arguments: 
%       samples - A sequence of samples corresponding to covariance
%       matrices
%=============================================
[dim, dim, num_blocks] = size(cov_matrices);
num_samples = num_blocks * block_length;
samples = zeros(num_samples, dim);
mu = zeros(1, dim);
for idx_blk=1:num_blocks
    samples((idx_blk - 1)*block_length + 1:idx_blk*block_length, :) = mvnrnd(mu, cov_matrices(:,:,idx_blk), block_length);
end

end

