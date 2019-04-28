function accuracy = test_statistic_lasso( dimension, block_length, num_blocks, lambda, node_idx)
%=============================================
%TEST_STATISTIC
%   Return 1 if the test statistic yields the true neighborhoods,
%   otherwise, return 0
%=============================================

% chain graph only for now
graph_type = 1;

% paper requires min eigenvalues to be at least 1.0
TARGET_MIN_EIG = 1.1;

% controls the factor at which eigenvalues of first and last block differ
betaRatio = 1.5;

% this paramets controls scale of entries in the precision matrix
c = 2;

% generate precision and covariance matrices
K =  zeros(dimension, dimension, num_blocks);
C = zeros(dimension, dimension, num_blocks);

precision_matrix = prec_mat_generator(dimension, graph_type,c);
cov_matrix = inv(precision_matrix);
cov_matrix = TARGET_MIN_EIG * cov_matrix / min(eig(cov_matrix)); % scale covariance matrix

for i = 1:1:num_blocks
    C(:,:,i) = (1 + (betaRatio-1)*(i-1)/num_blocks) * cov_matrix;
    precision_matrix = inv(C(:,:,i));
    precision_matrix(precision_matrix < 1e-10) = 0; % due to the numerical errors this may be slightly greater than zero
    K(:,:,i) = precision_matrix;
end

% create samples from covariance matrices
samples = sample_generator(C, block_length); 

% find its true neighbors
true_neighbors = neighbor_set(node_idx, graph_type, dimension); 

% find LASSO estimate
estimated_neigbours = lasso(node_idx, samples(:,node_idx),samples(:,[1:node_idx-1 node_idx+1:end]), block_length, lambda); % estimate neighbours using group lasso

% compare estimated neigbours to true ones
accuracy = 0;
if( ~isempty(estimated_neigbours))
    if(all(estimated_neigbours == true_neighbors))
        accuracy = 1;
    end
end

end