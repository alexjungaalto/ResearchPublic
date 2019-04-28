function [accuracyProj] = test_statistic( dimension, block_length, num_blocks, graph_type)
%=============================================
%TEST_STATISTIC
%   Return 1 if the test statistic yields the true neighborhoods,
%   otherwise, return 0
%=============================================

% paper requires min eigenvalues to be at least 1.0
TARGET_MIN_EIG = 1.1;

betaRatio = 1.5; % controls the factor at which eigenvalues of first and last block differ
sparsity = 2; % for chain graph

% generate precision and covariance matrices
K =  zeros(dimension, dimension, num_blocks);
C = zeros(dimension, dimension, num_blocks);

% this paramets controls scale of entries in the precision matrix
c = 2;

precision_matrix = prec_mat_generator(dimension, graph_type, c);
cov_matrix = inv(precision_matrix);
cov_matrix = TARGET_MIN_EIG*cov_matrix / min(eig(cov_matrix)); % scale covariance matrix
for i = 1:1:num_blocks
    C(:,:,i) = (1 + (betaRatio-1)*(i-1)/num_blocks) * cov_matrix;
    precision_matrix = inv(C(:,:,i));
    precision_matrix(precision_matrix < 1e-10) = 0; % due to the numerical errors this may be slightly greater than zero
    K(:,:,i) = precision_matrix;
end

samples = sample_generator(C, block_length); % create samples from covariance matrices

% get rho_min and compute reguralizer
rho_min = get_rho_min(K);
lambda = rho_min/6;

%node_i = randsample(2:dimension-1,1); % sample non-border node
%node_i = randsample([1 dimension],1); % sample border node
node_i = 1
neighbors = neighbor_set(node_i, graph_type, dimension); % find its true neighbors

% enumerate all the candidates
sets_T = {};
for s = 1:1:sparsity
    sets = combnk([1:node_i-1 node_i+1:dimension],s);
    for i = 1:1:size(sets,1)
        sets_T{end+1} = sets(i,:);
    end
end

% allocate some memory
sets_T = sets_T';
num_set_T = size(sets_T,1);
Z_T = zeros(num_set_T, 1);

% calcualte statistics
for i=1:num_set_T
    setT = sets_T{i};
    card = numel(setT);
    Z_T(i) = calc_stats_T(samples, C, node_i, setT, block_length) + lambda*card;
end

% find the smallest score and corresponding neigbour set
[~, minIdx] = min(Z_T);
estimatedNeigbours = sets_T(minIdx,:);
estimatedNeigbours = estimatedNeigbours{1};

% compare estiamted neigbours to true neighbours
accuracyProj = 0;
if(isequal(estimatedNeigbours,neighbors))
    accuracyProj = 1;
end

end