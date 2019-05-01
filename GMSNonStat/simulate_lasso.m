close all;
clear all;

%% Parameters
num_blocks = 4;
block_length = 100;
p = 32; % dimensionality
lambda = 200; % regularization for lasso

num_tests = 100; % number of simulation runs

% node index to find neighbours
% should be from 1 to p
node_idx = 16;

% generate precision and covariance matrices
K =  zeros(dimension, dimension, num_blocks);
C = zeros(dimension, dimension, num_blocks);

% generate precision and covariance matrices
 K =  zeros(dimension, dimension, num_blocks);
 C = zeros(dimension, dimension, num_blocks);

 C_0 = (1/3)*rho_min_fac_vals(iter_param)*(ones(dimension,dimension)-eye(dimension)) + eye(dimension);
    
 adj_matrix = adj_chain_graph(dimension) ;
    
 K_0 = C_0.*adj_matrix;
    
 %precision_matrix = prec_mat_generator_Ju(dimension, graph_type, c);
 %cov_matrix = inv(precision_matrix);
 %cov_matrix = TARGET_MIN_EIG*cov_matrix / min(eig(cov_matrix)); % scale covariance matrix
 
 for iter_blk_inner = 1:1:num_blocks
  %  C(:,:,i) = (1 + (betaRatio-1)*(i-1)/num_blocks) * cov_matrix;
    prec_tmp =  K_0 ; %inv(C(:,:,i));
    prec_tmp  (iter_blk_inner,iter_blk_inner+1)=0 ; % remove particular edge from the chain for each block 
    prec_tmp (iter_blk_inner+1,iter_blk_inner)=0 ; 
    prec_tmp(prec_tmp < 1e-10) = 0; % due to the numerical errors this may be slightly greater than zero
    K(:,:,iter_blk_inner) = prec_tmp;
    C(:,:,iter_blk_inner) = inv(prec_tmp) ; 
 end
 
 neighbors = neighbor_set(node_i, graph_type, dimension); % find its true neighbors
 true_N = zeros(dimension,1); 
% true_N(neighbors) = 1; 
 
% get rho_min and compute reguralizer
 rho_min = get_rho_min(K);
% lambda = rho_min/6;
 

%precision_matrix = prec_mat_generator(dimension, graph_type,c);
%cov_matrix = inv(precision_matrix);
%cov_matrix = TARGET_MIN_EIG * cov_matrix / min(eig(cov_matrix)); % scale covariance matrix

% for i = 1:1:num_blocks
%    % C(:,:,i) = (1 + (betaRatio-1)*(i-1)/num_blocks) * cov_matrix;
%     precision_matrix = inv(C(:,:,i));
%     precision_matrix(precision_matrix < 1e-10) = 0; % due to the numerical errors this may be slightly greater than zero
%     K(:,:,i) = precision_matrix;
% end

% create samples from covariance matrices

success_count = 0;
for idx=1:num_tests
  samples = sample_generator(C, block_length); 
  % find its true neighbors
  true_neighbors = neighbor_set(node_idx, graph_type, dimension); 
  % find LASSO estimate
  estimated_neigbours = lasso(node_idx, samples(:,node_idx),samples(:,[1:node_idx-1 node_idx+1:end]), block_length, lambda); % estimate neighbours using group lasso
  % compare estimated neigbours to true ones
  tmp = 0;
  if( ~isempty(estimated_neigbours))
   if(all(estimated_neigbours == true_neighbors))
    tmp = 1;
   end
  end

 % success = test_statistic_lasso_Ju(p, block_length, num_blocks, lambda, node_idx);
  success_count  = success_count + tmp;
end

accuracy = success_count/num_tests



