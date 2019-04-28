close all;
clear all;


task = 0 ; 
dim_vals = [64 128 256 512] ;
%dim_vals = [64 128 256 512] ;
graph_on = 1;
num_tests = 300; % number of simulations for each point
num_tests  = 100; 
rho_min_fac_vals = [1.2 0.9 0.8 1.2] ; 
graph_type=1 ; 


%% Parameters
num_blocks = 4;


lambda = 20; % regularization for lasso
blklenvals = ceil(10:10:100); 
block_length = 40;

num_tests = 100; % number of simulation runs

for iter_param=1:length(dim_vals) 
 p = dim_vals(iter_param); % dimensionality
 dimension=p; 

 % node index to find neighbours
% should be from 1 to p
node_idx = 2;
node_i=2;
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

 true_neighbors = neighbor_set(node_idx, graph_type, dimension);
accuracy = zeros(length(blklenvals),1); 

for iter_blklen=1:length(blklenvals)
    block_length=blklenvals(iter_blklen) 
success_count = 0;
for idx=1:num_tests
  samples = sample_generator(C, block_length); 
  % find its true neighbors
 
  % find LASSO estimate
  estimated_neigbours = lasso_Ju(node_idx, samples(:,node_idx),samples(:,[1:node_idx-1 node_idx+1:end]), block_length, lambda,rho_min); % estimate neighbours using group lasso
  % compare estimated neigbours to true ones
   hat_N = zeros(dimension,1); 
 hat_N(true_neighbors) = 1; 
 hat_N(estimated_neigbours)=hat_N(estimated_neigbours)-1; 
 
  tmp = 0;
   if(norm(hat_N) < 1/100)
    tmp = 1;
   end
  

 % success = test_statistic_lasso_Ju(p, block_length, num_blocks, lambda, node_idx);
  success_count  = success_count + tmp;
end

accuracy(iter_blklen) = success_count/num_tests  ; 
end


    x_vals = blklenvals*num_blocks; 
    y_vals = 1-accuracy(:) ;

   filename = sprintf('SampleGLasso_p%02d_task%02d_runs%02d.csv',dim_vals(iter_param),task,num_tests) ; 
   T = array2table([x_vals' y_vals],'VariableNames',{'x','y'});
   writetable(T,filename); 
     
   x_vals= blklenvals*num_blocks*((rho_min_fac_vals(iter_param))^2)/log(dim_vals(iter_param)); 

  
    y_vals = 1-accuracy(:) ;

   filename = sprintf('ScaledSampleGLasso_p%02d_task%02d_runs%02d.csv',dim_vals(iter_param),task,num_tests) ; 
   T = array2table([x_vals' y_vals],'VariableNames',{'x','y'});
   writetable(T,filename); 
   
end



