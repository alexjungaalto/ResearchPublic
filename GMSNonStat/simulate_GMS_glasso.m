function out= simulate_GMS_glasso(task,dim_vals,num_tests,rho_min_fac_vals,graph_on)

switch nargin
    case 0
        task = 0 ;
        dim_vals = [64 128 256 512] ;
        graph_on = 1;
        num_tests = 100; % number of simulations for each point
        rho_min_fac_vals = [1 0.9 0.8 1.2] ;
end

rng('shuffle') ;

%% Parameters
graph_type = 1; % type 1: chain graph
gamma = 20; % group lasso regularizer

num_blocks = 4;
sparsity = 2; % for chain graph
node_i = 2;

% specify range here
blklenvals = ceil(10:10:100);

glasso_results = zeros(length(blklenvals),length(dim_vals));

%% repeat over each dimension value
for iter_param = 1:length(dim_vals)
    
    dimension = dim_vals(iter_param) ; 
    
    % generate precision and covariance matrices
    K =  zeros(dimension, dimension, num_blocks);
    C = zeros(dimension, dimension, num_blocks);
    
    C_0 = (1/3)*rho_min_fac_vals(iter_param)*(ones(dimension,dimension)-eye(dimension)) + eye(dimension);
    
    adj_matrix = adj_chain_graph(dimension) ;
    
    K_0 = C_0.*adj_matrix;
    
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
    true_N(neighbors) = 1;
    
    % get rho_min and compute reguralizer
    rho_min = get_rho_min(K);
    lambda = rho_min/6;
    
    %for each block length 
    for  iter_blklen =  1:length(blklenvals)
        
        block_length=blklenvals(iter_blklen)
        gls_accuracy = 0;
              
        %%% loop over i.i.d. runs
        for idx=1:num_tests
            
            % [dmy, dim, num_blocks] = size(cov_matrices);
            num_samples = num_blocks * block_length;
            samples = zeros(num_samples, dimension);
            mu = zeros(1, dimension);
            for idx_blk=1:num_blocks
                
                dmy_samples =  mvnrnd(mu, C(:,:,idx_blk), block_length);
                %    dmy_samples = filter([1 0 0],[1],dmy_samples) ;
                %    dmy_samples = fft(dmy_samples)/sqrt(block_length) ;
                samples((idx_blk - 1)*block_length + 1:idx_blk*block_length, :) = dmy_samples ; %mvnrnd(mu, C(:,:,idx_blk), block_length);
                
            end
            % estimate neighbours using group lasso
            glasso_estimator = lasso_Ju(node_i, samples(:,node_i),samples(:,[1:node_i-1 node_i+1:end]), block_length, gamma, rho_min)
            glasso_N = zeros(dimension,1);
            glasso_N(glasso_estimator) = 1;
            
            if (norm(glasso_N-true_N)<0.1)
                gls_accuracy = gls_accuracy + 1;
            end
            
        end
        
        glasso_results(iter_blklen,iter_param) = gls_accuracy/num_tests;
        
    end
end

% plot the results
if graph_on==1
    figure
    hold on
end

for iter_param=1:length(dim_vals)
    x_vals = blklenvals*num_blocks;
    y_vals_gls = 1-glasso_results(:,iter_param) ;
    
    if graph_on==1
        plot(x_vals,y_vals_gls,'-o');
    end

    gfilename = sprintf('glsRawSample_p%02d_task%02d_runs%02d.csv',dim_vals(iter_param),task,num_tests) ;
    gT = array2table([x_vals' y_vals_gls],'VariableNames',{'x','y'});
    writetable(gT,gfilename);
    
end
if graph_on==1
    hold off
    figure
    hold on
end
for iter_param=1:length(dim_vals)
    x_vals = blklenvals*num_blocks*((rho_min_fac_vals(iter_param))^2)/log(dim_vals(iter_param));
    y_vals_gls = 1-glasso_results(:,iter_param) ;
    
    if graph_on==1
        plot(x_vals,y_vals_gls,'-o');
    end
    
    gfilename = sprintf('glsScaledSample_p%02d_task%02d_runs%02d.csv',dim_vals(iter_param),task,num_tests) ;
    gT = array2table([x_vals' y_vals_gls],'VariableNames',{'x','y'});
    writetable(gT,gfilename);
    
end