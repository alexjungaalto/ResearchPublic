function out=simulate_GMS_Ju(task,dim_vals,num_tests,rho_min_fac_vals,graph_on) 

switch nargin
    case 0
        task = 0 ; 
        dim_vals = [64 128 256 512] ; 
        graph_on = 1;
        num_tests = 300; % number of simulations for each point
        num_tests  = 100; 
        rho_min_fac_vals = [1 0.9 0.8 1.2] ; 
%         opt2 = 17;
%         opt3 = @magic;
%     case 3
%         opt2 = 17;
%         opt3 = @magic;
%     case 4
%         opt3 = @magic;
end

rng('shuffle') ; 

if graph_on ==1
close all;
%clear all;
end

%% Parameters
graph_type = 1; % type 1: chain graph

% this paramets controls scale of entries in the precision matrix
c = 2;

num_blocks = 4;



% paper requires min eigenvalues to be at least 1.0
TARGET_MIN_EIG = 1.1;

betaRatio = 1.5; % controls the factor at which eigenvalues of first and last block differ
sparsity = 2; % for chain graph


node_i = 2;

% specify range here
blklenvals = ceil(((2*sparsity):10:100));
blklenvals = ceil(10:10:100);





Results = zeros(length(blklenvals),length(dim_vals));
Results_naive = zeros(size(Results)) ; 

for iter_param = 1:length(dim_vals) 
    
   dimension = dim_vals(iter_param) ;  %  [8 16 64]

 % enumerate all candidate sets for the correct neighborhoud 
 
 sets_T = [];
 num_set_T = 0 ; 
 vec_s = []; 
 for s = 1:1:sparsity
    sets = combnk([(1:(node_i-1)) ((node_i+1):dimension)],s);
    num_set_s = size(sets,1); 
    current_sets= zeros(num_set_s,sparsity) ; 
    num_set_T = num_set_T+num_set_s; 
    vec_s = [vec_s; ones(num_set_s,1)*s] ; 
    for iter_sets = 1:1:num_set_s
        %sets_T =[sets_T;sets(iter_sets,:)]; 
        
        current_sets(iter_sets,1:s) = sets(iter_sets,:) ; 
        
        %sets_T{end+1} = sets(i,:);
    end
    sets_T = [sets_T; current_sets]; 
 end
 
 %sets_T = sets_T';
 %num_set_T = size(sets_T,1);
 Z_T = zeros(num_set_T, 1);

 stat = [];
 
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
 true_N(neighbors) = 1; 
 
% get rho_min and compute reguralizer
 rho_min = get_rho_min(K);
 lambda = rho_min/6;
 
 for  iter_blklen =  1:length(blklenvals)
     
  block_length=blklenvals(iter_blklen) 
  totalAccuracy = 0;
  totalAccuracy_naive=0; 
  
  
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
    
    hat_C = zeros(dimension,dimension,num_blocks) ; 
    
    for idx_blk=1:num_blocks
          block_range = ((idx_blk - 1)*block_length)+ (1:block_length);
          tmp = samples(block_range, :) ; 
              
          hat_C(:,:,idx_blk) = tmp'*tmp; 
          
          %X_i_b = samples(block_range, node_i);
         % X_T_b = samples(block_range, setT);
        %  proj_T_b = (X_T_b / (X_T_b' * X_T_b)) * X_T_b';
        %  projec_stat = projec_stat + norm(X_i_b - proj_T_b * X_i_b,2)^2;
     end
        

   %samples = sample_generator(C, block_length); % create samples from covariance matrices
    %accuracy = test_statistic(dimension, block_length, num_blocks, graph_type);
   
  Z_T = zeros(num_set_T,num_blocks); 
   % card = numel(setT); 
  %projec_stat = 0;
  N_hat_naive = zeros(dimension,1) ; 

  for idx_blk=1:num_blocks
        iter_T = 0 ; 
    for iter_s=1:sparsity 
       dmy = combnk([(1:(node_i-1)) ((node_i+1):dimension)],iter_s);
       num_set_s = size(dmy,1); 
       for iter_j=1:num_set_s
       
       
        %setT = sets_T(iter_T,1:vec_s(iter_T));
        setT = dmy(iter_j,:);
       
        card = iter_s; %vec_s(iter_T) ; 
       
       %block_range = ((idx_blk - 1)*block_length)+ (1:block_length);
       %X_i_b = samples(block_range, node_i);
       %X_T_b = samples(block_range, setT);
       %proj_T_b = (X_T_b / (X_T_b' * X_T_b)) * X_T_b';
        new_err = hat_C(node_i,node_i,idx_blk) - hat_C(setT,node_i,idx_blk)'*inv(hat_C(setT,setT,idx_blk))*hat_C(setT,node_i,idx_blk); 
       %projec_stat = projec_stat + new_err;
       %projec_stat = projec_stat + norm(X_i_b - proj_T_b * X_i_b,2)^2;
        Z_T(iter_T+iter_j,idx_blk) = new_err + block_length*lambda*card;
       end
        iter_T = iter_T +num_set_s ; 
     end
   [~, minIdx] = min(Z_T(:,idx_blk));
   %N_hat = zeros(dimension,1) ; 
   N_hat_naive(sets_T(minIdx,1:vec_s(minIdx))) = 1 ; 
  end
% find the smallest score and corresponding neigbour set
   [~, minIdx] = min(sum(Z_T,2));
   N_hat = zeros(dimension,1) ; 
   N_hat(sets_T(minIdx,1:vec_s(minIdx))) = 1 ; 
   
 
  % estimatedNeigbours = sets_T(minIdx,:);
  % estimatedNeigbours = estimatedNeigbours{1};

% compare estiamted neigbours to true neighbours
  accuracyProj = 0;
  if (norm(N_hat-true_N)<0.1) 
     accuracyProj = 1;
  end
  totalAccuracy  = totalAccuracy + accuracyProj;
  
  accuracyProj = 0;
  if (norm(N_hat_naive-true_N)<0.1) 
     accuracyProj = 1;
  end
  totalAccuracy_naive  = totalAccuracy_naive + accuracyProj;
  
  end
  %stat = [stat; num_blocks*block_length totalAccuracy];
  Results(iter_blklen,iter_param) = totalAccuracy/num_tests; 
  Results_naive(iter_blklen,iter_param) = totalAccuracy_naive/num_tests; 
 end
 %stat(:,2) = stat(:,2)/num_tests;

 %Results{end+1} = {dimension stat}; 
end




% plot the results
if graph_on==1
figure
hold on 
end

for iter_param=1:length(dim_vals)
    x_vals = blklenvals*num_blocks; 
    y_vals = 1-Results(:,iter_param) ;
    y_vals_naiv = 1-Results_naive(:,iter_param) ;
%for c = Results
   % dim = c{1}{1};
   % data = c{1}{2};
   % mtx = [mtx (1-data(:,2))]; 
   filename = sprintf('RawSample_p%02d_task%02d_runs%02d.csv',dim_vals(iter_param),task,num_tests) ; 
   if graph_on==1
   plot(x_vals,y_vals,'-o');
   plot(x_vals,y_vals_naiv,'-x') ; 
   end
   T = array2table([x_vals' y_vals],'VariableNames',{'x','y'});
   writetable(T,filename); 
end
if graph_on==1
hold off 
figure
hold on
end
for iter_param=1:length(dim_vals)
    x_vals = blklenvals*num_blocks*((rho_min_fac_vals(iter_param))^2)/log(dim_vals(iter_param)); 
    y_vals = 1-Results(:,iter_param) ; 
    y_vals_naiv = 1-Results_naive(:,iter_param) ;
%for c = Results
   % dim = c{1}{1};
   % data = c{1}{2};
   % mtx = [mtx (1-data(:,2))]; 
   
   if graph_on==1 
    plot(x_vals,y_vals,'-o');
    plot(x_vals,y_vals_naiv,'-x') ; 
   end
   filename = sprintf('ScaledSample_p%02d_task%02d_runs%02d.csv',dim_vals(iter_param),task,num_tests) ; 
   T = array2table([x_vals' y_vals],'VariableNames',{'x','y'});
   writetable(T,filename); 
   filename = sprintf('ScaledSampleNaive_p%02d_task%02d_runs%02d.csv',dim_vals(iter_param),task,num_tests) ; 
   T = array2table([x_vals' y_vals_naiv],'VariableNames',{'x','y'});
   writetable(T,filename); 
end
%hold off
%xlabel('Sample size');
%ylabel('Error rate')
%legend
%grid on

% x_vals = data(:,1); 
% x_vals = x_vals/max(x_vals) ; 
% mtx=[x_vals mtx]; 
% T = array2table(mtx,'VariableNames',{'x','y','v','w'});
% %csvwrite('hingelosswoheader.csv',mtx);
% writetable(T,'ErrorVsSampeSize.csv');
% 
% figure
% hold on
% for c = Results
%     dim = c{1}{1};
%     data = c{1}{2};
%     plot(data(:,1)*rho_min_fac_vals(dim)/log(dim),1-data(:,2),'-o','DisplayName',num2str(dim));
% end
% hold off
% xlabel('Scaled sample size, N/log(p)');
% ylabel('Error rate')
% legend
% grid on
% 
out = 1 ;
end