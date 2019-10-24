function out = nLassoLinModel(G_SBM,MAX_ITER,lambda,sigma2,X_mtx,samplingset,true_sig) 

[d,N] = size(X_mtx); 
Adjac = triu(G_SBM,1) ; 
A_undirected = Adjac+Adjac' ; 
degrees = sum(A_undirected,1); 
inv_degrees = 1./degrees';

%%%% create weighted incidence matrix 
G = digraph(triu(G_SBM,1)) ;
D = sparse(incidence(G)') ;
D_block= kron(sparse(D),sparse(eye(d,d))) ; 

y = zeros(1,N); 
y(samplingset)= true_sig(samplingset); 

[M, N] = size(D); 

Lambda = diag(sparse(1./(sum(abs(D),2)))) ; 
Lambda_block = kron(Lambda,eye(d,d)) ;
Gamma_vec=(sparse(1./(sum(abs(D),1))))' ;
Gamma = diag(sparse(Gamma_vec));  
Gamma_block = sparse(kron(Gamma,eye(d,d))) ; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Algorithm Initialisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


hatx = zeros(N*d,1); 
running_average =zeros(N*d,1);
%haty = ((1:(N-1))/(N-1))'; 
haty = zeros(M*d,1); 

dmy = length(Gamma_vec) ;
mtx_A_block = zeros(N*d,N*d); 
mtx_B_block = zeros(N*d,N*d); 
for iter_node=1:N 
    msk_dmy = zeros(N,N) ; 
    msk_dmy(iter_node,iter_node) = 1; 
    tilde_tau = length(samplingset)/(2*Gamma_vec(iter_node)) ; 
    
    mtx_A = tilde_tau*eye(d,d)*inv((1/sigma2)*X_mtx(:,iter_node)*X_mtx(:,iter_node)'+tilde_tau*eye(d,d)); 
    mtx_B = inv((1/sigma2)*X_mtx(:,iter_node)*X_mtx(:,iter_node)'+tilde_tau*eye(d,d)); 
    mtx_A_block = sparse(mtx_A_block) + sparse(kron(msk_dmy,mtx_A)) ; 
     mtx_B_block = sparse(mtx_B_block) + sparse(kron(msk_dmy,mtx_B)) ; 
end

%tilde_tau = length(samplingset)*(1./(2*diag(Gamma_block))) ; 
%mtx_A = diag(ones(dmy,1)./(ones(dmy,1)+2*Gamma_vec/length(samplingset))) ; 

%mtx_B = diag(ones(dmy,1)./(ones(dmy,1)+tilde_tau)) ; 
vec_B = (1/sigma2)*mtx_B_block*reshape(X_mtx*diag(y),N*d,1) ; 

%%% now compute nLasso 

for iterk=1:MAX_ITER
    
  % LP iteration 
  %  hatxLP = inv_degrees.*(A_undirected*hatxLP); 
  %  hatxLP(samplingset) = graphsig(samplingset); 
    
    
    newx = hatx - 0.9*Gamma_block*(D_block'*haty) ; 
    
    %%%% update for least absoluate deviation
    %%%% newx = block_thresholding(newx,samplingset,y,X_mtx,d,N,Gamma_vec) ; 
    
    %%%% update for least squared linear regression 
    old = newx ; 
    newx = update_x_linreg(newx,samplingset,d,N,mtx_A_block,vec_B) ; 
    %newx(samplingset) = vec_B(samplingset); 
    
    %%%% update for sparse label propagation 
    %newx(samplingset) = graphsig(samplingset) ;
    
    tildex = 2*newx-hatx; 
    newy = haty + Lambda_block*(D_block*tildex); 
    haty = block_clipping(newy,d,M,lambda)  ;
  % haty=newy;
    hatx = newx; 
    running_average = (running_average*(iterk-1) +hatx)/iterk; 
    
    
    %running_averagey = (running_average*(iterk-1) +haty)/iterk;+
    %  dual = sign(D*running_average); 
    % dual(iterk:(N-1)) = 0 ;
    % log_conv(iterk) = sum(abs(D*running_average)); 
    % log_bound(iterk) =(1/(2*iterk))*(primSLP'*inv(Gamma)*primSLP)+((dualSLP-dual)'*inv(Lambda)*(dualSLP-dual)) ;
end

out = running_average ; 

end



function weights_out = block_clipping (weights_in,feature_dim,nr_edges,lambda) 
%%% input: weights_in vector of length featuredim*nr_datapoints, scalar
%%% lambda 
 mtx = reshape(weights_in,feature_dim,nr_edges);
 weights_out = mtx; 
 x_norm = sqrt(sum(mtx.^2,1)) ; 

 
 idx_exceed = find(x_norm>lambda) ; 
 factor = ones(1,nr_edges); 
 
 for iter_idx=1:length(idx_exceed) 
     idx = idx_exceed(iter_idx) ; 
     tmp = weights_out(:,idx); 
     weights_out(:,idx) = tmp*lambda./x_norm(idx);
    % factor(idx) = lambda./x_norm(idx); 
 end
 
 
 %factor = max([x_norm;lambda*ones(1,nr_edges)],[],1) ;  
 %weights_out= mtx*diag(factor);
 
 weights_out = reshape(weights_out,feature_dim*nr_edges,1) ; 
 


end 

function weights_out = block_thresholding (weights_in,sampling_set,y,feature_mtx,feature_dim,nr_nodes,tau_vec) 
%%% input: weights_in vector of length featuredim*nr_datapoints, scalar
%%% lambda 
 
 mask_sampling_set = zeros(1,nr_nodes); 
 mask_sampling_set(sampling_set) = 1 ; 
 mask_sampling_set= ones(feature_dim,1)*mask_sampling_set; 
 
 mask_mtx =kron(eye(nr_nodes,nr_nodes),ones(feature_dim,feature_dim)) ;
 
 feature_norm_squared_vec = sum(feature_mtx.^2,1) ; 
 
 norm_features_2 = kron(diag((feature_norm_squared_vec)),eye(feature_dim,feature_dim)) ; 
 
 feature_mtx_vec= reshape(feature_mtx,nr_nodes*feature_dim,1); 
 
 proj_feature = ((feature_mtx_vec*feature_mtx_vec').*mask_mtx)*inv(norm_features_2);
 
 %%%% project input weight vector on orthogonal complement of feature
 %%%% vector space 
 
 out_of_feature = (eye(nr_nodes*feature_dim,nr_nodes*feature_dim)-proj_feature)*weights_in ; 
 
 %%%% compute coefficient for input weights on feature direction input
 
 weights_in_row = reshape(weights_in,feature_dim,nr_nodes); 
 
 weights_in_coeff = sum(feature_mtx.*weights_in_row,1)./feature_norm_squared_vec ; 
 
 coeffs = zeros(nr_nodes,1); 
 
 %coeffs(sampling_set) = y(sampling_set)./norm_sqared_features(sampling_set) ;
 
 for iter_node=1:length(sampling_set) 
    node_idx = sampling_set(iter_node); 
    dmy = weights_in_coeff(node_idx) - (y(node_idx)/feature_norm_squared_vec(node_idx)) ; 
    dmy = wthresh(dmy,'s',tau_vec(node_idx)) ; 
    coeffs (node_idx) = (y(node_idx)/feature_norm_squared_vec(node_idx))+dmy; 
 end  
 
 
 %mtx_w = reshape(weights_in,feature_dim,nr_nodes);
 
 
%  component1 = zeros(size(mtx_w)) ; 
%  component2 = zeros(size(mtx_w)) ; 
%  
%  mtx_x = feature_mtx; 
%  
%  tau_vec = reshape(tau_vec,nr_nodes,1); 
%  tau_mtx = ones(feature_dim,1)*tau_vec'; 
 
%  y_hat = sum(mtx_w.*mtx_x,1) ; 
%  y     = reshape(y,1,nr_nodes); 
%  
%  
%  x_norm = sum(mtx_x.^2,1).*tau_vec' ;
%  
%  index_vec = find(y > (y_hat + x_norm)); 
%  
%  dmy = (mtx_x.*tau_mtx).*mask_sampling_set ; 
%  
%  component1(:,index_vec) = dmy(:,index_vec) ; 
%  
%  index_vec = find(y < (y_hat + x_norm)); 
%  
%  dmy = - (mtx_x.*tau_mtx).*mask_sampling_set ; 
%  component2(:,index_vec) = dmy(:,index_vec) ; 
%  
%  mtx_w = mtx_w + component1+component2; 
%  
 
 tmp = feature_mtx*diag(coeffs) + reshape(out_of_feature,feature_dim,nr_nodes); 

 weights_out = reshape(weights_in,feature_dim,nr_nodes) ; 
 weights_out(:,sampling_set) = tmp(:,sampling_set); 
 weights_out = reshape(weights_out,feature_dim*nr_nodes,1) ; 
 
end 


function weights_out = update_x_linreg (weights_in,sampling_set,feature_dim,nr_nodes,mtx_A,vec_B) 
%%% input: weights_in vector of length featuredim*nr_datapoints, scalar
%%% lambda 
 %DiagTAU = length(sampling_set).*inv(2*diag(tau_vec));
 

 tmp = mtx_A*weights_in + vec_B ; %((feature_mtx*diag(y))+w_in*DiagTAU).*inv(eye(length(tau_vec))  %  + reshape(out_of_feature,feature_dim,nr_nodes); 
% tmp = mtx_A_block*weights_in + vec_B ; %((feature_mtx*diag(y))+w_in*DiagTAU).*inv(eye(length(tau_vec))  %  + reshape(out_of_feature,feature_dim,nr_nodes); 
 
 tmp = reshape(tmp,feature_dim,nr_nodes); 
 weights_out = reshape(weights_in,feature_dim,nr_nodes); 
 weights_out(:,sampling_set) = tmp(:,sampling_set); 
 weights_out = reshape(weights_out,feature_dim*nr_nodes,1) ; 
 
end 

