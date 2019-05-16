clear all
close all
%%%% Generate two separate ER graphs of size N/2 with edge probability p
%%%% connect nodes between two cluter ERs using edge probability q
restoredefaultpath
rehash toolboxcache

[pathtothismfile,name,ext] = fileparts(mfilename('fullpath')) ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graph Initialisation SBM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Intra-edge probability of a cluster:
p = .1;
% Inter-edge probability between clusters
%q = 0.99;



% Number of nodes per cluster
N1 = 40;
N2 = 40;
% total number of nodes 
N=N1+N2;
RUNS = 10; 

% numer of features at each data point 

d=2; 


%G = zeros(N,N) ;




avg_degree=10; 
boundary_size_values = [5,6,7,8,9,10,12,14,15,16,18,20,22,25,27,30,35,40,50,60,70,80]';
boundary_size_values = [(2:2:10),20,40,60]';


mse_over_boundary = zeros(length(boundary_size_values),1); 

ratio_bound_flow = zeros(length(boundary_size_values),RUNS); 

for iter_boundary=1:length(boundary_size_values)
    boundary_edges = boundary_size_values(iter_boundary); 

    %boundary_edges = 12; 
for iter_RUNS=1:RUNS
    
%%%%%%%%%%%%%
% generate empirical graph by sparsely connected two ER graph 
%%%%%%%%%%%%%%
G_SBM= twocluster(N1,N2,avg_degree,boundary_edges,1); 

%%%%% 
% generate feature vectors for each data point 
%%%%%

X_mtx = randn(d,N) ; 
X_mtx = X_mtx*inv(diag(sqrt(diag(X_mtx'*X_mtx)))) ; 
w1 = ones(d,1) ; %randn(d,1) ;
w1 = w1/norm(w1);
w2 = -ones(d,1) ; 
%w2 = randn(d,1) ; 
w2 = w2/norm(w2) ; 
barw = (1/100)*[w1*ones(1,N1),w2*ones(1,N2)]; % matrix barw hold in its colums the true weight vector for each node 
barw_vec = reshape(barw,d*N,1); %vectorized form of barw

%triu(G_SBM,1) ;

% define training set and true graph signal 
%graphsig = 0.1*[ones(N1,1);-1*ones(N2,1)] ;

y   = sum(X_mtx.*barw,1) ; 

samplingset = [1:3 (N-2):N]; 
%samplingset = [1:N]; 

%% determine flow from sampled node 1 to boundary between first 
%% N1 nodes and second N2 nodes 

G1 = G_SBM; 
G1(1:N1,(N1+1):N) = 2*G1(1:N1,(N1+1):N) ; 
G1((N1+1):N,(N1+1):N) = (10^3) *ones(N2,N2) ;
%G1(1:N1,1:N1) = (10^10) *ones(N1,N1) ;
G1 = triu(G1,1); 
G1 = G1+G1' ; 

%%flowgraph1(((N1+1):N):((N1+1):N)) = 10^10*ones( ; 
 flowgraph = digraph(G1); 
 flowsamplednodes=maxflow(flowgraph,1,N) ;
 flow_boundary=sum(sum(G_SBM(1:N1,(N1+1):N))) ; 
 
 if (flow_boundary < 1) 
     flow_boundary = 1 ; 
 end
 
 ratio_bound_flow(iter_boundary,iter_RUNS) = flowsamplednodes/flow_boundary ; 
 
% %samplingset = [samplingnode1 samplingnode2]; 




%plot(nodes(samplingset,1),nodes(samplingset,2),'cx','Color',[1,0,0])


% for iter_node=1:N 
%     for iter_node1=1:iter_node 
%         
%         if G(iter_node,iter_node1) >0
%             plot([nodes(iter_node,1) nodes(iter_node1,1)],[nodes(iter_node,2) nodes(iter_node1,2)],'ro')
% 
%            hline = gline; % Connect circles
%             set(hline,'Color','r')
%         end
%     end
% end

%D = triu(G); 
Adjac = triu(G_SBM,1) ; 
A_undirected = Adjac+Adjac' ; 
degrees = sum(A_undirected,1); 
inv_degrees = 1./degrees';

%%%% create weighted incidence matrix 
G = digraph(triu(G_SBM,1)) ;
D = sparse(incidence(G)') ;
D_block= kron(D,eye(d,d)) ; 

[M, N] = size(D); 
edge_weights = zeros(M,1); 
%for iter_edge=1:M
%    [s,t] = findedge(G,iter_edge); %finds the source and target nodes of the edges specified by idx.
%     edge_weights(iter_edge) = sqrt(A_undirected(s,t)) ; 
%end
%D = diag(edge_weights)*D ; 

%%%%% some visio

%scatter(nodes(:,1),nodes(:,2)) ; 
%figure(1);
%plot(G);   
%hold on 

Lambda = diag((1./(sum(abs(D),2)))) ; 
Lambda_block = kron(Lambda,eye(d,d)) ;
Gamma_vec=(1./(sum(abs(D),1)))' ;
Gamma = diag(Gamma_vec);  
Gamma_block = kron(Gamma,eye(d,d)) ; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Algorithm Initialisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


primSLP = ones(N,1); 
primSLp(N) = 0 ; 
%running_average;
dualSLP =(1:(N-1))'/(N-1) ; %running_averagey; 
dualSLP = zeros(M,1); 


hatx = zeros(N,1); 


lambda =1/10 ; 
lambda = 1/9; 

hatx = zeros(N*d,1); 
running_average =zeros(N*d,1);
%haty = ((1:(N-1))/(N-1))'; 
haty = zeros(M*d,1); 
running_averagey = zeros(M*d,1); 
hatxLP = zeros(N*d,1); 
hatw = zeros(d*N,1) ; 


%log_conv=zeros(N,1); 
%log_bound=zeros(N,1); 

for iterk=1:1000
    
  % LP iteration 
  %  hatxLP = inv_degrees.*(A_undirected*hatxLP); 
  %  hatxLP(samplingset) = graphsig(samplingset); 
    
    % SLP iteration
    newx = hatx - 0.5*Gamma_block*(D_block'*haty) ; 
    newx = block_thresholding(newx,samplingset,y,X_mtx,d,N,Gamma_vec) ; 
    %newx(samplingset) = graphsig(samplingset) ; 
    tildex = 2*newx-hatx; 
    newy = haty + 0.5*Lambda_block*(D_block*tildex); 
    haty = block_clipping(newy,d,M,lambda)  ;
    hatx = newx; 
    running_average = (running_average*(iterk-1) +hatx)/iterk; 
    
    
    %running_averagey = (running_average*(iterk-1) +haty)/iterk;+
    %  dual = sign(D*running_average); 
    % dual(iterk:(N-1)) = 0 ;
    % log_conv(iterk) = sum(abs(D*running_average)); 
    % log_bound(iterk) =(1/(2*iterk))*(primSLP'*inv(Gamma)*primSLP)+((dualSLP-dual)'*inv(Lambda)*(dualSLP-dual)) ;
end

mse_over_boundary(iter_boundary)=mse_over_boundary(iter_boundary)+norm(barw_vec-running_average)^2;
end
end
%figure(1); 
%stem(primSLP);
%title('primal SLP')
%figure(2); 
%stem(dualSLP); 
%title('dual SLP')
figure(2); 
norm_w_hat=reshape(running_average,d,N) ; 
norm_w_hat=sum(norm_w_hat.^2,1); 

stem(norm_w_hat./sum(barw.^2,1));
title('output nLasso');
figure(4); 
x_vals = boundary_size_values/(avg_degree*length(samplingset)); 
x_vals = sum(ratio_bound_flow,2)/RUNS ; 
mselog = mse_over_boundary/(RUNS) ; 
mselog=mselog/norm(barw_vec)^2 ; %norm(graphsig)^2 ; 
stem(x_vals,mselog); 

mtx=[x_vals mselog]; 
%mtx = flipup(mtx); 
T = array2table(mtx,'VariableNames',{'a','b'});
%csvwrite('hingelosswoheader.csv',mtx);
filename = sprintf('MSEoverBoundary_%s.csv',date) ;
writetable(T,fullfile(pathtothismfile,filename));


%figure(4); 
%stem(dual);
%bound =log_bound+(1./(2*(1:K))*(hatx'*inv(Lambda)*hatx) ; 
%plot(1:N,log([log_conv log_bound])); 



function G_SBM=twocluster(N1,N2,avg_degree,boundary_edges,intra_weight) 

%%%% generates two cluster network; each nodes has avg_degree neighbours in
%%%% same cluster with weight intra_weight 
%% number of edges between clusters is boundary edges (with weight 1) 

p_in = avg_degree/N1; 

A = [[(rand(N1)<p_in) (rand(N1,N2)<(boundary_edges/(N1*N2)))];[zeros(N2,N1) (rand(N2)<p_in)]]; 
G_SBM = triu(A,1);
G_SBM = G_SBM + G_SBM' ; 

check_degrees=sum(G_SBM,2); 
idx_sing = find(check_degrees <0.5) ; 
for iter_i=1:length(idx_sing) 
    node_a = idx_sing(iter_i); 
    if node_a > N1
        
        G_SBM(node_a,N1) = 1 ; 
        G_SBM(N1,node_a) = 1 ; 
    else 
        G_SBM(node_a,N1+1) = 1 ; 
        G_SBM(N1+1,node_a) = 1 ; 
    end
    
end

G_SBM(1:N1,1:N1) = intra_weight*G_SBM(1:N1,1:N1)  ; 
G_SBM((N1+1):(N1+N2),(N1+1):(N1+N2)) = intra_weight*G_SBM((N1+1):(N1+N2),(N1+1):(N1+N2)) ; 

end


%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
% generate random geometric graph on unit square [0,1] x [0,1] %%% 
%%%%%%%%%%%%%%%%

% function [G,samplingset]=RandGeomGraph(N1,N2)
% 
% mean1 =[0,0]; 
% mean2 = [4,0]; 
% cluster1 = randn(N1,2)*diag([1,0.01])+ones(N1,2)*diag(mean1); % sample 50 nodes from Guassian with mean [0.2,0.2] 
% cluster2 = randn(N2,2)*diag([1,0.01])+ones(N2,2)*diag(mean2); % sample 50 nodes from Guassian with mean [0.2,0.2] 
% cluster1(1,:) = [mean1]; 
% cluster2(N2,:) = [mean2]; 
% nodes = [cluster1;cluster2]; 
% diff_mtx_x = repmat(nodes(:,1),1,N) ; 
% diff_mtx_x = diff_mtx_x - diff_mtx_x'; 
% diff_mtx_x = diff_mtx_x.^2 ; 
% nodes = [cluster1;cluster2]; 
% diff_mtx_y = repmat(nodes(:,2),1,N) ; 
% diff_mtx_y = diff_mtx_y - diff_mtx_y'; 
% diff_mtx_y = diff_mtx_y.^2 ; 
% diff_mtx = diff_mtx_x+diff_mtx_y ; 
% 
% [dmy samplingnode1] = min(nodes(1:N1,1)); 
% [dmy samplingnode2] = max(nodes((N1+1):N,1)); 
% samplingnode2=samplingnode2+N1;
% samplingset = [samplingset1 samplingset2]; 
% 
% 
% %%%% nearest neighbour graph
% radius = 0.3 ; 
% G(diff_mtx<radius^2)= 1 ; 
% %%%%%%%%%%%
% 
% %%%% gaussian filter graph
% %sigma = 0.4; 
% %G = exp(-diff_mtx/(2*sigma^2)); 
% 
% %%%%%%%%%%%
% end


function weights_out = block_clipping (weights_in,feature_dim,nr_edges,lambda) 
%%% input: weights_in vector of length featuredim*nr_datapoints, scalar
%%% lambda 
 mtx = reshape(weights_in,feature_dim,nr_edges);
 x_norm = sqrt(sum(mtx.^2,1)) ; 

 
 idx_exceed = find(x_norm>lambda) ; 
 factor = ones(1,nr_edges); 
 
 for iter_idx=1:length(idx_exceed) 
     idx = idx_exceed(iter_idx) ; 
     factor(idx) = lambda./x_norm(idx); 
 end
 
 
 %factor = max([x_norm;lambda*ones(1,nr_edges)],[],1) ;  
 weights_out= mtx*diag(factor);
 
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


