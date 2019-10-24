%%% source code for the experiment discussed in Sec. VI-A of the paper 
%%% "Semi-supervised Learning in Network-Structured Data via 
%%%   Total Variation Minimization"

clear all
close all

restoredefaultpath
rehash toolboxcache

[pathtothismfile,name,ext] = fileparts(mfilename('fullpath')) ; 

% Number of nodes per cluster
N1 = 100;
N2 = 100;
% total number of nodes 
N=N1+N2;
RUNS = 100; 
RUNS = 10; 
RUNS = 1 ; 

lambda_vals = [1/10000,1/1000,1/100] ; % values for nLasso lambda
nLasso_out = zeros(length(lambda_vals),N); % resutls of nLasso for different lambdas are stored here

%G = zeros(N,N) ;

d=1 ; % parameter for networked linear model (used for nLasso)
sigma2 = 1 ; % noise variance in networked linear model 

X_mtx = ones(d,N);   % set node features equal to 1 so that linear Gaussian model reduces to scalar signal in noise model 

avg_degree=20; 
boundary_size_values = [1,5,6,7,8,9,10,12,14,15,16,18,20,22,25,27,30,35,40,50,60,70,80]';
boundary_size_values = [80]';
mse_over_boundary = zeros(length(boundary_size_values),1); 
mse_over_boundary_nLasso = 0*mse_over_boundary; 
ratio_bound_flow = zeros(length(boundary_size_values),RUNS); 

for iter_boundary=1:length(boundary_size_values)
    boundary_edges = boundary_size_values(iter_boundary); 

    %boundary_edges = 12; 
for iter_RUNS=1:RUNS
    
%%%%%%%%%%%%%
% generate graph by sparsely connected two ER graph 
%%%%%%%%%%%%%%
G_SBM= twocluster(N1,N2,avg_degree,boundary_edges); 
        
    

%triu(G_SBM,1) ;

% define training set and true graph signal 
graphsig = 0.1*[ones(N1,1);-1*ones(N2,1)] ;
samplingset = [1 N]; 

%% determine flow from sampled node 1 to boundary between first 
%% N1 nodes and second N2 nodes 

G1 = G_SBM; 
G1(1:N1,(N1+1):N) = 3*G1(1:N1,(N1+1):N) ; 
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

Lambda = ((1./(sum(abs(D),2)))) ; 
Gamma = ((1./(sum(abs(D),1))))'; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Algorithm Initialisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


primSLP = ones(N,1); 
primSLp(N) = 0 ; 
%running_average;
dualSLP =(1:(N-1))'/(N-1) ; %running_averagey; 
dualSLP = zeros(M,1); 


hatx = zeros(N,1); 


running_average = 0*hatx; 
running_averagey = 0*hatx; 
hatx = zeros(N,1); 
haty = ((1:(N-1))/(N-1))'; 
haty = zeros(M,1); 
hatxLP = zeros(N,1); 

log_conv=zeros(N,1); 
log_bound=zeros(N,1); 

for iterk=1:1000
    
    % LP iteration 
    hatxLP = inv_degrees.*(A_undirected*hatxLP); 
    hatxLP(samplingset) = graphsig(samplingset); 
    
    % SLP iteration
    newx = hatx - 0.9*Gamma.*(D'*haty) ; 
    newx(samplingset) = graphsig(samplingset) ; 
    tildex = 2*newx-hatx; 
    newy = haty + 0.99*Lambda.*(D*tildex); 
    haty = newy./max([abs(newy),ones(M,1)],[],2) ; 
    hatx = newx; 
    running_average = (running_average*(iterk-1) +hatx)/iterk; 
    
    
    %running_averagey = (running_average*(iterk-1) +haty)/iterk;+
    %  dual = sign(D*running_average); 
    % dual(iterk:(N-1)) = 0 ;
    % log_conv(iterk) = sum(abs(D*running_average)); 
    % log_bound(iterk) =(1/(2*iterk))*(primSLP'*inv(Gamma)*primSLP)+((dualSLP-dual)'*inv(Lambda)*(dualSLP-dual)) ;
end
if iter_boundary==length(boundary_size_values) 
    for iter_lambda=1:length(lambda_vals)
       
      out_nLasso = nLassoLinModel(G_SBM,1000,lambda_vals(iter_lambda),1,X_mtx,samplingset,graphsig); 
      nLasso_out(iter_lambda,:) = out_nLasso; 
    end
    
end
mse_over_boundary_nLasso(iter_boundary) =mse_over_boundary_nLasso(iter_boundary)+norm(out_nLasso-graphsig)^2;
mse_over_boundary(iter_boundary)=mse_over_boundary(iter_boundary)+norm(running_average-graphsig)^2;
end
end
%figure(1); 
%stem(primSLP);
%title('primal SLP')
%figure(2); 
%stem(dualSLP); 
%title('dual SLP')
figure(2); 
stem(running_average);
title('output SLP');
figure(3); 
stem(out_nLasso); 
title('Network Lasso'); 
figure(4); 
x_vals = boundary_size_values/(avg_degree*length(samplingset)); 
x_vals = sum(ratio_bound_flow,2)/RUNS ; 
mselog = mse_over_boundary/(RUNS) ; 
mselog=mselog/norm(graphsig)^2 ; 
stem(x_vals,mselog); 

mtx=[x_vals mselog]; 
%mtx = flipup(mtx); 
T = array2table(mtx,'VariableNames',{'a','b'});
%csvwrite('hingelosswoheader.csv',mtx);

filename = sprintf('MSEoverBoundary_%s.csv',datetime(now,'ConvertFrom','datenum')) ; 

writetable(T,fullfile(pathtothismfile,filename));


mtx=[(1:N)' nLasso_out' running_average graphsig]; 
idx = [1,10:10:190,200] ; 
mtx = mtx(idx,:) ; 

%mtx = flipup(mtx); 
emptyCell = cell(length(lambda_vals)+3,1) ; 
for iter_lambda=1:length(lambda_vals)  
    emptyCell{iter_lambda+1} = sprintf('lambdaval%c',96+iter_lambda) ; 
end
emptyCell{1} = 'Node' ; 
emptyCell{length(lambda_vals)+2} = 'TVMin' ; 
emptyCell{length(lambda_vals)+3} = 'true' ; 


T = array2table(mtx,'VariableNames',emptyCell);
filename = sprintf('TVMinNLasso_%s.csv',datetime(now,'ConvertFrom','datenum')) ; 
writetable(T,fullfile(pathtothismfile,filename));


%figure(4); 
%stem(dual);
%bound =log_bound+(1./(2*(1:K))*(hatx'*inv(Lambda)*hatx) ; 
%plot(1:N,log([log_conv log_bound])); 



function G_SBM=twocluster(N1,N2,avg_degree,boundary_edges) 

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

end


%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
% generate random geometric graph on unit square [0,1] x [0,1] %%% 
%%%%%%%%%%%%%%%%

function [G,samplingset]=RandGeomGraph(N1,N2)

mean1 =[0,0]; 
mean2 = [4,0]; 
cluster1 = randn(N1,2)*diag([1,0.01])+ones(N1,2)*diag(mean1); % sample 50 nodes from Guassian with mean [0.2,0.2] 
cluster2 = randn(N2,2)*diag([1,0.01])+ones(N2,2)*diag(mean2); % sample 50 nodes from Guassian with mean [0.2,0.2] 
cluster1(1,:) = [mean1]; 
cluster2(N2,:) = [mean2]; 
nodes = [cluster1;cluster2]; 
diff_mtx_x = repmat(nodes(:,1),1,N) ; 
diff_mtx_x = diff_mtx_x - diff_mtx_x'; 
diff_mtx_x = diff_mtx_x.^2 ; 
nodes = [cluster1;cluster2]; 
diff_mtx_y = repmat(nodes(:,2),1,N) ; 
diff_mtx_y = diff_mtx_y - diff_mtx_y'; 
diff_mtx_y = diff_mtx_y.^2 ; 
diff_mtx = diff_mtx_x+diff_mtx_y ; 

[dmy samplingnode1] = min(nodes(1:N1,1)); 
[dmy samplingnode2] = max(nodes((N1+1):N,1)); 
samplingnode2=samplingnode2+N1;
samplingset = [samplingset1 samplingset2]; 


%%%% nearest neighbour graph
radius = 0.3 ; 
G(diff_mtx<radius^2)= 1 ; 
%%%%%%%%%%%

%%%% gaussian filter graph
%sigma = 0.4; 
%G = exp(-diff_mtx/(2*sigma^2)); 

%%%%%%%%%%%
end


