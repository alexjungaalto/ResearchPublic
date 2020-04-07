%%% source code for the experiment discussed in the paper 
%%% "Clustering in Partially Labeled Stochastic Block Models 
%%% via Total Variation Minimization" 

clear all
close all

restoredefaultpath
rehash toolboxcache

[pathtothismfile,name,ext] = fileparts(mfilename('fullpath')) ; 


RUNS = 100; 
RUNS = 10; 



p_out_vals = [0.05, 0.1, 0.15,0.2] ; 
p_out_vals = 1./linspace(1/10,1/0.07,20); 
%p_out_vals=1;

p_in_vals = 0.6*ones(length(p_out_vals),1) ;
p_in_vals = 0.9*ones(length(p_out_vals),1) ;

alpha_vals = [0.1,0.2,0.3]  ; %% ratio of labeled nodes with known cluster assignments 


mse_log = zeros(length(p_in_vals),length(alpha_vals)); 
accuracy = zeros(length(p_in_vals),length(alpha_vals)); 
for iter_alpha=1:length(alpha_vals) 
    alpha = alpha_vals(iter_alpha) ; 
    
for iter_param=1:length(p_in_vals)

p_in = p_in_vals(iter_param);      % edge probability within cluster
p_out = p_out_vals(iter_param) ;    % edge probability between clsuters
nodes_in_cluster = [50;50] ; % cluster sizes 

first_node = cumsum(nodes_in_cluster); 
first_node = [1;first_node(1:(length(first_node)-1))+1] ; 


    for iter_RUNS=1:RUNS
    





[nr_clusters dmy]= size(nodes_in_cluster) ; 
nr_nodes = sum(nodes_in_cluster) ; 

G = rand(nr_nodes,nr_nodes) ;
for iter_cluster=1:nr_clusters 
    for iter_cluster_1=1:nr_clusters 
        idx = first_node(iter_cluster):(first_node(iter_cluster)+nodes_in_cluster(iter_cluster)-1); 
        idx1 = first_node(iter_cluster_1):(first_node(iter_cluster_1)+nodes_in_cluster(iter_cluster_1)-1);
        tmp = G(idx,idx1); 
        threshold = p_out; 
        if iter_cluster== iter_cluster_1 
            threshold = p_in ; 
        end
        
        G(idx,idx1) = tmp < threshold; 
    end
end

G_SBM = G ; 
        
    
%graphsig = 0.1*[ones(N1,1);-1*ones(N2,1)] ;

graphsig = zeros(nr_nodes,1); 

for iter_cluster=1:nr_clusters 
    idx = first_node(iter_cluster):(first_node(iter_cluster)+nodes_in_cluster(iter_cluster)-1); 
    graphsig(idx) = iter_cluster*ones(length(idx),1); 
end


samplingset = [first_node];%first_node+2;first_node+3;first_node+4];  
samplingset = [1:ceil(alpha*nodes_in_cluster(1)),(nodes_in_cluster(1)+((nodes_in_cluster(2)-ceil(alpha*nodes_in_cluster(2))):nodes_in_cluster(2)))];
%samplingset = [first_node;first_node+1;first_node+2;first_node+3;first_node+4];  


Adjac = triu(G_SBM,1) ; 
A_undirected = Adjac+Adjac' ; 
degrees = sparse(sum(A_undirected,1)); 
inv_degrees = 1./degrees';

%%%% create weighted incidence matrix 
G = digraph(triu(G_SBM,1)) ;
D = sparse(incidence(G)') ; 
[M, N] = size(D); 
edge_weights = zeros(M,1); 

Lambda = sparse((1./sparse(sum(abs(D),2)))) ; 
Gamma = sparse((1./sparse(sum(abs(D),1))))'; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Algorithm Initialisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dualSLP = zeros(M,1); 
hatx = zeros(N,1); 
haty = zeros(M,1); 
hatxLP = zeros(N,1); 
running_average = 0*hatx; 

for iterk=1:1000
    
    % LP iteration 
    hatxLP = inv_degrees.*(A_undirected*hatxLP); 
    hatxLP(samplingset) = graphsig(samplingset); 
    
    % SLP iteration
    newx = hatx - Gamma.*(D'*haty) ; 
    newx(samplingset) = graphsig(samplingset) ; 
    tildex = 2*newx-hatx; 
    newy = haty + 0.9*Lambda.*(D*tildex); 
    haty = newy./max([abs(newy),ones(M,1)],[],2) ; 
    hatx = newx; 
    running_average = (running_average*(iterk-1) +hatx)/iterk; 
    
    
    %running_averagey = (running_average*(iterk-1) +haty)/iterk;+
    %  dual = sign(D*running_average); 
    % dual(iterk:(N-1)) = 0 ;
    % log_conv(iterk) = sum(abs(D*running_average)); 
    % log_bound(iterk) =(1/(2*iterk))*(primSLP'*inv(Gamma)*primSLP)+((dualSLP-dual)'*inv(Lambda)*(dualSLP-dual)) ;
end
tmp = running_average-graphsig ; 
accuracy(iter_param,iter_alpha) = accuracy(iter_param,iter_alpha) + length(find(abs(tmp)<1/2))-length(samplingset) ; 
mse_log(iter_param,iter_alpha)=mse_log(iter_param,iter_alpha)+norm(running_average-graphsig)^2;
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
stem(hatxLP); 
title('Label Propagation'); 
figure(4); 

x_vals = p_in_vals./p_out_vals' ; 
mselog = mse_log/(RUNS) ; 
accuracy(:,iter_alpha) = accuracy(:,iter_alpha)/(RUNS*(sum(nodes_in_cluster)-length(samplingset))); 
mselog(:,iter_alpha)=mselog(:,iter_alpha)/norm(graphsig)^2 ; 
%figure(iter_alpha)


%accuracy_vs_alpha(:,iter_alpha) = accuracy; 
%mse_vs_alpha(:,iter_alpha) = mselog; 
end

for iter=1:length(alpha_vals)
    alpha = alpha_vals(iter); 
    S=ceil(alpha*nodes_in_cluster(1)); 
    %stem(x_vals,mselog); 
    mtx=[S*x_vals accuracy(:,iter)]; 
    T = array2table(mtx,'VariableNames',{'a','b'});
  %  filename = sprintf('ACCoverSBMParam_%02d.csv',datetime(now,'ConvertFrom','datenum'),S) ; 
    filename = sprintf('ACCoverSBMParam_%02d.csv',S) ; 
 
    writetable(T,fullfile(pathtothismfile,filename));
end





