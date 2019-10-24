%%% source code for the experiment discussed in Sec. VI-B of the paper 
%%% "Semi-supervised Learning in Network-Structured Data via 
%%%   Total Variation Minimization"

clear all
close all

restoredefaultpath
rehash toolboxcache

[pathtothismfile,name,ext] = fileparts(mfilename('fullpath')) ; 


RUNS = 100; 
RUNS = 100; 



p_out_vals = [0.05, 0.1, 0.15,0.2] ; 
p_out_vals = 1./linspace(1/10,1/0.07,20); 
p_in_vals = 0.6*ones(length(p_out_vals),1) ;
p_in_vals = 0.9*ones(length(p_out_vals),1) ;

mse_log = zeros(length(p_in_vals),1); 

for iter_param=1:length(p_in_vals)

p_in = p_in_vals(iter_param);      % edge probability within cluster
p_out = p_out_vals(iter_param) ;    % edge probability between clsuters
nodes_in_cluster = [10;10;10] ; % cluster sizes 

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
    graphsig(idx) = (1/100)*((-1)^iter_cluster)*ones(length(idx),1); 
end


samplingset = [first_node;first_node+1;first_node+2;first_node+3;first_node+4];  


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

mse_log(iter_param)=mse_log(iter_param)+norm(running_average-graphsig)^2;
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
mselog=mselog/norm(graphsig)^2 ; 
stem(x_vals,mselog); 

mtx=[x_vals mselog]; 
%mtx = flipup(mtx); 
T = array2table(mtx,'VariableNames',{'a','b'});
%csvwrite('hingelosswoheader.csv',mtx);

filename = sprintf('MSEoverSBMParam_%s.csv',datetime(now,'ConvertFrom','datenum')) ; 

writetable(T,fullfile(pathtothismfile,filename));




