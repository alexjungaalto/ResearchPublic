clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%% DATA PREPROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%

% this is block length
L = 3500;
data = readtable("../TurkuPedestrian.csv");
all_ids = data{:,1};
unique_ids = unique(all_ids);

nodes_samples = [];
for id = unique_ids(:)'
    indices = find(id == all_ids);
    
    % discard nodes which have too few samples
    if (numel(data{indices,5}) < L)
        continue
    end
    
    % differentiate time series to turn data into stationary process:
    % inspired by this:
    % https://www.researchgate.net/post/How_to_convert_non-stationery_time_series_into_stationery
    % https://www.researchgate.net/post/How_can_I_make_a_time-series_stationary
    samples = diff(data{indices,5});
    nodes_samples = [nodes_samples samples(1:L)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% According to the geography of the experiment (see explanation slides), the
% following nodes should have a link between each other:
% 1 and 2
% 3 and 4
% ..
% i and i+1, for i = 1,3,5.....


% specify node id to find neigbours
node_i = 4;
% set sparsity here, for example 3
s = 3;

%%%%%%%%%%%%%%%%%%%%%%% OUR METHOD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = size(nodes_samples,2);

% estimate rho_min and lambda
% maybe there is better way?
C = cov(nodes_samples);
K = inv(C);
rho_min = get_rho_min(inv(C));

% enumerate all the candidates
sets_T = {};
for s = 1:1:s
    sets = combnk([1:node_i-1 node_i+1:p],s);
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
    Z_T(i) = calc_stats_T(nodes_samples, C, node_i, setT, L) + (rho_min/6)*card;
end

% find the smallest score and corresponding neigbour set
[~, minIdx] = min(Z_T);
estimated_neigbours_ours = sets_T(minIdx,:);
estimated_neigbours_ours = estimated_neigbours_ours{1}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%% LASSO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set reguralizer
lambda = 300;
estimated_neigbours_lasso = lasso(node_i, nodes_samples(:,node_i),nodes_samples(:,[1:node_i-1 node_i+1:end]), L, lambda)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%