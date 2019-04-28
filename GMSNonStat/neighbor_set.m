function [ neighbors ] = neighbor_set( node_i, graph_type, dim )
%=============================================
%NEIGHBOR_SET 
%   Return the set of neighbors of node_i if we know the graph_type
%   graph_type = 1 for chain graph and 2 for grid graph.
%=============================================
if graph_type == 1
    if node_i == 1
        neighbors = [node_i + 1];
    elseif node_i == dim
        neighbors = [node_i - 1];
    else
        neighbors = [node_i - 1, node_i + 1];
    end
elseif graph_type == 2
    grid_size = round(sqrt(dim));
    if node_i == 1
        neighbors = [node_i + 1, node_i + grid_size];
    elseif node_i == grid_size
        neighbors = [node_i - 1, node_i + grid_size];
    elseif node_i == dim
        neighbors = [node_i - 1, node_i - grid_size];
    elseif node_i == dim - grid_size + 1
        neighbors = [node_i + 1, node_i - grid_size];
    elseif node_i < grid_size
        neighbors = [node_i + 1, node_i - 1, node_i + grid_size];
    elseif node_i > dim - grid_size
        neighbors = [node_i + 1, node_i - 1, node_i - grid_size];
    elseif mod(node_i, grid_size) ==  1
        neighbors = [node_i - grid_size, node_i + 1, node_i + grid_size];
    elseif mod(node_i, grid_size) ==  0
        neighbors = [node_i - grid_size, node_i - 1, node_i + grid_size];
    else 
        neighbors = [node_i - grid_size, node_i - 1, node_i + 1, node_i + grid_size];
    end
end

end

