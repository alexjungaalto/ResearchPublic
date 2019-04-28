function [ adj_chain ] = adj_chain_graph( dim )
%=============================================
%ADJ_CHAIN_GRAPH 
%   Return the adjacent matrix for chain graph
%=============================================
adj_chain =  eye(dim);
for idx=1:dim-1
    adj_chain(idx,idx+1) = 1;
    adj_chain(idx + 1,idx) = 1;
end
end

