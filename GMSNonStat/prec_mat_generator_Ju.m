function [ precision_matrix ] = prec_mat_generator(dimension, graph_type, c)
%=============================================
%PREC_MAT_GENERATOR Summary of this function goes here
%   Create a precision matrix for predfined graph
%=============================================
% c allows to control size of partial correlations

%while 1
    %Cov_mtx = randn(dimension,dimension);
    %C_0 = Cov_mtx*Cov_mtx';   
    C_0 = (1/3)*ones(dimension,dimension) + eye(dimension);
    
    if graph_type == 1
        adj_matrix = adj_chain_graph(dimension) ;
    else
        adj_matrix = adj_grid_graph(dimension) ;
    end
    
    precision_matrix = C_0.*adj_matrix;
   % if(all(eigs(precision_matrix) >= 0)) %ensure that matrix is a PSD one
   %     break;
   % end
%end

end

