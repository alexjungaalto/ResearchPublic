close all;
clear all;

% paper requires min eigenvalues to be at least 1.0
TARGET_MIN_EIG = 1.1;
betaRatio = 1.5; % controls the factor at which eigenvalues of first and last block differ
c = 1.1; % this paramets controls scale of entries in the precision matrix
dimension = 10;
num_blocks = 4;
num_sims = 10;

%
sparsity = dimension-1;

total_err = [];
for block_length = ceil(linspace(100,500,10))

    success_count = 0;
    total_count = 0;
    for sim = 1:num_sims
        
        % generate precision and covariance matrices
        K =  zeros(dimension, dimension, num_blocks);
        C = zeros(dimension, dimension, num_blocks);

        while 1
            % generate random adjacency matrix which corresponds to random
            % graph
            A = randi(2,dimension,dimension) - 1;
            A = A - tril(A,-1) + triu(A,1)';
            A(1:1+size(A,1):end) = 1;
            
            C_0 = ones(dimension,dimension) + c*eye(dimension);
            precision_matrix = C_0.*A;
            
            % ensure that matrix is a PSD one since
            % not every random matrix is PSD one
            if(all(eigs(precision_matrix,dimension) >= 0)) 
                break;
            end
        end

        cov_matrix = inv(precision_matrix);
        cov_matrix = TARGET_MIN_EIG*cov_matrix / min(eig(cov_matrix)); % scale covariance matrix
        for i = 1:1:num_blocks
            C(:,:,i) = (1 + (betaRatio-1)*(i-1)/num_blocks) * cov_matrix;
            precision_matrix = inv(C(:,:,i));
            %eigs(precision_matrix)
            precision_matrix(precision_matrix < 1e-10) = 0; % due to the numerical errors this may be slightly greater than zero
            K(:,:,i) = precision_matrix;
        end

        samples = sample_generator(C, block_length); % create samples from covariance matrices

        % get rho_min and compute reguralizer
        rho_min = get_rho_min(K);
        lambda = rho_min/6;

        for node_i = 1:dimension
            % find true neigbours
            row = precision_matrix(node_i,:);
            row(node_i) = 0;
            neighbors = find(row);

            % enumerate all the candidates
            sets_T = {};
            for s = 1:1:sparsity
                sets = combnk([1:node_i-1 node_i+1:dimension],s);
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
                Z_T(i) = calc_stats_T(samples, node_i, setT, block_length) + lambda*card;
            end

            % find the smallest score and corresponding neigbour set
            [~, minIdx] = min(Z_T);
            estimatedNeigbours = sets_T(minIdx,:);
            estimatedNeigbours = estimatedNeigbours{1};

            % compare estiamted neigbours to true neighbours
            success = isequal(estimatedNeigbours,neighbors);
            success_count = success_count + success;
            total_count  = total_count + 1;
        end

    end
total_err = [total_err; block_length * num_blocks 1 - success_count / total_count];
end



