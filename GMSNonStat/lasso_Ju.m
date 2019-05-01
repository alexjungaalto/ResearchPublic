function neigbours = lasso_Ju(node_i ,x_i, samples, L, lambda,rho_min)
% node_i - index of the node we want to approximate
% x_i - samples of the node we want to approximate
% samples - samples of the nodes, using which we want to
% approximate x_i
% L - block length

p = size(samples,2);
B = size(samples,1)/L;

X_b = cell(B,1);
for i = 1:1:B
    X_b{i} = sparse(samples(1+(i-1)*L : i*L,:));
end
X_1 = blkdiag(X_b{:});

X = cell(p,1);
for j = 1:1:p
    X{j} =  X_1(:,j:p:end);
end

X_V  =[X{:}];

% set regularization and run lasso
partition = B*ones(p, 1);
% lasso outputs regression coeffitients
[coeff, ~] = group_lasso(X_V, x_i,lambda, partition, 1, 1);

% for finding regularization proper lambda
% you may want to plot regression coefficents
% plot(coeff)

% calculate L2-norms of obtained coefficients on the block level
% if block is forced to zero, then its norm would be zero too
scores = zeros(1,p);
for j = 1:1:p
   scores(j) = norm(coeff(1+(j-1)*B: j*B))^2;
end

% estimated neigbours correspond to nodes with non-zero norms
neigbours = find(scores>rho_min/2);

% since we extracted original node, apply index correction
neigbours(neigbours >= node_i) = neigbours(neigbours >= node_i)+1;
end

