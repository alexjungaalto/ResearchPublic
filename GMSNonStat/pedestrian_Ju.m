% source code for the simulations discussed in Sec. 4-B "Pedestrian Coutns"
% of the paper "On the Sample Complexity of Graphical Model Selection from
% Non-Stationary Samples"
% Last Mod.: 09.10.2019

clear all;
close all;

%restoredefaultpath
%rehash toolbox

%%%%%%%%%%%%%%%%%%%%%%% DATA PREPROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%

% this is block length
%L = 3500;

num_blocks = 129;
block_length = 24 ; 

%data = readtable("../TurkuPedestrian.csv");
%all_ids = data{:,2};
%unique_ids = unique(all_ids);

%data = readtable("../TurkuPedestrian.csv",'TextType','string'); 

data = readtable("selectedrows.csv"); 

%all_ids =data(:,1) ;
%unique_ids = unique(all_ids);

nr_counters = 10 ; 

[nr_data_points dmy] = size(data); 

time_series = [data.c73159,data.c73160,data.c73163,data.c73164,data.c73165,data.c73166,data.c73167,data.c73168,data.c73169,data.c73170];  %zeros(nr_data_points,nr_counters); 
%time_offset = size(time_series); 
%valid_series = 0; 
time_series = [time_series 1000*randn(nr_data_points,(nr_counters-5)*2)];



% series_idx= []; 
% 
% for iter_lines=1:nr_counting_lines
%     rows = data.Id==unique_ids{iter_lines,1} ;
%     signal = data.Var5(rows) ;
%     
%     if (length(signal) >= (2*num_blocks*block_length+1)) 
%         valid_series= valid_series+1;
%         series_idx = [series_idx iter_lines]; 
%     end
%     
%     zeitpunkt =datetime(data.To(rows),'InputFormat','yyyy-MM-dd HH:mm:ss') - datetime(2018,6,1); 
%     time_offset(1:length(signal),iter_lines) = seconds(zeitpunkt); 
%     time_series(1:length(signal),iter_lines) = signal;  
% end
% 

%count_lines = unique_ids{series_idx,1};

%time_series = time_series(1:(2*num_blocks*block_length+1),series_idx) ; 
%time_offset = time_offset(1:(2*num_blocks*block_length+1),series_idx) ;

% for each counting station, combine the counts for two directions 
% to an average count representing both directions

sum_directions=zeros(2*nr_counters,nr_counters); 

for iter_i=1:nr_counters 
    
sum_directions(2*(iter_i-1)+1,iter_i) = 1; 
sum_directions(2*iter_i,iter_i) = 1; 
end

time_series = time_series*sum_directions/2; 








% differentiate time series to turn data into stationary process:
% inspired by this:
% https://www.researchgate.net/post/How_to_convert_non-stationery_time_series_into_stationery
% https://www.researchgate.net/post/How_can_I_make_a_time-series_stationary
% samples = diff(time_series);

% remove the 24 hour seasonal component 

 time_steps = (1:(num_blocks*block_length))'; 
 samples = time_series(time_steps,:) ; 
 
 [dmy,dimension] = size(time_series); 
 
% samples = time_series((num_blocks*block_length)+time_steps,:) ; 
% Mdl = arima(2,0,1); 
% EstMdl = estimate(Mdl,train(:,1)) ; 
% 
% hat_y = zeros(length(time_steps),1); 
% 
% for iter_time=1:length(time_steps)
%     hat_y(iter_time) = forecast(EstMdl,1,time_series((iter_time-1)+time_steps,1)); 
% end



% 
% samples = samples - ones(length(time_steps),1)*sum(samples,1)/length(time_steps); 
% for period=23:25
%  sinus_func = exp(-1i*2*pi*(time_steps-1)/period) ; 
%  sinus_func = sinus_func/norm(sinus_func,2);
%  samples = (eye(length(time_steps)) - sinus_func * sinus_func')*samples; 
%  sinus_func = conj(sinus_func); 
%  samples = (eye(length(time_steps)) - sinus_func * sinus_func')*samples; 
% end

Prec_matrix = zeros(dimension,dimension,num_blocks); 

% for iter_block=1:num_blocks 
%     idx_block = ((iter_block-1)*block_length+1):(iter_block*block_length);
%     dmy_samples = fft(samples(idx_block,:));
%     samples(idx_block,:) = dmy_samples ;
%     C = cov(dmy_samples);
%     Prec_matrix(:,:,iter_block) = inv(C);
% end


[rows,dimension]= size(samples); 

new_samples=zeros(rows-24,dimension); 
diff_sig = zeros(rows-24,dimension); 
for iter_line=1:dimension 
    
    sigvec=samples(:,iter_line) ; 
    sigvec=diff((reshape(sigvec,[],rows/24))');   %%% difference at lag 24 
    dmy = sigvec'; 
    dmy = reshape(dmy,[],1) ; 
    diff_sig(:,iter_line) = dmy; 
    
   % sigvec=[sigvec;zeros(1,24)];
    
    dmy = reshape(dmy,[],((num_blocks-1) *block_length)/(128*24)) ; 
    dmy = fft(dmy) ; 
    sigvec = dmy; 
  %  sigvec = fft(sigvec'); 
  %  sigvec=sigvec';
    sigvec=reshape(sigvec,[],1); 
    new_samples(:,iter_line)=sigvec; 
    
 
end

for iter_line=1:dimension
    
 filename = sprintf('Raw_Counts_%d.csv',iter_line) ; 
 x_vals=1:(24*10); 
 T = array2table([x_vals' time_series(x_vals,iter_line)],'VariableNames',{'x','y'});
 writetable(T,filename); 

 filename = sprintf('Diff_Counts_%d.csv',iter_line) ; 
 x_vals=1:(24*10); 
 T = array2table([x_vals' diff_sig(x_vals,iter_line)],'VariableNames',{'x','y'});
 writetable(T,filename); 
end

%num_blocks = num_blocks-1; 
%new_samples = samples(1:24:rows,:); 


[rows,cols] = size(new_samples); 
%new_samples = new_samples - ones(rows,1)*sum(new_samples,1)/rows; 

%new_samples = new_samples*diag(1./sqrt(sum(abs(new_samples).^2,1))) ; 


dimension=cols; 
samples = new_samples; 
%%%%%%%%%%

  % test stationarity by commputing average power over short intervals

% tmp_power = zeros(floor(rows/96),cols); 
% 
% for iter_line=1:length(count_lines) 
%     
%    
%     sigvec=samples(:,iter_line); 
%     
%   
%     sigvec=sum(reshape(sigvec,96,[]).^2,1);
%     tmp_power(:,iter_line) = reshape(sigvec,[],1); 
%     
% end





%%% test_variance 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% According to the geography of the experiment (see explanation slides), the
% following nodes should have a link between each other:
% 1 and 2
% 3 and 4
% ..
% i and i+1, for i = 1,3,5.....


% specify node id to find neigbours
node_i = 5;
% set sparsity here, for example 3
sparsity = 3;

[sample_size dmy] = size(samples); 

block_length=4*sparsity ; 
num_blocks=floor(sample_size/block_length) ; 


%%%%%%%%%%%%%%%%%%%%%%% OUR METHOD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = dimension; 

 for idx_blk=1:num_blocks
            block_range = (idx_blk - 1)*block_length + (1:block_length);
            dmy = samples(block_range, :);
            dmy = dmy - ones(block_length,1)*sum(dmy,1)/block_length; 
            samples(block_range, :) = dmy; 
        %    X_T_b = samples(block_range, setT);
        %    proj_T_b = (X_T_b / (X_T_b' * X_T_b)) * X_T_b';
        %    Z_T_tmp = Z_T_tmp + power(norm(X_i_b - proj_T_b * X_i_b),2);
 end
        
 
rho_min = get_rho_min(Prec_matrix);

% enumerate all the candidates
sets_T = {};
min_val = zeros(sparsity,1); 
min_val_sanity = zeros(sparsity,1); 

for s = 1:1:sparsity
    sets_T = combnk([1:node_i-1 node_i+1:p],s);
    num_set_T = size(sets_T,1);
    Z_T = zeros(num_set_T, 1);
    for iter_set = 1:1:num_set_T
        setT = sets_T(iter_set,:);
     %   Z_T(iter_set) = calc_stats_T(samples, node_i, setT, block_length) ;
        Z_T_tmp = 0;
        for idx_blk=1:num_blocks
            block_range = (idx_blk - 1)*block_length + (1:block_length);
            X_i_b = samples(block_range, node_i);
            X_T_b = samples(block_range, setT);
            proj_T_b = (X_T_b / (X_T_b' * X_T_b)) * X_T_b';
            Z_T_tmp = Z_T_tmp + power(norm(X_i_b - proj_T_b * X_i_b),2);
        end

        Z_T(iter_set) = Z_T_tmp/(block_length*num_blocks);
    end
    [min_val(s), minIdx] = min(Z_T);
  % estimated_neigbours_ours = count_lines(2*sets_T(minIdx,:))
    estimated_neigbours_ours = sets_T(minIdx,:)
    
    setT = sets_T(minIdx,:);
    setT = [setT,7]; 
    %%%% sanity check by adding completely random fake measurement
    Z_T_tmp = 0;
    for idx_blk=1:num_blocks
       block_range = (idx_blk - 1)*block_length + (1:block_length);
       X_i_b = samples(block_range, node_i);
       X_T_b = samples(block_range, setT);
       proj_T_b = (X_T_b / (X_T_b' * X_T_b)) * X_T_b';
       Z_T_tmp = Z_T_tmp + power(norm(X_i_b - proj_T_b * X_i_b),2);
     end

       
     min_val_sanity(s) =  Z_T_tmp/(block_length*num_blocks);
 %  estimated_neigbours_ours = estimated_neigbours_ours{1}     
    
end

%%% projection error when using to regressor at all (results in variance of
%%% tiem series)
min_val = [sum(abs(samples(:,node_i)).^2)/(block_length*num_blocks); min_val];

setT = [7]; 
    %%%% sanity check by adding completely random fake measurement
     Z_T_tmp = 0;
        for idx_blk=1:num_blocks
            block_range = (idx_blk - 1)*block_length + (1:block_length);
            X_i_b = samples(block_range, node_i);
            X_T_b = samples(block_range, setT);
            proj_T_b = (X_T_b / (X_T_b' * X_T_b)) * X_T_b';
            Z_T_tmp = Z_T_tmp + power(norm(X_i_b - proj_T_b * X_i_b),2);
        end

min_val_sanity=[Z_T_tmp/(block_length*num_blocks);min_val_sanity];

diff_min_val = diff(min_val); 
tmp =(min_val-min_val_sanity); 
significance_test = -diff_min_val./tmp(1:length(diff_min_val)) ; 

filename = sprintf('min_val_i_%d.csv',node_i) ; 
T = array2table([(1:length(min_val))' min_val/max(min_val)],'VariableNames',{'x','y'});
writetable(T,filename);
% % allocate some memory
% sets_T = sets_T';
% num_set_T = size(sets_T,1);
% Z_T = zeros(num_set_T, 1);
% % calcualte statistics
% for i=1:num_set_T
%     setT = sets_T{i};
%     card = numel(setT);
%     Z_T(i) = calc_stats_T(samples, node_i, setT, L) + (rho_min/6)*card;
% end

% find the smallest score and corresponding neigbour set
%[~, minIdx] = min(Z_T);
%estimated_neigbours_ours = sets_T(minIdx,:);
%estimated_neigbours_ours = estimated_neigbours_ours{1}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%% LASSO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set reguralizer
%lambda = 300;
%estimated_neigbours_lasso = lasso(node_i, samples(:,node_i),samples(:,[1:node_i-1 node_i+1:end]), block_length, lambda)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%