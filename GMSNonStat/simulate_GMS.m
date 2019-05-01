close all;
clear all;

%% Parameters
graph_type = 1; % type 1: chain graph
num_blocks = 4;
num_tests = 300; % number of simulations for each point
num_tests  = 10; 

% specify range here
length = ceil((2:20:442));

C = {};
for dimension = [8 16 32]
    stat = [];
    for  block_length =  length
        totalAccuracy = 0;
        for idx=1:num_tests
            accuracy = test_statistic(dimension, block_length, num_blocks, graph_type);
            totalAccuracy  = totalAccuracy + accuracy;
        end
        stat = [stat; num_blocks*block_length totalAccuracy];
    end
    stat(:,2) = stat(:,2)/num_tests;
    C{end+1} = {dimension stat}
end




% plot the results

figure
hold on
mtx =[]; 
for c = C
    dim = c{1}{1};
    data = c{1}{2};
    mtx = [mtx (1-data(:,2))]; 
    plot(data(:,1),1-data(:,2),'-o','DisplayName',num2str(dim));
end
hold off
xlabel('Sample size');
ylabel('Error rate')
legend
grid on

x_vals = data(:,1); 
x_vals = x_vals/max(x_vals) ; 
mtx=[x_vals mtx]; 
T = array2table(mtx,'VariableNames',{'x','y1','y2','y3'});
%csvwrite('hingelosswoheader.csv',mtx);
writetable(T,'ErrorVsSampeSize.csv');

figure
hold on
for c = C
    dim = c{1}{1};
    data = c{1}{2};
    plot(data(:,1)/log(dim),1-data(:,2),'-o','DisplayName',num2str(dim));
end
hold off
xlabel('Scaled sample size, N/log(p)');
ylabel('Error rate')
legend
grid on

