%%%% This is the source code for num. experiments discussed in the paper 
%%%% "Localized Graph Clustering with Network Lasso" 
%%%% A. Jung, April 2020 

clear all; 
[pathtothismfile,name,ext] = fileparts(mfilename('fullpath')) ; 

N = 10;  % nr of nodes in chain 
K = 100; %nr of iterations used for nLasso


B = diag(ones(N,1),0) - diag(ones(N-1,1),1) ; % incidence matrix 
B = B(1:(N-1),:) ; 


%kappa = 0.5 ; 



Lambda = 1*diag(1./(sum(abs(B),2))) ; 
Gamma_vec = (1./(sum(abs(B),1)))' ; 
Gamma = diag(Gamma_vec) ; 

cluster1 = [ones(N/2,1);zeros(N/2,1)]; 
cluster2 = [zeros(N/2,1);ones(N/2,1)]; 
c1 = 1 ; c2 = 0; 
graphsig = c1*cluster1 + c2*cluster2; 



%% here we can loop over different values for boundary weight 
%%% here we can loop over different values of nLasso lambda 

weight = eye(N-1,N-1); 
%weight(N/2,N/2) = 1/4 ; 
weight_vec = 1./(1:(N-1)) ; 
weight_vec = (5/4)*ones(N-1,1); 
eta = 1/4; 
eta=1; 
weight_vec(2) = eta ; 
weight =  diag(weight_vec) ; 
lambda_nLasso = 1/3 ;  %% nLasso parameter
lambda_nLasso = 1/3 ; 

anal_sig = (c1 - lambda_nLasso*eta)*cluster1 + (c2 + lambda_nLasso*eta)*cluster2 ; 


D = weight*B;

hatx = zeros(N,1); 
haty = ((1:(N-1))/(N-1))'; 
samplingset = [1 N]; 
running_average = 0*hatx; 
running_averagey = 0*haty; 

% 
% for iterk=1:K
%    
%     newx = hatx - Gamma*D'*haty ; 
%     newx(samplingset) = graphsig(samplingset) ; 
%     tildex = 2*newx-hatx; 
%     newy = haty + Lambda*D*tildex; 
%     haty = newy./max([abs(newy),ones(N-1,1)],[],2) ; 
%     hatx = newx; 
%     running_average = (running_average*(iterk-1) +hatx)/iterk; 
%     running_averagey = (running_averagey*(iterk-1) +haty)/iterk; 
%  %   log_conv(iterk) = sum(abs(D*running_average)); 
%  %   log_bound(iterk) =1/(2*iterk)*(hatx'*inv(Gamma)*hatx)+(haty'*inv(Lambda)*haty) ;
% end

primSLP = ones(N,1); 
primSLp(N) = 0 ; 
%running_average;
dualSLP =(1:(N-1))'/(N-1) ; %running_averagey; 


hatx = zeros(N,1); 

%% create sampling set = set of seed nodes 
samplingset = [1]; 

%% create indicator signal for nodes not being seed nodes  
seednodesindicator= zeros(N,1); 
seednodesindicator(samplingset) = 1; 
noseednodeindicator = ones(N,1); 
noseednodeindicator(samplingset)=0; 

running_average = 0*hatx; 
running_averagey = 0*hatx; 


log_conv=zeros(N,1); 
log_bound=zeros(N,1); 
newx = 0*hatx; 

hist_y = zeros(K,N-1); 
hist_x = zeros(K,N) ; 

alpha = 1/10; 
fac_alpha = 1./(Gamma_vec*alpha+1) ; 


for iterk=1:K
    % update (????)
    tildex = 2*newx-hatx;  
    % update (???) 
    newy = haty + (1/2)*B*tildex; 
    % update (???) 
    haty = newy./max([abs(newy)./(lambda_nLasso*weight_vec),ones(N-1,1)],[],2) ; 
    % update (???) 
    newx = hatx - Gamma_vec.*(B'*haty); 
    % update (???)
    for dmy=1:length(samplingset)
        idx_dmy = samplingset(dmy); 
        newx(idx_dmy) = (newx(idx_dmy) +Gamma_vec(idx_dmy)*graphsig(idx_dmy))/(1+Gamma_vec(idx_dmy));
    end
    newx = seednodesindicator.*newx + noseednodeindicator.*(newx.*fac_alpha); 
    hatx = newx; 
    
   % newx(samplingset) = graphsig(samplingset) ; 
 
  %  haty = newy./max([abs(newy),ones(N-1,1)],[],2) ; 
  
    
    hist_y(iterk,:) = haty; 
    hist_x(iterk,:) = hatx; 
    
    running_average = (running_average*(iterk-1) +hatx)/iterk; 
    %running_averagey = (running_average*(iterk-1) +haty)/iterk;+
    dual = sign(B*running_average); 
    dual(iterk:(N-1)) = 0 ;
    log_conv(iterk) = sum(abs(B*running_average)); 
    log_bound(iterk) =(1/(2*iterk))*(primSLP'*inv(Gamma)*primSLP)+((dualSLP-dual)'*inv(Lambda)*(dualSLP-dual)) ;
end

close all; 
figure(1); 
stem(hatx);
title('primal iterate')
figure(2); 
stem(haty); 
title('dual iterate')
figure(3); 
stem(running_average);
title('ouput SLP'); 
anal_est =  norm(anal_sig-graphsig)^2/norm(graphsig)^2 
err = norm(running_average-graphsig)^2/norm(graphsig)^2 


mtx=[(1:N)' running_average]; 
%mtx = flipup(mtx); 
T = array2table(mtx,'VariableNames',{'i','x'});
%csvwrite('hingelosswoheader.csv',mtx);


filename = sprintf('NumExpChainPrimal%s.csv',datetime(now,'ConvertFrom','datenum')) ; 
%filename = 'NumExpChainPrimal.csv' ; 

writetable(T,fullfile(pathtothismfile,filename));

mtx=[(1:(N-1))' haty]; 
%mtx = flipup(mtx); 
T = array2table(mtx,'VariableNames',{'i','y'});
%csvwrite('hingelosswoheader.csv',mtx);

filename = sprintf('NumExpChainDual%s.csv',datetime(now,'ConvertFrom','datenum')) ; 

writetable(T,fullfile(pathtothismfile,filename));

mtx=[(1:(N))' graphsig]; 
%mtx = flipup(mtx); 
T = array2table(mtx,'VariableNames',{'i','x'});
%csvwrite('hingelosswoheader.csv',mtx);

filename = sprintf('NumExpChainSig%s.csv',datetime(now,'ConvertFrom','datenum')) ; 

writetable(T,fullfile(pathtothismfile,filename));


%figure(4); 
%stem(dual);
%bound =log_bound+(1./(2*(1:K))*(hatx'*inv(Lambda)*hatx) ; 
%plot(1:N,log([log_conv log_bound])); 