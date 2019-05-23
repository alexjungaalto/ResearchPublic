%function [Z, history] = covsel(S_func, lambda, rho, alpha)



clear all;
close all ;



alpha = 1  ;
p = 64;
s =4 ;
rho_min = 0.4/s;
%N = 512 ;
MAX_ITER = 10;
close all ;

t_start = tic;

%% Global constants and defaults

QUIET    = 0;

ABSTOL   = 1e-4;
RELTOL   = 1e-2;

valsN=200:50:600 ;
valsN=128;
dec_fac = valsN(1)/4 ;
Pdvals = zeros(length(valsN),1) ;
rho = 100;
RUNS = 100 ;

lambdavals= linspace(0.5,3.0,10)*100*rho_min/s  ;
%lambdavals = linspace(0.5,1.5,10)*100*rho_min/s ;

%lambdavals=0 ;

PdROC = zeros(length(lambdavals),1) ;
PfaROC = zeros(length(lambdavals),1) ;


maxeigvaloptvar = 1 ;

%lambda = 3/s;
%lambda = 1.2/s;

k_max=2 ; %impulse response length of innovations filter

prec_mtx = 0.5*eye(p,p) ;
idx = randperm(p-1) ;
prec_mtx(1+idx(1:s),1) = rho_min ;
prec_mtx(1,1+idx(1:s)) = rho_min;
CIG = zeros(p,p) ;
CIG(1+idx(1:s),1) = 1 ;
CIG(1,1+idx(1:s)) = 1 ;

cov_matrix = inv(prec_mtx);

nr_correct_det = zeros(length(lambdavals),1) ;
nr_false_alarm = zeros(length(lambdavals),1) ;
success = zeros(length(lambdavals),1) ;


for iter_RUNS=1:RUNS
    
    %impulse_resp = randn(1,floor(k_max/2)+1)+3 ;
    %impulse_resp(1) = 3*sum(abs(impulse_resp(2:length(impulse_resp)))) ;
    %impulse_resp = impulse_resp/impulse_resp(1) ;
    
    impulse_resp = zeros(2,1) ;
    impulse_resp(1) = 2 ;
    impulse_resp(2) = min(abs(randn(1,1)*0.2),1 ) ;
    
    %impulse_resp = impulse_resp/norm(impulse_resp) ;
    
    for iterlambda=1:length(lambdavals) ;
        
        lambda = lambdavals(iterlambda)/s ;
        for iterN=1:length(valsN)
            N = valsN(iterN) ;
            nrfreqs=N/dec_fac;
            
            
            X = sqrtm(cov_matrix)*randn(p,N) ;
            
            emp_cov = (1/N)*(X*X') ;
            
            S_func = zeros(p,p,nrfreqs);
            
            inv_emp_cov = inv(emp_cov) ;
            
            mineigval = min(eig(cov_matrix)) ;
            %cov_matrix = 2*cov_matrix/mineigval;
            prec_matrix = inv(cov_matrix) ;
            rho_min = max(abs(prec_matrix(1,2:p))) ; % max can be used only if all non-zero values have same magnitude !!!
            
            
            for freqbin=1:nrfreqs
                
                S_func(:,:,freqbin) = emp_cov ;
                
            end
            process_in = sqrtm(cov_matrix)*randn(p,N+2*length(impulse_resp));
            process = zeros(p,N) ;
            
            for iter_p = 1:p
                %   in_dummy = [process(iter_p,:) zeros(1,2*length(impulse_resp))] ;
                out_dummy = conv(process_in(iter_p,:),impulse_resp) ;
                process(iter_p,:) = out_dummy((length(impulse_resp)+1):(N+length(impulse_resp))) ;
            end
            
            
            
            padded_observation = [process.'; zeros(k_max,p)] ;
            shifted_samples = padded_observation;
            
            ACF_hat = [] ;
            %%%% estimate the autocovariance function ACF_{x}[m] of x[n] up to m =
            %%%% k_max
            len_filter=length(impulse_resp) ;
            autoconvfun = xcorr(impulse_resp) ;
            
            k_max = 1 ;
            ACF = zeros (p,p,k_max+1) ;
            
            ACF(:,:,1) = padded_observation'*shifted_samples ;
            ACF(:,:,1) = cov_matrix*autoconvfun(len_filter) ;
            
            
            for iter_lag = 1:k_max
                shifted_samples = circshift(shifted_samples,[1 0]) ;
                ACF(:,:,iter_lag+1) = padded_observation'*shifted_samples ;
                ACF(:,:,iter_lag+1) = cov_matrix*autoconvfun(len_filter+iter_lag) ; %padded_observation'*shifted_samples ;
            end
            
            %% compute estimate of SDM via plain old periodogramm ---> yields an estimate which is psd (cf. Brillinger (proof of) Thm. 2.5.1. p. 393)
            
            SDM = zeros(p,p,N) ;
            hat_SDM = zeros(p,p,N) ;
            periog_dummy = fft([process zeros(p,N)].') ;
            for iter_freq=1:N
                hat_SDM(:,:,iter_freq) = (1/(N))*conj(periog_dummy(iter_freq,:)'*periog_dummy(iter_freq,:));
            end
            
            %    smoothing_lag_window = exp(-((1:(2*N))'-floor((2*N)/2)).^2/(2*((2*N)/k_max)^2)) ;
            %    smoothing_lag_window = exp(-((1:(2*N))'-floor((2*N)/2)).^2/(10*((2*N)/k_max)^2)) ;
            smoothing_lag_window = abs(fft(exp(-((-N+1):N).^2)))';
            
            %  periog_dummy = sum(smoothing_lag_window);
            
            for iter_p=1:p
                for iter_p_prime=1:p
                    smoothed_dummy =  cconv(hat_SDM(iter_p,iter_p_prime,:) ,smoothing_lag_window,2*N) ;
                    hat_SDM(iter_p,iter_p_prime,:) = smoothed_dummy(1:2:2*N) ;
                end
            end
            
            hat_SDM = hat_SDM / sum (smoothing_lag_window);
            
            for iter_freq=1:nrfreqs
                freqbin=(iter_freq-1)*dec_fac+1 ;
                
                S_func(:,:,iter_freq) = hat_SDM(:,:,freqbin) ;
                
                
                %  S_func(:,:,iter_freq) = cov_matrix ; %hat_SDM(:,:,freqbin) ;
            end
            
            
            
            %% Data preprocessing
            
            
            
            dimes = size(S_func) ;
            process_dim = dimes(1);
            n = process_dim;
            %nrfreqs = dimes(3) ;
            
            %% ADMM solver
            
            X_func = zeros(dimes);
            Z_func = zeros(dimes);
            U_func = zeros(dimes);
            
            X = zeros(n);
            Z = zeros(n);
            U = zeros(n);
            
            if ~QUIET
                fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', ...
                    'r norm', 'eps pri', 's norm', 'eps dual', 'objective');
            end
            
            for k = 1:MAX_ITER
                
                
                % x-update
                
                for iter_freq=1:nrfreqs
                    Z = Z_func(:,:,iter_freq) ;
                    U = U_func(:,:,iter_freq) ;
                    S = S_func(:,:,iter_freq) ;
                    
                    [Q,L] = eig(rho*(Z - U) - S);
                    es = diag(L);
                    xi = (es + sqrt(es.^2 + 4*rho))./(2*rho);
                    
                    xi( xi >= maxeigvaloptvar) = maxeigvaloptvar;
                    
                    X_func(:,:,iter_freq) = Q*diag(xi)*Q';
                end
                % z-update with relaxation
                
                Zold_func = Z_func;
                X_hat_func = alpha*X_func + (1 - alpha)*Zold_func;
                for iter_col=1:process_dim
                    for iter_row =1:process_dim
                        % if (iter_col==iter_row)
                        %    Z_func(iter_row,iter_col,:) =    X_hat_func(iter_row,iter_col,:) ;
                        % else
                        
                        Z_func(iter_row,iter_col,:) = max(1-sqrt(nrfreqs)*lambda/(rho*norm(squeeze(X_hat_func(iter_row,iter_col,:) + U_func(iter_row,iter_col,:)))),0)*(X_hat_func(iter_row,iter_col,:) + U_func(iter_row,iter_col,:)) ;
                        % end
                    end
                end
                
                U_func = U_func + (X_hat_func - Z_func);
                
                % diagnostics, reporting, termination checks
                
                history.objval(k)  = 0;% objective(S, X, Z, lambda);
                
                history.r_norm(k)  = norm(squeeze(reshape(X_func - Z_func,dimes(1)*dimes(2)*nrfreqs,1,1)), 'fro');
                history.s_norm(k)  = norm(-rho*squeeze(reshape(Zold_func - Z_func,dimes(1)*dimes(2)*nrfreqs,1,1)),'fro');
                
                history.eps_pri(k) = 0;%sqrt(n*n)*ABSTOL + RELTOL*max(norm(X_func,'fro'), norm(Z_func,'fro'));
                history.eps_dual(k)= 0;%sqrt(n*n)*ABSTOL + RELTOL*norm(rho*U_func,'fro');
                
                
                if ~QUIET
                    fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n', k, ...
                        history.r_norm(k), history.eps_pri(k), ...
                        history.s_norm(k), history.eps_dual(k), history.objval(k));
                end
                
                if (history.r_norm(k) < history.eps_pri(k) && ...
                        history.s_norm(k) < history.eps_dual(k))
                    break;
                end
            end
            
            
            out_mtx = zeros(p,p) ;
            for iter_i=1:p
                for iter_j=1:p
                    % if (iter_i<=iter_j)
                    out_mtx(iter_i,iter_j) = norm(squeeze(X_func(iter_i,iter_j,:)), 'fro') ;
                    % end
                end
            end
            out_mtx = out_mtx/(sqrt(nrfreqs)*norm(impulse_resp)) ;
            
            out_offdiag = out_mtx.*(ones(p,p)-eye(p)) ;
            max_val = max(max(out_offdiag)) ;
            out_graph = zeros(p,p) ;
            thresholdidx = find(out_offdiag>=rho_min/100) ;
            out_graph (thresholdidx) = 1 ;
            error_pattern = CIG+2*out_graph ;
            correct_det =  length(find(abs(error_pattern-3)<0.1)) ;
            false_al = length(find(abs(error_pattern-2)<0.1)) ;
            nr_correct_det(iterlambda) = nr_correct_det(iterlambda)+ correct_det;
            nr_false_alarm(iterlambda) = nr_false_alarm(iterlambda)+ false_al ;
            if correct_det==2*s
                if false_al<1
                    success(iterlambda) = success(iterlambda)+1 ;
                end
            end
            
        end
    end
end
nr_correct_det=nr_correct_det/RUNS ;
nr_false_alarm= nr_false_alarm/RUNS;


figure(1) ;
subplot(1,3,1) ;
imagesc(abs(prec_mtx)) ;
subplot(1,3,2) ;
imagesc(out_mtx) ;
subplot(1,3,3) ;
%plot([out_mtx(1,:)' abs(inv_emp_cov(:,1))]) ;
plot([out_mtx(1,:)']) ;
hold on;
stem(1+idx(1:s),ones(s,1)*0.3,'r.','MarkerSize',20)
legend('gLASSO','true edges') ;
datumuzeit = clock;
filename= sprintf('CSout%4d-%02d-%02d-%02d_%02d_N_%3d-p_%3d', datumuzeit(1), datumuzeit(2), datumuzeit(3),datumuzeit(4),datumuzeit(5),N,p);
save(filename) ;

if ~QUIET
    toc(t_start);
end

%end


