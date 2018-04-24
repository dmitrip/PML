function F0 = ML_unseen_symbols_uniform(n, T)
% maximum likelihood estimate for number of unseen symbols given empirical
% support n samples, set size T
%
% created by Dmitri S. Pavlichin on October 2, 2017
%
% Matlab version: R2015a
%
% Args:
%     * n (int) - sample size
%     * T (int) - empirical support set size
%
% Returns:
%     * F0 (int) - maximum likelihood estimator for number of unseen symbols

% check for F0 == 0
if n >= T*harmonic_number(T)
    F0 = 0;
else
    % check for F0 == Inf
    if T >= n
        F0 = Inf;
    else
        fun = @(f0) gammaln(f0+T+1)-gammaln(f0+1)-n*log(f0+T);
        F0 = maximize_on_interval_int(fun, 0, ceil((T^2-n)/(n-T))+1);
    end
end

