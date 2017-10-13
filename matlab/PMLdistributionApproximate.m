function [p_approx,F0,V_approx,have_continuous_part] = PMLdistributionApproximate(p_empirical,K)
% approximates pattern maximum likelihood (PML) distribution p_approx
% where p_approx \approx \arg \max_p P_p(unordered p_empirical)
% Approximation sums over all permutations that mix within blocks of
% constant p_approx
%
% created by Dmitri S. Pavlichin on May 14, 2017
%
% Matlab version: R2015a
%
% Args:
%     * p_empirical - (integer-valued vector) empirical histogram, entries sum to sample size
% Optional args:
%     * K - (integer) assumed support set size, must have K >= sum(p > 0).  
%         If K is not provided, then we attempt to estimate the support set size
%
% Returns:
%     * p_approx - (column vector) approximate PML distribution, sorted in
%         descending order
%     * F0 - (integer) estimate for number of symbols seen 0 times in
%         p_empirical.  If K is provided as an argument, then F0 = max(K - sum(p > 0), 0)
%     * V_approx - (double) approximate value of log(P_p(unordered p_empirical))
%         achieved by p = p_approx
%
% Example:
%     Optimize support set size:
%         >> p_approx = PMLdistributionApproximate([9 3 2 1 1])
%         yields p_approx == [9/16 7/80 7/80 7/80 7/80 7/80]'
%     OK to have zeros in input histogram, zeros are ignored:
%         >> p_approx = PMLdistributionApproximate([9 3 2 1 1 0 0 0 0 0 0]')
%         yields p_approx == [9/16 7/80 7/80 7/80 7/80 7/80]'
%     Assumed support set size:
%         >> p_approx = PMLdistributionApproximate([9 3 2 1 1]',5)
%         yields p_approx == [9/16 7/64 77/6/64 7/64 7/64]'

p_empirical = p_empirical(p_empirical > 0); % remove 0-counts

fing = int_hist(p_empirical(:));
fing_supp = find(fing); % sorted in ascending order
fing_values = fing(fing_supp);

B = length(fing_supp); % multibins, F_+
n = sum(fing_supp.*fing_values); % bin counts times fing_values, sum(bm) = sum(p)

% early detect case all symbols observed once?

% get counts and fing_values cumsum matrix
% n[i,j] = sum_{k = i to j} n[k]
counts = fing_supp(:).*fing_values(:);
counts_mat = cumsum_mat(counts);
mults_mat = cumsum_mat(fing_values);
prob_mat = counts_mat./(n.*mults_mat); % probability within each block

% C[i,j] = log reward assign multibins from i to j to same prob
log_reward_mat = gammaln(mults_mat+1) + counts_mat.*log(prob_mat);
log_reward_mat(isnan(log_reward_mat)) = 0.0; % handles 0 log 0 = 0 case, sets below diagonal entries to 0

% dynamic programming
% V[i] = best log reward for first i bins
V = zeros(B,1); % V[B+1] = 0, values matrix
backpointers = zeros(B,1)-1; % backpointers to get best log reward
for i = B:-1:1
    V(i) = log_reward_mat(i,B);
    backpointers(i) = B;
    for j = i:(B-1)
        proposal_reward = log_reward_mat(i,j) + V(j+1);
        if proposal_reward > V(i)
            V(i) = proposal_reward;
            backpointers(i) = j;
        end
    end
end

%% merge 0s with rest of symbols
V0_opt = -Inf; % optimal value
F0_opt = -1;
i0_opt = 0; % index of optimal fingerprint bin up to which to merge 0s
for i=1:B
    T = sum(fing_values(1:i));
    N = sum(fing_supp(1:i).*fing_values(1:i));
    if nargin == 2 % assumed alphabet size
        F0 = K - length(p_empirical);
    else
        % optimize alphabet size
        % estimate number of unseen symbols for uniform distribution, N samples, empirical support set size T
        F0 = ML_unseen_symbols_uniform(T, N);
    end
    % get value
    V0 = N*log(N/n);
    if i < B
        V0 = V0 + V(i+1); % V[B+1] = 0
    end
    if ~isinf(F0) % \lim_{F0 -> \infty} f(F0,T,n=T) = 0
        V0 = V0 + gammaln(F0+T+1) - gammaln(F0+1) - N*log(F0+T);
    end
    % check if best so far
    if V0 > V0_opt
        V0_opt = V0;
        F0_opt = F0;
        i0_opt = i;
    end
end
F0 = F0_opt;
V_approx = V0_opt;

%% get assignment
p_approx = p_empirical(:).*0;
if F0 > 0
    backpointers(1) = i0_opt; % merge 1s with same symbols as 0s, since guaranteed to merge 0s at least with 1s
end
i = 1;
while i <= B
    j = backpointers(i);
    if i == 1
        x1 = 1;
    else
        x1 = mults_mat(1,i-1)+1;
    end
    x2 = mults_mat(1,j);
    p_approx(x1:x2) = prob_mat(i,j);
    i = j + 1;
end
continuous_part = 0; % probability mass of continuous part
% add 0s in front
if F0 > 0
    x = mults_mat(1,i0_opt); % number of seen symbols in same level set as unseen symbols
    if ~isinf(F0)
        p_approx = [zeros(F0,1) ; p_approx];
        p_approx(1:(x+F0)) = counts_mat(1,i0_opt)/(n*(x+F0));
    else
        continuous_part = counts_mat(1,i0_opt)/n;
        p_approx(1:x) = [];
    end
end

%% get log probability of fingerprint
V_approx = V_approx + gammaln(n+1) - sum(fing_values(:).*gammaln(fing_supp(:)+1)) - sum(gammaln(fing_values+1));

%% sort PML approximation in descending order
p_approx = flipud(p_approx);

have_continuous_part = isinf(F0);

% % sort PML approximation in same order as input distribution
% [~, sortperm] = sort(p,'ascend');
% [~,sortinvperm] = sort(sortperm,'ascend');
% p_approx = p_approx(sortinvperm);
%
% p_approx = p_approx(:);

end

function C = cumsum_mat(x)
% computes cumultative sum matrix
% C[i,j] = sum(x(i:j))
%
% created by Dmitri S. Pavlichin on May 14, 2017
%
% Matlab version: R2015a

C = zeros(length(x),length(x));
for i = 1:length(x),
    C(i,i) = x(i);
    for j = i+1:length(x),
        C(i,j) = C(i,j-1) + x(j);
    end
end

end
