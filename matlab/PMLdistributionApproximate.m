function [p_approx,F0_est,V_approx,reached_max_F0] = PMLdistributionApproximate(p_empirical,K)
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
% %     * K_upper_bound - (integer) assumed maximum value for support set
%         size, use when estimating support set size.  If K_upper_bound is
%         not provided, then we set K_upper_bound = sum(p == 1)^2
%
% Returns:
%     * p_approx - (column vector) approximate PML distribution, sorted in
%         descending order
%     * F0_est - (integer) estimate for number of symbols seen 0 times in
%         p_empirical.  If K is provided as an argument, then F0_est = max(K - sum(p > 0), 0)
%     * V_approx - (double) approximate value of log(P_p(unordered p_empirical))
%         achieved by p = p_approx
%     * reached_max_F0 - (boolean) true iff we did not find an upper bound
%        when estimating support set size, hit upper bound sum(p > 0) + sum(p == 1)^2
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

fing = int_hist(p_empirical(:)+1);
fing = fing(2:end);
bin_values = flipud(find(fing)); % sort in descending order
multiplicities = fing(bin_values);
num_bins = sum(multiplicities);

B = length(bin_values); % multibins
n = sum(bin_values.*multiplicities); % bin counts times multiplicities, sum(bm) = sum(p)
n_1 = sum(p_empirical == 1); % number of symbols observed once

% early detect case n_1 == n?

% get counts and multiplicities cumsum matrix
% n[i,j] = sum_{k = i to j} n[k]
counts = bin_values(:).*multiplicities(:);
counts_mat = cumsum_mat(counts);
mults_mat = cumsum_mat(multiplicities);
prob_mat = counts_mat./(n.*mults_mat); % probability within each block

% C[i,j] = log reward assign multibins from i to j to same prob
log_reward_mat = gammaln(mults_mat+1) + counts_mat.*log(prob_mat);
log_reward_mat(isnan(log_reward_mat)) = 0.0; % handles 0 log 0 = 0 case, sets below diagonal entries to 0

% V[i] = best log reward for first i bins
V = zeros(B,1); % V[0] = 0, values matrix
backpointers = zeros(B,1); % backpointers to get best log reward
for i = 1:B
    V(i) = log_reward_mat(1,i);
    for j = 1:i-1
        proposal_reward = V(j) + log_reward_mat(j+1,i);
        if proposal_reward > V(i)
            V(i) = proposal_reward;
            backpointers(i) = j;
        end
    end
end

%% optimize alphabet size
% get upper bound on F0
if nargin == 1
    F0 = 0; % extra 0 entries
    step_size = 10;
    F0_max = num_bins^2;
    V0 = f0(V, counts, multiplicities, 0, n, num_bins, B) - gammaln(0+1);
    
    reached_max_F0 = false;
    done = false;
    while ~done,
        F0 = F0 + step_size;
        
        V0_prev = V0;
        V0 = f0(V, counts, multiplicities, F0, n, num_bins, B) - gammaln(F0+1);
        
        if V0_prev >= V0,
            done = true;
        else
            step_size = step_size*2;
            if F0 > F0_max,
                %             warning('inferred support size exceeds max of length(p)^2');
%                 warning('inferred support size exceeds max of 100*length(p)');
                reached_max_F0 = true;
                done = true;
            end
        end
    end
    F0_max = F0;
    
    % get optimum value of F0
    if ~reached_max_F0
        F0_lower = 0;
        F0_upper = F0_max;
        while F0_lower ~= F0_upper
            F1 = floor((F0_upper - F0_lower)/3)+F0_lower;
            F2 = 2*floor((F0_upper - F0_lower)/3)+F0_lower+1;
            V1 = f0(V, counts, multiplicities, F1, n, num_bins, B) - gammaln(F1+1);
            V2 = f0(V, counts, multiplicities, F2, n, num_bins, B) - gammaln(F2+1);
            if V1 - V2 <= -1e-10
                F0_lower = F1+1;
            else
                F0_upper = F2-1;
            end
        end
        F0 = F0_lower;
    else
        F0 = F0_max;
    end
    % test to see if have continuous part
    if bin_values(end) == 1
        V_non_cont = f0(V, counts, multiplicities, F0, n, num_bins, B) - gammaln(F0+1);
        if length(bin_values) > 1
            V_cont = V(end-1) + n_1*log(n_1/n);
        else
            V_cont = 0.0;
        end
        if V_cont > V_non_cont
            warning('optimal distribution has continuous part, F0 = Inf');
            F0 = Inf;
        end
    end
else
    F0 = K - num_bins;
    if F0 < 0
        F0 = 0;
    end
end
F0_est = F0;

% keyboard

if ~isinf(F0)
    [Vmax,imax] = f0(V, counts, multiplicities, F0, n, num_bins, B);
    backpointers(end+1) = imax;
    bin_values(end+1) = 0;
    multiplicities(end+1) = F0;
    B = B + 1;
    num_bins = num_bins + F0;
    counts = bin_values(:).*multiplicities(:);
    counts_mat = cumsum_mat(counts);
    mults_mat = cumsum_mat(multiplicities);
    prob_mat = counts_mat./(n.*mults_mat); % probability within each block
else
    Vmax = V_cont;
    backpointers(end) = [];
    bin_values(end) = [];
    multiplicities(end) = [];
    B = B - 1;
    num_bins = num_bins - n_1;
    counts = bin_values(:).*multiplicities(:);
    counts_mat = cumsum_mat(counts);
    mults_mat = cumsum_mat(multiplicities);
    prob_mat = counts_mat./(n.*mults_mat); % probability within each block
end

%% get assignment
p_approx = zeros(num_bins,1);

current_multibin = B;
ix = num_bins;

while ix > 0
    t_prev = current_multibin;
    t = backpointers(current_multibin);
    while current_multibin > t
        for m = 1:multiplicities(current_multibin)
            p_approx(ix) = prob_mat(t+1,t_prev);
            ix = ix - 1;
        end
        current_multibin = current_multibin - 1;
    end
end

%%

V_approx = Vmax + gammaln(n+1) - sum(multiplicities(:).*gammaln(bin_values(:)+1)) - sum(gammaln(multiplicities+1));

% % sort PML approximation in same order as input distribution
% [~, sortperm] = sort(p,'ascend');
% [~,sortinvperm] = sort(sortperm,'ascend');
% p_approx = p_approx(sortinvperm);
%
% p_approx = p_approx(:);

end

function [Vmax,imax] = f0(V, counts, multiplicities, F0, n, num_bins, B)
Vmax = gammaln(num_bins+F0+1) + n*log(1/(num_bins+F0));
imax = 0;
for i = 1:B
    c = sum(counts(i+1:end));
    m = sum(multiplicities(i+1:end));
    
    if c > 0
        proposal_reward = V(i) + gammaln(m+F0+1) + c*log((c/n)/(m+F0));
    else
        proposal_reward = V(i) + gammaln(m+F0+1);
    end
    
    if proposal_reward > Vmax
        Vmax = proposal_reward;
        imax = i;
    end
end

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
