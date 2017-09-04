function K_est = estimate_support_from_histogram_PML_approximate(hist_vec, min_value)
% estimates support set size of distribution given sample samp
%
% approximates pattern maximum likelihood (PML) distribution p
% where p \approx \arg \max_p' P(unordered histogram)
% Approximation sums over all permutations that mix within blocks of
% constant p
%
% Uses approximate PML distribution as plug-in for support functional
%
% created by Dmitri S. Pavlichin on June 8, 2017
%
% Matlab version: R2015a
%
% Args:
%     * hist_vec - (vector) histogram, vector of nonnegative integers
%     * min_value - (float) assumed minimum value of distribution on support
%
% Returns:
%     * K_est - (integer) estimate for support set size
%
% Example:

% estimate support set size
[p_approx,~,~,not_have_valid_K_est] = PMLdistributionApproximate(hist_vec);
% K_est = sum(hist_vec > 0) + F0_est;
K_est = length(p_approx);

% check for min element, adjust F0_est until min element matches
if not_have_valid_K_est
    %     warning('do not have valid estimate for support set size')
    % add 0 counts until value assigned to them is min_value
    cont_mass = 1 - sum(p_approx); % mass in continuous part of distribution
    K_est = length(p_approx) + cont_mass/min_value;
    %     if cont_mass < 1
    % %         keyboard
    %     end
    %     if min(p_approx) < min_value
    %         disp('there')
    %         keyboard
    %     end
elseif min(p_approx) < min_value
    mu = min(p_approx);
    r_mu = find(p_approx == mu);
    min_mass = sum(p_approx(r_mu)); % mass assigned to min
%     keyboard
    K_est = K_est - length(r_mu) + min_mass / min_value;
    
    %     % check if second-loweest value lower than min_value
    %     disp('here')
    %     mu = min(p_approx);
    %     if min(p_approx(p_approx > mu)) < min_value
    %         keyboard
    %     end
end
