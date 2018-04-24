function q = DrawFromMultinomial(p,n)
% draws type from pmf p, n samples
%
% created by Dmitri S. Pavlichin on May 14, 2017
%
% Matlab version: R2015a

p = p./sum(p);
q = p.*0; % current type
z = sum(q); % sum(q), overshoot, undershoot
num_trials = 0;
alpha = 0.999;
while z < n,
    q_ = random('poisson',alpha.*(p.*(n-z)));
    z_ = sum(q_);
    if z + z_ <= n, % accept
        q = q + q_;
        z = z + z_;
    end
    num_trials = num_trials + 1;
end