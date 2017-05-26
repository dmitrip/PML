function H = entropyOfDistribution(vec,base)
% computes Shannon entropy of vector vec, default base 2

if nargin == 1,
    base = 2;
end

support = find(vec~=0);

vec = vec./sum(vec);

H = -sum(vec(support).*log(vec(support)))/log(base);