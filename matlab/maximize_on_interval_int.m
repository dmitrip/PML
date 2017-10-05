function [x, val] = maximize_on_interval_int(fun, x_min, x_max, tol)
% maximizes semiconcave function fun, integer arguments, on interval [x_min, x_max]
%
% created by Dmitri S. Pavlichin on May 14, 2017
%
% Matlab version: R2015a
%
% Args:
%     * fun - function to maximize
%     * x_min - interval min
%     * x_max - interval max
%
% Optional:
%     * tol (float) - tolerance optimize
%
% Returns:
%     * x (int) - maximizing value of x
%     * val - fun(x)

if nargin == 3
    tol = 1e-14;
end

n_iter = 0;
while floor(x_max) > floor(x_min)
    n_iter = n_iter + 1;
    x1 = ((x_max - x_min)/3)+x_min;
    x2 = 2*((x_max - x_min)/3)+x_min;
    V1 = fun(x1);
    V2 = fun(x2);
    if V1 - V2 <= -abs(tol)
        x_min = x1;
    else
        x_max = x2;
    end
    %disp([x_min x_max])
end
%disp(['n_iter = ' num2str(n_iter)])

x = x_min;
if fun(floor(x)) >= fun(ceil(x))
    x = floor(x_min);
else
    x = ceil(x_min);
end

val = fun(x);

