function h = int_hist(x)
% int_hist(x) is a histogram of all integer values 1:max(max(x))
% this file created by Dmitri S. Pavlichin on May 23rd, 2017

wid_x = size(x,2);
large = max(max(x));
h = zeros(large,wid_x);
for iter = 1:wid_x
    for j = 1:size(x,1)
        h(x(j,iter),iter) = h(x(j,iter),iter) + 1;
    end
end
