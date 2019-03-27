function out = nuclearOptimal(in, tau)

[U, S, V] = svd(in, 'econ');
index = find(S>0);
S(index) = softshrink(S(index), tau);
out = U*S*V';

end