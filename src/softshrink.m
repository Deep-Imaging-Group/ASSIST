function out= softshrink(in, tau)
out = sign(in).*max(abs(in)-tau, 0);
end