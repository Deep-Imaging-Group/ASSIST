function out = TensorNuclear(in, tau)
N = size(in);
[in1, in2, in3] = unfold(in);
out = (fold1(nuclearOptimal(in1,tau),N) + fold2(nuclearOptimal(in2,tau),N)+fold3(nuclearOptimal(in3,tau),N))/3;
end

function y = fold1(X, N)
y = zeros(N);
for i = 1:N(2)
    for j = 1:N(3)
        y(:, i, j) = X(:, ((j-1)*N(2) + i));
    end
end
end

function y = fold2(X,N)
    y = zeros(N);
    for i = 1:N(1)
        for j = 1:N(3)
            y(i, :, j) = X(:, ((j-1)*N(1) + i)); 
        end
    end
end

function y = fold3(X,N)
    y = zeros(N);
    for i = 1:N(1)
        for j = 1:N(2)
            y(i, j, :) = X(:, ((j-1)*N(1) + i));
        end
    end
end



function [ Y1, Y2, Y3 ] = unfold( Y_t )
%UNFOLD Summary of this function goes here
%   Turn tensor Y to matrix Y1, Y2, Y3 as the unfold rules.

Ysiz = size(Y_t);
Y1 = zeros(Ysiz(1), Ysiz(2) *Ysiz(3));
Y2 = zeros(Ysiz(2), Ysiz(1) *Ysiz(3));
Y3 = zeros(Ysiz(3), Ysiz(2) *Ysiz(1));

for i = 1:Ysiz(2)
    for j = 1:Ysiz(3)
        Y1(:, ((j-1)*Ysiz(2) + i)) = Y_t(:, i, j);
    end
end

for i = 1:Ysiz(1)
    for j = 1:Ysiz(3)
        Y2(:, ((j-1)*Ysiz(1) + i)) = Y_t(i, :, j);        
    end
end

for i = 1:Ysiz(1)
    for j = 1:Ysiz(2)
        Y3(:, ((j-1)*Ysiz(1) + i)) = Y_t(i, j, :); 
    end
end

end