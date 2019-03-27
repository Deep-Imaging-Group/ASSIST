function res = chambolleAlgorithm(in, alpha, step, iter)
N = size(in);
px = zeros(N);
py = zeros(N);
for i = 1:iter
    temp = div(px, py)-in/alpha;
    dx = gradx(temp);
    dy = grady(temp);
    tv = sqrt(dx.^2+dy.^2);
    px = (px+step*dx)./(1+step*tv);
    py = (py+step*dy)./(1+step*tv);
end
res = in-alpha*div(px, py);
end

function res = gradx(in)
[m,~] = size(in);
res = in([2:m m],:)-in;
end

function res = grady(in)
[~,n] = size(in);
res = in(:,[2:n n])-in;
end

function res = div(inx, iny)
[m,n] = size(inx);
dx = inx - inx([1 1:m-1],:);
dx(1,:) = inx(1,:);
dx(m,:) = -inx(m-1,:);
dy = iny - iny(:,[1 1:n-1]);
dy(:,1) = iny(:,1);
dy(:,n) = -iny(:,n-1);
res = dx+dy;
end