function img = patchReturn(patch, N, groupMat)
M = size(patch);
patchSize = sqrt(M(1));
img = zeros(N);
patchWin = [N(1)-patchSize+1, N(2)-patchSize+1, N(3)];
temp = zeros(M(1), patchWin(1)*patchWin(2), N(3));
for i = 1:M(4)
    temp(:,groupMat(:,i),:) = temp(:,groupMat(:,i),:) + patch(:,:,:,i);
end
% temp(:,curMat,:) = permute(reshape(sum(patch,2),[M(1),M(3),M(4)]), [1,3,2]);

for i = 1:patchSize
    for j = 1:patchSize
        img(j:j+patchWin(1)-1, i:i+patchWin(2)-1,:) = img(j:j+patchWin(1)-1, i:i+patchWin(2)-1,:)+ reshape(temp((i-1)*patchSize+j,:,:), patchWin);
    end
end
% img = img./M(2);
% img = bsxfun(@rdivide, img, coef);

end