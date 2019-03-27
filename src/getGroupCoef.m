function coefMat = getGroupCoef(index, patchSize, groupMat)
patchWin = size(index);
coefMat = zeros(patchWin(1)+patchSize-1,patchWin(2)+patchSize-1);
mask = ones(patchSize);
gpM = zeros(patchWin(1)*patchWin(2),1);
for i = 1:size(groupMat,2)
    gpM(groupMat(:,i)) = gpM(groupMat(:,i))+1;
end
gpM = reshape(gpM,patchWin);
for i = 1:patchWin(1)
    for j = 1:patchWin(2)       
        coefMat(i:i+patchSize-1,j:j+patchSize-1) = coefMat(i:i+patchSize-1,j:j+patchSize-1)+mask*gpM(i,j);
    end
end

end