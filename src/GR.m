function groupPatch = GR(patchSet, groupMat)
N = size(patchSet);
[m, n] = size(groupMat);
groupPatch = zeros(N(1), m, N(3), n);

for i = 1:n    
    groupPatch(:,:,:,i) = patchSet(:,groupMat(:,i),:);    
end
end