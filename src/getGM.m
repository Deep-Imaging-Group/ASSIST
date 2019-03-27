function groupMat = getGM(patchSet, index, curMat, searchWin, groupNum, comBin)
N = size(patchSet);
M = size(index);
[m, n] = size(curMat);
WinSize = size(index);
groupMat = zeros(groupNum, m*n);

for i = 1:(m*n)
    cur = curMat(i); 
    curCol = ceil(cur/M(1));
    curRow = cur - M(1)*(curCol-1);
    rmin = max(1, curRow-searchWin);
    rmax = min(WinSize(1), curRow+searchWin);
    cmin = max(1, curCol-searchWin);
    cmax = min(WinSize(2), curCol+searchWin);
    winIndex = index(rmin:rmax,cmin:cmax);
    winIndex = winIndex(:);
    t = randperm(length(winIndex));
    winIndex = winIndex(t);
    error = bsxfun(@minus, patchSet(:,winIndex,comBin),patchSet(:,cur,comBin));
    error = sum(error.^2);
    [~, picIndex] = sort(error);
    groupMat(:,i) = winIndex(picIndex(1:groupNum)); 
    groupMat(1,i) = cur;
end
end