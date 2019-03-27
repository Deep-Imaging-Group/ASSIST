clear all;
clc;

% The size of data need to be [Nrow, Ncol, Nchannel]
load ImgData;
oringImg = X;
[option,params]=initialParams();
% Set parallel pools, will greatly increaed the speed,
% Recommended setting is the number of channels.
parpoolnum = 2;

Xsz = size(oringImg);

% Generate the system matrix and projection data
sysmat = GetSysMat(option);
params.sysmat = sysmat;
projData = zeros(option.angNum, option.detecNum, Xsz(3));
for i = 1:Xsz(3)
    projData(:,:,i) = projGen(oringImg(:,:,i),sysmat);
end
delete(gcp('nocreate'));
reconImg = reconstruct(oringImg, projData, params,parpoolnum);

%The display window
showWin = [(-160+500+1024)/3000,(240+500+1024)/3000];
for i = 1:Xsz(3)
    figure;imshow(reconImg(:,:,i),showWin);
end
delete(gcp('nocreate'));
