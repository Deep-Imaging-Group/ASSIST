function [option, params]=initialParams()
% The parameters of the geometry
option.angNum = 65;                 % The number of angles
option.angIncrement = 5.6;          % The increment of angle
option.imgSize = 256;               % The size of image, need to be square
option.Source2Center = 5;           % The distance from source to center
option.Source2Detector = 10;        % The distance from source to detector
option.detecNum = 512;              % The number of detector elements
option.detecLen = 4.5;              % The length of the detetcor
option.imgPysicSize = 2;            % The pysical size of the image
option.mode = 'fan';                % The projection mode

% The regular parameters
% The model is:
% min     lambda1 * |AX-Y|^2 
%       + lambda2 * |L|_* 
%       + lambda3 * |S|_TV
%       + rho * |R-L-S|^2
% Details can be seen in our paper
params.lambda1 = 1e8;
params.lambda2 = 10000;
params.lambda3 = 10000;
params.rho = 1e5;

% The parameters of group
params.patchSize = 6;
params.stride = 3;
params.groupNum = 12;
params.searchWin = 5;

% The number of iteration
params.iter = 300;