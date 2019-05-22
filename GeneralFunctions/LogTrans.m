function result = LogTrans(M, forward, numToAdd)
% LogTrans
%   Transforms data to log scale
% Input:
%   M               A matrix to transform
%   forward         True if transforming to log scale, false if from log
%                   scale.
%   numToAdd        (optional) The number to add to all data before log.
%                   Defaults to 1.
% Usage: result = LogTrans(M, true);
%
% Johan Gustafsson, 2019-05-20
%
if nargin < 3
    numToAdd = 1;
end
if (forward)
    result = log(M+numToAdd);
else
    result = exp(1).^M - numToAdd;
end