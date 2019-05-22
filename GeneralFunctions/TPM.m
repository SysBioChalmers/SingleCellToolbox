function outObj = TPM(obj)
% TPM
%   Converts an object to TPM/CPM
% Input:
%   obj             Could either be a matrix or an object with a .data
%                   member
% Usage: Two possibilities:
%   1. ds.data = TPM(ds.data)
%   2. ds = TPM(ds)
%
% Johan Gustafsson, 2019-05-21
%
    if isobject(obj)
        outObj = obj;
        currSum = sum(obj.data);
        outObj.data = outObj.data * 1000000 ./ currSum;
        outObj.data(isnan(outObj.data)) = 0;%problems with division by zero
    else
        outObj = obj;
        currSum = sum(obj);
        outObj = outObj * 1000000 ./ currSum;
        outObj(isnan(outObj)) = 0;
    end
end