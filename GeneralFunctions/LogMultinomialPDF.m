function ll = LogMultinomialPDF(obs, prob, specCalcLimit)
% LogMultinomialPDF
%   Calculates the probability density function at a certain value for a 
%   multinomial distribution. The value is a vector with one count per bin.
%   The idea with this function is that it is not easy to calculate the PDF of
%   a multinomial if n is large, since we cannot practically calculate n! of a
%   large number. The logarithm of the pdf is fine to calculate however, since 
%   log(n!) = sum(log(n) + log(n-1) + ... + log(1)). So we cannot use the
%   built-in function for the pdf and afterwards log it, but need to implement
%   the function ourselves instead.
% Input:
%   obs             Vertical vector of observed values
%   prob            Vertical vector of probabilities
%   specCalcLimit   Any n! where n is larger than this number will be
%                   calculated using a specialized method. Defaults to 150.
%                   This can be set lower from the test case, that is the
%                   only purpose of this parameter.
% Usage: ll = LogMultinomialPDF(obs, prob);
%
% Johan Gustafsson, 2019-05-20
%

if nargin < 3
    specCalcLimit = 150;%this can be specified to a low number for testing, that's the only purpose of not hardcoding this value
end
n = sum(obs,1);

%formula:
%x = obs;
%p = prob;
% PDF = n!/(x1! ... xn!) * p1^x1 * .... * pn^xn
% =>
% log(PDF) = log(n!) - log(x1!) ... - log(xn!) + log(p1^x1) + .... + log(pn*xn)


%generate list of n, n-1, ..1, cannot calculate n! for large n otherwise
facs = n:-1:1;

s1 = sum(log(facs),2);
s2 = log(factorial(obs));%don't sum yet
s3 = nansum(log(prob).*obs);

% Some of the s2s may have become inf due to that matlab cannot calculate
% n! for large n (> 170, we used 150 to stay on the safe side). Replace those values by values calculated by the
% n method
sel = obs > specCalcLimit; 
ind = (1:size(obs,1)).';
selInd = ind(sel);
for i = 1:size(selInd)
    facs = obs(selInd(i)):-1:1;
    s2(selInd(i)) = sum(log(facs),2);
end

s2 = -sum(s2,1);

ll = s1 + s2 + s3;

end


