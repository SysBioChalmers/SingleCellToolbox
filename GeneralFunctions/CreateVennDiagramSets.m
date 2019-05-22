function [s1,s2,s3,s1s2,s1s3,s2s3,s1s2s3] = CreateVennDiagramSets(s1in,s2in,s3in)
% CreateVennDiagramSets
%   This function takes 3 string sets and returns the intersection sets for
%   a Venn diagram. Make sure there are no duplicates in the sets that are
%   expected to be cell arrays of char arrays
% Input:
%   s1in,s2in,s3in  The three string sets.
% Output:
%   s1,s2,s3,s1s2,s1s3,s2s3,s1s2s3 String sets where for example s1s2
%   contains all string present in s1in and s2in, but not in s3in.
% Usage: [s1,s2,s3,s1s2,s1s3,s2s3,s1s2s3] = CreateVennDiagramSets(s1in,s2in,s3in);
%
% Johan Gustafsson, 2019-05-20
%

s1s2 = intersect(s1in,s2in);
s1s3 = intersect(s1in,s3in);
s2s3 = intersect(s2in,s3in);

s1 = setdiff(s1in,s1s2);
s1 = setdiff(s1,s1s3);
s2 = setdiff(s2in,s1s2);
s2 = setdiff(s2,s2s3);
s3 = setdiff(s3in,s1s3);
s3 = setdiff(s3,s2s3);

s1s2s3 = intersect(s1s2,s1s3);

s1s2 = setdiff(s1s2,s1s2s3);
s1s3 = setdiff(s1s3,s1s2s3);
s2s3 = setdiff(s2s3,s1s2s3);

end