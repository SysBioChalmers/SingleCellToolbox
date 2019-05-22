function s = ImportTabSepSamples(filename)
% ImportTabSepSamples
%   Reads a text file into a samples object
% Input:
%   filename  Filename
% Usage: s = ImportTabSepSamples(filename);
%
% Johan Gustafsson, 2019-05-20
%
    f = readtable(filename, 'ReadVariableNames',true, 'ReadRowNames', true, 'Delimiter', '\t');
    s = Samples;
    s.data = table2array(f);
    s.genes = f.Properties.RowNames;
    s.sampleIds = f.Properties.VariableNames;
end


