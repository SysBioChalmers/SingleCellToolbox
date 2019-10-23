%this does not work for the format where R does not place an empty string
%followed by a tab above the row names, since the column names row will be
%one element shorter than the rest of the rows. I currently don't know how
%to solve this in a good way.
function s = ImportSamplesFromRText(filename)
    
    f = readtable(filename, 'ReadVariableNames',true, 'ReadRowNames', true, 'Delimiter', '\t');
    s = Samples;
    s.data = table2array(f).';
    s.sampleIds = f.Properties.VariableNames.';
    s.genes = f.Properties.RowNames.';
end


