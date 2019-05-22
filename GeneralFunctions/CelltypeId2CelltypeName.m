%

function res = CelltypeId2CelltypeName(x)
% CelltypeId2CelltypeName
%   Converts an numeric value of cell type into a readable text string
% Input:
%   x               The input constant from the class Celltype. Can also be
%                   a scalar.
% Output:
%   res             A single char array or a cell array of such, depending on input.
% Usage: res = CelltypeId2CelltypeName(ds.cellType);
%
% Johan Gustafsson, 2019-05-20
%
[r,c] = size(x);
if (r > 1) || (c > 1)
    xscalar = 0;
else
    xscalar = 1;
end

if ~xscalar
    %[r,c] = size(x);
    %res = cell(r,c);
    %this can probably be done smarter, but it works
    res = num2cell(x);
    res = cellfun(@(b) ConvCT(b), res, 'UniformOutput', false);
else
   res = ConvCT(x); 
end
    
    function str = ConvCT(a)
        switch a
            case Celltype.Invalid
                str = 'Invalid';
            case Celltype.Unknown
                str = 'Unknown';
            case Celltype.Malignant
                str = 'Malignant';
            case Celltype.Stromal
                str = 'Stromal';
            case Celltype.Myeloid
                str = 'Myeloid';
            case Celltype.TCell
                str = 'T cell';
            case Celltype.BCell
                str = 'B cell';
            case Celltype.NKCell
                str = 'NK cell';
            case Celltype.Fibroblast
                str = 'Fibroblast';
            case Celltype.Macrophage
                str = 'Macrophage';
            case Celltype.Endothelial
                str = 'Endothelial';
            case Celltype.Myocyte
                str = 'Myocyte';
            case Celltype.Dendritic
                str = 'Dendritic';
            case Celltype.Mast
                str = 'Mast';
            case Celltype.Epithelial
                str = 'Epithelial';
            case Celltype.MacrophageOrMonocyte
                str = 'Macrophage/Monocyte';
            case Celltype.Monocyte
                str = 'Monocyte';
            case Celltype.TCellCD8Pos
                str = 'T cell CD8+';
            case Celltype.TCellCD4Pos
                str = 'T cell CD4+';
            case Celltype.TCellReg
                str = 'Reg T cell';
            case Celltype.OvarianCarcinoma
                str = 'Ovarian carcinoma';
            case Celltype.Melanoma
                str = 'Melanoma';
            case Celltype.Neutrophil
                str = 'Neutrophil';
            case Celltype.Alveolar
                str = 'Alveolar';
            case Celltype.Erythroblast
                str = 'Erythroblast';
            case Celltype.Langerhans
                str = 'Langerhans';
            case Celltype.Granulocyte
                str = 'Granulocyte';
            case Celltype.Megakaryocyte
                str = 'Megakaryocyte';
            case Celltype.HematopeticStemOrProgenitor
                str = 'Hematopetic Stem Or Progenitor';
            case Celltype.TCellCD4Memory
                str = 'T cell CD4+ Memory';
            case Celltype.Erythrocyte
                str = 'Erythrocyte';
            otherwise
                str = strcat('non-supported: ', num2str(a));
        end
    end
end
