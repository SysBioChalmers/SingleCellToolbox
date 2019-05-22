
classdef Celltype
% Celltype
%   Enumeration for human cell types. We avoided using an enum due to some 
%   irritating properties, so now we just use constants in a class.
%   Also includes a cell type hierarchy tree.
%
% Johan Gustafsson, 2019-05-20
%

    properties( Constant = true )
 		Invalid = 0
		Unknown = 1
		Malignant = 2
		Stromal = 3
		Myeloid = 4
		TCell = 5
		BCell = 6
		NKCell = 7
		Fibroblast = 8
		Macrophage = 9
		Endothelial = 10
		Myocyte = 11
		Dendritic = 12
		Mast = 13
        Epithelial = 14
        MacrophageOrMonocyte = 15
        Monocyte = 16
        TCellCD8Pos = 17
        TCellCD4Pos = 18
        TCellReg = 19
        OvarianCarcinoma = 20
        Melanoma = 21
        Neutrophil = 22
        Alveolar = 23
        Erythroblast = 24
        Langerhans = 25
        Granulocyte = 26
        Megakaryocyte = 27
        HematopeticStemOrProgenitor = 28
        TCellCD4Memory = 29 % CD4+/CD45RO+
        Erythrocyte = 30
        %TCellCD4Naive = 23
        %TCellCD8Naive = 24
        %NKTCell = 25
        
        
        %Cell hierarchy. Celltype, parent. No entries for top level types.
        cellHierarchy = [Celltype.TCellCD4Pos Celltype.TCell; ...
                         Celltype.TCellCD8Pos Celltype.TCell; ...
                         Celltype.TCellReg Celltype.TCellCD4Pos; ...
                         Celltype.TCellReg Celltype.TCellCD4Pos; ...
                         Celltype.Macrophage Celltype.MacrophageOrMonocyte; ...
                         Celltype.Monocyte Celltype.MacrophageOrMonocyte; ...
                         Celltype.OvarianCarcinoma Celltype.Malignant; ...
                         Celltype.Melanoma Celltype.Malignant
                         Celltype.TCellCD4Memory Celltype.TCellCD4Pos]
    end
    
end