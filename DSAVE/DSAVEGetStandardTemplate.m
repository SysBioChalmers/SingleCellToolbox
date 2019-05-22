function ret = DSAVEGetStandardTemplate()
% DSAVEGetStandardTemplate
%   Gets the standard DSAVE template. The template is generated if needed,
%   but cached both in memory and on disc.
%
% Usage: templInfo = DSAVEGetStandardTemplate();
%
% Johan Gustafsson, 2019-05-21
%

SCDep.init();
persistent v;
if isempty(v) 
    disp('reading DSAVE standard template...');
    prevDir = SCDep.setPathToSource();
    filename = '../../TempData/DSAVE_std_template.mat';
    if(~exist(filename,'file'))    
        disp('no .mat file found, regenerating template');
        ovm = SCDep.scd_ovasc.cellSubset(SCDep.scd_ovasc.cellType == Celltype.MacrophageOrMonocyte);
        bc2t = SCDep.scd_bc2.cellSubset(SCDep.scd_bc2.cellType == Celltype.TCellCD4Pos | SCDep.scd_bc2.cellType == Celltype.TCellCD8Pos | SCDep.scd_bc2.cellType == Celltype.TCellReg);
        bc2t_bc4tumor = bc2t.cellSubset(strcmp(bc2t.sampleIds, 'BC4_TUMOR'));
        b10000 = SCDep.scd_pbmcb10000;
        [scd_GSE112845_pat_a,scd_GSE112845_pat_b,scd_GSE112845_cd8] = SCDep.scd_GSE112845;

        datasets = {ovm,bc2t_bc4tumor, b10000, scd_GSE112845_cd8};
        v = DSAVEGenerateTemplateInfo(bc2t_bc4tumor, datasets, 2000, 750, 0.025, 0.025);
        save(filename, 'v');
    else
        a = load(filename);
        v = a.v;
    end
    SCDep.restoreDir(prevDir);
end
ret = v;

end
