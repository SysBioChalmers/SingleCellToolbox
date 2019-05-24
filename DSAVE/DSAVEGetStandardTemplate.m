function ret = DSAVEGetStandardTemplate()
% DSAVEGetStandardTemplate
%   Gets the standard DSAVE template. The template is generated if needed,
%   but cached both in memory and on disc.
%
% Usage: templInfo = DSAVEGetStandardTemplate();
%
% Johan Gustafsson, 2019-05-21
%

DsHelper.init(); % to set current directory relative to source
persistent v;
if isempty(v) 
    disp('reading DSAVE standard template...');
    prevDir = DsHelper.setPathToSource();
    filename = '../../TempData/DSAVE_std_template.mat';
    if(~exist(filename,'file'))    
        disp('no .mat file found, regenerating template');
        ovm = DsHelper.scd_ovasc.cellSubset(DsHelper.scd_ovasc.cellType == Celltype.MacrophageOrMonocyte);
        bc2t = DsHelper.scd_bc2.cellSubset(DsHelper.scd_bc2.cellType == Celltype.TCellCD4Pos | DsHelper.scd_bc2.cellType == Celltype.TCellCD8Pos | DsHelper.scd_bc2.cellType == Celltype.TCellReg);
        bc2t_bc4tumor = bc2t.cellSubset(strcmp(bc2t.sampleIds, 'BC4_TUMOR'));
        b10000 = DsHelper.scd_pbmcb10000;
        [scd_GSE112845_pat_a,scd_GSE112845_pat_b,scd_GSE112845_cd8] = DsHelper.scd_GSE112845;

        datasets = {ovm,bc2t_bc4tumor, b10000, scd_GSE112845_cd8};
        v = DSAVEGenerateTemplateInfo(bc2t_bc4tumor, datasets, 2000, 750, 0.025, 0.025);
        save(filename, 'v');
    else
        a = load(filename);
        v = a.v;
    end
    DsHelper.restoreDir(prevDir);
end
ret = v;

end
