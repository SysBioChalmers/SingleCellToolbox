% DSAVE Test cases
% Johan Gustafsson, 2019-05-21
% See VerificationMatrix.csv for how test cases map to code

%% T0001: SCDataset: Test geneSubset and cellSubset
ds1 = SCDataset;
ds1.genes = {'A';'B';'C';'D'};
ds1.data = [1 2 3 4; 5 6 7 8; 9 10 11 12; 13 14 15 16];
ds1.cellIds = {'1_1','1_2','1_3','1_4'};
ds1.sampleIds = {'1_11','1_21','1_31','1_41'};
ds1.name = 'ds1';
ds1.cellType = [1 2 3 4];
ds1.subCellType = [11 12 13 14];
ds1.extraCellInfo = {'A','C','F','H'};

gss1 = ds1.geneSubset(strcmp(ds1.genes,'A') | strcmp(ds1.genes,'C'));
gss2 = ds1.geneSubset({'A';'C'});
gssExpData = [1 2 3 4; 9 10 11 12];
gssExpGenes = {'A';'C'};

css = ds1.cellSubset(~strcmp(ds1.cellIds,'1_3'));
cssExpData = [1 2 4; 5 6 8; 9 10 12; 13 14 16];
cssExpCellIds = {'1_1','1_2','1_4'};
cssExpSampleIds = {'1_11','1_21','1_41'};
cssExpCellType = [1 2 4];
cssExpSubCellType = [11 12 14];
cssExpExtraCellInfo = {'A','C','H'};

%verify
if ~all(strcmp(gss1.genes, gssExpGenes))
    error('T0001: SCDataset.geneSubset: genes incorrect');
elseif ~all(strcmp(gss2.genes, gssExpGenes))
    error('T0001: SCDataset.geneSubset: genes incorrect');
elseif ~isequal(gss1.data, gssExpData)
    error('T0001: SCDataset.geneSubset: data incorrect');
elseif ~isequal(gss2.data, gssExpData)
    error('T0001: SCDataset.geneSubset: data incorrect');
elseif ~all(strcmp(gss1.cellIds, ds1.cellIds))
    error('T0001: SCDataset.geneSubset: cellIds incorrect');
elseif ~all(strcmp(gss2.cellIds, ds1.cellIds))
    error('T0001: SCDataset.geneSubset: cellIds incorrect');
elseif ~all(strcmp(gss1.sampleIds, ds1.sampleIds))
    error('T0001: SCDataset.geneSubset: sampleIds incorrect');
elseif ~all(strcmp(gss2.sampleIds, ds1.sampleIds))
    error('T0001: SCDataset.geneSubset: sampleIds incorrect');
elseif ~isequal(gss1.cellType, ds1.cellType)
    error('T0001: SCDataset.geneSubset: cellType incorrect');
elseif ~isequal(gss2.cellType, ds1.cellType)
    error('T0001: SCDataset.geneSubset: cellType incorrect');
elseif ~isequal(gss1.subCellType, ds1.subCellType)
    error('T0001: SCDataset.geneSubset: subCellType incorrect');
elseif ~isequal(gss2.subCellType, ds1.subCellType)
    error('T0001: SCDataset.geneSubset: subCellType incorrect');
elseif ~all(strcmp(gss1.extraCellInfo, ds1.extraCellInfo))
    error('T0001: SCDataset.geneSubset: extraCellInfo incorrect');
elseif ~all(strcmp(gss2.extraCellInfo, ds1.extraCellInfo))
    error('T0001: SCDataset.geneSubset: extraCellInfo incorrect');
    
elseif ~all(strcmp(css.genes, ds1.genes))
    error('T0001: SCDataset.cellSubset: genes incorrect');
elseif ~isequal(css.data, cssExpData)
    error('T0001: SCDataset.cellSubset: data incorrect');
elseif ~all(strcmp(css.cellIds, cssExpCellIds))
    error('T0001: SCDataset.cellSubset: cellIds incorrect');
elseif ~all(strcmp(css.sampleIds, cssExpSampleIds))
    error('T0001: SCDataset.cellSubset: sampleIds incorrect');
elseif ~isequal(css.cellType, cssExpCellType)
    error('T0001: SCDataset.cellSubset: cellType incorrect');
elseif ~isequal(css.subCellType, cssExpSubCellType)
    error('T0001: SCDataset.cellSubset: subCellType incorrect');
elseif ~all(strcmp(css.extraCellInfo, cssExpExtraCellInfo))
    error('T0001: SCDataset.cellSubset: extraCellInfo incorrect');
else
    disp('T0001: SCDataset.geneSubset and SCDataset.cellSubset: ok');
end
   

%% T0014: Samples: Test geneSubset and sampleSubset
s1 = Samples;
s1.genes = {'A';'B';'C';'D'};
s1.data = [1 2 3 4; 5 6 7 8; 9 10 11 12; 13 14 15 16];
s1.sampleIds = {'1_11','1_21','1_31','1_41'};

gss1 = s1.geneSubset(strcmp(s1.genes,'A') | strcmp(s1.genes,'C'));
gss2 = s1.geneSubset({'A';'C'});
gssExpData = [1 2 3 4; 9 10 11 12];
gssExpGenes = {'A';'C'};

sss = s1.sampleSubset(~strcmp(s1.sampleIds,'1_31'));
sssExpData = [1 2 4; 5 6 8; 9 10 12; 13 14 16];
sssExpSampleIds = {'1_11','1_21','1_41'};

%verify
if ~all(strcmp(gss1.genes, gssExpGenes))
    error('T-0014: Samples.geneSubset: genes incorrect');
elseif ~all(strcmp(gss2.genes, gssExpGenes))
    error('T-0014: Samples.geneSubset: genes incorrect');
elseif ~isequal(gss1.data, gssExpData)
    error('T-0014: Samples.geneSubset: data incorrect');
elseif ~isequal(gss2.data, gssExpData)
    error('T-0014: Samples.geneSubset: data incorrect');
elseif ~all(strcmp(gss1.sampleIds, s1.sampleIds))
    error('T-0014: Samples.geneSubset: sampleIds incorrect');
elseif ~all(strcmp(gss2.sampleIds, s1.sampleIds))
    error('T-0014: Samples.geneSubset: sampleIds incorrect');
    
elseif ~all(strcmp(sss.genes, s1.genes))
    error('T-0014: Samples.sampleSubset: genes incorrect');
elseif ~isequal(sss.data, sssExpData)
    error('T-0014: Samples.sampleSubset: data incorrect');
elseif ~all(strcmp(sss.sampleIds, sssExpSampleIds))
    error('T-0014: Samples.sampleSubset: sampleIds incorrect');
else
    disp('T-0014: Samples.geneSubset and SCDataset.cellSubset: ok');
end
   


%% T0002: Test dataset joins
joinds1 = SCDataset;
joinds2 = SCDataset;

joinds1.genes = {'A';'B';'C';'D'};
joinds2.genes = {'C';'A';'E'};
joinds1.cellIds = {'1_1','1_2','1_3','1_4'};
joinds2.cellIds = {'2_1','2_2'};
joinds1.data = [1 2 3 4; 5 6 7 8; 9 10 11 12; 13 14 15 16];
joinds2.data = [20 21; 22 23; 24 25];
joinds1.name = 'ds1';
joinds2.name = 'ds2';
joinds1.cellType = [1 2 3 4];
joinds2.cellType = [5 6];
joinds1.subCellType = [11 12 13 14];
joinds2.subCellType = [15 16];
joinds1.extraCellInfo = {'A','C','F','H'};
joinds2.extraCellInfo = {'J','L'};
expInnerJoinGenes = {'A';'C'};
expOuterJoinGenes = {'A';'B';'C';'D';'E'};
expJoinCellIds = {'1_1','1_2','1_3','1_4','2_1','2_2'};
expInnerJoinData = [1 2 3 4 22 23; 9 10 11 12 20 21];
expOuterJoinData = [1 2 3 4 22 23; 5 6 7 8 0 0; 9 10 11 12 20 21; 13 14 15 16 0 0; 0 0 0 0 24 25];
expInnerJoinName = 'inner join (ds1,ds2)';
expOuterJoinName = 'full outer join (ds1,ds2)';
expJoinCellType = [1 2 3 4 5 6];
expJoinSubCellType = [11 12 13 14 15 16];
expJoinExtraCellInfo = {'A','C','F','H','J','L'};

innerds = joinds1.innerJoin(joinds2);
outerds = joinds1.fullOuterJoin(joinds2);

%verify
if ~all(strcmp(innerds.genes, expInnerJoinGenes))
    error('T0002:SCDataset.innerJoin: genes incorrect');
elseif ~all(strcmp(outerds.genes, expOuterJoinGenes))
    error('T0002: SCDataset.fullOuterJoin: genes incorrect');
elseif ~all(strcmp(innerds.cellIds, expJoinCellIds))
    error('T0002: SCDataset.innerJoin: cell ids incorrect');
elseif ~all(strcmp(outerds.cellIds, expJoinCellIds))
    error('T0002: SCDataset.fullOuterJoin: cell ids incorrect');
elseif ~isequal(innerds.data, expInnerJoinData)
    error('T0002: SCDataset.innerJoin: data incorrect');
elseif ~isequal(outerds.data, expOuterJoinData)
    error('T0002: SCDataset.fullOuterJoin: data incorrect');
elseif ~isequal(innerds.cellType, expJoinCellType)
    error('T0002: SCDataset.innerJoin: cellType incorrect');
elseif ~isequal(outerds.cellType, expJoinCellType)
    error('T0002: SCDataset.fullOuterJoin: cellType incorrect');
elseif ~isequal(innerds.subCellType, expJoinSubCellType)
    error('T0002: SCDataset.innerJoin: custClust incorrect');
elseif ~isequal(outerds.subCellType, expJoinSubCellType)
    error('T0002: SCDataset.fullOuterJoin: custClust incorrect');
elseif ~all(strcmp(innerds.extraCellInfo, expJoinExtraCellInfo))
    error('T0002: SCDataset.innerJoin: extra cell info incorrect');
elseif ~all(strcmp(outerds.extraCellInfo, expJoinExtraCellInfo))
    error('T0002: SCDataset.fullOuterJoin: extra cell info incorrect');
elseif ~strcmp(innerds.name, expInnerJoinName)
    error('T0002: SCDataset.innerJoin: name incorrect');
elseif ~strcmp(outerds.name, expOuterJoinName)
    error('T0002: SCDataset.fullOuterJoin: name incorrect');
else
    disp('T0002: SCDataset.innerJoin and SCDataset.fullOuterJoin: ok');
end

%% T0003: SCDataset.FillEmpties
ds1 = SCDataset;
ds1.genes = {'A';'B';'C';'D'};
ds1.data = [1 2 3 4; 5 6 7 8; 9 10 11 12; 13 14 15 16];
ds2 = ds1.fillEmpties();
expCellIds = {'_1','_2','_3','_4'};
expSampleIds = {'Unknown','Unknown','Unknown','Unknown'};
expCellType = [Celltype.Unknown Celltype.Unknown Celltype.Unknown Celltype.Unknown];
expSubCellType = [0 0 0 0];
cssExpExtraCellInfo = {'','','',''};
%verify
if ~all(strcmp(ds1.genes, ds2.genes))
    error('T0003: SCDataset.fillEmpties: genes incorrect');
elseif ~isequal(ds1.data, ds2.data)
    error('T0003: SCDataset.fillEmpties: data incorrect');
elseif ~all(strcmp(ds2.cellIds, expCellIds))
    error('T0003: SCDataset.fillEmpties: cellIds incorrect');
elseif ~all(strcmp(ds2.sampleIds, expSampleIds))
    error('T0003: SCDataset.fillEmpties: sampleIds incorrect');
elseif ~isequal(ds2.cellType, expCellType)
    error('T0003: SCDataset.fillEmpties: cellType incorrect');
elseif ~isequal(ds2.subCellType, expSubCellType)
    error('T0003: SCDataset.fillEmpties: subCellType incorrect');
elseif ~all(strcmp(ds2.extraCellInfo, cssExpExtraCellInfo))
    error('T0003: SCDataset.fillEmpties: extraCellInfo incorrect');
else
    disp('T0003: SCDataset.fillEmpties: ok');
end


%% T0004: SynchronizeGenes 
ds = SCDataset();
samp = Samples();

ds.genes = {'A';'B';'C';'E'};
samp.genes = {'C';'A';'D'};
ds.data = [1 2 3 4; 5 6 7 8; 9 10 11 12; 13 14 15 16];
samp.data = [20 21; 22 23; 24 25];

[resDsDisc,resSampDisc] = SynchronizeGenes(ds, samp, 1);
[resDsKeep,resSampKeep] = SynchronizeGenes(ds, samp, 0);

expDiscGenes = {'A';'C'};
expKeepGenes = {'A';'B';'C';'D';'E'};
expDsDiscData = [1 2 3 4; 9 10 11 12];
expDsKeepData = [1 2 3 4 ; 5 6 7 8; 9 10 11 12; 0 0 0 0; 13 14 15 16];
expSampDiscData = [22 23; 20 21];
expSampKeepData = [22 23; 0 0;20 21; 24 25; 0 0];
if ~all(strcmp(resDsDisc.genes, expDiscGenes))
    error('T0004: SynchronizeGenes: genes incorrect (discard)');
elseif ~all(strcmp(resSampDisc.genes, expDiscGenes))
    error('T0004: SynchronizeGenes: genes incorrect (discard)');
elseif ~all(strcmp(resDsKeep.genes, expKeepGenes))
    error('T0004: SynchronizeGenes: genes incorrect (keep)');
elseif ~all(strcmp(resSampKeep.genes, expKeepGenes))
    error('T0004: SynchronizeGenes: genes incorrect (keep)');
elseif ~isequal(resDsDisc.data, expDsDiscData)
    error('T0004: SynchronizeGenes: data incorrect (discard)');
elseif ~isequal(resDsKeep.data, expDsKeepData)
    error('T0004: SynchronizeGenes: data incorrect (keep)');
elseif ~isequal(resSampDisc.data, expSampDiscData)
    error('T0004: SynchronizeGenes: data incorrect (discard)');
elseif ~isequal(resSampKeep.data, expSampKeepData)
    error('T0004: SynchronizeGenes: data incorrect (keep)');
else
    disp('T0004: SynchronizeGenes: ok');
end

%% T0015 TPM
s1 = Samples;
s1.genes = {'A';'B';'C';'D'};
s1.data = [1 2 3 4; 5 6 7 8; 9 10 0 12; 5 32 0 76];
s1.sampleIds = {'1_11','1_21','1_31','1_41'};

s2 = TPM(s1);
a3 = TPM(s1.data);

if a3(2,1) ~= 250000 || a3(3,4) ~= 120000
    error('T0015: TPM function incorrect');
elseif ~isequal(s2.data, a3)
    error('T0015: different results between matrix and object');
else
    disp('T0015: TPM ok')
end

%T0005: Multiple testing
unadj = [0.001, 0.14, 0.00454, 0.000565, 0.23, 0.0041, 0.0056];
adj = AdjustPval(unadj,'benjamini');
expAdj = [0.0035,0.1633333,0.00784,0.0035,0.23,0.00784,0.00784];%calculated using the online calculator https://www.sdmproject.com/utilities/?show=FDR
if max(abs(adj-expAdj)) < 10^-7%this reflects the number of decimals produced by the tool
    disp('T0005: Multiple testing Benjamini: ok');
else
    error('T0005: Multiple testing Benjamini failed');
end

%% T0025 LogTrans
logTestData = [0 9 99 999; 9999 99999 999999 9999999];
expLogRes = [0 1 2 3; 4 5 6 7];
logRes = LogTrans(logTestData, true);
restoredLogData = LogTrans(logRes, false); 
if ~(isequal(logTestData, restoredLogData) && isequal(logRes, expLogRes))
    disp('T0025: LogTrans failed');
else
    disp('T0025: LogTrans: ok');
end


%T0016: LogMultinomialPDF
%tests that this function gives the same
%results as the built-in matlab function (which can't be used for large
%number of bins due to 64-bit double overflow)
obs = [1;1;4;2;0;1];
prob = [0.2;0.2;0.1;0.3;0.1;0.1];
y = mnpdf(obs,prob);
a = log(y);
ll = LogMultinomialPDF(obs, prob, 3);
if abs(a - ll) > 0.00001 %check that they are the same, with some reasonable round-off error
    error('T0016: LogMultinomialPDF failed');
else
    disp('T0016: LogMultinomialPDF: ok');
end


%T0017 CreateVennDiagramSets
s1in = {'A','B','C','E','G','H','I','J','K'};
s2in = {'B','C','D','E','G','H','J','L','M'};
s3in = {'A','B','C','F','G','H','I','J','L','N','O','P'};
[s1,s2,s3,s1s2,s1s3,s2s3,s1s2s3] = CreateVennDiagramSets(s1in,s2in,s3in);
exps1 = {'K'};
exps2 = {'D','M'};
exps3 = {'F','N','O','P'};
exps1s2 = {'E'};
exps1s3 = {'A','I'};
exps2s3 = {'L'};
exps1s2s3 = {'B','C','G','H','J'};
if ~(all(strcmp(s1, exps1)) && all(strcmp(s2, exps2)) && all(strcmp(s3, exps3)) && all(strcmp(s1s2, exps1s2)) && all(strcmp(s1s3, exps1s3)) && all(strcmp(s2s3, exps2s3)) && all(strcmp(s1s2s3, exps1s2s3)))
    error('T0017: CreateVennDiagramSets failed');
else
    disp('T0017: CreateVennDiagramSets ok');
end




%% T0006: LC dataset
[lc,hlc] = SCDep.scd_lc;
%cell AAATCCCTCTAGAC_1, gene NOC2L should be 2
geneSel = strcmp(lc.genes, 'NOC2L');
cellSel = strcmp(lc.cellIds,'AAATCCCTCTAGAC_1');
if lc.data(geneSel,cellSel) ~= 2
    error('T0006: LC dataset data incorrect');
elseif sum(strcmp('AAACATACTCAGAC_6',lc.cellIds)) ~= 0 || sum(strcmp('AAACATACTCAGAC_6',hlc.cellIds)) ~= 1
    error('T0006: LC dataset tumor/healthy incorrect');
elseif lc.cellType(strcmp(lc.cellIds,'TTTGGTTTCGAGAGCA_15')) ~= Celltype.BCell
    error('T0006: LC dataset cell type incorrect');
elseif ~strcmp(lc.sampleIds(strcmp(lc.cellIds,'TTTGGTTTCGAGAGCA_15')), '3')
    error('T0006: LC dataset sample id incorrect');
else
    disp('T0006: LC dataset ok')
end

%% T0007: BC2 dataset
%cell 332, gene A2M should be 3
bc = SCDep.scd_bc2;
geneSel = strcmp(bc.genes, 'A2M');
cellSel = strcmp(bc.cellIds,'332');
if bc.data(geneSel,cellSel) ~= 3
    error('T0007: BC2 dataset data incorrect');
elseif bc.cellType(strcmp(bc.cellIds,'211')) ~= Celltype.TCellCD8Pos
    error('T0007: BC dataset cell type incorrect');
elseif ~strcmp(bc.sampleIds(strcmp('8073',bc.cellIds)), 'BC2_LYMPHNODE')
    error('T0007: BC dataset tissue incorrect');
else
    disp('T0007: BC2 dataset ok')
end

%% T0008: GSE112845 dataset, test CD8 only
%cell AAACCTGAGCCTTGAT-1, gene RPL22 should be 9
[xPat,yPat,cd8] = SCDep.scd_GSE112845;
geneSel = strcmp(cd8.genes, 'RPL22');
cellSel = strcmp(cd8.cellIds,'AAACCTGAGCCTTGAT-1');
if cd8.data(geneSel,cellSel) ~= 9
    error('T0008: GSE112845 dataset data incorrect');
elseif cd8.cellType(strcmp(cd8.cellIds,'AAACCTGAGCCTTGAT-1')) ~= Celltype.TCellCD8Pos
    error('T0008: GSE112845 dataset cell type incorrect');
else
    disp('T0008: GSE112845, CD8 dataset ok')
end

%% T0009: livt dataset
%cell PTC143-0205, gene ADA should be 556
ds = SCDep.scd_livt;
geneSel = strcmp(ds.genes, 'ADA');
cellSel = strcmp(ds.cellIds,'PTC143-0205');
if ds.data(geneSel,cellSel) ~= 556
    error('T0009: livt dataset data incorrect');
elseif ds.cellType(strcmp(ds.cellIds,'PTC143-0205')) ~= Celltype.TCell
    error('T0009: livt dataset cell type incorrect');
elseif ~strcmp(ds.sampleIds(strcmp(ds.cellIds,'PTC143-0205')), '0205')
    error('T0009: livt dataset sample ids incorrect');
else
    disp('T0009: livt dataset ok')
end

%% T0010: ovasc dataset
%cell ovasc_2, gene AAAS should be 11
ds = SCDep.scd_ovasc;
geneSel = strcmp(ds.genes, 'AAAS');
cellSel = strcmp(ds.cellIds,'ovasc_2');
if ds.data(geneSel,cellSel) ~= 11
    error('T0010: ovasc dataset data incorrect');
elseif ds.cellType(strcmp(ds.cellIds,'ovasc_3069')) ~= Celltype.MacrophageOrMonocyte
    error('T0010: ovasc dataset cell type incorrect');
elseif ~strcmp(ds.sampleIds(strcmp(ds.cellIds,'ovasc_3097')), '7990M2')
    error('T0010: ovasc dataset sample ids incorrect');
else
    disp('T0010: ovasc dataset ok')
end

%% T0011: PBMC68k
%cell AAACATACACCCAA-1, gene RPL3 should be 20
ds = SCDep.scd_pbmc68000;
geneSel = strcmp(ds.genes, 'RPL3');
cellSel = strcmp(ds.cellIds,'AAACATACACCCAA-1');
if ds.data(geneSel,cellSel) ~= 20
    error('T0011: pbmc68k dataset data incorrect');
elseif ds.cellType(strcmp(ds.cellIds,'AAACGCTGTTATCC-1')) ~= Celltype.TCellReg
    error('T0011: pbmc68k dataset cell type incorrect');
else
    disp('T0011: pbmc68k dataset ok')
end

%% T0012: B10k
%cell AAACATACAATGCC-1, gene MT-CO1 should be 17
ds = SCDep.scd_pbmcb10000;
geneSel = strcmp(ds.genes, 'MT-CO1');
cellSel = strcmp(ds.cellIds,'AAACATACAATGCC-1');
if ds.data(geneSel,cellSel) ~= 17
    error('T0012: b10k dataset data incorrect');
elseif ds.cellType(strcmp(ds.cellIds,'AAACATACAATGCC-1')) ~= Celltype.BCell
    error('T0012: b10k dataset cell type incorrect');
else
    disp('T0012: b10k dataset ok')
end

%% T0013: CD4MEM
%cell AAACATACACCTAG-1, gene RPL13A should be 29
ds = SCDep.scd_pbmctcd4mem10000;
geneSel = strcmp(ds.genes, 'RPL13A');
cellSel = strcmp(ds.cellIds,'AAACATACACCTAG-1');
if ds.data(geneSel,cellSel) ~= 29
    error('T0013: cd4mem dataset data incorrect');
elseif ds.cellType(strcmp(ds.cellIds,'AAACATACACCTAG-1')) ~= Celltype.TCellCD4Memory
    error('T0013: cd4mem dataset cell type incorrect');
else
    disp('T0013: cd4mem dataset ok')
end

% T0018: AlignDataset
% This is not a full test of the function, it tests some basics. A
% complementary test is the result shown in Figure 2 B-C, which shows that
% the sampling noise of all aligned datasets is virtually the same
templInfo = DSAVEGetStandardTemplate();
ds = SCDep.scd_hca_cb;
ds = ds.cellSubset(1:2001);
ds2 = DSAVEAlignDataset(ds, templInfo);
UMIDistr = sum(ds2.data,1);
if ~isequal(sort(UMIDistr), sort(templInfo.UMIDistr))
    error('T0018: DSAVEAlignDataset: UMIDistr not ok');%the HCA dataset is able to provide identical UMIDistr
elseif ~all(strcmp(sort(ds2.genes), sort(templInfo.geneSet)))
    error('T0018: : DSAVEAlignDataset: genes not ok');
else
    disp('T0018: DSAVEAlignDataset: ok');
end

%% T0019: DSAVEGenerateTemplateInfo
% This is not a full test of the function, it tests some basics. 
% The code is implicitly tested in Figure 2 B-C, which shows that
% the sampling noise of all aligned datasets is virtually the same
templInfo = DSAVEGetStandardTemplate();

genes = SCDep.scd_ovasc.genes;
genes = intersect(genes,SCDep.scd_bc2.genes);
genes = intersect(genes,SCDep.scd_pbmcb10000.genes);
[scd_GSE112845_pat_a,scd_GSE112845_pat_b,scd_GSE112845_cd8] = SCDep.scd_GSE112845;
genes = intersect(genes,scd_GSE112845_cd8.genes);

if mean(templInfo.UMIDistr) ~= 750
    error('T0019: DSAVEGenerateTemplateInfo: UMIDistr not ok');
elseif ~all(strcmp(sort(genes), sort(templInfo.geneSet)))
    error('T0019: DSAVEGenerateTemplateInfo: genes not ok');
else
    disp('T0019: DSAVEGenerateTemplateInfo: ok');
end

% T0020: DSAVEGenerateSNODataset
% This is not a full test of the function, it tests some basics. 
% The code is implicitly verified since Figure 3A in 
% supplementary information shows that different datasets with the same 
% noise get similar scores, and that this function is used for generating them.
ds = SCDep.scd_ovasc;
ds2 = DSAVEGenerateSNODataset(ds);
UMIs1 = sum(ds.data,1);
UMIs2 = sum(ds2.data,1);

if ~isequal(UMIs1, UMIs2)
    error('T0020: DSAVEGenerateSNODataset: UMIDistr not ok');
else
    disp('T0020: DSAVEGenerateSNODataset: ok');
end

%% T0021: DSAVEGetSingleCellDivergence
ds = SCDataset;
ds.data = [0 1 2 1;3 3 1 2];
ds.genes = {'A';'B'};
ds = ds.fillEmpties();
lls = DSAVEGetSingleCellDivergence(ds, 3);
ds2 = ds;
ds2.data = [0 14 2 1;17 0 1 2];
lls2 = DSAVEGetSingleCellDivergence(ds2, 3);
ds3 = ds;
ds3.data = [0 3 2 1;3 0 1 2];
lls3 = DSAVEGetSingleCellDivergence(ds3, 3);


if ~((lls(1,4) > lls(1,3)) && (lls(1,4) > lls(1,1)) && (lls(1,4) >= lls(1,2))) % 2 could become the same or worse than 4, depending on which reads gets discarded
    error('T0021: DSAVEGenerateSNODataset: UMIDistr not ok');
elseif ~isequal(lls2, lls3)
    error('T0021: DSAVEGenerateSNODataset: down-sampling not ok');
else
    disp('T0020: DSAVEGenerateSNODataset: ok');
end

%% T0022: DSAVEGetTotalVariationFromBulk
% We can ignore that the data should be somewhat TPM here.
s = Samples;
s.data = [1 1 1 1 2 2 2 2;300 300 300 300 300 300 300 300; 0 0 0 0 0 0 0 0];
s.genes = {'A';'B';'C'};
s = s.fillEmpties();

%1 vs 1:
t1vs1 = DSAVEGetTotalVariationFromBulk(s, false, 250, 0.5);
%calculate (each value will be compared with 3 samples of the same value and 4 samples with the double value)
exp1vs1 = log(2.05/1.05)*4/7;

%4 vs 4
t4vs4 = DSAVEGetTotalVariationFromBulk(s, true, 250, 0.5);
%calculate this differently than in the function, with loops :)
%calculate the distribution over number of ones (vs twos) in the first of
%the two sets for each combination
numones = zeros(1,5);%represents 0 1 2 3 4 ones, i.e. the first the number of combinations with 0 ones, etc
for i = 1:5
   for j = i+1:6
       for k = j+1:7
           for m = k+1:8
               index = (i<5)+(j<5)+(k<5)+(m<5)+1;%represents number of ones in the combination
               numones(1,index) = numones(1,index) + 1;
           end
       end
   end
end
a = numones./sum(numones);%find number to scale each comb type

exp4vs4 = log(2.05/1.05)*a(1,1)*2 + log(1.80/1.30)*a(1,2)*2;

if abs(exp1vs1 - t1vs1) > 10^-10 % compensate for round off effects
    error('T0022: DSAVEGetTotalVariationFromBulk: 1 vs 1 not ok');
elseif abs(exp4vs4 - t4vs4) > 10^-10 % compensate for round off effects
    error('T0022: DSAVEGetTotalVariationFromBulk: 4 vs 4 not ok');
else
    disp('T0020: DSAVEGetTotalVariationFromBulk: ok');
end

%% T0023: DSAVEGetTotalVariationVsPoolSize
% Create a very simple dataset with only two cells and one valid gene, this should be
% deterministic with the same value for all points! This does not fully
% test the function, but at least partly.
ds = SCDataset;
ds.data = [1 2;999998.9 999997.9; 0.1 0.1];
ds.genes = {'A';'B';'C'};
ds = ds.fillEmpties();
vals = DSAVEGetTotalVariationVsPoolSize(ds, 1, 100, 0.5);
vals = vals(2,:);%throw away the x:es, we don't need them here
expVal = log(2.05/1.05);
expVals = repmat(expVal,1,100);
% compensate for round off effects:
expVals = round(expVals,10);
vals = round(vals,10);

if ~isequal(expVals,vals)
    error('T0023: DSAVEGetTotalVariationVsPoolSize: not ok');
else
    disp('T0023: DSAVEGetTotalVariationVsPoolSize: ok');
end

% T0024 is the progress bar test case


%% T0026: DSAVEGetGeneVariation
%Create a dataset where we can easily calculate the p value
%and see that it gets somewhat right. 
%There is a risk that this test will have problems with stability, since
%there are random numbers involved
dsTest = SCDataset;
dsTest.genes = {'A';'B';'C'};
dsTest.cellIds = {'1', '2', '3', '4', '5'};
dsTest.data = [0 4 0 0 0; 0 0 0 0 0; 4 0 4 4 4];
%the number of combinations in total are 5^4 = 625
%the number of combinations with 4 in the same pile is 5
%the p value should thus be 5/625 = 0.0080
[genesTest, logCVTest,pvalsTest] = DSAVEGetGeneVariation(dsTest,0,100000,10000);
expGenes = {'A';'C'};
%calculate the expected CV
CVDS = log(std(dsTest.data(1,:) .* 250000, [], 2) ./ (200000 + 0.05) + 1);
%calculate the expected SNO CV mean:

%generate all combinations (not all are unique!) using a loop; if all combinations are there once
%it will be the true mean
avgSNO = zeros(625,5);
%the loop variables represents the cell index where each of the four counts should be
%placed
for i = 1:5
    for j = 1:5
        for k = 1:5
            for m = 1:5
                index = (i-1)*125 + (j-1)*25 + (k-1)*5 + m;
                avgSNO(index,i) = avgSNO(index,i) + 1;
                avgSNO(index,j) = avgSNO(index,j) + 1;
                avgSNO(index,k) = avgSNO(index,k) + 1;
                avgSNO(index,m) = avgSNO(index,m) + 1;
            end
        end
    end
end

CVSNO = log(std(avgSNO .* 250000, [], 2) ./ (200000 + 0.05) + 1);
CVDiffExp = CVDS - mean(CVSNO,1);
%now test with different number of UMIs
%here, the result should be 1/1000 = 0.001
dsTest2 = SCDataset;
dsTest2.genes = {'A';'B';'C'};
dsTest2.cellIds = {'1', '2', '3', '4', '5'};
dsTest2.data = [0 3 0 0 0; 0 0 0 0 0; 6 0 6 6 9];
[~, ~,pvalsTest2] = DSAVEGetGeneVariation(dsTest2,0,100000,10000);

if abs(pvalsTest(1,1) - 0.008) > 0.002
    error('T0026: DSAVEGetGeneVariation: not ok, pVal not right');
elseif abs(pvalsTest2(1,1) - 0.001) > 0.0007
    error('T0026: DSAVEGetGeneVariation: not ok, pVal2 not right');
elseif ~all(strcmp(genesTest, expGenes))
    error('T0026: DSAVEGetGeneVariation: not ok, genes not correct');
elseif abs(logCVTest - CVDiffExp) > 0.01
    error('T0026: DSAVEGetGeneVariation: not ok, pVal2 not right');
else
    disp('T0026: DSAVEGetTotalVariationVsPoolSize: ok');
end


