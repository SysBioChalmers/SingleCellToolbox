SingleCellToolbox:

SingleCellToolbox is a matlab toolbox for single-cell RNA-Seq data.

Installation:
The only installation needed is to clone the git repo and set the path. Recommended file structure:
SingleCellToolbox
	SingleCellToolbox - this is the repo
	ImportableData - This is where all the existing public dataset files go
	TempData - Will be created by the toolbox and used for .mat files for caching purposes

Setting the path and cleaning up once you don't want this package anymore can be done by using the commands:
SingleCellToolboxInstall.install
SingleCellToolboxInstall.uninstall


The toolbox includes code for importing many existing public single-cell datasets into matlab. For this to work, the 
ImportableData folder needs to be filled with files from public datasets. Where to copy this from is yet to be worked out.