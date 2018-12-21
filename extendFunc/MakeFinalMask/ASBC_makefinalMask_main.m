function ASBC_makefinalMask_main(CalPara)
Outdir = CalPara.Outdir;
Indat = CalPara.Indat;
Inroi = CalPara.Inroi;
stnum = 0;
for i = 1:size(Indat,1)
    [v,dat,namelist] = Dynamic_read_dir_NIFTI(Indat{i});
    V_DAT{i} = v;
    if i==1
        DAT_DAT = sum(logical(dat),2);
    else
        try
            DAT_DAT = DAT_DAT+sum(logical(dat),2);
        catch
            error(['different Input Data dims: Group ',num2str(i)]);
        end
    end
    stnum = stnum+size(dat,2);
    NAMELIST_namelist{i} = namelist;
end
for i = 1:size(Inroi,1)
    [vroi,datroi] = Dynamic_read_dir_NIFTI(Inroi{i});
    V_ROI{i} = vroi;
    try
        ROI_DAT(:,i) = datroi;
    catch
        error(['different Input ROI dims: ROI ',num2str(i)]);
    end
end
save(fullfile(Outdir,'Datainfor.mat'),'V_DAT','NAMELIST_namelist','V_ROI');
if size(DAT_DAT,1)~=size(ROI_DAT,1)
%     (PI,PO,NewVoxSize,hld,TargetSpace)
%     dynamicBC_Reslice(maskfiles,Maskimg0,vinput.dim,0,vinput.fname);
    mkdir([Outdir,filesep,'ReslicedROI']);
    clear ROI_DAT V_ROI;
    for i = 1:size(Inroi,1)
        nametemp = Inroi{i};
        [pat,nam,ext] = fileparts(nametemp);
        nametempO = fullfile([Outdir,filesep,'ReslicedROI'],[nam,ext]);
        dynamicBC_Reslice(nametemp,nametempO,V_DAT{1}.dim,0,V_DAT{1}.fname);
        [vroi,datroi] = Dynamic_read_dir_NIFTI(nametempO);
        V_ROI{i} = vroi;
        ROI_DAT(:,i) = datroi;
    end
    save(fullfile(Outdir,'Datainfor.mat'),'V_DAT','NAMELIST_namelist','V_ROI');
end
NumAll = stnum;
DAT_Allsum = DAT_DAT;
DAT_AllC = (DAT_Allsum==NumAll);
DAT_All99 = (DAT_Allsum>=(NumAll*0.99));
DAT_All95 = (DAT_Allsum>=(NumAll*0.95));
DAT_All90 = (DAT_Allsum>=(NumAll*0.90));
save(fullfile(Outdir,'DATsuit.mat'),'DAT_Allsum');
mkdir([Outdir,filesep,'Percent100']);
mkdir([Outdir,filesep,'Percent99']);
mkdir([Outdir,filesep,'Percent95']);
mkdir([Outdir,filesep,'Percent90']);
DAT_AllCre = reshape(DAT_AllC,V_DAT{1}.dim(1),V_DAT{1}.dim(2),V_DAT{1}.dim(3));
DynamicBC_write_NIFTI(DAT_AllCre,V_DAT{1},fullfile([Outdir,filesep,'Percent100'],'DAT_All.nii'));
for i = 1:size(ROI_DAT,2)
    ROI_DATtemp = ROI_DAT(:,i).*DAT_AllC;
    ROI_Allre = reshape(ROI_DATtemp,V_DAT{1}.dim(1),V_DAT{1}.dim(2),V_DAT{1}.dim(3));
    nametemp = Inroi{i};
    [pat,nam,ext] = fileparts(nametemp);
    DynamicBC_write_NIFTI(ROI_Allre,V_DAT{1},fullfile([Outdir,filesep,'Percent100'],[nam,ext]));
end
DAT_AllCre99 = reshape(DAT_All99,V_DAT{1}.dim(1),V_DAT{1}.dim(2),V_DAT{1}.dim(3));
DynamicBC_write_NIFTI(DAT_AllCre99,V_DAT{1},fullfile([Outdir,filesep,'Percent99'],'DAT_All99.nii'));
for i = 1:size(ROI_DAT,2)
    ROI_DATtemp = ROI_DAT(:,i).*DAT_All99;
    ROI_Allre = reshape(ROI_DATtemp,V_DAT{1}.dim(1),V_DAT{1}.dim(2),V_DAT{1}.dim(3));
    nametemp = Inroi{i};
    [pat,nam,ext] = fileparts(nametemp);
    DynamicBC_write_NIFTI(ROI_Allre,V_DAT{1},fullfile([Outdir,filesep,'Percent99'],[nam,ext]));
end
DAT_AllCre95 = reshape(DAT_All95,V_DAT{1}.dim(1),V_DAT{1}.dim(2),V_DAT{1}.dim(3));
DynamicBC_write_NIFTI(DAT_AllCre95,V_DAT{1},fullfile([Outdir,filesep,'Percent95'],'DAT_All95.nii'));
for i = 1:size(ROI_DAT,2)
    ROI_DATtemp = ROI_DAT(:,i).*DAT_All95;
    ROI_Allre = reshape(ROI_DATtemp,V_DAT{1}.dim(1),V_DAT{1}.dim(2),V_DAT{1}.dim(3));
    nametemp = Inroi{i};
    [pat,nam,ext] = fileparts(nametemp);
    DynamicBC_write_NIFTI(ROI_Allre,V_DAT{1},fullfile([Outdir,filesep,'Percent95'],[nam,ext]));
end
DAT_AllCre90 = reshape(DAT_All90,V_DAT{1}.dim(1),V_DAT{1}.dim(2),V_DAT{1}.dim(3));
DynamicBC_write_NIFTI(DAT_AllCre90,V_DAT{1},fullfile([Outdir,filesep,'Percent90'],'DAT_All90.nii'));
for i = 1:size(ROI_DAT,2)
    ROI_DATtemp = ROI_DAT(:,i).*DAT_All90;
    ROI_Allre = reshape(ROI_DATtemp,V_DAT{1}.dim(1),V_DAT{1}.dim(2),V_DAT{1}.dim(3));
    nametemp = Inroi{i};
    [pat,nam,ext] = fileparts(nametemp);
    DynamicBC_write_NIFTI(ROI_Allre,V_DAT{1},fullfile([Outdir,filesep,'Percent90'],[nam,ext]));
end
end