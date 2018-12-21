function AS_GCA_computemain(Parameter)
pathfiles = which('AS_GCA_computemain.m');
[path nam ext] = fileparts(pathfiles);
Caltype = Parameter.Caltype;
Imgord = Parameter.Imgord;

if Caltype % map mode
    Outputdir = Parameter.Outputdir;
    RealCompPara.Outputdir = Outputdir;
    Inputdir = Parameter.Inputdir;
    RealCompPara.Inputdir = Inputdir;
    SeedROItype = Parameter.SeedROItype;
    RealCompPara.SeedROItype = SeedROItype;
    
    [vinput,datainput] = Dynamic_read_dir_NIFTI(Inputdir);
    dims = vinput.dim;
    Ttrans = vinput.mat;
    VOXELSIZE = abs([vinput.mat(1,1),vinput.mat(2,2),vinput.mat(3,3)]);
    if strcmp(Parameter.masks,'Defaults')
        maskfiles = fullfile(path,'grey.nii');
        [vmask,datamask] = Dynamic_read_dir_NIFTI(maskfiles);
        datamask = datamask>0.2;
    else
        maskfiles = Parameter.masks;
        [vmask,datamask] = Dynamic_read_dir_NIFTI(maskfiles);
    end
    if any(vmask.dim-vinput.dim)
        Maskimg0 = fullfile(Outputdir,'mask.nii');
        %     (PI,PO,NewVoxSize,hld,TargetSpace)
        dynamicBC_Reslice(maskfiles,Maskimg0,vinput.dim,0,vinput.fname);
        clear vmask datamask
        [vmask,datamask] = Dynamic_read_dir_NIFTI(Maskimg0);
        if strcmp(Parameter.masks,'Defaults')
            datamask = datamask>0.2;
        end
        MaskimgU = fullfile(Outputdir,'maskused.nii');
        DATMASK = reshape(datamask,dims);
        DynamicBC_write_NIFTI(DATMASK,vmask,MaskimgU);
    end
    DATMASK = reshape(datamask,dims);
    switch SeedROItype
        case 1 % MNI condition
            mnis = Parameter.SeedROI_mni;
            radius = Parameter.SeedROI_radius;
            indexs = mniroi(mnis,radius,Ttrans,dims,VOXELSIZE,DATMASK);
            ROIsignals = mean(datainput(indexs,:),1)';
            DATROI = zeros(dims);
            DATROI(indexs) = 1;
            ROIname = fullfile(Outputdir,'roi.nii');
            DynamicBC_write_NIFTI(DATROI,vmask,ROIname);
            Signalname = fullfile(Outputdir,'ROIsignal.mat');
            save(Signalname,'ROIsignals');
        case 2
            rois = Parameter.SeedROI_nifti;
            [vroi,dataroi] = Dynamic_read_dir_NIFTI(rois);
            if any(vroi.dim-vinput.dim)
                Roisimg0 = fullfile(Outputdir,'roi.nii');
                dynamicBC_Reslice(rois,Roisimg0,vinput.dim,0,vinput.fname);
                [vroi,dataroi] = Dynamic_read_dir_NIFTI(Roisimg0);
            end
            dat00 = unique(dataroi);
            nrois = length(unique(dataroi))-1;
            for i = 1:nrois
                indexs = find(dataroi==dat00(i+1));
                ROIsignals(:,i) = mean(datainput(indexs,:),1)';
            end
            Signalname = fullfile(Outputdir,'ROIsignal.mat');
            save(Signalname,'ROIsignals');
        case 3
            mats = Parameter.SeedROI_mat;
            temps = load(Parameter.SeedROI_mat);
            varname = Parameter.SeedROI_varname;
            eval(['ROIsignals = temps.',varname,';']);
            Signalname = fullfile(Outputdir,'ROIsignal.mat');
            save(Signalname,'ROIsignals');
        case 4
            txts = Parameter.SeedROI_txt;
            ROIsignals = load(txts);
            Signalname = fullfile(Outputdir,'ROIsignal.mat');
            save(Signalname,'ROIsignals');
        otherwise
            disp('error ROI defined')
    end    
    COVcond = Parameter.covs;
    RealCompPara.COVcond = COVcond;
    COVdir = Parameter.COVtext;
    RealCompPara.COVdir = COVdir;
    if COVcond==1
        COV_S = load(COVdir);
    end
    Npoint = size(ROIsignals,1);
    if strcmp(Imgord,'Defaults')
        IOrds = 1:Npoint;
    else
        IOrdst = load(Imgord);
        [ix,IOrds] = sort(IOrdst);
    end
    maskedSignal = datainput(find(DATMASK),:)';
    maskedsigname = fullfile(Outputdir,'maskedSignal.mat');
    save(maskedsigname,'maskedSignal','DATMASK');
    calmethod = Parameter.Calmethod;
    Order = Parameter.GCAORD;
    if calmethod % residual
        for i = 1:size(ROIsignals,2)
            if i<10
                Outfilenametemp1 = fullfile(Outputdir,['GCAx2y_ROI00000',num2str(i),'.nii']);
                Outfilenametemp2 = fullfile(Outputdir,['GCAy2x_ROI00000',num2str(i),'.nii']);
                PmapOutfilenametemp1 = fullfile(Outputdir,['Pmap_GCAx2y_ROI00000',num2str(i),'.nii']);
                PmapOutfilenametemp2 = fullfile(Outputdir,['Pmap_GCAy2x_ROI00000',num2str(i),'.nii']);
                Outfilenametemp3 = fullfile(Outputdir,['GCAx2y_transformed_ROI00000',num2str(i),'.nii']);
                Outfilenametemp4 = fullfile(Outputdir,['GCAy2x_transformed_ROI00000',num2str(i),'.nii']);
                Outfilenametemp5 = fullfile(Outputdir,['NetFx2y_ROI00000',num2str(i),'.nii']);
            elseif i<100
                Outfilenametemp1 = fullfile(Outputdir,['GCAx2y_ROI0000',num2str(i),'.nii']);
                Outfilenametemp2 = fullfile(Outputdir,['GCAy2x_ROI0000',num2str(i),'.nii']);
                PmapOutfilenametemp1 = fullfile(Outputdir,['Pmap_GCAx2y_ROI0000',num2str(i),'.nii']);
                PmapOutfilenametemp2 = fullfile(Outputdir,['Pmap_GCAy2x_ROI0000',num2str(i),'.nii']);
                Outfilenametemp3 = fullfile(Outputdir,['GCAx2y_transformed_ROI0000',num2str(i),'.nii']);
                Outfilenametemp4 = fullfile(Outputdir,['GCAy2x_transformed_ROI0000',num2str(i),'.nii']);
                Outfilenametemp5 = fullfile(Outputdir,['NetFx2y_ROI0000',num2str(i),'.nii']);
            elseif i<1000
                Outfilenametemp1 = fullfile(Outputdir,['GCAx2y_ROI000',num2str(i),'.nii']);
                Outfilenametemp2 = fullfile(Outputdir,['GCAy2x_ROI000',num2str(i),'.nii']);
                PmapOutfilenametemp1 = fullfile(Outputdir,['Pmap_GCAx2y_ROI000',num2str(i),'.nii']);
                PmapOutfilenametemp2 = fullfile(Outputdir,['Pmap_GCAy2x_ROI000',num2str(i),'.nii']);
                Outfilenametemp3 = fullfile(Outputdir,['GCAx2y_transformed_ROI000',num2str(i),'.nii']);
                Outfilenametemp4 = fullfile(Outputdir,['GCAy2x_transformed_ROI000',num2str(i),'.nii']);
                Outfilenametemp5 = fullfile(Outputdir,['NetFx2y_ROI000',num2str(i),'.nii']);
            else
                Outfilenametemp1 = fullfile(Outputdir,['GCAx2y_ROI00',num2str(i),'.nii']);
                Outfilenametemp2 = fullfile(Outputdir,['GCAy2x_ROI00',num2str(i),'.nii']);
                PmapOutfilenametemp1 = fullfile(Outputdir,['Pmap_GCAx2y_ROI00',num2str(i),'.nii']);
                PmapOutfilenametemp2 = fullfile(Outputdir,['Pmap_GCAy2x_ROI00',num2str(i),'.nii']);
                Outfilenametemp3 = fullfile(Outputdir,['GCAx2y_transformed_ROI00',num2str(i),'.nii']);
                Outfilenametemp4 = fullfile(Outputdir,['GCAy2x_transformed_ROI00',num2str(i),'.nii']);
                Outfilenametemp5 = fullfile(Outputdir,['NetFx2y_ROI00',num2str(i),'.nii']);
            end
            if COVcond==1
                theCovariables = COV_S;
            else
                theCovariables = ones(size(ROIsignals,1),1);
            end
%             save(fullfile(Outputdir,)
            [ResultMap1,ResultMap2,ResultMap3,ResultMap4,ResultMap5] = restgca_residual(ROIsignals(IOrds,i), maskedSignal(IOrds,:),Order,theCovariables(IOrds,:));
            PvalResultMap1 = wgr_pwGC_F(ResultMap1,size(ROIsignals,1),Order);
            PvalResultMap2 = wgr_pwGC_F(ResultMap2,size(ROIsignals,1),Order);
            DataOut = zeros(dims);
            DataOut(find(DATMASK)) = ResultMap1;
            DynamicBC_write_NIFTI(DataOut,vmask,Outfilenametemp1);
            DataOut = zeros(dims);
            DataOut(find(DATMASK)) = ResultMap2;
            DynamicBC_write_NIFTI(DataOut,vmask,Outfilenametemp2);
            DataOut = zeros(dims);
            DataOut(find(DATMASK)) = ResultMap3;
            DynamicBC_write_NIFTI(DataOut,vmask,Outfilenametemp3);
            DataOut = zeros(dims);
            DataOut(find(DATMASK)) = ResultMap4;
            DynamicBC_write_NIFTI(DataOut,vmask,Outfilenametemp4);
            DataOut = zeros(dims);
            DataOut(find(DATMASK)) = ResultMap5;
            DynamicBC_write_NIFTI(DataOut,vmask,Outfilenametemp5);
            
            
            DataOut = zeros(dims);
            DataOut(find(DATMASK)) = PvalResultMap1;
            DynamicBC_write_NIFTI(DataOut,vmask,PmapOutfilenametemp1);
            DataOut = zeros(dims);
            DataOut(find(DATMASK)) = PvalResultMap2;
            DynamicBC_write_NIFTI(DataOut,vmask,PmapOutfilenametemp2);
        end
    else % coef
        for i = 1:size(ROIsignals,2)
            if i<10
                Outfilenametemp1 = ['Coef_GCAx2y_ROI00000',num2str(i),'.nii'];
                Outfilenametemp2 = ['Coef_GCAy2x_ROI00000',num2str(i),'.nii'];
                Outfilenametemp3 = ['Coef_GCAx2y_AR_ROI00000',num2str(i),'.nii'];
                Outfilenametemp4 = ['Coef_GCAy2x_AR_ROI00000',num2str(i),'.nii'];
            elseif i<100
                Outfilenametemp1 = ['Coef_GCAx2y_ROI0000',num2str(i),'.nii'];
                Outfilenametemp2 = ['Coef_GCAy2x_ROI0000',num2str(i),'.nii'];
                Outfilenametemp3 = ['Coef_GCAx2y_AR_ROI0000',num2str(i),'.nii'];
                Outfilenametemp4 = ['Coef_GCAy2x_AR_ROI0000',num2str(i),'.nii'];
            elseif i<1000
                Outfilenametemp1 = ['Coef_GCAx2y_ROI000',num2str(i),'.nii'];
                Outfilenametemp2 = ['Coef_GCAy2x_ROI000',num2str(i),'.nii'];
                Outfilenametemp3 = ['Coef_GCAx2y_AR_ROI000',num2str(i),'.nii'];
                Outfilenametemp4 = ['Coef_GCAy2x_AR_ROI000',num2str(i),'.nii'];
            else
                Outfilenametemp1 = ['Coef_GCAx2y_ROI00',num2str(i),'.nii'];
                Outfilenametemp2 = ['Coef_GCAy2x_ROI00',num2str(i),'.nii'];
                Outfilenametemp3 = ['Coef_GCAx2y_AR_ROI00',num2str(i),'.nii'];
                Outfilenametemp4 = ['Coef_GCAy2x_AR_ROI00',num2str(i),'.nii'];
            end
            
            if COVcond==1
                theCovariables = COV_S;
            else
                theCovariables = ones(size(ROIsignals,1),1);
            end
            [ResultMap1,ResultMap2,ResultMap3,ResultMap4]=restgca_coefficient(ROIsignals(IOrds,i), maskedSignal(IOrds,:),Order,theCovariables(IOrds,:));
            for j = 1:Order
                DataOut = zeros(dims);
                DataOut(find(DATMASK)) = ResultMap1{j};
                Outfiletemp = fullfile(Outputdir,['Order',num2str(j),'_',Outfilenametemp1]);
                DynamicBC_write_NIFTI(DataOut,vmask,Outfiletemp);
                DataOut = zeros(dims);
                DataOut(find(DATMASK)) = ResultMap2{j};
                Outfiletemp = fullfile(Outputdir,['Order',num2str(j),'_',Outfilenametemp2]);
                DynamicBC_write_NIFTI(DataOut,vmask,Outfiletemp);
                DataOut = zeros(dims);
                DataOut(find(DATMASK)) = ResultMap3{j};
                Outfiletemp = fullfile(Outputdir,['Order',num2str(j),'_',Outfilenametemp3]);
                DynamicBC_write_NIFTI(DataOut,vmask,Outfiletemp);
                DataOut = zeros(dims);
                DataOut(find(DATMASK)) = ResultMap4{j};
                Outfiletemp = fullfile(Outputdir,['Order',num2str(j),'_',Outfilenametemp4]);
                DynamicBC_write_NIFTI(DataOut,vmask,Outfiletemp);
            end
        end
    end
else    % matrix mode
    Outputdir = Parameter.Outputdir;
    RealCompPara.Outputdir = Outputdir;
    Inputdir = Parameter.Inputdir;
    RealCompPara.Inputdir = Inputdir;
    SeedROItype = Parameter.SeedROItype;
    RealCompPara.SeedROItype = SeedROItype;
    
    [vinput,datainput] = Dynamic_read_dir_NIFTI(Inputdir);
    dims = vinput.dim;
    Ttrans = vinput.mat;
    VOXELSIZE = abs([vinput.mat(1,1),vinput.mat(2,2),vinput.mat(3,3)]);
    switch SeedROItype
        case 1 % MNI condition
            roipath = Parameter.SeedROI_seproi;
            [vroi,dataroi,roinamelist] = Dynamic_read_dir_NIFTI(roipath);
            if any(vroi.dim-vinput.dim)
                error('wrong dimension for the ROI defination!');
            end
            nrois = size(dataroi,2);
            for i = 1:nrois
                indexs = find(dataroi(:,i));
                ROIsignals(:,i) = mean(datainput(indexs,:),1)';
            end
            Signalname = fullfile(Outputdir,'ROIsignal.mat');
            Signalnamelist = fullfile(Outputdir,'roinamelist.mat');
            save(Signalname,'ROIsignals');
            save(Signalnamelist,'roinamelist');
        case 2
            rois = Parameter.SeedROI_nifti;
            [vroi,dataroi] = Dynamic_read_dir_NIFTI(rois);
            if any(vroi.dim-vinput.dim)
                Roisimg0 = fullfile(Outputdir,'roi.nii');
                dynamicBC_Reslice(rois,Roisimg0,vinput.dim,0,vinput.fname);
                [vroi,dataroi] = Dynamic_read_dir_NIFTI(Roisimg0);
            end
            dat00 = unique(dataroi);
            nrois = length(unique(dataroi))-1;
            for i = 1:nrois
                indexs = find(dataroi==dat00(i+1));
                ROIsignals(:,i) = mean(datainput(indexs,:),1)';
            end
            Signalname = fullfile(Outputdir,'ROIsignal.mat');
            save(Signalname,'ROIsignals');
        case 3
            mats = Parameter.SeedROI_mat;
            temps = load(Parameter.SeedROI_mat);
            varname = Parameter.SeedROI_varname;
            eval(['ROIsignals = temps.',varname,';']);
            Signalname = fullfile(Outputdir,'ROIsignal.mat');
            save(Signalname,'ROIsignals');
        case 4
            txts = Parameter.SeedROI_txt;
            ROIsignals = load(txts);
            Signalname = fullfile(Outputdir,'ROIsignal.mat');
            save(Signalname,'ROIsignals');
        otherwise
            disp('error ROI defined')
    end
    
    COVcond = Parameter.covs;
    RealCompPara.COVcond = COVcond;
    COVdir = Parameter.COVtext;
    RealCompPara.COVdir = COVdir;
    calmethod = Parameter.Calmethod;
    Order = Parameter.GCAORD;
    
    
    Npoint = size(ROIsignals,1);
    if strcmp(Imgord,'Defaults')
        IOrds = 1:Npoint;
    else
        IOrdst = load(Imgord);
        [ix,IOrds] = sort(IOrdst);
    end
    save(Signalname,'ROIsignals','IOrds');
    if COVcond==1
        COV_S = load(COVdir);
    end
    if calmethod % residual
        Outfilename1 =  fullfile(Outputdir,'GCAx2y.mat');
        Outfilename2 =  fullfile(Outputdir,'GCAy2x.mat');
        pOutfilename1 =  fullfile(Outputdir,'PvalGCAx2y.mat');
        pOutfilename2 =  fullfile(Outputdir,'PvalGCAy2x.mat');
        Outfilename3 =  fullfile(Outputdir,'GCAx2y_transformed.mat');
        Outfilename4 =  fullfile(Outputdir,'GCAy2x_transformed.mat');
        Outfilename5 =  fullfile(Outputdir,'NetFx2y.mat');
        if COVcond==1
            theCovariables = COV_S;
        else
            theCovariables = ones(size(ROIsignals,1),1);
        end
        for i = 1:size(ROIsignals,2)-1
            for j = i+1:size(ROIsignals,2)
                ROIsignalst = [ROIsignals(:,i),ROIsignals(:,j)];
                [ResultMap1,ResultMap2,ResultMap3,ResultMap4,ResultMap5] = restgca_FROI(ROIsignalst(IOrds,:),Order,theCovariables(IOrds,:));
                X2Y(i,j) =  ResultMap1;%X2Y(j,i) =  ResultMap1;
                Y2X(i,j) =  ResultMap2;%Y2X(j,i) =  ResultMap2;
                Pval_X2Y(i,j) = wgr_pwGC_F(X2Y(i,j),size(ROIsignals,1),Order);
                Pval_Y2X(i,j) = wgr_pwGC_F(Y2X(i,j),size(ROIsignals,1),Order);
                X2Y_transform(i,j) =  ResultMap3;%X2Y_transform(j,i) =  ResultMap3;
                Y2X_transform(i,j) =  ResultMap4;%Y2X_transform(j,i) =  ResultMap4;
                NetFx2y(i,j) = ResultMap5; NetFx2y(j,i) = ResultMap5;
            end
        end
        save(Outfilename1,'X2Y');
        save(Outfilename2,'Y2X');
        save(Outfilename3,'X2Y_transform');
        save(Outfilename4,'Y2X_transform');
        save(Outfilename5,'NetFx2y');
        save(pOutfilename1,'Pval_X2Y');
        save(pOutfilename2,'Pval_Y2X');
    else % coef
        Outfilename1 = fullfile(Outputdir,'Coef_GCAx2y.mat');
        Outfilename2 = fullfile(Outputdir,'Coef_GCAy2x.mat');
        Outfilename3 = fullfile(Outputdir,'Coef_GCAx2x.mat');
        Outfilename4 = fullfile(Outputdir,'Coef_GCAy2y.mat');
        
        if COVcond==1
            theCovariables = COV_S;
        else
            theCovariables = ones(size(ROIsignals,1),1);
        end
        
        
        for i = 1:size(ROIsignals,2)-1
            for j = i+1:size(ROIsignals,2)
                ROIsignalst = [ROIsignals(:,i),ROIsignals(:,j)];
                [Result_X2Y,Result_Y2X,ROI_sequence] = restgca_CROI_Bivariate(ROIsignalst(IOrds,:),Order,theCovariables(IOrds,:));
                X2Y(i,j,:) =  Result_X2Y(1:Order);X2X(i,j,:) =  Result_X2Y(Order+1:Order*2);
                Y2X(i,j,:) =  Result_Y2X(1:Order);Y2Y(i,j,:) =  Result_Y2X(Order+1:Order*2);
%                 X2Y_transform(i,j) =  ResultMap3;X2Y_transform(j,i) =  ResultMap3;
%                 Y2X_transform(i,j) =  ResultMap4;Y2X_transform(j,i) =  ResultMap4;
%                 NetFx2y(i,j) = ResultMap5; NetFx2y(j,i) = ResultMap5;
            end
        end
%         [Result_X2Y,Result_Y2X,ROI_sequence] = restgca_CROI_Bivariate(ROIsignals,Order,theCovariables);
        save(Outfilename1,'X2Y');
        save(Outfilename2,'Y2X');
        save(Outfilename3,'X2X');
        save(Outfilename4,'Y2Y');
    end
end
disp('Compute Finished!')
uiwait(msgbox('Compute Finished!'));
end

function indexs = mniroi(mnis,radius,Ttrans,dims,VOXELSIZE,DATAMASK)
coordinate = mni2cor(mnis, Ttrans);
sizes = round(radius./VOXELSIZE);
MATS = zeros(dims);
if coordinate(1)-sizes(1)>0
    Xrange(1) = coordinate(1)-sizes(1);
else
    Xrange(1) = 1;
end
if coordinate(1)+sizes(1)<=dims(1)
    Xrange(2) = coordinate(1)+sizes(1);
else
    Xrange(2) = dims(1);
end
%
if coordinate(2)-sizes(2)>0
    Yrange(1) = coordinate(2)-sizes(2);
else
    Yrange(1) = 1;
end
if coordinate(2)+sizes(2)<=dims(2)
    Yrange(2) = coordinate(2)+sizes(2);
else
    Yrange(2) = dims(2);
end
%
if coordinate(3)-sizes(3)>0
    Zrange(1) = coordinate(3)-sizes(3);
else
    Zrange(1) = 1;
end
if coordinate(3)+sizes(3)<=dims(3)
    Zrange(2) = coordinate(3)+sizes(3);
else
    Zrange(2) = dims(3);
end
MATS(Xrange(1):Xrange(2),Yrange(1):Yrange(2),Zrange(1):Zrange(2)) = 1;
MATS = MATS.*DATAMASK;
ind = find(MATS~=0);
[ix iy iz] = ind2sub(dims,ind);
exportnum = 1;
indexs = [];
for i = 1:length(ind)
    if sqrt(((ix(i)-coordinate(1))*VOXELSIZE(1))^2+...
            ((iy(i)-coordinate(2))*VOXELSIZE(2))^2+...
            ((iz(i)-coordinate(3))*VOXELSIZE(3))^2)<radius
        indexs(exportnum) = ind(i);
        exportnum = exportnum+1;
    end
end
if isempty(indexs)
   error('wrong defination of ROI(MNI)') 
end
end