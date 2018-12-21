function AS_GCA_pvalcomp_perm(Parameter)
Permnum = Parameter.Permnum;
Outdir = Parameter.Outputdir;
Caltype = Parameter.Caltype;
Order = Parameter.GCAORD;
if Caltype % map
    ROIdir = fullfile(Outdir,'ROIsignal.mat');
    load(ROIdir);
    COVcond = Parameter.covs;
    COVdir = Parameter.COVtext;
    if COVcond==1
        COV_S = load(COVdir);
    end
    maskedsignaldir = fullfile(Outdir,'maskedSignal.mat');
    load(maskedsignaldir)
    if COVcond==1
        theCovariables = COV_S;
    else
        theCovariables = ones(size(ROIsignals,1),1);
    end
    Calmethod = Parameter.Calmethod;
    if Calmethod % residual
        for i = 1:size(ROIsignals,2)
            Outputdir = Outdir;
            if i<10
                Outfilenametemp1 = fullfile(Outputdir,['GCAx2y_ROI00000',num2str(i),'.nii']);
                Outfilenametemp2 = fullfile(Outputdir,['GCAy2x_ROI00000',num2str(i),'.nii']);
                Outfilenametemp3 = fullfile(Outputdir,['GCAx2y_transformed_ROI00000',num2str(i),'.nii']);
                Outfilenametemp4 = fullfile(Outputdir,['GCAy2x_transformed_ROI00000',num2str(i),'.nii']);
                Outfilenametemp5 = fullfile(Outputdir,['NetFx2y_ROI00000',num2str(i),'.nii']);
                POutfilenametemp1 = fullfile(Outputdir,['Perm_GCAx2y_ROI00000',num2str(i),'.nii']);
                POutfilenametemp2 = fullfile(Outputdir,['Perm_GCAy2x_ROI00000',num2str(i),'.nii']);
                POutfilenametemp3 = fullfile(Outputdir,['Perm_GCAx2y_transformed_ROI00000',num2str(i),'.nii']);
                POutfilenametemp4 = fullfile(Outputdir,['Perm_GCAy2x_transformed_ROI00000',num2str(i),'.nii']);
                POutfilenametemp5 = fullfile(Outputdir,['Perm_NetFx2y_ROI00000',num2str(i),'.nii']);
            elseif i<100
                Outfilenametemp1 = fullfile(Outputdir,['GCAx2y_ROI0000',num2str(i),'.nii']);
                Outfilenametemp2 = fullfile(Outputdir,['GCAy2x_ROI0000',num2str(i),'.nii']);
                Outfilenametemp3 = fullfile(Outputdir,['GCAx2y_transformed_ROI0000',num2str(i),'.nii']);
                Outfilenametemp4 = fullfile(Outputdir,['GCAy2x_transformed_ROI0000',num2str(i),'.nii']);
                Outfilenametemp5 = fullfile(Outputdir,['NetFx2y_ROI0000',num2str(i),'.nii']);
                POutfilenametemp1 = fullfile(Outputdir,['Perm_GCAx2y_ROI0000',num2str(i),'.nii']);
                POutfilenametemp2 = fullfile(Outputdir,['Perm_GCAy2x_ROI0000',num2str(i),'.nii']);
                POutfilenametemp3 = fullfile(Outputdir,['Perm_GCAx2y_transformed_ROI0000',num2str(i),'.nii']);
                POutfilenametemp4 = fullfile(Outputdir,['Perm_GCAy2x_transformed_ROI0000',num2str(i),'.nii']);
                POutfilenametemp5 = fullfile(Outputdir,['Perm_NetFx2y_ROI0000',num2str(i),'.nii']);
            elseif i<1000
                Outfilenametemp1 = fullfile(Outputdir,['GCAx2y_ROI000',num2str(i),'.nii']);
                Outfilenametemp2 = fullfile(Outputdir,['GCAy2x_ROI000',num2str(i),'.nii']);
                Outfilenametemp3 = fullfile(Outputdir,['GCAx2y_transformed_ROI000',num2str(i),'.nii']);
                Outfilenametemp4 = fullfile(Outputdir,['GCAy2x_transformed_ROI000',num2str(i),'.nii']);
                Outfilenametemp5 = fullfile(Outputdir,['NetFx2y_ROI000',num2str(i),'.nii']);
                POutfilenametemp1 = fullfile(Outputdir,['Perm_GCAx2y_ROI000',num2str(i),'.nii']);
                POutfilenametemp2 = fullfile(Outputdir,['Perm_GCAy2x_ROI000',num2str(i),'.nii']);
                POutfilenametemp3 = fullfile(Outputdir,['Perm_GCAx2y_transformed_ROI000',num2str(i),'.nii']);
                POutfilenametemp4 = fullfile(Outputdir,['Perm_GCAy2x_transformed_ROI000',num2str(i),'.nii']);
                POutfilenametemp5 = fullfile(Outputdir,['Perm_NetFx2y_ROI000',num2str(i),'.nii']);
            else
                Outfilenametemp1 = fullfile(Outputdir,['GCAx2y_ROI00',num2str(i),'.nii']);
                Outfilenametemp2 = fullfile(Outputdir,['GCAy2x_ROI00',num2str(i),'.nii']);
                Outfilenametemp3 = fullfile(Outputdir,['GCAx2y_transformed_ROI00',num2str(i),'.nii']);
                Outfilenametemp4 = fullfile(Outputdir,['GCAy2x_transformed_ROI00',num2str(i),'.nii']);
                Outfilenametemp5 = fullfile(Outputdir,['NetFx2y_ROI00',num2str(i),'.nii']);
                POutfilenametemp1 = fullfile(Outputdir,['Perm_GCAx2y_ROI00',num2str(i),'.nii']);
                POutfilenametemp2 = fullfile(Outputdir,['Perm_GCAy2x_ROI00',num2str(i),'.nii']);
                POutfilenametemp3 = fullfile(Outputdir,['Perm_GCAx2y_transformed_ROI00',num2str(i),'.nii']);
                POutfilenametemp4 = fullfile(Outputdir,['Perm_GCAy2x_transformed_ROI00',num2str(i),'.nii']);
                POutfilenametemp5 = fullfile(Outputdir,['Perm_NetFx2y_ROI00',num2str(i),'.nii']);
            end
            [v1 map1] = Dynamic_read_dir_NIFTI(Outfilenametemp1);
            [v2 map2] = Dynamic_read_dir_NIFTI(Outfilenametemp2);
            [v3 map3] = Dynamic_read_dir_NIFTI(Outfilenametemp3);
            [v4 map4] = Dynamic_read_dir_NIFTI(Outfilenametemp4);
            [v5 map5] = Dynamic_read_dir_NIFTI(Outfilenametemp5);
            INDMAP = find(DATMASK);
            try
                ResultMap1 = zeros(Permnum,size(maskedSignal,2));
                ResultMap2 = zeros(Permnum,size(maskedSignal,2));
                ResultMap3 = zeros(Permnum,size(maskedSignal,2));
                ResultMap4 = zeros(Permnum,size(maskedSignal,2));
                ResultMap5 = zeros(Permnum,size(maskedSignal,2));
                for iperm = 1:Permnum
                    IOrds = randperm(size(ROIsignals,1));
                    [ResultMap1(iperm,:),ResultMap2(iperm,:),ResultMap3(iperm,:),ResultMap4(iperm,:),ResultMap5(iperm,:)] = restgca_residual(ROIsignals(IOrds,i), maskedSignal(IOrds,:),Order,theCovariables(IOrds,:));
                end
                
                [mu,sig,~,sigci] = normfit(ResultMap1);
                P_map1 = normcdf(map1(INDMAP)',mu,sig);
                DataOut = zeros(size(DATMASK));
                DataOut(INDMAP) = P_map1;
                DynamicBC_write_NIFTI(DataOut,v1,POutfilenametemp1);
                
                [mu,sig,~,sigci] = normfit(ResultMap2);
                P_map2 = normcdf(map2(INDMAP)',mu,sig);
                DataOut = zeros(size(DATMASK));
                DataOut(INDMAP) = P_map2;
                DynamicBC_write_NIFTI(DataOut,v1,POutfilenametemp2);
                
                [mu,sig,~,sigci] = normfit(ResultMap3);
                P_map3 = normcdf(map3(INDMAP)',mu,sig);
                DataOut = zeros(size(DATMASK));
                DataOut(INDMAP) = P_map3;
                DynamicBC_write_NIFTI(DataOut,v1,POutfilenametemp3);
                
                [mu,sig,~,sigci] = normfit(ResultMap4);
                P_map4 = normcdf(map4(INDMAP)',mu,sig);
                DataOut = zeros(size(DATMASK));
                DataOut(INDMAP) = P_map4;
                DynamicBC_write_NIFTI(DataOut,v1,POutfilenametemp4);
                
                [mu,sig,~,sigci] = normfit(ResultMap5);
                P_map5 = normcdf(map5(INDMAP)',mu,sig);
                DataOut = zeros(size(DATMASK));
                DataOut(INDMAP) = P_map5;
                DynamicBC_write_NIFTI(DataOut,v1,POutfilenametemp5);
            catch
                totallen = size(maskedSignal,2);
                piece = 1000;
                Nlen = ceil(totallen/piece);
                DataOut1 = zeros(size(DATMASK));
                DataOut2 = zeros(size(DATMASK));
                DataOut3 = zeros(size(DATMASK));
                DataOut4 = zeros(size(DATMASK));
                DataOut5 = zeros(size(DATMASK));
                for ilen = 1:Nlen
                    if ilen~=Nlen
                        PIECEORD = 1+piece*(ilen-1):piece*ilen;
                    else
                        PIECEORD = 1+piece*(ilen-1):totallen;
                    end
                    clear ResultMap1 ResultMap2 ResultMap3 ResultMap4 ResultMap5
                    for iperm = 1:Permnum
                        IOrds = randperm(size(ROIsignals,1));
                        [ResultMap1(iperm,:),ResultMap2(iperm,:),ResultMap3(iperm,:),ResultMap4(iperm,:),ResultMap5(iperm,:)] = restgca_residual(ROIsignals(IOrds,i), maskedSignal(IOrds,PIECEORD),Order,theCovariables(IOrds,:));
                    end
                    
                    [mu,sig,~,sigci] = normfit(ResultMap1);
                    P_map1 = normcdf(map1(INDMAP(PIECEORD))',mu,sig);
                    DataOut1(INDMAP(PIECEORD)) = P_map1;
                    
                    [mu,sig,~,sigci] = normfit(ResultMap2);
                    P_map2 = normcdf(map2(INDMAP(PIECEORD))',mu,sig);
                    DataOut2(INDMAP(PIECEORD)) = P_map2;
                    
                    [mu,sig,~,sigci] = normfit(ResultMap3);
                    P_map3 = normcdf(map3(INDMAP(PIECEORD))',mu,sig);
                    DataOut3(INDMAP(PIECEORD)) = P_map3;
                    
                    [mu,sig,~,sigci] = normfit(ResultMap4);
                    P_map4 = normcdf(map4(INDMAP(PIECEORD))',mu,sig);
                    DataOut4(INDMAP(PIECEORD)) = P_map4;
                    
                    [mu,sig,~,sigci] = normfit(ResultMap5);
                    P_map5 = normcdf(map5(INDMAP(PIECEORD))',mu,sig);
                    DataOut5(INDMAP(PIECEORD)) = P_map5;
                end
                DynamicBC_write_NIFTI(DataOut1,v1,POutfilenametemp1);
                DynamicBC_write_NIFTI(DataOut2,v1,POutfilenametemp2);
                DynamicBC_write_NIFTI(DataOut3,v1,POutfilenametemp3);
                DynamicBC_write_NIFTI(DataOut4,v1,POutfilenametemp4);
                DynamicBC_write_NIFTI(DataOut5,v1,POutfilenametemp5);
            end
        end
    else % coef
        INDMAP = find(DATMASK);
        for i = 1:size(ROIsignals,2)
            if i<10
                Outfilenametemp1 = ['Coef_GCAx2y_ROI00000',num2str(i),'.nii'];
                Outfilenametemp2 = ['Coef_GCAy2x_ROI00000',num2str(i),'.nii'];
                Outfilenametemp3 = ['Coef_GCAx2y_AR_ROI00000',num2str(i),'.nii'];
                Outfilenametemp4 = ['Coef_GCAy2x_AR_ROI00000',num2str(i),'.nii'];
                POutfilenametemp1 = ['Perm_Coef_GCAx2y_ROI00000',num2str(i),'.nii'];
                POutfilenametemp2 = ['Perm_Coef_GCAy2x_ROI00000',num2str(i),'.nii'];
                POutfilenametemp3 = ['Perm_Coef_GCAx2y_AR_ROI00000',num2str(i),'.nii'];
                POutfilenametemp4 = ['Perm_Coef_GCAy2x_AR_ROI00000',num2str(i),'.nii'];
            elseif i<100
                Outfilenametemp1 = ['Coef_GCAx2y_ROI0000',num2str(i),'.nii'];
                Outfilenametemp2 = ['Coef_GCAy2x_ROI0000',num2str(i),'.nii'];
                Outfilenametemp3 = ['Coef_GCAx2y_AR_ROI0000',num2str(i),'.nii'];
                Outfilenametemp4 = ['Coef_GCAy2x_AR_ROI0000',num2str(i),'.nii'];
                POutfilenametemp1 = ['Perm_Coef_GCAx2y_ROI0000',num2str(i),'.nii'];
                POutfilenametemp2 = ['Perm_Coef_GCAy2x_ROI0000',num2str(i),'.nii'];
                POutfilenametemp3 = ['Perm_Coef_GCAx2y_AR_ROI0000',num2str(i),'.nii'];
                POutfilenametemp4 = ['Perm_Coef_GCAy2x_AR_ROI0000',num2str(i),'.nii'];
            elseif i<1000
                Outfilenametemp1 = ['Coef_GCAx2y_ROI000',num2str(i),'.nii'];
                Outfilenametemp2 = ['Coef_GCAy2x_ROI000',num2str(i),'.nii'];
                Outfilenametemp3 = ['Coef_GCAx2y_AR_ROI000',num2str(i),'.nii'];
                Outfilenametemp4 = ['Coef_GCAy2x_AR_ROI000',num2str(i),'.nii'];
                POutfilenametemp1 = ['Perm_Coef_GCAx2y_ROI000',num2str(i),'.nii'];
                POutfilenametemp2 = ['Perm_Coef_GCAy2x_ROI000',num2str(i),'.nii'];
                POutfilenametemp3 = ['Perm_Coef_GCAx2y_AR_ROI000',num2str(i),'.nii'];
                POutfilenametemp4 = ['Perm_Coef_GCAy2x_AR_ROI000',num2str(i),'.nii'];
            else
                Outfilenametemp1 = ['Coef_GCAx2y_ROI00',num2str(i),'.nii'];
                Outfilenametemp2 = ['Coef_GCAy2x_ROI00',num2str(i),'.nii'];
                Outfilenametemp3 = ['Coef_GCAx2y_AR_ROI00',num2str(i),'.nii'];
                Outfilenametemp4 = ['Coef_GCAy2x_AR_ROI00',num2str(i),'.nii'];
                POutfilenametemp1 = ['Perm_Coef_GCAx2y_ROI00',num2str(i),'.nii'];
                POutfilenametemp2 = ['Perm_Coef_GCAy2x_ROI00',num2str(i),'.nii'];
                POutfilenametemp3 = ['Perm_Coef_GCAx2y_AR_ROI00',num2str(i),'.nii'];
                POutfilenametemp4 = ['Perm_Coef_GCAy2x_AR_ROI00',num2str(i),'.nii'];
            end
            for iord = 1:Order
                Outfiletemp = fullfile(Outdir,['Order',num2str(iord),'_',Outfilenametemp1]);
                [v1 map1] = Dynamic_read_dir_NIFTI(Outfiletemp);
                MAP1{iord} = map1;
                Outfiletemp = fullfile(Outdir,['Order',num2str(iord),'_',Outfilenametemp2]);
                [v2 map2] = Dynamic_read_dir_NIFTI(Outfiletemp);
                MAP2{iord} = map2;
                Outfiletemp = fullfile(Outdir,['Order',num2str(iord),'_',Outfilenametemp3]);
                [v3 map3] = Dynamic_read_dir_NIFTI(Outfiletemp);
                MAP3{iord} = map3;
                Outfiletemp = fullfile(Outdir,['Order',num2str(iord),'_',Outfilenametemp4]);
                [v4 map4] = Dynamic_read_dir_NIFTI(Outfiletemp);
                MAP4{iord} = map4;
            end
            totallen = size(maskedSignal,2);
            piece = 1000;
            Nlen = ceil(totallen/piece);
            for ilen = 1:Nlen
                if ilen~=Nlen
                    PIECEORD = 1+piece*(ilen-1):piece*ilen;
                else
                    PIECEORD = 1+piece*(ilen-1):totallen;
                end
                for iord = 1:Order
                    eval(['clear ResultMap1_',num2str(iord),' ResultMap2_',num2str(iord),' ResultMap3_',num2str(iord),' ResultMap4_',num2str(iord),';']);
                end
                for iperm = 1:Permnum
                    IOrds = randperm(size(ROIsignals,1));
                    clear ResultMap1 ResultMap2 ResultMap3 ResultMap4
%                     eval(['clear ResultMap1_',num2str(iord),' ResultMap2_',num2str(iord),' ResultMap3_',num2str(iord),' ResultMap4_',num2str(iord),';']);
                    [ResultMap1,ResultMap2,ResultMap3,ResultMap4]=restgca_coefficient(ROIsignals(IOrds,i), maskedSignal(IOrds,PIECEORD),Order,theCovariables(IOrds,:));
                    
                    for iord = 1:Order
                        eval(['ResultMap1_',num2str(iord),'(iperm,:) = ResultMap1{',num2str(iord),'};']);
                        eval(['ResultMap2_',num2str(iord),'(iperm,:) = ResultMap2{',num2str(iord),'};']);
                        eval(['ResultMap3_',num2str(iord),'(iperm,:) = ResultMap3{',num2str(iord),'};']);
                        eval(['ResultMap4_',num2str(iord),'(iperm,:) = ResultMap4{',num2str(iord),'};']);
                    end
                end
                for iord = 1:Order
                    eval(['[mu,sig,~,sigci] = normfit(ResultMap1_',num2str(iord),');']);
                    map1_t = MAP1{iord};
                    P_map1(PIECEORD) = normcdf(map1_t(INDMAP(PIECEORD))',mu,sig);
                    PMAP1{iord} = P_map1;
                    
                    eval(['[mu,sig,~,sigci] = normfit(ResultMap2_',num2str(iord),');']);
                    map2_t = MAP2{iord};
                    P_map2(PIECEORD) = normcdf(map2_t(INDMAP(PIECEORD))',mu,sig);
                    PMAP2{iord} = P_map2;
                    
                    eval(['[mu,sig,~,sigci] = normfit(ResultMap3_',num2str(iord),');']);
                    map3_t = MAP3{iord};
                    P_map3(PIECEORD) = normcdf(map3_t(INDMAP(PIECEORD))',mu,sig);
                    PMAP3{iord} = P_map3;
                    
                    eval(['[mu,sig,~,sigci] = normfit(ResultMap4_',num2str(iord),');']);
                    map4_t = MAP4{iord};
                    P_map4(PIECEORD) = normcdf(map4_t(INDMAP(PIECEORD))',mu,sig);
                    PMAP4{iord} = P_map4;
                end
            end
            for iord = 1:Order
                Outfiletemp = fullfile(Outdir,['Order',num2str(iord),'_',POutfilenametemp1]);
                DataOut = zeros(size(DATMASK));
                DataOut(find(DATMASK)) = PMAP1{iord};
                DynamicBC_write_NIFTI(DataOut,v1,Outfiletemp);
                
                Outfiletemp = fullfile(Outdir,['Order',num2str(iord),'_',POutfilenametemp2]);
                DataOut = zeros(size(DATMASK));
                DataOut(find(DATMASK)) = PMAP2{iord};
                DynamicBC_write_NIFTI(DataOut,v1,Outfiletemp);
                
                Outfiletemp = fullfile(Outdir,['Order',num2str(iord),'_',POutfilenametemp3]);
                DataOut = zeros(size(DATMASK));
                DataOut(find(DATMASK)) = PMAP3{iord};
                DynamicBC_write_NIFTI(DataOut,v1,Outfiletemp);
                
                Outfiletemp = fullfile(Outdir,['Order',num2str(iord),'_',POutfilenametemp4]);
                DataOut = zeros(size(DATMASK));
                DataOut(find(DATMASK)) = PMAP4{iord};
                DynamicBC_write_NIFTI(DataOut,v1,Outfiletemp);                
            end
        end
        
    end
else % matrix
    ROIdir = fullfile(Outdir,'ROIsignal.mat');
    load(ROIdir);
    Calmethod = Parameter.Calmethod;
    COVcond = Parameter.covs;
    COVdir = Parameter.COVtext;
    if COVcond==1
        COV_S = load(COVdir);
    end
    if COVcond==1
        theCovariables = COV_S;
    else
        theCovariables = ones(size(ROIsignals,1),1);
    end
    
    if Calmethod % Residual
        for iperm = 1:Permnum
            IOrds = randperm(size(ROIsignals,1));
            for i = 1:size(ROIsignals,2)-1
                for j = i+1:size(ROIsignals,2)
                    ROIsignalst = [ROIsignals(:,i),ROIsignals(:,j)];
                    [ResultMap1,ResultMap2,ResultMap3,ResultMap4,ResultMap5] = restgca_FROI(ROIsignalst(IOrds,:),Order,theCovariables(IOrds,:));
                    X2Y(i,j,iperm) =  ResultMap1;%X2Y(j,i) =  ResultMap1;
                    Y2X(i,j,iperm) =  ResultMap2;%Y2X(j,i) =  ResultMap2;
                    X2Y_transform(i,j,iperm) =  ResultMap3;%X2Y_transform(j,i) =  ResultMap3;
                    Y2X_transform(i,j,iperm) =  ResultMap4;%Y2X_transform(j,i) =  ResultMap4;
                    NetFx2y(i,j,iperm) = ResultMap5; NetFx2y(j,i,iperm) = ResultMap5;
                end
            end
        end
        X2Yori = load(fullfile(Outdir,'GCAx2y.mat'));
        Y2Xori = load(fullfile(Outdir,'GCAy2x.mat'));
        X2Ytransori = load(fullfile(Outdir,'GCAx2y_transformed.mat'));
        Y2Xtransori = load(fullfile(Outdir,'GCAy2x_transformed.mat'));
        NetFx2yori = load(fullfile(Outdir,'NetFx2y.mat'));
        for i = 1:size(ROIsignals,2)-1
            for j = i+1:size(ROIsignals,2)
                [mu,sig,~,sigci] = normfit(squeeze(X2Y(i,j,:)));
                P_X2Y(i,j) = normcdf(X2Yori.X2Y(i,j),mu,sig);
                [mu,sig,~,sigci] = normfit(squeeze(Y2X(i,j,:)));
                P_Y2X(i,j) = normcdf(Y2Xori.Y2X(i,j),mu,sig);
                [mu,sig,~,sigci] = normfit(squeeze(X2Y_transform(i,j,:)));
                P_X2Y_transformed(i,j) = normcdf(X2Ytransori.X2Y_transform(i,j),mu,sig);
                [mu,sig,~,sigci] = normfit(squeeze(Y2X_transform(i,j,:)));
                P_Y2X_transformed(i,j) = normcdf(Y2Xtransori.Y2X_transform(i,j),mu,sig);
                [mu,sig,~,sigci] = normfit(squeeze(NetFx2y(i,j,:)));
                P_NetFx2y(i,j) = normcdf(NetFx2yori.NetFx2y(i,j),mu,sig);
            end
        end
        X2Yoridir = fullfile(Outdir,'PermPval_GCAx2y.mat');
        Y2Xoridir = fullfile(Outdir,'PermPval_GCAy2x.mat');
        X2Ytransoridir = fullfile(Outdir,'PermPval_GCAx2y_transformed.mat');
        Y2Xtransoridir = fullfile(Outdir,'PermPval_GCAy2x_transformed.mat');
        NetFx2yoridir = fullfile(Outdir,'PermPval_NetFx2y.mat');
        save(X2Yoridir,'P_X2Y');
        save(Y2Xoridir,'P_Y2X');
        save(X2Ytransoridir,'P_X2Y_transformed');
        save(Y2Xtransoridir,'P_Y2X_transformed');
        save(NetFx2yoridir,'P_NetFx2y');
    else % coef
        for iperm = 1:Permnum
            IOrds = randperm(size(ROIsignals,1));
            for i = 1:size(ROIsignals,2)-1
                for j = i+1:size(ROIsignals,2)
                    ROIsignalst = [ROIsignals(:,i),ROIsignals(:,j)];
                    [Result_X2Y,Result_Y2X,ROI_sequence] = restgca_CROI_Bivariate(ROIsignalst(IOrds,:),Order,theCovariables(IOrds,:));
                    X2Y(i,j,:) =  Result_X2Y(1:Order);X2X(i,j,:) =  Result_X2Y(Order+1:Order*2);
                    Y2X(i,j,:) =  Result_Y2X(1:Order);Y2Y(i,j,:) =  Result_Y2X(Order+1:Order*2);
                end
            end
            X2Y_t{iperm} = X2Y;
            X2X_t{iperm} = X2X;
            Y2X_t{iperm} = Y2X;
            Y2Y_t{iperm} = Y2Y;
        end
        Ori_x2x = load(fullfile(Outdir,'Coef_GCAx2x.mat'));
        Ori_x2y = load(fullfile(Outdir,'Coef_GCAx2y.mat'));
        Ori_y2x = load(fullfile(Outdir,'Coef_GCAy2x.mat'));
        Ori_y2y = load(fullfile(Outdir,'Coef_GCAy2y.mat'));
        for k = 1:Order
            for i = 1:Permnum
                temps = X2Y_t{i};
                X2Y_T(:,:,i) = temps(:,:,k);
                temps = X2X_t{i};
                X2X_T(:,:,i) = temps(:,:,k);
                temps = Y2X_t{i};
                Y2X_T(:,:,i) = temps(:,:,k);
                temps = Y2Y_t{i};
                Y2Y_T(:,:,i) = temps(:,:,k);
            end
            for i = 1:size(ROIsignals,2)-1
                for j = i+1:size(ROIsignals,2)
                    [mu,sig,~,sigci] = normfit(squeeze(X2Y_T(i,j,:)));
                    P_X2Y(i,j,k) = normcdf(Ori_x2y.X2Y(i,j),mu,sig);
                    [mu,sig,~,sigci] = normfit(squeeze(X2X_T(i,j,:)));
                    P_X2X(i,j,k) = normcdf(Ori_x2x.X2X(i,j),mu,sig);
                    [mu,sig,~,sigci] = normfit(squeeze(Y2X_T(i,j,:)));
                    P_Y2X(i,j,k) = normcdf(Ori_y2x.Y2X(i,j),mu,sig);
                    [mu,sig,~,sigci] = normfit(squeeze(Y2Y_T(i,j,:)));
                    P_Y2Y(i,j,k) = normcdf(Ori_y2y.Y2Y(i,j),mu,sig);
                end
            end
        end        
        dirOri_x2x = fullfile(Outdir,'PermPval_Coef_GCAx2x.mat');
        dirOri_x2y = fullfile(Outdir,'PermPval_Coef_GCAx2y.mat');
        dirOri_y2x = fullfile(Outdir,'PermPval_Coef_GCAy2x.mat');
        dirOri_y2y = fullfile(Outdir,'PermPval_Coef_GCAy2y.mat');
        save(dirOri_x2x,'P_X2X');
        save(dirOri_x2y,'P_X2Y');
        save(dirOri_y2x,'P_Y2X');
        save(dirOri_y2y,'P_Y2Y');
    end
end