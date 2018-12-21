function AS_WTA_computemain_master(Parameter)
Outputdir = Parameter.Outputdir; % ���·��
RealCompPara.Outputdir = Outputdir;
Inputdir = Parameter.Inputdir; % ����·��
RealCompPara.Inputdir = Inputdir;
Methods = Parameter.methodused; % ����1 Ƥ��ѷ��أ�����2 ƫ���
RealCompPara.Methods = Methods;
COVcond = Parameter.covs; % �Ƿ�ʹ��Э����
RealCompPara.COVcond = COVcond;
COVdir = Parameter.COVtext; % Э��������·��
RealCompPara.COVdir = COVdir;
Seeddir = Parameter.SeedROI; % Ƥ��ROI
RealCompPara.Seeddir = Seeddir;
Targetdir = Parameter.TargetROI; % Ŀ��ROI
RealCompPara.Targetdir = Targetdir;
[vseed,dataseed] = Dynamic_read_dir_NIFTI(Seeddir); % ��ȡ����
[vtarget,datatarget] = Dynamic_read_dir_NIFTI(Targetdir);
[vinput,datainput] = Dynamic_read_dir_NIFTI(Inputdir);
%% deal with the wrong slice condition.
% if any(vseed.dim-vinput.dim)|any(vseed.mat(:)-vseed.mat(:))
if any(vseed.dim-vinput.dim) % ��Ƥ��ROI���������ݲ�һ�µ�ʱ��
    Maskimg = vseed.fname;
    [pat nam ext] = fileparts(Maskimg);
    Maskimg0 = fullfile(Outputdir,['seed_',nam,'_reslice',ext]);
    dynamicBC_Reslice(Maskimg,Maskimg0,vinput.dim,0,vinput.fname);
    Seeddir = Maskimg0;
end
clear vseed dataseed
[vseed,dataseed] = Dynamic_read_dir_NIFTI(Seeddir);
RealCompPara.Seeddir = Seeddir;
% if any(vtarget.dim-vinput.dim)|any(vtarget.mat(:)-vinput.mat(:))
if any(vtarget.dim-vinput.dim) % ��Ƥ����ROI���������ݲ�һ�µ�ʱ�� % �˴����ǵ����ݵĲ������ԣ���׼���ܴ�������
    Maskimg = vtarget.fname;
    [pat nam ext] = fileparts(Maskimg);
    Maskimg0 = fullfile(Outputdir,['Target_',nam,'_reslice',ext]);
    dynamicBC_Reslice(Maskimg,Maskimg0,vinput.dim,0,vinput.fname);
    Targetdir = Maskimg0;
end
clear vtarget datatarget
[vtarget,datatarget] = Dynamic_read_dir_NIFTI(Targetdir);
RealCompPara.Targetdir = Targetdir;
%
CompMat = fullfile(Outputdir,'RealCompute.mat');
save(CompMat,'RealCompPara');
%%
% indexseed = find(dataseed);
datatarget(isnan(datatarget)) = 0;  % ����Ƥ��ROI����NaN�����
dataseed(isnan(dataseed)) = 0;  % ����Ƥ��ROI����NaN�����
indexstarget = find(datatarget); % Ѱ��Ƥ��ROI��index
seedorders = unique(dataseed); % Ѱ��Ƥ��ROI������ֵ
seednum = length(seedorders)-1; % ����Ƥ��ROI��Ŀ
for i = 1:seednum % ��ȡƤ��ROI�ź�
    datSeed(i,:) = mean(datainput(dataseed==seedorders(i+1),:),1);
end
datTarget = datainput(indexstarget,:); % ��ȡƤ��ROI���������ź�
if COVcond==1 % ��ȡCOV����ʽTXT
    COV = load(COVdir); 
else
    COV = []; % ����COV�ǣ�����Ϊ��
end
if Methods==1 %pearson
    if COVcond==1
        [r p] = partialcorr(datSeed',datTarget',COV); % ��Э������ƫ��ط���
    else
        [r p] = corr(datSeed',datTarget'); % ��ط���
    end
else % patial correlation
    for i = 1:seednum  % ��Ҫÿ��ROI��μ���
        COVind = 1:seednum; % �����ȥ��i��ROI��ʣ���ź�
        COVind(i) = [];
        COVadd = datSeed(COVind,:)'; % ��ȡ��ȥ��i��ROI���ʣ���ź���ΪЭ����
        if COVcond==1 % ����������Э������ʱ��
            COVall = [COVadd,COV]; % ƴ��Э����
            [r(i,:) p(i,:)] = partialcorr(datSeed(i,:)',datTarget',COVall); % ƫ��ط���
        else
            [r(i,:) p(i,:)] = partialcorr(datSeed(i,:)',datTarget',COVadd); % ƫ��ط���
        end
    end
end
Signalsavename = fullfile(Outputdir,'Signals.mat'); 
save(Signalsavename,'datSeed','datTarget','COV'); % ��Ƥ��ROI��Ŀ����ROI��ԭʼ�źű���ΪSignals.mat
Rmatsavename = fullfile(Outputdir,'computeval.mat'); %
IND0 = sum(datTarget,2);
IND1 = find(IND0~=0);
save(Rmatsavename,'r','p','IND1','vtarget','indexstarget'); % ����¼��Ϣ����Ϊcomputeval.mat������IND1Ϊ�����
%%
for i = 1:seednum %��r��pд��map
    rval = r(i,:);
    rval(isnan(rval)) = 0;
    pval = p(i,:);
    pval(isnan(rval)) = 0;
    if i<10
        outnamer = fullfile(Outputdir,['R_ROI000',num2str(i),'.nii']);
        outnamep = fullfile(Outputdir,['P_ROI000',num2str(i),'.nii']);
    elseif i<100
        outnamer = fullfile(Outputdir,['R_ROI00',num2str(i),'.nii']);
        outnamep = fullfile(Outputdir,['P_ROI00',num2str(i),'.nii']);
    elseif i<1000
        outnamer = fullfile(Outputdir,['R_ROI0',num2str(i),'.nii']);
        outnamep = fullfile(Outputdir,['P_ROI0',num2str(i),'.nii']);
    else
        outnamer = fullfile(Outputdir,['R_ROI',num2str(i),'.nii']);
        outnamep = fullfile(Outputdir,['P_ROI',num2str(i),'.nii']);
    end
    Rmap = zeros(vtarget.dim);
    Rmap(indexstarget) = rval;
    Pmap = zeros(vtarget.dim);
    Pmap(indexstarget) = pval;
    DynamicBC_write_NIFTI(Rmap,vtarget,outnamer);
    DynamicBC_write_NIFTI(Pmap,vtarget,outnamep);
end
% �������ֵ����Ӧ��Ƥ��ROI
[maxval maxind] = max(r);
MAXmap = zeros(vtarget.dim);
MAXidmap = zeros(vtarget.dim);
MAXmap(indexstarget) = maxval;
MAXidmap(indexstarget) = maxind;
DynamicBC_write_NIFTI(MAXmap,vtarget,fullfile(Outputdir,'mixedMaxVal.nii'));
DynamicBC_write_NIFTI(MAXidmap,vtarget,fullfile(Outputdir,'mixedMaxID.nii'));
%%
[maxvalA maxindA] = max(abs(r));
MAXmapA = zeros(vtarget.dim);
MAXidmapA = zeros(vtarget.dim);
for i = 1:length(maxindA)
    MAXmapA(indexstarget(i)) = r(maxindA(i),i);
end
MAXidmapA(indexstarget) = maxindA;
DynamicBC_write_NIFTI(MAXmapA,vtarget,fullfile(Outputdir,'ABSmixedMaxVal.nii'));
DynamicBC_write_NIFTI(MAXidmapA,vtarget,fullfile(Outputdir,'ABSmixedMaxID.nii'));
%%
for i = 1:seednum
    SEDID = zeros(vtarget.dim);
    SEDID(MAXidmap==i) = 1;
    SEDIDA = zeros(vtarget.dim);
    SEDIDA(MAXidmapA==i) = 1;
    if i<10
        outnameS = fullfile(Outputdir,['SepID_ROI000',num2str(i),'.nii']);
        outnameSA = fullfile(Outputdir,['ABSSepID_ROI000',num2str(i),'.nii']);
    elseif i<100
        outnameS = fullfile(Outputdir,['SepID_ROI00',num2str(i),'.nii']);
        outnameSA = fullfile(Outputdir,['ABSSepID_ROI00',num2str(i),'.nii']);
    elseif i<1000
        outnameS = fullfile(Outputdir,['SepID_ROI0',num2str(i),'.nii']);
        outnameSA = fullfile(Outputdir,['ABSSepID_ROI0',num2str(i),'.nii']);
    else
        outnameS = fullfile(Outputdir,['SepID_ROI',num2str(i),'.nii']);
        outnameSA = fullfile(Outputdir,['ABSSepID_ROI',num2str(i),'.nii']);
    end
    DynamicBC_write_NIFTI(SEDID,vtarget,outnameS);
    DynamicBC_write_NIFTI(SEDIDA,vtarget,outnameSA);
end
save(fullfile(Outputdir,'LabedVal.mat'),'maxval','maxind','maxvalA','maxindA','seednum');
%%
%% SCATTER
if strcmp(Parameter.cal_types,'Yes')
%     Methods = 1;
    if Methods==1 %pearson 
        if COVcond==0 
            for i = 1:size(datSeed,1)
                [R P labs] = scatcorr(datTarget',datSeed(i,:)'); % scat corr
                r(i,:) = R;
                p(i,:) = P;
                LABS{i} = labs;
            end
        else % withcov
            DATCOMB = [datSeed',datTarget'];
            COVz = [ones(size(COV,1),1),COV];
            resid = DATCOMB - COVz*(COVz \ DATCOMB); % ȥ��Э����
            datSeedN = resid(:,1:size(datSeed,1)); % ǰN��ROIΪSeed
            datTargetN = resid(:,size(datSeed,1)+1:size(resid,2)); % N+1���ΪTarget
            sig_pearson_covsavename = fullfile(Outputdir,'Scatter_PearsonCOVsignal.mat'); 
            save(sig_pearson_covsavename,'datSeedN','datTargetN'); % ����ȥ��Э�������seed��target�ź�
            
            for i = 1:size(datSeedN,2)
                [R P labs] = scatcorr(datTargetN,datSeedN(:,i));
                r(i,:) = R;
                p(i,:) = P;
                LABS{i} = labs;
            end
        end
    else % partial correlation
        for i = 1:seednum
            COVind = 1:seednum;
            COVind(i) = [];
            COVadd = datSeed(COVind,:)';
            DATCOMB = [datSeed(i,:)',datTarget'];
            if COVcond==0
                
                COVz = [ones(size(COVadd,1),1),COVadd];
                resid = DATCOMB - COVz*(COVz \ DATCOMB);
                datSeedN = resid(:,1);
                datTargetN = resid(:,2:size(resid,2));
                sig_pearson_covsavename = fullfile(Outputdir,['Scatter_ROI',num2str(i),'_Partial_NoCOVsignal.mat']);
                save(sig_pearson_covsavename,'datSeedN','datTargetN');
                [R P labs] = scatcorr(datTargetN,datSeedN);
                r(i,:) = R;
                p(i,:) = P;
                LABS{i} = labs;
            else
                COVz = [ones(size(COVadd,1),1),COVadd,COV];
                resid = DATCOMB - COVz*(COVz \ DATCOMB);
                datSeedN = resid(:,1);
                datTargetN = resid(:,2:size(resid,2));
                sig_pearson_covsavename = fullfile(Outputdir,['Scatter_ROI',num2str(i),'_Partial_WithCOVsignal.mat']);
                save(sig_pearson_covsavename,'datSeedN','datTargetN');
                [R P labs] = scatcorr(datTargetN,datSeedN);
                r(i,:) = R;
                p(i,:) = P;
                LABS{i} = labs;
            end
        end
    end
    Rmatsavename = fullfile(Outputdir,'Scatter_computeval.mat');
    IND0 = sum(datTarget,2);
    IND1 = find(IND0~=0);
    save(Rmatsavename,'r','p','IND1','vtarget','LABS','indexstarget');
    %%
    for i = 1:seednum
        rval = r(i,:);
        rval(isnan(rval)) = 0;
        pval = p(i,:);
        pval(isnan(rval)) = 0;
        if i<10
            outnamer = fullfile(Outputdir,['Scatter_R_ROI000',num2str(i),'.nii']);
            outnamep = fullfile(Outputdir,['Scatter_P_ROI000',num2str(i),'.nii']);
        elseif i<100
            outnamer = fullfile(Outputdir,['Scatter_R_ROI00',num2str(i),'.nii']);
            outnamep = fullfile(Outputdir,['Scatter_P_ROI00',num2str(i),'.nii']);
        elseif i<1000
            outnamer = fullfile(Outputdir,['Scatter_R_ROI0',num2str(i),'.nii']);
            outnamep = fullfile(Outputdir,['Scatter_P_ROI0',num2str(i),'.nii']);
        else
            outnamer = fullfile(Outputdir,['Scatter_R_ROI',num2str(i),'.nii']);
            outnamep = fullfile(Outputdir,['Scatter_P_ROI',num2str(i),'.nii']);
        end
        Rmap = zeros(vtarget.dim);
        Rmap(indexstarget) = rval;
        Pmap = zeros(vtarget.dim);
        Pmap(indexstarget) = pval;
        DynamicBC_write_NIFTI(Rmap,vtarget,outnamer);
        DynamicBC_write_NIFTI(Pmap,vtarget,outnamep);
    end
    [maxval maxind] = max(r);
    MAXmap = zeros(vtarget.dim);
    MAXidmap = zeros(vtarget.dim);
    MAXmap(indexstarget) = maxval;
    MAXidmap(indexstarget) = maxind;
    DynamicBC_write_NIFTI(MAXmap,vtarget,fullfile(Outputdir,'Scatter_mixedMaxVal.nii'));
    DynamicBC_write_NIFTI(MAXidmap,vtarget,fullfile(Outputdir,'Scatter_mixedMaxID.nii'));
    
    %%
    [maxvalA maxindA] = max(abs(r));
    MAXmapA = zeros(vtarget.dim);
    MAXidmapA = zeros(vtarget.dim);
    for i = 1:length(maxindA)
        MAXmapA(indexstarget(i)) = r(maxindA(i),i);
    end
    MAXidmapA(indexstarget) = maxindA;
    DynamicBC_write_NIFTI(MAXmapA,vtarget,fullfile(Outputdir,'Scatter_ABSmixedMaxVal.nii'));
    DynamicBC_write_NIFTI(MAXidmapA,vtarget,fullfile(Outputdir,'Scatter_ABSmixedMaxID.nii'));
    %%
    for i = 1:seednum
        SEDID = zeros(vtarget.dim);
        SEDID(MAXidmap==i) = 1;
        SEDIDA = zeros(vtarget.dim);
        SEDIDA(MAXidmapA==i) = 1;
        if i<10
            outnameS = fullfile(Outputdir,['Scatter_SepID_ROI000',num2str(i),'.nii']);
            outnameSA = fullfile(Outputdir,['Scatter_ABSSepID_ROI000',num2str(i),'.nii']);
        elseif i<100
            outnameS = fullfile(Outputdir,['Scatter_SepID_ROI00',num2str(i),'.nii']);
            outnameSA = fullfile(Outputdir,['Scatter_ABSSepID_ROI00',num2str(i),'.nii']);
        elseif i<1000
            outnameS = fullfile(Outputdir,['Scatter_SepID_ROI0',num2str(i),'.nii']);
            outnameSA = fullfile(Outputdir,['Scatter_ABSSepID_ROI0',num2str(i),'.nii']);
        else
            outnameS = fullfile(Outputdir,['Scatter_SepID_ROI',num2str(i),'.nii']);
            outnameSA = fullfile(Outputdir,['Scatter_ABSSepID_ROI',num2str(i),'.nii']);
        end
        DynamicBC_write_NIFTI(SEDID,vtarget,outnameS);
        DynamicBC_write_NIFTI(SEDIDA,vtarget,outnameSA);
    end
    save(fullfile(Outputdir,'Scatter_LabedVal.mat'),'maxval','maxind','maxvalA','maxindA','seednum');
end
%%

uiwait(msgbox('Compute Finished!'));
end