function AS_WTA_Stat_main_master(Parameter)
Outdir = Parameter.Outdir;
Indir1 = Parameter.Input1;
Indir2 = Parameter.Input2;
Incalmat1 = load(fullfile(Indir1,'SetUpparameter.mat')); 
Incalmat2 = load(fullfile(Indir2,'SetUpparameter.mat'));
COVcond1 = Incalmat1.Parameter.covs; % g1 Э�������
COVcond2 = Incalmat2.Parameter.covs; % g2 Э�������
Calmethod1 = Incalmat1.Parameter.methodused; % G1 ���㷽�� 1 pearson 2 partialcorr
Calmethod2 = Incalmat2.Parameter.methodused; % G2 ���㷽�� 1 pearson 2 partialcorr 
if COVcond1~=COVcond2
    error('The two groups used different cov conditions');
end
if Calmethod1~=Calmethod2
    error('The tow groups used different calculate methods');
end
COVCOND = COVcond1;
PermMethod = Parameter.permtype; % ͳ�Ʒ����� 1 permutation 2 interaction 

computeval1 = load(fullfile(Indir1,'computeval.mat'));
computeval2 = load(fullfile(Indir2,'computeval.mat'));
OrigR1 = computeval1.r; % G1 no scat��rֵ
OrigR2 = computeval2.r; % G2 no scat��rֵ
DELTR = OrigR1-OrigR2; % ����ֵ����Ҫ����permutation
SIG1 = load(fullfile(Indir1,'Signals.mat')); % ��ȡԭʼ�ź�
SIG2 = load(fullfile(Indir2,'Signals.mat')); % ��ȡԭʼ�ź�
datTarget1 = SIG1.datTarget; % G1 target�źţ�ԭʼ��
datTarget2 = SIG2.datTarget; % G2 target�źţ�ԭʼ��
datSeed1 = SIG1.datSeed; % G1 seed�źţ�ԭʼ��
datSeed2 = SIG2.datSeed; % G2 seed�źţ�ԭʼ��
NROI = size(datSeed1,1);
if PermMethod==1 % permutation 
    Permnum = Parameter.permnum; % permutation��Ŀ
    N1 = size(datTarget1,2); % G1 ��������ʱ��㳤�ȣ�
    N2 = size(datTarget2,2); % G2 ��������ʱ��㳤�ȣ�
    Ntotal = N1+N2; % ���屻����Ŀ��ʱ��㳤�ȣ�
    dattarget = [datTarget1,datTarget2]'; % ƴ�ӵ�target����
    if COVCOND==0 % no cov 
        if Calmethod1==1 % pearson & nocov
            for i = 1:NROI
                sigseedt = [datSeed1(i,:),datSeed2(i,:)]'; % ��i��ROI��seed�ź�ƴ��
                deltr = zeros(Permnum,size(datTarget1,1)); % Ԥ��ֵdelta-pֵ����
                parfor iperm = 1:Permnum % ����û�
                    Nrandp = randperm(Ntotal); % ����
                    seedsig1 = sigseedt(Nrandp(1:N1));  % ǰN1��Seed
                    seedsig2 = sigseedt(Nrandp(N1+1:Ntotal)); % ��N2��Seed
                    targetsig1 = dattarget(Nrandp(1:N1),:);  % ǰN1��target
                    targetsig2 = dattarget(Nrandp(N1+1:Ntotal),:); % ��N2��target
                    [r1 p1] = corr(seedsig1,targetsig1); %��ط���
                    [r2 p2] = corr(seedsig2,targetsig2); %��ط���
                    deltr(iperm,:) = r1-r2; % delta-pֵ¼��
                end
                [mu,sig,~,sigci] = normfit(deltr); %��̬�ֲ�
                P = normcdf(DELTR(i,:),mu,sig); % ����pֵ��<0.05ΪG1<G2��>0.95ΪG1>G2)
                Pval(i,:) = P; % ¼��pֵ
            end
            INDEXS = computeval1.indexstarget;
            vmat = computeval1.vtarget;
            for i = 1:NROI
                if i<10
                    OutPname = fullfile(Outdir,['Pval_ROI0000',num2str(i),'.nii']);
                elseif i<100
                    OutPname = fullfile(Outdir,['Pval_ROI000',num2str(i),'.nii']);
                elseif i<1000
                    OutPname = fullfile(Outdir,['Pval_ROI00',num2str(i),'.nii']);
                else
                    error('too many ROIs');
                end
                Ptemp = zeros(vmat.dim);
                Ptemp(INDEXS) = Pval(i,:);
                DynamicBC_write_NIFTI(Ptemp,vmat,OutPname);
            end
        else % particalcorr & nocov
            for i = 1:NROI % ���ROI
                sigseedt = [datSeed1(i,:),datSeed2(i,:)]';  % ƴ�����ӵ�
                deltr = zeros(Permnum,size(datTarget1,1));  % ����r�������
                Rtemp1 = zeros(Permnum,size(datTarget1,1));  
                Rtemp2 = zeros(Permnum,size(datTarget1,1));  
                INDUSED = 1:NROI;INDUSED(i) = []; % ����Э����
                covt = [datSeed1(INDUSED,:),datSeed2(INDUSED,:)]'; % Э����ƴ��
                parfor iperm = 1:Permnum % ����û�
                    Nrandp = randperm(Ntotal); % �������
                    seedsig1 = sigseedt(Nrandp(1:N1)); % ǰN1��Seed
                    seedsig2 = sigseedt(Nrandp(N1+1:Ntotal)); % ��N2��seed
                    targetsig1 = dattarget(Nrandp(1:N1),:); % ǰN1��target
                    targetsig2 = dattarget(Nrandp(N1+1:Ntotal),:); % ��N2��target
                    covt1 = covt(Nrandp(1:N1),:); % ǰN1��COV
                    covt2 = covt(Nrandp(N1+1:Ntotal),:); % ��N2��COV
                    [r1 p1] = partialcorr(seedsig1,targetsig1,covt1); % ƫ���
                    [r2 p2] = partialcorr(seedsig2,targetsig2,covt2); % ƫ���
                    deltr(iperm,:) = r1-r2; % rֵ����¼��
                    Rtemp1(iperm,:) = r1; 
                    Rtemp2(iperm,:) = r2;
                end
                save([Outdir,filesep,'ExtRval_',num2str(i),'.mat'],'Rtemp1','Rtemp2');
                [mu,sig,~,sigci] = normfit(deltr); % ��̬�ֲ�
                P = normcdf(DELTR(i,:),mu,sig); % pֵ����
                Pval(i,:) = P; %pֵ¼��
            end
            INDEXS = computeval1.indexstarget;
            vmat = computeval1.vtarget;
            for i = 1:NROI
                if i<10
                    OutPname = fullfile(Outdir,['Pval_ROI0000',num2str(i),'.nii']);
                elseif i<100
                    OutPname = fullfile(Outdir,['Pval_ROI000',num2str(i),'.nii']);
                elseif i<1000
                    OutPname = fullfile(Outdir,['Pval_ROI00',num2str(i),'.nii']);
                else
                    error('too many ROIs');
                end
                Ptemp = zeros(vmat.dim);
                Ptemp(INDEXS) = Pval(i,:);
                DynamicBC_write_NIFTI(Ptemp,vmat,OutPname);
            end
        end
    else % with cov
        if Calmethod1==1 % pearson & withcov
            for i = 1:NROI
                sigseedt = [datSeed1(i,:),datSeed2(i,:)]'; % ��i��ROI
                deltr = zeros(Permnum,size(datTarget1,1)); % ����delta-r
                Rtemp1 = zeros(Permnum,size(datTarget1,1));  
                Rtemp2 = zeros(Permnum,size(datTarget1,1));  
                covt = [SIG1.COV;SIG2.COV]; % Э����
                parfor iperm = 1:Permnum % ����û�
                    Nrandp = randperm(Ntotal);
                    seedsig1 = sigseedt(Nrandp(1:N1));
                    seedsig2 = sigseedt(Nrandp(N1+1:Ntotal));
                    targetsig1 = dattarget(Nrandp(1:N1),:);
                    targetsig2 = dattarget(Nrandp(N1+1:Ntotal),:);
                    covt1 = covt(Nrandp(1:N1),:);
                    covt2 = covt(Nrandp(N1+1:Ntotal),:);
                    [r1 p1] = partialcorr(seedsig1,targetsig1,covt1);
                    [r2 p2] = partialcorr(seedsig2,targetsig2,covt2);
                    deltr(iperm,:) = r1-r2;
                    
                    Rtemp1(iperm,:) = r1; 
                    Rtemp2(iperm,:) = r2;
                end
                save([Outdir,filesep,'ExtRval_',num2str(i),'.mat'],'Rtemp1','Rtemp2');
                [mu,sig,~,sigci] = normfit(deltr);
                P = normcdf(DELTR(i,:),mu,sig);
                Pval(i,:) = P;
            end
            INDEXS = computeval1.indexstarget;
            vmat = computeval1.vtarget;
            for i = 1:NROI
                if i<10
                    OutPname = fullfile(Outdir,['Pval_ROI0000',num2str(i),'.nii']);
                elseif i<100
                    OutPname = fullfile(Outdir,['Pval_ROI000',num2str(i),'.nii']);
                elseif i<1000
                    OutPname = fullfile(Outdir,['Pval_ROI00',num2str(i),'.nii']);
                else
                    error('too many ROIs');
                end
                Ptemp = zeros(vmat.dim);
                Ptemp(INDEXS) = Pval(i,:);
                DynamicBC_write_NIFTI(Ptemp,vmat,OutPname);
            end
        else % partialcorr with cov
            for i = 1:NROI
                sigseedt = [datSeed1(i,:),datSeed2(i,:)]';
                deltr = zeros(Permnum,size(datTarget1,1));
                Rtemp1 = zeros(Permnum,size(datTarget1,1));  
                Rtemp2 = zeros(Permnum,size(datTarget1,1));  
                covtCOV = [SIG1.COV;SIG2.COV];
                INDUSED = 1:NROI;INDUSED(i) = [];
                covt = [covtCOV,[datSeed1(INDUSED,:),datSeed2(INDUSED,:)]']; % cov��COV������seed�źŹ���
                parfor iperm = 1:Permnum
                    Nrandp = randperm(Ntotal);
                    seedsig1 = sigseedt(Nrandp(1:N1));
                    seedsig2 = sigseedt(Nrandp(N1+1:Ntotal));
                    targetsig1 = dattarget(Nrandp(1:N1),:);
                    targetsig2 = dattarget(Nrandp(N1+1:Ntotal),:);
                    covt1 = covt(Nrandp(1:N1),:);
                    covt2 = covt(Nrandp(N1+1:Ntotal),:);
                    [r1 p1] = partialcorr(seedsig1,targetsig1,covt1);
                    [r2 p2] = partialcorr(seedsig2,targetsig2,covt2);
                    deltr(iperm,:) = r1-r2;
                    Rtemp1(iperm,:) = r1; 
                    Rtemp2(iperm,:) = r2;
                end
                save([Outdir,filesep,'ExtRval_',num2str(i),'.mat'],'Rtemp1','Rtemp2');
                [mu,sig,~,sigci] = normfit(deltr);
                P = normcdf(DELTR(i,:),mu,sig);
                Pval(i,:) = P;
            end
            INDEXS = computeval1.indexstarget;
            vmat = computeval1.vtarget;
            for i = 1:NROI
                if i<10
                    OutPname = fullfile(Outdir,['Pval_ROI0000',num2str(i),'.nii']);
                elseif i<100
                    OutPname = fullfile(Outdir,['Pval_ROI000',num2str(i),'.nii']);
                elseif i<1000
                    OutPname = fullfile(Outdir,['Pval_ROI00',num2str(i),'.nii']);
                else
                    error('too many ROIs');
                end
                Ptemp = zeros(vmat.dim);
                Ptemp(INDEXS) = Pval(i,:);
                DynamicBC_write_NIFTI(Ptemp,vmat,OutPname);
            end
        end
    end
else % interaction
    
end
%% Scatter!!
Scatcond1 = dir(fullfile(Indir1,'Scatter_computeval.mat'));
Scatcond2 = dir(fullfile(Indir2,'Scatter_computeval.mat'));
if isempty(Scatcond1)||isempty(Scatcond2)
    return;
end
Scat_comp1 = load(fullfile(Indir1,'Scatter_computeval.mat'));
Scat_comp2 = load(fullfile(Indir2,'Scatter_computeval.mat'));
ScatR1 = Scat_comp1.r;
ScatR2 = Scat_comp2.r;
clear deltr
scatDELTR = ScatR1-ScatR2;
if PermMethod==1
    Permnum = Parameter.permnum; % permutation��Ŀ
    if COVCOND==0 % no cov
        if Calmethod1==1 % pearson & no cov % ��ʱ�ź�Ϊԭʼ�źţ�û���κ�����������Ҫ����
            for i = 1:NROI % ���ROI
                LABS1 = Scat_comp1.LABS{i}; % G1 ��i��ROI��Ӧ�ı�ǩ ������Ŀ����target��Ŀ
                LABS2 = Scat_comp2.LABS{i}; % G2 ��i��ROI��Ӧ�ı�ǩ
                for ivox = 1:size(datTarget1,1) % ÿ��target����
                    ind1 = ~LABS1(:,ivox); % G1 ��ivox�����أ�ʣ��������
                    ind2 = ~LABS2(:,ivox); % G2 ��ivox�����أ�ʣ��������
                    N1 = nnz(ind1); % G1 ʣ��㣨���ԣ���Ŀ
                    N2 = nnz(ind2); % G2 ʣ��㣨���ԣ���Ŀ
                    Ntotal = N1+N2; % G1+G2����Ŀ
                    sigseedt = [datSeed1(i,ind1),datSeed2(i,ind2)]'; % ƴ�ӵ�i��ROI��seed
                    dattarget = [datTarget1(ivox,ind1),datTarget2(ivox,ind2)]'; % ƴ�ӵ�ivox�����ص�ʱ������
                    parfor iperm = 1:Permnum % �û�
                        Nrandp = randperm(Ntotal);
                        seedsig1 = sigseedt(Nrandp(1:N1));
                        seedsig2 = sigseedt(Nrandp(N1+1:Ntotal));
                        targetsig1 = dattarget(Nrandp(1:N1),:);
                        targetsig2 = dattarget(Nrandp(N1+1:Ntotal),:);
                        [r1 p1] = corr(seedsig1,targetsig1);
                        [r2 p2] = corr(seedsig2,targetsig2);
                        deltr(iperm) = r1-r2;
                    end
                    [mu,sig,~,sigci] = normfit(deltr);
                    P = normcdf(scatDELTR(i,ivox),mu,sig);
                    Pval(i,ivox) = P;
                end
            end
            INDEXS = computeval1.indexstarget;
            vmat = computeval1.vtarget;
            for i = 1:NROI
                if i<10
                    OutPname = fullfile(Outdir,['Scatt_Pval_ROI0000',num2str(i),'.nii']);
                elseif i<100
                    OutPname = fullfile(Outdir,['Scatt_Pval_ROI000',num2str(i),'.nii']);
                elseif i<1000
                    OutPname = fullfile(Outdir,['Scatt_Pval_ROI00',num2str(i),'.nii']);
                else
                    error('too many ROIs');
                end
                Ptemp = zeros(vmat.dim);
                Ptemp(INDEXS) = Pval(i,:);
                DynamicBC_write_NIFTI(Ptemp,vmat,OutPname);
            end
        else % partialcorr & nocov           
            for i = 1:NROI
                LABS1 = Scat_comp1.LABS{i};
                LABS2 = Scat_comp2.LABS{i};
                ScatDat1 = load(fullfile(Indir1,['Scatter_ROI',num2str(i),'_Partial_NoCOVsignal.mat'])); %ÿ��ROI����Э������ͬ����Ҫ��������
                ScatDat2 = load(fullfile(Indir2,['Scatter_ROI',num2str(i),'_Partial_NoCOVsignal.mat'])); %
                datSeed1 = ScatDat1.datSeedN; % seed����
                datSeed2 = ScatDat2.datSeedN; %
                datTarget1 = ScatDat1.datTargetN; % target����
                datTarget2 = ScatDat2.datTargetN;
                for ivox = 1:size(datTarget1,1)
                    ind1 = ~LABS1(:,ivox);
                    ind2 = ~LABS2(:,ivox);
                    N1 = nnz(ind1);
                    N2 = nnz(ind2);
                    Ntotal = N1+N2;
                    sigseedt = [datSeed1(ind1);datSeed2(ind2)];
                    dattarget = [datTarget1(ind1,ivox);datTarget2(ind2,ivox)];
                    parfor iperm = 1:Permnum
                        Nrandp = randperm(Ntotal);
                        seedsig1 = sigseedt(Nrandp(1:N1));
                        seedsig2 = sigseedt(Nrandp(N1+1:Ntotal));
                        targetsig1 = dattarget(Nrandp(1:N1),:);
                        targetsig2 = dattarget(Nrandp(N1+1:Ntotal),:);
                        [r1 p1] = corr(seedsig1,targetsig1); % �����Ѿ�ȥ����Э��������˴˴�ֻ��ʹ��corr
                        [r2 p2] = corr(seedsig2,targetsig2);
                        deltr(iperm) = r1-r2;
                    end
                    [mu,sig,~,sigci] = normfit(deltr);
                    P = normcdf(scatDELTR(i,ivox),mu,sig);
                    Pval(i,ivox) = P;
                end
            end
            INDEXS = computeval1.indexstarget;
            vmat = computeval1.vtarget;
            for i = 1:NROI
                if i<10
                    OutPname = fullfile(Outdir,['Scatt_Pval_ROI0000',num2str(i),'.nii']);
                elseif i<100
                    OutPname = fullfile(Outdir,['Scatt_Pval_ROI000',num2str(i),'.nii']);
                elseif i<1000
                    OutPname = fullfile(Outdir,['Scatt_Pval_ROI00',num2str(i),'.nii']);
                else
                    error('too many ROIs');
                end
                Ptemp = zeros(vmat.dim);
                Ptemp(INDEXS) = Pval(i,:);
                DynamicBC_write_NIFTI(Ptemp,vmat,OutPname);
            end
        end
    else % withcov & 
        if Calmethod1==1 % pearson & withcov
            scatt_signal1 = load(fullfile(Indir1,'Scatter_PearsonCOVsignal.mat')); % ����ȥ��ͬ����Э����������źŲ��õ�������
            scatt_signal2 = load(fullfile(Indir2,'Scatter_PearsonCOVsignal.mat')); % 
            Scatt_datTarget1 = scatt_signal1.datTargetN; % ��ȡtarget
            Scatt_datTarget2 = scatt_signal2.datTargetN; 
            Scatt_datSeed1 = scatt_signal1.datSeedN;
            Scatt_datSeed2 = scatt_signal2.datSeedN;
            for i = 1:NROI
                LABS1 = Scat_comp1.LABS{i};
                LABS2 = Scat_comp2.LABS{i};
                for ivox = 1:size(datTarget1,1)
                    ind1 = ~LABS1(:,ivox);
                    ind2 = ~LABS2(:,ivox);
                    N1 = nnz(ind1);
                    N2 = nnz(ind2);
                    Ntotal = N1+N2;
                    sigseedt = [Scatt_datSeed1(ind1,i);Scatt_datSeed2(ind2,i)];
                    dattarget = [Scatt_datTarget1(ind1,ivox);Scatt_datTarget2(ind2,ivox)];
                    parfor iperm = 1:Permnum
                        Nrandp = randperm(Ntotal);
                        seedsig1 = sigseedt(Nrandp(1:N1));
                        seedsig2 = sigseedt(Nrandp(N1+1:Ntotal));
                        targetsig1 = dattarget(Nrandp(1:N1),:);
                        targetsig2 = dattarget(Nrandp(N1+1:Ntotal),:);
                        [r1 p1] = corr(seedsig1,targetsig1);
                        [r2 p2] = corr(seedsig2,targetsig2);
                        deltr(iperm) = r1-r2;
                    end
                    [mu,sig,~,sigci] = normfit(deltr);
                    P = normcdf(scatDELTR(i,ivox),mu,sig);
                    Pval(i,ivox) = P;
                end
            end
            INDEXS = computeval1.indexstarget;
            vmat = computeval1.vtarget;
            for i = 1:NROI
                if i<10
                    OutPname = fullfile(Outdir,['Scatt_Pval_ROI0000',num2str(i),'.nii']);
                elseif i<100
                    OutPname = fullfile(Outdir,['Scatt_Pval_ROI000',num2str(i),'.nii']);
                elseif i<1000
                    OutPname = fullfile(Outdir,['Scatt_Pval_ROI00',num2str(i),'.nii']);
                else
                    error('too many ROIs');
                end
                Ptemp = zeros(vmat.dim);
                Ptemp(INDEXS) = Pval(i,:);
                DynamicBC_write_NIFTI(Ptemp,vmat,OutPname);
            end
        else % partialcorr & withcov
            for i = 1:NROI
                LABS1 = Scat_comp1.LABS{i};
                LABS2 = Scat_comp2.LABS{i};
                ScatDat1 = load(fullfile(Indir1,['Scatter_ROI',num2str(i),'_Partial_WithCOVsignal.mat'])); % ÿ��ROI�в�ͬ��Э����
                ScatDat2 = load(fullfile(Indir2,['Scatter_ROI',num2str(i),'_Partial_WithCOVsignal.mat'])); 
                datSeed1 = ScatDat1.datSeedN;
                datSeed2 = ScatDat2.datSeedN;
                datTarget1 = ScatDat1.datTargetN;
                datTarget2 = ScatDat2.datTargetN;
                for ivox = 1:size(datTarget1,1)
                    ind1 = ~LABS1(:,ivox);
                    ind2 = ~LABS2(:,ivox);
                    N1 = nnz(ind1);
                    N2 = nnz(ind2);
                    Ntotal = N1+N2;
                    sigseedt = [datSeed1(ind1);datSeed2(ind2)];
                    dattarget = [datTarget1(ind1,ivox);datTarget2(ind2,ivox)];
                    parfor iperm = 1:Permnum
                        Nrandp = randperm(Ntotal);
                        seedsig1 = sigseedt(Nrandp(1:N1));
                        seedsig2 = sigseedt(Nrandp(N1+1:Ntotal));
                        targetsig1 = dattarget(Nrandp(1:N1),:);
                        targetsig2 = dattarget(Nrandp(N1+1:Ntotal),:);
                        [r1 p1] = corr(seedsig1,targetsig1);
                        [r2 p2] = corr(seedsig2,targetsig2);
                        deltr(iperm) = r1-r2;
                    end
                    [mu,sig,~,sigci] = normfit(deltr);
                    P = normcdf(scatDELTR(i,ivox),mu,sig);
                    Pval(i,ivox) = P;
                end
            end
            INDEXS = computeval1.indexstarget;
            vmat = computeval1.vtarget;
            for i = 1:NROI
                if i<10
                    OutPname = fullfile(Outdir,['Scatt_Pval_ROI0000',num2str(i),'.nii']);
                elseif i<100
                    OutPname = fullfile(Outdir,['Scatt_Pval_ROI000',num2str(i),'.nii']);
                elseif i<1000
                    OutPname = fullfile(Outdir,['Scatt_Pval_ROI00',num2str(i),'.nii']);
                else
                    error('too many ROIs');
                end
                Ptemp = zeros(vmat.dim);
                Ptemp(INDEXS) = Pval(i,:);
                DynamicBC_write_NIFTI(Ptemp,vmat,OutPname);
            end
        end
    end
else % interaction
    
end
end