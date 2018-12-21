function [DOUTSHOWaxial,DOUTSHOWcornoral,DOUTSHOWsagittal,DOUTSHOWaxialPOS,DOUTSHOWcornoralPOS,DOUTSHOWsagittalPOS,DOUTSHOWaxialNEG,DOUTSHOWcornoralNEG,DOUTSHOWsagittalNEG,dshownum]=AS_MapFCshow3sliceUP(outdir,inmixdir,backgroundmap,pthrv,colmaxt,CLUSNUM)
uiwait(msgbox('This is only for the permutation Pvalue!'));
[pat,nam.ext] = fileparts(which('AS_MapFCshow.m'));
% backgroundmap = fullfile(pat,'mni_icbm152_t1_tal_nlin_asym_09a.nii');
[vbg dbg] = Dynamic_read_dir_NIFTI(backgroundmap);
dbgre = reshape(dbg,vbg.dim(1),vbg.dim(2),vbg.dim(3));
[vmid dmid] = Dynamic_read_dir_NIFTI(inmixdir);
% if vbg.mat(1)<0
%     dbgre = flipdim(dbgre,1);
% end
% if vmid.mat(1)<0
%     dmid2 = reshape(dmid,vmid.dim(1),vmid.dim(2),vmid.dim(3));
%     dmid2 = flipdim(dmid2,1);
%     dmid = reshape(dmid2,vmid.dim(1)*vmid.dim(2)*vmid.dim(3),1);
% end
dmid(dmid<0.5&dmid>0) = dmid(dmid<0.5&dmid>0)-1;
dmid(dmid==0.5) = 0;
% dmid(dmid>0.5) = 1-dmid(dmid>0.5);
% maxvshow = max(dmid);
% minvshow = min(dmid);
pthr = pthrv;
% dof2 = dof20;
% labmark =labmarkinput;
% dof = dof1;
% [Z, P] = AS_TFRtoZ(dmid,labmark,dof,dof2);
% [Zout,out] = AS_PtoTRF(pthr,labmark,dof,dof2);
% changed later
% colmaxt = max(colmaxt);
%
Valthr = colmaxt;
out = min(colmaxt);
colmaxt = max(colmaxt);
Hsize = get(0,'ScreenSize');
OutshowRes = reshape(dmid,vmid.dim(1),vmid.dim(2),vmid.dim(3));

load(fullfile(pat,'graycolmap.mat'));
colormapshow0 = AFNICOLORMAP(64);
colormapshow(1:32,:) = colormapshow0(1:32,:);
colormapshow(33:96,:) = colgray;
colormapshow(97:128,:) = colormapshow0(33:64,:);
dshownum = Valthr(2)-Valthr(1);


OUTSHOWRESNew = OutshowRes;
OUTSHOWRESNew(abs(OutshowRes)<=out) = 0;
[L NUM] = bwlabeln(OUTSHOWRESNew,18); % 18: spm connect type.
for i = 1:NUM
    Clusize(i,1) = length(find(L==i));
end
indcsize = find(Clusize>=CLUSNUM);
OUTSHOWRESNew2 = zeros(size(OUTSHOWRESNew));
for i = 1:length(indcsize)
    OUTSHOWRESNew2(L==(indcsize(i))) = 1;
end
OutshowResNEW = OutshowRes.*OUTSHOWRESNew2;
indnew = find(abs(OutshowResNEW)>out);
[indnewx, indnewy, indnewz] = ind2sub(vmid.dim,indnew);
% save temp0706
indnewsub = [indnewx, indnewy, indnewz];
indnewsubT = trans2solution(indnewsub,vmid.mat,vbg.mat);
diffx = abs(vmid.mat(1,1)/vbg.mat(1,1))/2*abs(vbg.mat(1,1));
diffy = abs(vmid.mat(2,2)/vbg.mat(2,2))/2*abs(vbg.mat(2,2));
diffz = abs(vmid.mat(3,3)/vbg.mat(3,3))/2*abs(vbg.mat(3,3));
MATOUTT = zeros(vbg.dim);
for i = 1:size(indnewsub,1)
    indu = indnewsubT(i,:);
    induextendx = round(indu(1)-diffx):round(indu(1)+diffx);
    induextendy = round(indu(2)-diffy):round(indu(2)+diffy);
    induextendz = round(indu(3)-diffz):round(indu(3)+diffz);
    MATOUTT(induextendx,induextendy,induextendz) = OutshowResNEW(indnew(i));
end


for i = 1:size(dbgre,3)
% for i = 70:90
    Doutshowtemp = squeeze(dbgre(:,:,i));
    Doutshowtemp(Doutshowtemp<30*0.9) = 0;
    Doutshowtemp(Doutshowtemp>0&Doutshowtemp<=30*0.9) = 30+(dshownum)/(85-30);
    Doutshowtemp(Doutshowtemp>=85) = 85-dshownum/55;
    Doutshowtemp(Doutshowtemp>0) = (Doutshowtemp(Doutshowtemp>0)-30)/(85-30)*(dshownum*0.9)+(dshownum*0.05);
%     Doutshowtemp(Doutshowtemp>0) = (Doutshowtemp(Doutshowtemp>0)-(max(Doutshowtemp(Doutshowtemp>0))+min(Doutshowtemp(Doutshowtemp>0)))/2)*2;
%     Doutshowtemp = (Doutshowtemp-(max(Doutshowtemp(Doutshowtemp>0))+min(Doutshowtemp(Doutshowtemp>0)))/2)*2;
    Doutshowtemp = (Doutshowtemp-dshownum/2)*2;
    DoutshowtempPOS = Doutshowtemp;
    DoutshowtempNEG = Doutshowtemp;
%     figure;imagesc(Doutshowtemp)
    Doutmatout = squeeze(MATOUTT(:,:,i));
    ind_postemp = find(Doutmatout>0);
    ind_negtemp = find(Doutmatout<0);
    Doutshowtemp(ind_postemp) = Doutmatout(ind_postemp)+dshownum-out;
    DoutshowtempPOS(ind_postemp) = Doutmatout(ind_postemp)+dshownum-out;
    Doutshowtemp(ind_negtemp) = Doutmatout(ind_negtemp)+out-dshownum;
    DoutshowtempNEG(ind_negtemp) = Doutmatout(ind_negtemp)+out-dshownum;
    Doutshowtemp1 = rot90(Doutshowtemp);
    DOUTSHOWaxial(:,:,i) = Doutshowtemp1;
    DOUTSHOWaxialPOS(:,:,i) = rot90(DoutshowtempPOS);
    DOUTSHOWaxialNEG(:,:,i) = rot90(DoutshowtempNEG);
end
for i = 1:size(dbgre,2)
    Doutshowtemp = squeeze(dbgre(:,i,:));
    Doutshowtemp(Doutshowtemp<30*0.9) = 0;
    Doutshowtemp(Doutshowtemp>0&Doutshowtemp<=30*0.9) = 30+(dshownum)/(85-30);
    Doutshowtemp(Doutshowtemp>=85) = 85-dshownum/55;
    Doutshowtemp(Doutshowtemp>0) = (Doutshowtemp(Doutshowtemp>0)-30)/(85-30)*(dshownum*0.9)+(dshownum*0.05);
    Doutshowtemp = (Doutshowtemp-dshownum/2)*2;
    DoutshowtempPOS = Doutshowtemp;
    DoutshowtempNEG = Doutshowtemp;
%     figure;imagesc(Doutshowtemp)
    Doutmatout = squeeze(MATOUTT(:,i,:));
    ind_postemp = find(Doutmatout>0);
    ind_negtemp = find(Doutmatout<0);
    Doutshowtemp(ind_postemp) = Doutmatout(ind_postemp)+dshownum-out;
    DoutshowtempPOS(ind_postemp) = Doutmatout(ind_postemp)+dshownum-out;
    Doutshowtemp(ind_negtemp) = Doutmatout(ind_negtemp)+out-dshownum;
    DoutshowtempNEG(ind_negtemp) = Doutmatout(ind_negtemp)+out-dshownum;
    
    Doutshowtemp1 = rot90(Doutshowtemp);
    DOUTSHOWcornoral(:,:,i) = Doutshowtemp1;
    DOUTSHOWcornoralPOS(:,:,i) = rot90(DoutshowtempPOS);
    DOUTSHOWcornoralNEG(:,:,i) = rot90(DoutshowtempNEG);
end
for i = 1:size(dbgre,1)
    Doutshowtemp = squeeze(dbgre(i,:,:));
    Doutshowtemp(Doutshowtemp<30*0.9) = 0;
    Doutshowtemp(Doutshowtemp>0&Doutshowtemp<=30*0.9) = 30+(dshownum)/(85-30);
    Doutshowtemp(Doutshowtemp>=85) = 85-dshownum/55;
    Doutshowtemp(Doutshowtemp>0) = (Doutshowtemp(Doutshowtemp>0)-30)/(85-30)*(dshownum*0.9)+(dshownum*0.05);
    Doutshowtemp = (Doutshowtemp-dshownum/2)*2;
    DoutshowtempPOS = Doutshowtemp;
    DoutshowtempNEG = Doutshowtemp;
%     figure;imagesc(Doutshowtemp)
    Doutmatout = squeeze(MATOUTT(i,:,:));
    ind_postemp = find(Doutmatout>0);
    ind_negtemp = find(Doutmatout<0);
    Doutshowtemp(ind_postemp) = Doutmatout(ind_postemp)+dshownum-out;
    DoutshowtempPOS(ind_postemp) = Doutmatout(ind_postemp)+dshownum-out;
    Doutshowtemp(ind_negtemp) = Doutmatout(ind_negtemp)+out-dshownum;
    DoutshowtempNEG(ind_negtemp) = Doutmatout(ind_negtemp)+out-dshownum;
    
    Doutshowtemp1 = rot90(Doutshowtemp);
    DOUTSHOWsagittal(:,:,i) = Doutshowtemp1;
    DOUTSHOWsagittalPOS(:,:,i) = rot90(DoutshowtempPOS);
    DOUTSHOWsagittalNEG(:,:,i) = rot90(DoutshowtempNEG);
end
save([outdir,filesep,'Showinfo.mat'],'DOUTSHOWsagittal','DOUTSHOWcornoral','DOUTSHOWaxial')
save([outdir,filesep,'ShowinfoPOS.mat'],'DOUTSHOWsagittalPOS','DOUTSHOWcornoralPOS','DOUTSHOWaxialPOS')
save([outdir,filesep,'ShowinfoNEG.mat'],'DOUTSHOWsagittalNEG','DOUTSHOWcornoralNEG','DOUTSHOWaxialNEG')