function AS_MapFCshow_used_Axial(outdir,inmixdir,OutnameLab,pthrv,labmarkinput,dof1,dof2,colmaxt)
[pat,nam.ext] = fileparts(which('AS_MapFCshow.m'));
backgroundmap = fullfile(pat,'mni_icbm152_t1_tal_nlin_asym_09a.nii');
% [patnam,indir,ext]= uigetfile({'*.nii';'*.img';'*.*'},'Resultmap for show');
% outdir = uigetdir(pwd,'Pic Out Dir');
[vbg dbg] = Dynamic_read_dir_NIFTI(backgroundmap);
dbgre = reshape(dbg,vbg.dim(1),vbg.dim(2),vbg.dim(3));

% inmixdir = fullfile(indir,patnam);
[vmid dmid] = Dynamic_read_dir_NIFTI(inmixdir);
maxvshow = max(dmid);
minvshow = min(dmid);
desipinfo = vmid.descrip;
inddes1 = find(desipinfo=='[');
inddes2 = find(desipinfo==']');
if nargin<3
    error('not enough input');
end    
if nargin<4
    pthrans = inputdlg('threshold','threshold',1,{'0.01'});
    pthr = str2num(pthrans{1});
else
    pthr = pthrv;
end
dof2 = [];

if nargin>=7
    labmark =labmarkinput;
    dof = dof1;
%     dof2 = [];
% elseif nargin==8
%     labmark =labmarkinput;
%     dof = dof1;
elseif isempty(inddes1)
    labmarkans = inputdlg('labels of stat','labels of stat',1,'T');
    labmark = labmarkans{1};
    if strcmpi(labmark,'F')
        dofans = inputdlg('dof of stat(F)','dof of stat(F)',1,{'[10,2]'});
        doftemp = str2num(dofans{1});
        dof = doftemp(1);
        dof2 = doftemp(2);
    else
        dofans = inputdlg('dof of stat','dof of stat',1,{'10'});
        dof = str2num(dofans{1});
        dof2 = [];
    end        
else
    dof = str2num(desipinfo(inddes1+1:inddes2-1));
    labmark = desipinfo(5);
end
[Z, P] = AS_TFRtoZ(dmid,labmark,dof,dof2);
% P1 = P(Z>0);
% P2 = P(Z<0);
% [h1 pi1] = fdr(P1, 0.05);
% [h2 pi2] = fdr(P2, 0.05);
% PIND1 = find(P<pthr);
% Texist = dmid(PIND1);
[Zout,out] = AS_PtoTRF(pthr,labmark,dof,[]);
% maxvabs = max(abs(dmid));
if max(dmid)>out+1&&min(dmid)<-out-1
    maxvabs = min(max(dmid),-min(dmid));
elseif max(dmid)<=out+1&&min(dmid)<-out-1
    maxvabs = -min(dmid);
elseif max(dmid)>out+1&&min(dmid)>=-out-1
    maxvabs = max(dmid);
else
    maxvabs = out+1;
end
    

if nargin<8
    Valthr = [out,maxvabs];
else
    Valthr = [out,colmaxt];
end
OutshowRes = reshape(dmid,vmid.dim(1),vmid.dim(2),vmid.dim(3));
%%
% load(fullfile(pat,'SEEDCOLOR.mat'));
load(fullfile(pat,'graycolmap.mat'));
colormapshow0 = AFNICOLORMAP(128);
% colormapshowO = AFNICOLORMAP(12);
colormapshow(1:64,:) = colormapshow0(1:64,:);
colormapshow(65:128,:) = colgray;
colormapshow(129:192,:) = colormapshow0(65:end,:);
dshownum = Valthr(2)-Valthr(1);
%%
indnew = find(abs(OutshowRes)>out);
[indnewx, indnewy, indnewz] = ind2sub(vmid.dim,indnew);

indnewsub = [indnewx, indnewy, indnewz];
indnewsubT = trans2solution(indnewsub,vmid.mat,vbg.mat);
diffx = abs(vmid.mat(1,1)/vbg.mat(1,1))/2*abs(vbg.mat(1,1));
diffy = abs(vmid.mat(2,2)/vbg.mat(2,2))/2*abs(vbg.mat(2,2));
diffz = abs(vmid.mat(3,3)/vbg.mat(3,3))/2*abs(vbg.mat(3,3));
% INDNEWPLOT = [];
MATOUTT = zeros(vbg.dim);
% indnewsubT = find(abs(OutshowRes)>out);
for i = 1:size(indnewsub,1)
    indu = indnewsubT(i,:);
    induextendx = round(indu(1)-diffx):round(indu(1)+diffx);
    induextendy = round(indu(2)-diffy):round(indu(2)+diffy);
    induextendz = round(indu(3)-diffz):round(indu(3)+diffz);
    MATOUTT(induextendx,induextendy,induextendz) = OutshowRes(indnew(i));
end
for i = 1:size(dbgre,3)
% for i = 100:100
    Doutshowtemp = squeeze(dbgre(:,:,i));
    Doutshowtemp(Doutshowtemp<30*0.9) = 0;
    Doutshowtemp(Doutshowtemp>0&Doutshowtemp<=30*0.9) = 30+(dshownum)/(85-30);
    Doutshowtemp(Doutshowtemp>=85) = 85-dshownum/55;
    Doutshowtemp(Doutshowtemp>0) = (Doutshowtemp(Doutshowtemp>0)-30)/(85-30)*(dshownum*0.9)+(dshownum*0.05);
%     figure;imagesc(Doutshowtemp)
    Doutmatout = squeeze(MATOUTT(:,:,i));
    ind_postemp = find(Doutmatout>0);
    ind_negtemp = find(Doutmatout<0);
    Doutshowtemp(ind_postemp) = Doutmatout(ind_postemp)+dshownum-out;
    Doutshowtemp(ind_negtemp) = Doutmatout(ind_negtemp)+out;
    
    Doutshowtemp1 = rot90(Doutshowtemp);
    H = figure('pos',[100,100,size(Doutshowtemp1,2)*2,size(Doutshowtemp1,1)*2]);
    imagesc(Doutshowtemp1,[-dshownum,dshownum*2]);colormap(colormapshow);
    axis off;
    saveas(H,[outdir,filesep,[OutnameLab,'_z',num2str(i),'.fig']])
    set(H,'PaperPositionMode','manual');
    set(H,'PaperUnits','inch')
    XSIZE = size(Doutshowtemp1,2);
    YSIZE = size(Doutshowtemp1,1);
    factor = 1:100;
    XSIZEnew = XSIZE*factor;
    YSIZEnew = YSIZE*factor;
    FACTORS1 = find(XSIZEnew>400);
    FACTORS2 = find(YSIZEnew>400);
    FACTORS = max(FACTORS1(1),FACTORS2(1));
    XSIZEU = XSIZEnew(FACTORS);
    YSIZEU = YSIZEnew(FACTORS);    
    set(H,'Paperposition',[1,1,XSIZEU/300,YSIZEU/300]);
    print(H,[outdir,filesep,[OutnameLab,'_z',num2str(i),'.tif']],'-dtiff','-r300')
    close(H)
end
% for i = 1:size(dbgre,2)
%     Doutshowtemp = squeeze(dbgre(:,i,:));
%     Doutshowtemp(Doutshowtemp<30*0.9) = 0;
%     Doutshowtemp(Doutshowtemp>0&Doutshowtemp<=30*0.9) = 30+(dshownum)/(85-30);
%     Doutshowtemp(Doutshowtemp>=85) = 85-dshownum/55;
%     Doutshowtemp(Doutshowtemp>0) = (Doutshowtemp(Doutshowtemp>0)-30)/(85-30)*(dshownum*0.9)+(dshownum*0.05);
% %     figure;imagesc(Doutshowtemp)
%     Doutmatout = squeeze(MATOUTT(:,i,:));
%     ind_postemp = find(Doutmatout>0);
%     ind_negtemp = find(Doutmatout<0);
%     Doutshowtemp(ind_postemp) = Doutmatout(ind_postemp)+dshownum-out;
%     Doutshowtemp(ind_negtemp) = Doutmatout(ind_negtemp)+out;
%     
%     Doutshowtemp1 = rot90(Doutshowtemp);
%     H = figure('pos',[100,100,size(Doutshowtemp1,2)*2,size(Doutshowtemp1,1)*2]);
%     imagesc(Doutshowtemp1,[-dshownum,dshownum*2]);colormap(colormapshow);
%     axis off;
%     saveas(H,[outdir,filesep,[OutnameLab,'_y',num2str(i),'.fig']])
%     set(H,'PaperPositionMode','manual');
%     set(H,'PaperUnits','inch')
%     XSIZE = size(Doutshowtemp1,2);
%     YSIZE = size(Doutshowtemp1,1);
%     factor = 1:100;
%     XSIZEnew = XSIZE*factor;
%     YSIZEnew = YSIZE*factor;
%     FACTORS1 = find(XSIZEnew>400);
%     FACTORS2 = find(YSIZEnew>400);
%     FACTORS = max(FACTORS1(1),FACTORS2(1));
%     XSIZEU = XSIZEnew(FACTORS);
%     YSIZEU = YSIZEnew(FACTORS);    
%     set(H,'Paperposition',[1,1,XSIZEU/300,YSIZEU/300]);
%     print(H,[outdir,filesep,[OutnameLab,'_y',num2str(i),'.tif']],'-dtiff','-r300')
%     close(H)
% end
% for i = 1:size(dbgre,1)
%     Doutshowtemp = squeeze(dbgre(i,:,:));
%     Doutshowtemp(Doutshowtemp<30*0.9) = 0;
%     Doutshowtemp(Doutshowtemp>0&Doutshowtemp<=30*0.9) = 30+(dshownum)/(85-30);
%     Doutshowtemp(Doutshowtemp>=85) = 85-dshownum/55;
%     Doutshowtemp(Doutshowtemp>0) = (Doutshowtemp(Doutshowtemp>0)-30)/(85-30)*(dshownum*0.9)+(dshownum*0.05);
% %     figure;imagesc(Doutshowtemp)
%     Doutmatout = squeeze(MATOUTT(i,:,:));
%     ind_postemp = find(Doutmatout>0);
%     ind_negtemp = find(Doutmatout<0);
%     Doutshowtemp(ind_postemp) = Doutmatout(ind_postemp)+dshownum-out;
%     Doutshowtemp(ind_negtemp) = Doutmatout(ind_negtemp)+out;
%     
%     Doutshowtemp1 = rot90(Doutshowtemp);
%     H = figure('pos',[100,100,size(Doutshowtemp1,2)*2,size(Doutshowtemp1,1)*2]);
%     imagesc(Doutshowtemp1,[-dshownum,dshownum*2]);colormap(colormapshow);
%     axis off;
%     saveas(H,[outdir,filesep,[OutnameLab,'_x',num2str(i),'.fig']])
%     set(H,'PaperPositionMode','manual');
%     set(H,'PaperUnits','inch')
%     XSIZE = size(Doutshowtemp1,2);
%     YSIZE = size(Doutshowtemp1,1);
%     factor = 1:100;
%     XSIZEnew = XSIZE*factor;
%     YSIZEnew = YSIZE*factor;
%     FACTORS1 = find(XSIZEnew>400);
%     FACTORS2 = find(YSIZEnew>400);
%     FACTORS = max(FACTORS1(1),FACTORS2(1));
%     XSIZEU = XSIZEnew(FACTORS);
%     YSIZEU = YSIZEnew(FACTORS);    
%     set(H,'Paperposition',[1,1,XSIZEU/300,YSIZEU/300]);
%     print(H,[outdir,filesep,[OutnameLab,'_x',num2str(i),'.tif']],'-dtiff','-r300')
%     close(H)
% end

