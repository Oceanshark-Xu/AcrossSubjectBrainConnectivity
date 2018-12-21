function AS_MapFCshow_used_simpleU(outdir,inmixdir,backgroundmap,OutnameLab,pthrv,labmarkinput,dof1,dof20,colmaxt,CLUSNUM,ShowtypeTry)
[pat,nam.ext] = fileparts(which('AS_MapFCshow.m'));
% backgroundmap = fullfile(pat,'mni_icbm152_t1_tal_nlin_asym_09a.nii');
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

pthr = pthrv;
dof2 = dof20;
labmark =labmarkinput;
dof = dof1;
[Z, P] = AS_TFRtoZ(dmid,labmark,dof,dof2);
% P1 = P(Z>0);
% P2 = P(Z<0);
% [h1 pi1] = fdr(P1, 0.05);
% [h2 pi2] = fdr(P2, 0.05);
% PIND1 = find(P<pthr);
% Texist = dmid(PIND1);
[Zout,out] = AS_PtoTRF(pthr,labmark,dof,[]);
% maxvabs = max(abs(dmid));
colmint = min(colmaxt);
colmaxt = max(colmaxt);
Valthr = [out,colmaxt];

OutshowRes = reshape(dmid,vmid.dim(1),vmid.dim(2),vmid.dim(3));
%%
% load(fullfile(pat,'SEEDCOLOR.mat'));
load(fullfile(pat,'graycolmap.mat'));
colormapshow0 = AFNICOLORMAP(64);
% colormapshowO = AFNICOLORMAP(12);
colormapshow(1:32,:) = colormapshow0(1:32,:);
colormapshow(33:96,:) = colgray;
colormapshow(97:128,:) = colormapshow0(33:64,:);
dshownum = Valthr(2)-Valthr(1);
%%
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
% indnew = find(abs(OutshowRes)>out);
% [indnewx, indnewy, indnewz] = ind2sub(vmid.dim,indnew);
indnew = find(abs(OutshowResNEW)>out);
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
%     MATOUTT(induextendx,induextendy,induextendz) = OutshowRes(indnew(i));
    MATOUTT(induextendx,induextendy,induextendz) = OutshowResNEW(indnew(i));
end
if ShowtypeTry(1)
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
        %     figure;imagesc(Doutshowtemp)
        Doutmatout = squeeze(MATOUTT(:,:,i));
        ind_postemp = find(Doutmatout>0);
        ind_negtemp = find(Doutmatout<0);
        Doutshowtemp(ind_postemp) = Doutmatout(ind_postemp)+dshownum-out;
        Doutshowtemp(ind_negtemp) = Doutmatout(ind_negtemp)+out-dshownum;
        
        Doutshowtemp1 = rot90(Doutshowtemp);
        H = figure('pos',[100,100,size(Doutshowtemp1,2)*2,size(Doutshowtemp1,1)*2]);
        imagesc(Doutshowtemp1,[-dshownum*2,dshownum*2]);colormap(colormapshow);
        axis off;
        saveas(H,[outdir,filesep,[OutnameLab,'_cor_z',num2str(i),'.fig']])
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
        print(H,[outdir,filesep,[OutnameLab,'_cor_z',num2str(i),'.tif']],'-dtiff','-r300')
        close(H)
    end
    for i = 1:size(dbgre,2)
        Doutshowtemp = squeeze(dbgre(:,i,:));
        Doutshowtemp(Doutshowtemp<30*0.9) = 0;
        Doutshowtemp(Doutshowtemp>0&Doutshowtemp<=30*0.9) = 30+(dshownum)/(85-30);
        Doutshowtemp(Doutshowtemp>=85) = 85-dshownum/55;
        Doutshowtemp(Doutshowtemp>0) = (Doutshowtemp(Doutshowtemp>0)-30)/(85-30)*(dshownum*0.9)+(dshownum*0.05);
        Doutshowtemp = (Doutshowtemp-dshownum/2)*2;
        %     figure;imagesc(Doutshowtemp)
        Doutmatout = squeeze(MATOUTT(:,i,:));
        ind_postemp = find(Doutmatout>0);
        ind_negtemp = find(Doutmatout<0);
        Doutshowtemp(ind_postemp) = Doutmatout(ind_postemp)+dshownum-out;
        Doutshowtemp(ind_negtemp) = Doutmatout(ind_negtemp)+out;
        
        Doutshowtemp1 = rot90(Doutshowtemp);
        H = figure('pos',[100,100,size(Doutshowtemp1,2)*2,size(Doutshowtemp1,1)*2]);
        imagesc(Doutshowtemp1,[-dshownum*2,dshownum*2]);colormap(colormapshow);
        axis off;
        saveas(H,[outdir,filesep,[OutnameLab,'_cor_y',num2str(i),'.fig']])
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
        print(H,[outdir,filesep,[OutnameLab,'_cor_y',num2str(i),'.tif']],'-dtiff','-r300')
        close(H)
    end
    for i = 1:size(dbgre,1)
        Doutshowtemp = squeeze(dbgre(i,:,:));
        Doutshowtemp(Doutshowtemp<30*0.9) = 0;
        Doutshowtemp(Doutshowtemp>0&Doutshowtemp<=30*0.9) = 30+(dshownum)/(85-30);
        Doutshowtemp(Doutshowtemp>=85) = 85-dshownum/55;
        Doutshowtemp(Doutshowtemp>0) = (Doutshowtemp(Doutshowtemp>0)-30)/(85-30)*(dshownum*0.9)+(dshownum*0.05);
        Doutshowtemp = (Doutshowtemp-dshownum/2)*2;
        %     figure;imagesc(Doutshowtemp)
        Doutmatout = squeeze(MATOUTT(i,:,:));
        ind_postemp = find(Doutmatout>0);
        ind_negtemp = find(Doutmatout<0);
        Doutshowtemp(ind_postemp) = Doutmatout(ind_postemp)+dshownum-out;
        Doutshowtemp(ind_negtemp) = Doutmatout(ind_negtemp)+out;
        
        Doutshowtemp1 = rot90(Doutshowtemp);
        H = figure('pos',[100,100,size(Doutshowtemp1,2)*2,size(Doutshowtemp1,1)*2]);
        imagesc(Doutshowtemp1,[-dshownum*2,dshownum*2]);colormap(colormapshow);
        axis off;
        saveas(H,[outdir,filesep,[OutnameLab,'_cor_x',num2str(i),'.fig']])
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
        print(H,[outdir,filesep,[OutnameLab,'_cor_x',num2str(i),'.tif']],'-dtiff','-r300')
        close(H)
    end
elseif ShowtypeTry(2)
    colormapshow(1:32,:) = [];
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
        %     figure;imagesc(Doutshowtemp)
        Doutmatout = squeeze(MATOUTT(:,:,i));
        ind_postemp = find(Doutmatout>0);
%         ind_negtemp = find(Doutmatout<0);
        Doutshowtemp(ind_postemp) = Doutmatout(ind_postemp)+dshownum-out;
%         Doutshowtemp(ind_negtemp) = Doutmatout(ind_negtemp)+out-dshownum;
        
        Doutshowtemp1 = rot90(Doutshowtemp);
%         Doutshowtemp1(Douttemp1<-dshownum) = -dshownum;
        H = figure('pos',[100,100,size(Doutshowtemp1,2)*2,size(Doutshowtemp1,1)*2]);
        imagesc(Doutshowtemp1,[-dshownum*1,dshownum*2]);colormap(colormapshow);
        axis off;
        saveas(H,[outdir,filesep,[OutnameLab,'_cor_z',num2str(i),'.fig']])
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
        print(H,[outdir,filesep,[OutnameLab,'_cor_z',num2str(i),'.tif']],'-dtiff','-r300')
        close(H)
    end
    for i = 1:size(dbgre,2)
        Doutshowtemp = squeeze(dbgre(:,i,:));
        Doutshowtemp(Doutshowtemp<30*0.9) = 0;
        Doutshowtemp(Doutshowtemp>0&Doutshowtemp<=30*0.9) = 30+(dshownum)/(85-30);
        Doutshowtemp(Doutshowtemp>=85) = 85-dshownum/55;
        Doutshowtemp(Doutshowtemp>0) = (Doutshowtemp(Doutshowtemp>0)-30)/(85-30)*(dshownum*0.9)+(dshownum*0.05);
        Doutshowtemp = (Doutshowtemp-dshownum/2)*2;
        %     figure;imagesc(Doutshowtemp)
        Doutmatout = squeeze(MATOUTT(:,i,:));
        ind_postemp = find(Doutmatout>0);
%         ind_negtemp = find(Doutmatout<0);
        Doutshowtemp(ind_postemp) = Doutmatout(ind_postemp)+dshownum-out;
%         Doutshowtemp(ind_negtemp) = Doutmatout(ind_negtemp)+out-dshownum;
        
        Doutshowtemp1 = rot90(Doutshowtemp);
%         Doutshowtemp1(Douttemp1<-dshownum) = -dshownum;
        H = figure('pos',[100,100,size(Doutshowtemp1,2)*2,size(Doutshowtemp1,1)*2]);
        imagesc(Doutshowtemp1,[-dshownum*1,dshownum*2]);colormap(colormapshow);
        axis off;
        saveas(H,[outdir,filesep,[OutnameLab,'_cor_y',num2str(i),'.fig']])
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
        print(H,[outdir,filesep,[OutnameLab,'_cor_y',num2str(i),'.tif']],'-dtiff','-r300')
        close(H)
    end
    for i = 1:size(dbgre,1)
        Doutshowtemp = squeeze(dbgre(i,:,:));
        Doutshowtemp(Doutshowtemp<30*0.9) = 0;
        Doutshowtemp(Doutshowtemp>0&Doutshowtemp<=30*0.9) = 30+(dshownum)/(85-30);
        Doutshowtemp(Doutshowtemp>=85) = 85-dshownum/55;
        Doutshowtemp(Doutshowtemp>0) = (Doutshowtemp(Doutshowtemp>0)-30)/(85-30)*(dshownum*0.9)+(dshownum*0.05);
        Doutshowtemp = (Doutshowtemp-dshownum/2)*2;
        %     figure;imagesc(Doutshowtemp)
        Doutmatout = squeeze(MATOUTT(i,:,:));
        ind_postemp = find(Doutmatout>0);
%         ind_negtemp = find(Doutmatout<0);
        Doutshowtemp(ind_postemp) = Doutmatout(ind_postemp)+dshownum-out;
%         Doutshowtemp(ind_negtemp) = Doutmatout(ind_negtemp)+out-dshownum;
        
        Doutshowtemp1 = rot90(Doutshowtemp);
%         Doutshowtemp1(Douttemp1<-dshownum) = -dshownum;
        H = figure('pos',[100,100,size(Doutshowtemp1,2)*2,size(Doutshowtemp1,1)*2]);
        imagesc(Doutshowtemp1,[-dshownum*1,dshownum*2]);colormap(colormapshow);
        axis off;
        saveas(H,[outdir,filesep,[OutnameLab,'_cor_x',num2str(i),'.fig']])
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
        print(H,[outdir,filesep,[OutnameLab,'_cor_x',num2str(i),'.tif']],'-dtiff','-r300')
        close(H)
    end
elseif Showtype(3)
    
    colormapshow(97:128,:) = [];
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
        %     figure;imagesc(Doutshowtemp)
        Doutmatout = squeeze(MATOUTT(:,:,i));
%         ind_postemp = find(Doutmatout>0);
        ind_negtemp = find(Doutmatout<0);
%         Doutshowtemp(ind_postemp) = Doutmatout(ind_postemp)+dshownum-out;
        Doutshowtemp(ind_negtemp) = Doutmatout(ind_negtemp)+out-dshownum;
        
        Doutshowtemp1 = rot90(Doutshowtemp);
%         Doutshowtemp1(Douttemp1>dshownum) = dshownum;
        H = figure('pos',[100,100,size(Doutshowtemp1,2)*2,size(Doutshowtemp1,1)*2]);
        imagesc(Doutshowtemp1,[-dshownum*2,dshownum*1]);colormap(colormapshow);
        axis off;
        saveas(H,[outdir,filesep,[OutnameLab,'_cor_z',num2str(i),'.fig']])
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
        print(H,[outdir,filesep,[OutnameLab,'_cor_z',num2str(i),'.tif']],'-dtiff','-r300')
        close(H)
    end
    for i = 1:size(dbgre,2)
        Doutshowtemp = squeeze(dbgre(:,i,:));
        Doutshowtemp(Doutshowtemp<30*0.9) = 0;
        Doutshowtemp(Doutshowtemp>0&Doutshowtemp<=30*0.9) = 30+(dshownum)/(85-30);
        Doutshowtemp(Doutshowtemp>=85) = 85-dshownum/55;
        Doutshowtemp(Doutshowtemp>0) = (Doutshowtemp(Doutshowtemp>0)-30)/(85-30)*(dshownum*0.9)+(dshownum*0.05);
        Doutshowtemp = (Doutshowtemp-dshownum/2)*2;
        %     figure;imagesc(Doutshowtemp)
        Doutmatout = squeeze(MATOUTT(:,i,:));
%         ind_postemp = find(Doutmatout>0);
        ind_negtemp = find(Doutmatout<0);
%         Doutshowtemp(ind_postemp) = Doutmatout(ind_postemp)+dshownum-out;
        Doutshowtemp(ind_negtemp) = Doutmatout(ind_negtemp)+out-dshownum;
        
        Doutshowtemp1 = rot90(Doutshowtemp);
%         Doutshowtemp1(Douttemp1>dshownum) = dshownum;
        H = figure('pos',[100,100,size(Doutshowtemp1,2)*2,size(Doutshowtemp1,1)*2]);
        imagesc(Doutshowtemp1,[-dshownum*2,dshownum*1]);colormap(colormapshow);
        axis off;
        saveas(H,[outdir,filesep,[OutnameLab,'_cor_y',num2str(i),'.fig']])
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
        print(H,[outdir,filesep,[OutnameLab,'_cor_y',num2str(i),'.tif']],'-dtiff','-r300')
        close(H)
    end
    for i = 1:size(dbgre,1)
        Doutshowtemp = squeeze(dbgre(i,:,:));
        Doutshowtemp(Doutshowtemp<30*0.9) = 0;
        Doutshowtemp(Doutshowtemp>0&Doutshowtemp<=30*0.9) = 30+(dshownum)/(85-30);
        Doutshowtemp(Doutshowtemp>=85) = 85-dshownum/55;
        Doutshowtemp(Doutshowtemp>0) = (Doutshowtemp(Doutshowtemp>0)-30)/(85-30)*(dshownum*0.9)+(dshownum*0.05);
        Doutshowtemp = (Doutshowtemp-dshownum/2)*2;
        %     figure;imagesc(Doutshowtemp)
        Doutmatout = squeeze(MATOUTT(i,:,:));
%         ind_postemp = find(Doutmatout>0);
        ind_negtemp = find(Doutmatout<0);
%         Doutshowtemp(ind_postemp) = Doutmatout(ind_postemp)+dshownum-out;
        Doutshowtemp(ind_negtemp) = Doutmatout(ind_negtemp)+out-dshownum;
        
        Doutshowtemp1 = rot90(Doutshowtemp);
%         Doutshowtemp1(Douttemp1>dshownum) = dshownum;
        H = figure('pos',[100,100,size(Doutshowtemp1,2)*2,size(Doutshowtemp1,1)*2]);
        imagesc(Doutshowtemp1,[-dshownum*2,dshownum*1]);colormap(colormapshow);
        axis off;
        saveas(H,[outdir,filesep,[OutnameLab,'_cor_x',num2str(i),'.fig']])
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
        print(H,[outdir,filesep,[OutnameLab,'_cor_x',num2str(i),'.tif']],'-dtiff','-r300')
        close(H)
    end
end
