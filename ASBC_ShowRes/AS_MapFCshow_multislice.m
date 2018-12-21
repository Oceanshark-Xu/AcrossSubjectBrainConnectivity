function AS_MapFCshow_multislice(outdir,inmixdir,pthrv,labmarkinput,dof1,dof20,colmaxt,CLUSNUM,sliceorder,showposlab,showlab)
[pat,nam.ext] = fileparts(which('AS_MapFCshow.m'));
backgroundmap = fullfile(pat,'mni_icbm152_t1_tal_nlin_asym_09a.nii');
[vbg dbg] = Dynamic_read_dir_NIFTI(backgroundmap);
dbgre = reshape(dbg,vbg.dim(1),vbg.dim(2),vbg.dim(3));
[vmid dmid] = Dynamic_read_dir_NIFTI(inmixdir);
maxvshow = max(dmid);
minvshow = min(dmid);
pthr = pthrv;
dof2 = dof20;
labmark =labmarkinput;
dof = dof1;
[Z, P] = AS_TFRtoZ(dmid,labmark,dof,dof2);
[Zout,out] = AS_PtoTRF(pthr,labmark,dof,dof2);
Valthr = [out,colmaxt];
Hsize = get(0,'ScreenSize');
OutshowRes = reshape(dmid,vmid.dim(1),vmid.dim(2),vmid.dim(3));
%%
load(fullfile(pat,'graycolmap.mat'));
colormapshow0 = AFNICOLORMAP(64);
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
indnew = find(abs(OutshowResNEW)>out);
[indnewx, indnewy, indnewz] = ind2sub(vmid.dim,indnew);

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

if size(sliceorder,1)<size(sliceorder,2)
    sliceorder = sliceorder';
end

NW = ceil(sqrt(length(sliceorder)));
for i = 1:length(sliceorder)
    if mod(i,NW)==0
        POSX(i,1) = i/NW;
        POSY(i,1) = NW;
    else
        POSX(i,1) = floor(i/NW)+1;
        POSY(i,1) = i-NW*(POSX(i,1)-1);
    end
end
POSXY = [POSX,POSY];
NL = length(unique(POSX));
switch showposlab
    case 1 % axiel
        if showlab==1 % mni condition
            Showlabtemp = [zeros(length(sliceorder),1),zeros(length(sliceorder),1),sliceorder];
            Showlab = mni2cor(Showlabtemp,vbg.mat);
            SliceorderCor = Showlab(:,3);
            SliceorderMNI = sliceorder;
        else
            SliceorderCor = sliceorder;
            Sliceordertemp = cor2mni(sliceorder,vbg.mat);
            SliceorderMNI = Sliceordertemp(:,3);
        end
        Doutshowtemp2 = zeros(size(dbgre,2)*NW,size(dbgre,1)*NL);
        for i = SliceorderCor
            Doutshowtemp = squeeze(dbgre(:,:,i));
            Doutshowtemp(Doutshowtemp<30*0.9) = 0;
            Doutshowtemp(Doutshowtemp>0&Doutshowtemp<=30*0.9) = 30+(dshownum)/(85-30);
            Doutshowtemp(Doutshowtemp>=85) = 85-dshownum/55;
            Doutshowtemp(Doutshowtemp>0) = (Doutshowtemp(Doutshowtemp>0)-30)/(85-30)*(dshownum*0.9)+(dshownum*0.05);
            Doutshowtemp = (Doutshowtemp-dshownum/2)*2;
            Doutmatout = squeeze(MATOUTT(:,:,i));
            ind_postemp = find(Doutmatout>0);
            ind_negtemp = find(Doutmatout<0);
            Doutshowtemp(ind_postemp) = Doutmatout(ind_postemp)+dshownum-out;
            Doutshowtemp(ind_negtemp) = Doutmatout(ind_negtemp)+out-dshownum;
            Doutshowtemp1 = rot90(Doutshowtemp);
            IX = i-SliceorderCor(1)+1;
            POSXT = POSX(IX);
            POSYT = POSY(IX);
            Doutshowtemp2(1+size(Doutshowtemp1,1)*(POSXT-1):size(Doutshowtemp1,1)*(POSXT),1+size(Doutshowtemp1,2)*(POSYT-1):size(Doutshowtemp1,2)*(POSYT)) = Doutshowtemp1;
        end
        Hexist1 = Hsize(3)-200;
        Hexist2 = Hsize(4)-200;
        Dsize = [size(Doutshowtemp2,2),size(Doutshowtemp2,1)];
        factorD = min(Dsize(1)/Hexist2,Dsize(2)/Hexist1);
        H = figure('pos',[100,100,Dsize(1)*factorD,Dsize(2)*factorD]);
        %             H = figure('pos',[100,100,size(Doutshowtemp1,2)*2,size(Doutshowtemp1,1)*2]);
        
        imagesc(Doutshowtemp2,[-dshownum*2,dshownum*2]);colormap(colormapshow);
        axis off;
        saveas(H,[outdir,filesep,[OutnameLab,'_z',num2str(i),'.fig']])
        set(H,'PaperPositionMode','manual');
        set(H,'PaperUnits','inch')
        XSIZE = size(Doutshowtemp2,2);
        YSIZE = size(Doutshowtemp2,1);
        factor = 1:100;
        XSIZEnew = XSIZE*factor;
        YSIZEnew = YSIZE*factor;
        FACTORS1 = find(XSIZEnew>Hsize(3));
        FACTORS2 = find(YSIZEnew>Hsize(4));
        FACTORS = max(FACTORS1(1),FACTORS2(1));
        XSIZEU = XSIZEnew(FACTORS);
        YSIZEU = YSIZEnew(FACTORS);
        set(H,'Paperposition',[1,1,XSIZEU/300,YSIZEU/300]);
        print(H,[outdir,filesep,[OutnameLab,'_z',num2str(i),'.tif']],'-dtiff','-r300')
        close(H)
    case 2 % cornoal
        if showlab==1 % mni condition
            Showlabtemp = [zeros(length(sliceorder),1),sliceorder,zeros(length(sliceorder),1)];
            Showlab = mni2cor(Showlabtemp,vbg.mat);
            SliceorderCor = Showlab(:,2);
            SliceorderMNI = sliceorder;
        else
            SliceorderCor = sliceorder;
            Sliceordertemp = cor2mni(sliceorder,vbg.mat);
            SliceorderMNI = Sliceordertemp(:,2);
        end
        
        Doutshowtemp2 = zeros(size(dbgre,3)*NW,size(dbgre,1)*NL);
        for i = SliceorderCor
            Doutshowtemp = squeeze(dbgre(:,i,:));
            Doutshowtemp(Doutshowtemp<30*0.9) = 0;
            Doutshowtemp(Doutshowtemp>0&Doutshowtemp<=30*0.9) = 30+(dshownum)/(85-30);
            Doutshowtemp(Doutshowtemp>=85) = 85-dshownum/55;
            Doutshowtemp(Doutshowtemp>0) = (Doutshowtemp(Doutshowtemp>0)-30)/(85-30)*(dshownum*0.9)+(dshownum*0.05);
            Doutshowtemp = (Doutshowtemp-dshownum/2)*2;
            Doutmatout = squeeze(MATOUTT(:,i,:));
            ind_postemp = find(Doutmatout>0);
            ind_negtemp = find(Doutmatout<0);
            Doutshowtemp(ind_postemp) = Doutmatout(ind_postemp)+dshownum-out;
            Doutshowtemp(ind_negtemp) = Doutmatout(ind_negtemp)+out-dshownum;
            Doutshowtemp1 = rot90(Doutshowtemp);
            IX = i-SliceorderCor(1)+1;
            POSXT = POSX(IX);
            POSYT = POSY(IX);
            Doutshowtemp2(1+size(Doutshowtemp1,1)*(POSXT-1):size(Doutshowtemp1,1)*(POSXT),1+size(Doutshowtemp1,2)*(POSYT-1):size(Doutshowtemp1,2)*(POSYT)) = Doutshowtemp1;
        end
        Hexist1 = Hsize(3)-200;
        Hexist2 = Hsize(4)-200;
        Dsize = [size(Doutshowtemp2,2),size(Doutshowtemp2,1)];
        factorD = min(Dsize(1)/Hexist2,Dsize(2)/Hexist1);
        H = figure('pos',[100,100,Dsize(1)*factorD,Dsize(2)*factorD]);
        %             H = figure('pos',[100,100,size(Doutshowtemp1,2)*2,size(Doutshowtemp1,1)*2]);
        
        imagesc(Doutshowtemp2,[-dshownum*2,dshownum*2]);colormap(colormapshow);
        axis off;
        saveas(H,[outdir,filesep,[OutnameLab,'_y',num2str(i),'.fig']])
        set(H,'PaperPositionMode','manual');
        set(H,'PaperUnits','inch')
        XSIZE = size(Doutshowtemp2,2);
        YSIZE = size(Doutshowtemp2,1);
        factor = 1:100;
        XSIZEnew = XSIZE*factor;
        YSIZEnew = YSIZE*factor;
        FACTORS1 = find(XSIZEnew>Hsize(3));
        FACTORS2 = find(YSIZEnew>Hsize(4));
        FACTORS = max(FACTORS1(1),FACTORS2(1));
        XSIZEU = XSIZEnew(FACTORS);
        YSIZEU = YSIZEnew(FACTORS);
        set(H,'Paperposition',[1,1,XSIZEU/300,YSIZEU/300]);
        print(H,[outdir,filesep,[OutnameLab,'_y',num2str(i),'.tif']],'-dtiff','-r300')
        close(H)
        
    case 3 % sagittal
        if showlab==1 % mni condition
            Showlabtemp = [sliceorder,zeros(length(sliceorder),1),zeros(length(sliceorder),1)];
            Showlab = mni2cor(Showlabtemp,vbg.mat);
            SliceorderCor = Showlab(:,1);
            SliceorderMNI = sliceorder;
        else
            SliceorderCor = sliceorder;
            Sliceordertemp = cor2mni(sliceorder,vbg.mat);
            SliceorderMNI = Sliceordertemp(:,1);
        end
        
        Doutshowtemp2 = zeros(size(dbgre,3)*NW,size(dbgre,2)*NL);
        for i = SliceorderCor
            Doutshowtemp = squeeze(dbgre(i,:,:));
            Doutshowtemp(Doutshowtemp<30*0.9) = 0;
            Doutshowtemp(Doutshowtemp>0&Doutshowtemp<=30*0.9) = 30+(dshownum)/(85-30);
            Doutshowtemp(Doutshowtemp>=85) = 85-dshownum/55;
            Doutshowtemp(Doutshowtemp>0) = (Doutshowtemp(Doutshowtemp>0)-30)/(85-30)*(dshownum*0.9)+(dshownum*0.05);
            Doutshowtemp = (Doutshowtemp-dshownum/2)*2;
            Doutmatout = squeeze(MATOUTT(i,:,:));
            ind_postemp = find(Doutmatout>0);
            ind_negtemp = find(Doutmatout<0);
            Doutshowtemp(ind_postemp) = Doutmatout(ind_postemp)+dshownum-out;
            Doutshowtemp(ind_negtemp) = Doutmatout(ind_negtemp)+out-dshownum;
            Doutshowtemp1 = rot90(Doutshowtemp);
            IX = i-SliceorderCor(1)+1;
            POSXT = POSX(IX);
            POSYT = POSY(IX);
            Doutshowtemp2(1+size(Doutshowtemp1,1)*(POSXT-1):size(Doutshowtemp1,1)*(POSXT),1+size(Doutshowtemp1,2)*(POSYT-1):size(Doutshowtemp1,2)*(POSYT)) = Doutshowtemp1;
        end
        Hexist1 = Hsize(3)-200;
        Hexist2 = Hsize(4)-200;
        Dsize = [size(Doutshowtemp2,2),size(Doutshowtemp2,1)];
        factorD = min(Dsize(1)/Hexist2,Dsize(2)/Hexist1);
        H = figure('pos',[100,100,Dsize(1)*factorD,Dsize(2)*factorD]);
        %             H = figure('pos',[100,100,size(Doutshowtemp1,2)*2,size(Doutshowtemp1,1)*2]);
        
        imagesc(Doutshowtemp2,[-dshownum*2,dshownum*2]);colormap(colormapshow);
        axis off;
        saveas(H,[outdir,filesep,[OutnameLab,'_x',num2str(i),'.fig']])
        set(H,'PaperPositionMode','manual');
        set(H,'PaperUnits','inch')
        XSIZE = size(Doutshowtemp2,2);
        YSIZE = size(Doutshowtemp2,1);
        factor = 1:100;
        XSIZEnew = XSIZE*factor;
        YSIZEnew = YSIZE*factor;
        FACTORS1 = find(XSIZEnew>Hsize(3));
        FACTORS2 = find(YSIZEnew>Hsize(4));
        FACTORS = max(FACTORS1(1),FACTORS2(1));
        XSIZEU = XSIZEnew(FACTORS);
        YSIZEU = YSIZEnew(FACTORS);
        set(H,'Paperposition',[1,1,XSIZEU/300,YSIZEU/300]);
        print(H,[outdir,filesep,[OutnameLab,'_x',num2str(i),'.tif']],'-dtiff','-r300')
        close(H)        
    otherwise
        error('wrong output')
end