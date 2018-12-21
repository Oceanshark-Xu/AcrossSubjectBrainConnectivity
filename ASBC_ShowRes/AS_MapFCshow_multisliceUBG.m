function AS_MapFCshow_multisliceUBG(pat0,showlabs,sliceorder,outdir,OutnameLab,colormapshow,climsuse,printlab)
[pat,nam.ext] = fileparts(which('AS_MapFCshow.m'));
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
% POSXY = [POSX,POSY];
NL = length(unique(POSX));
load([pat0,filesep,'TempForOrigShow',filesep,showlabs,num2str(1),'.mat']);

Doutshowtemp2 = zeros(size(Slice,2)*NW,size(Slice,1)*NL);
for i = 1:length(sliceorder)    
    load([pat0,filesep,'TempForOrigShow',filesep,showlabs,num2str(sliceorder(i)),'.mat']);
    Doutshowtemp1 = rot90(Slice);
    IX = i;
    POSXT = POSX(IX);
    POSYT = POSY(IX);
    Doutshowtemp2(1+size(Doutshowtemp1,1)*(POSXT-1):size(Doutshowtemp1,1)*(POSXT),1+size(Doutshowtemp1,2)*(POSYT-1):size(Doutshowtemp1,2)*(POSYT)) = Doutshowtemp1;
end
Hsize = get(0,'screensize');
Hexist1 = Hsize(3)-200;
Hexist2 = Hsize(4)-200;
Dsize = [size(Doutshowtemp2,2),size(Doutshowtemp2,1)];
factorD = min(Dsize(1)/Hexist2,Dsize(2)/Hexist1);
H = figure('pos',[100,100,Dsize(1)*factorD,Dsize(2)*factorD]);
%             H = figure('pos',[100,100,size(Doutshowtemp1,2)*2,size(Doutshowtemp1,1)*2]);

imagesc(Doutshowtemp2,climsuse);colormap(colormapshow);
axis off;
if printlab
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
end
end
