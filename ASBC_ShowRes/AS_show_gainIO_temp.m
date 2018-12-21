function AS_show_gainIO_temp(Hshow)
Filepath = get(Hshow.IO_tempputed,'string');
pat = Hshow.pat;
[Vbg,Dbg] = Dynamic_read_dir_NIFTI(Filepath);
% [Vbg,Dbg] = Dynamic_read_dir_NIFTI(templatepath);
DbgRe = reshape(Dbg,Vbg.dim(1),Vbg.dim(2),Vbg.dim(3));
mnicoord = [0 0 0];
corcoord = mni2cor(mnicoord,Vbg.mat);
maxv = max(Dbg);
minv = min(Dbg);
Htempchose.fig = figure('units','norm','pos',[0.2 0.3 0.6 0.4]);
Htempchose.ax1 = axes('parent',Htempchose.fig,'units','norm','pos',[0.05 0.2 0.3 0.7]);
Htempchose.ax2 = axes('parent',Htempchose.fig,'units','norm','pos',[0.35 0.2 0.3 0.7]);
Htempchose.ax3 = axes('parent',Htempchose.fig,'units','norm','pos',[0.65 0.2 0.3 0.7]);
Htempchose.Low = uicontrol('parent',Htempchose.fig,'units','norm','pos',[0.1 0.1 0.3 0.08],'style','slider','min',minv,'max',maxv/2,'sliderstep',[0.01 0.1],'value',minv);
Htempchose.High = uicontrol('parent',Htempchose.fig,'units','norm','pos',[0.45 0.1 0.3 0.08],'style','slider','min',maxv/2,'max',maxv,'sliderstep',[0.01 0.1],'value',maxv);
Htempchose.ok = uicontrol('parent',Htempchose.fig,'units','norm','pos',[0.8 0.1 0.1 0.08],'style','pushbutton','string','OK');
AxialDATnew = rot90(squeeze(DbgRe(:,:,corcoord(3))));
CornoralDATnew = rot90(squeeze(DbgRe(:,corcoord(2),:)));
SagittalDATnew = rot90(squeeze(DbgRe(corcoord(1),:,:)));
image(AxialDATnew,'parent',Htempchose.ax3,'CDataMapping','scaled');axis(Htempchose.ax3,'off');colormap(Htempchose.ax3,'gray')
image(CornoralDATnew,'parent',Htempchose.ax2,'CDataMapping','scaled');axis(Htempchose.ax2,'off');colormap(Htempchose.ax2,'gray')
image(SagittalDATnew,'parent',Htempchose.ax1,'CDataMapping','scaled');axis(Htempchose.ax1,'off');colormap(Htempchose.ax1,'gray')
set(Htempchose.ax3,'Clim',[minv maxv]);
set(Htempchose.ax2,'Clim',[minv maxv]);
set(Htempchose.ax1,'Clim',[minv maxv]);
Htempchose.pat = pat;

set(Htempchose.ok,'callback',{@OKpb,Htempchose,Vbg,DbgRe})
set(Htempchose.Low,'callback',{@minvsel,Htempchose,AxialDATnew,CornoralDATnew,SagittalDATnew})
set(Htempchose.High,'callback',{@maxvsel,Htempchose,AxialDATnew,CornoralDATnew,SagittalDATnew})
end

function minvsel(varargin)
Htempchose = varargin{3};
AxialDATnew = varargin{4};
CornoralDATnew = varargin{5};
SagittalDATnew = varargin{6};
vallow = get(Htempchose.Low,'val');
valhigh = get(Htempchose.High,'val');
AxialDATnew(AxialDATnew<vallow*0.9) = 0;
AxialDATnew(AxialDATnew>=valhigh)= valhigh;
CornoralDATnew(CornoralDATnew<vallow*0.9) = 0;
CornoralDATnew(CornoralDATnew>=valhigh)= valhigh;
SagittalDATnew(SagittalDATnew<vallow*0.9) = 0;
SagittalDATnew(SagittalDATnew>=valhigh)= valhigh;
image(AxialDATnew,'parent',Htempchose.ax3,'CDataMapping','scaled');axis(Htempchose.ax3,'off');colormap(Htempchose.ax3,'gray')
image(CornoralDATnew,'parent',Htempchose.ax2,'CDataMapping','scaled');axis(Htempchose.ax2,'off');colormap(Htempchose.ax2,'gray')
image(SagittalDATnew,'parent',Htempchose.ax1,'CDataMapping','scaled');axis(Htempchose.ax1,'off');colormap(Htempchose.ax1,'gray')
set(Htempchose.ax3,'Clim',[vallow valhigh]);
set(Htempchose.ax2,'Clim',[vallow valhigh]);
set(Htempchose.ax1,'Clim',[vallow valhigh]);
end
function maxvsel(varargin)
Htempchose = varargin{3};
AxialDATnew = varargin{4};
CornoralDATnew = varargin{5};
SagittalDATnew = varargin{6};
vallow = get(Htempchose.Low,'val');
valhigh = get(Htempchose.High,'val');
AxialDATnew(AxialDATnew<vallow*0.9) = 0;
AxialDATnew(AxialDATnew>=valhigh)= valhigh;
CornoralDATnew(CornoralDATnew<vallow*0.9) = 0;
CornoralDATnew(CornoralDATnew>=valhigh)= valhigh;
SagittalDATnew(SagittalDATnew<vallow*0.9) = 0;
SagittalDATnew(SagittalDATnew>=valhigh)= valhigh;
image(AxialDATnew,'parent',Htempchose.ax3,'CDataMapping','scaled');axis(Htempchose.ax3,'off');colormap(Htempchose.ax3,'gray')
image(CornoralDATnew,'parent',Htempchose.ax2,'CDataMapping','scaled');axis(Htempchose.ax2,'off');colormap(Htempchose.ax2,'gray')
image(SagittalDATnew,'parent',Htempchose.ax1,'CDataMapping','scaled');axis(Htempchose.ax1,'off');colormap(Htempchose.ax1,'gray')
set(Htempchose.ax3,'Clim',[vallow valhigh]);
set(Htempchose.ax2,'Clim',[vallow valhigh]);
set(Htempchose.ax1,'Clim',[vallow valhigh]);
end

function OKpb(varargin)
Htempchose = varargin{3};
Vbg = varargin{4};
DbgRe = varargin{5};
vallow = get(Htempchose.Low,'val');
valhigh = get(Htempchose.High,'val');
pat = Htempchose.pat;
%%
DbgRe(DbgRe<vallow*0.9) = 0;
DbgRe(DbgRe>=valhigh)= valhigh;
if isempty(dir([pat,filesep,'TempForOrigShow']))
    mkdir([pat,filesep,'TempForOrigShow']);
else
    rmdir([pat,filesep,'TempForOrigShow'],'s');
    mkdir([pat,filesep,'TempForOrigShow']);
end

for i = 1:size(DbgRe,1)
    Slice = squeeze(DbgRe(i,:,:));
    save([pat,filesep,'TempForOrigShow',filesep,'X',num2str(i),'.mat'],'Slice');
end
for i = 1:size(DbgRe,2)
    Slice = squeeze(DbgRe(:,i,:));
    save([pat,filesep,'TempForOrigShow',filesep,'Y',num2str(i),'.mat'],'Slice');
end
for i = 1:size(DbgRe,3)
    Slice = squeeze(DbgRe(:,:,i));
    save([pat,filesep,'TempForOrigShow',filesep,'Z',num2str(i),'.mat'],'Slice');
end
clear Slice;
save([pat,filesep,'TempForOrigShow',filesep,'Vbg.mat'],'Vbg')
close(Htempchose.fig)
end