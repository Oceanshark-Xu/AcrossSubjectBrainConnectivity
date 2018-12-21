function aswtashowinfo(has)
ASINFO.fig = figure('units','norm','pos',[0.3 0.4 0.4 0.2]);
ASINFO.input = uicontrol('parent',ASINFO.fig,'units','norm','pos',[0.1 0.7 0.1 0.15],'style','text','string','input');
ASINFO.inputed = uicontrol('parent',ASINFO.fig,'units','norm','pos',[0.25 0.7 0.5 0.15],'style','edit','string','null');
ASINFO.inputsel = uicontrol('parent',ASINFO.fig,'units','norm','pos',[0.8 0.7 0.1 0.15],'style','pushbutton','string','...');

ASINFO.output = uicontrol('parent',ASINFO.fig,'units','norm','pos',[0.1 0.5 0.1 0.15],'style','text','string','output');
ASINFO.outputed = uicontrol('parent',ASINFO.fig,'units','norm','pos',[0.25 0.5 0.5 0.15],'style','edit','string','null');
ASINFO.outputsel = uicontrol('parent',ASINFO.fig,'units','norm','pos',[0.8 0.5 0.1 0.15],'style','pushbutton','string','...');

ASINFO.target = uicontrol('parent',ASINFO.fig,'units','norm','pos',[0.1 0.3 0.1 0.15],'style','text','string','targetROI');
ASINFO.targeted = uicontrol('parent',ASINFO.fig,'units','norm','pos',[0.25 0.3 0.5 0.15],'style','edit','string','null');
ASINFO.targetsel = uicontrol('parent',ASINFO.fig,'units','norm','pos',[0.8 0.3 0.1 0.15],'style','pushbutton','string','...');

ASINFO.seedN = uicontrol('parent',ASINFO.fig,'units','norm','pos',[0.1 0.1 0.1 0.15],'style','text','string','seed number');
ASINFO.seedNed = uicontrol('parent',ASINFO.fig,'units','norm','pos',[0.25 0.1 0.1 0.15],'style','edit','string','5');
ASINFO.targetN = uicontrol('parent',ASINFO.fig,'units','norm','pos',[0.4 0.1 0.1 0.15],'style','text','string','Target number');
ASINFO.targetNed = uicontrol('parent',ASINFO.fig,'units','norm','pos',[0.55 0.1 0.1 0.15],'style','edit','string','1');
ASINFO.OK = uicontrol('parent',ASINFO.fig,'units','norm','pos',[0.7 0.1 0.2 0.15],'style','pushbutton','string','OK');

set(ASINFO.inputsel,'callback',{@inputsel,ASINFO});
set(ASINFO.outputsel,'callback',{@outputsel,ASINFO});
set(ASINFO.targetsel,'callback',{@targetsel,ASINFO});
set(ASINFO.OK,'callback',{@SetOK,ASINFO,has});

end
function inputsel(varargin)
ASINFO = varargin{3};
path = uigetdir(pwd,'input directory');
set(ASINFO.inputed,'string',path);
end
function outputsel(varargin)
ASINFO = varargin{3};
path = uigetdir(pwd,'out directory');
set(ASINFO.outputed,'string',path);
end
function targetsel(varargin)
ASINFO = varargin{3};
[name,path ext] = uigetfile({'*.nii';'*.img';'*.*'},'target ROIs');
set(ASINFO.targeted,'string',fullfile(path,name));
end
function SetOK(varargin)
ASINFO = varargin{3};
has = varargin{4};
Info.input = get(ASINFO.inputed,'string');
Info.output = get(ASINFO.outputed,'string');
Info.target = get(ASINFO.targeted,'string');
targetnum = get(ASINFO.targetNed,'string');
seednum = get(ASINFO.seedNed,'string');
Info.targetnum = str2num(targetnum);
Info.seednum = str2num(seednum);
mkdir([Info.output,filesep,'Infoforshow']);
save([Info.output,filesep,'Infoforshow',filesep,'Info.mat'],'Info');
set(has.infocoled,'string',['Info saved in:',Info.output,filesep,'Infoforshow',filesep,'Info.mat']);
close(ASINFO.fig);
clear ASINFO
end
