function AS_show_gainIO(Hshow)
inputdir = get(Hshow.IO_inputed,'string');
% [vol dol] = Dynamic_read_dir_NIFTI(inputdir);
Hsize = get(0,'screensize');

[vmid dmid] = Dynamic_read_dir_NIFTI(inputdir);
maxvshow = max(dmid);
minvshow = min(dmid);
% this is for spm stat
desipinfo = vmid.descrip;
inddes1 = find(desipinfo=='[');
inddes2 = find(desipinfo==']');

if ~isempty(inddes1)
    pat = Hshow.pat;
    labmark = desipinfo(5);
    styleofstat = labmark;
    dof = desipinfo(inddes1+1:inddes2-1);
    if strcmp(styleofstat,'T')
        dof1 = str2num(dof);
        dof2 = [];
    else
        doftemp = str2num(dof);
        dof1 = doftemp(1);
        dof2 = doftemp(2);
    end
    save(fullfile(pat,'InfoofStatForShow','DOFinfo.mat'),'dof1','dof2','styleofstat');
else
    H.fig = figure('pos',[Hsize(3)/2-100,Hsize(4)/2-50,200,100]);
    H.pat = Hshow.pat;
    H.h1 = uicontrol('parent',H.fig,'units','norm','pos',[0.05 0.5 0.2 0.3],'style','rad','string','T','val',1);
    H.h2 = uicontrol('parent',H.fig,'units','norm','pos',[0.275 0.5 0.2 0.3],'style','rad','string','F','val',0);
    H.h3 = uicontrol('parent',H.fig,'units','norm','pos',[0.525 0.5 0.2 0.3],'style','rad','string','R','val',0);
    H.h4 = uicontrol('parent',H.fig,'units','norm','pos',[0.75 0.5 0.2 0.3],'style','rad','string','P','val',0);
    H.ok = uicontrol('parent',H.fig,'units','norm','pos',[0.3 0.1 0.4 0.2],'style','pushbutton','string','OK');
    set(H.h1,'callback',{@SelT,H});
    set(H.h2,'callback',{@SelF,H});
    set(H.h3,'callback',{@SelR,H});
    set(H.h4,'callback',{@SelP,H});
    set(H.ok,'callback',{@OKpb,H});
end
%

end

function SelT(varargin)
H = varargin{3};
set(H.h1,'val',1);
set(H.h2,'val',0);
set(H.h3,'val',0);
set(H.h4,'val',0);
end
function SelF(varargin)
H = varargin{3};
set(H.h1,'val',0);
set(H.h2,'val',1);
set(H.h3,'val',0);
set(H.h4,'val',0);
end
function SelR(varargin)
H = varargin{3};
set(H.h1,'val',0);
set(H.h2,'val',0);
set(H.h3,'val',1);
set(H.h4,'val',0);
end
function SelP(varargin)
H = varargin{3};
set(H.h1,'val',0);
set(H.h2,'val',0);
set(H.h3,'val',0);
set(H.h4,'val',1);
end
function OKpb(varargin)
H = varargin{3};
VALSHOW = [get(H.h1,'val'),get(H.h2,'val'),get(H.h3,'val'),get(H.h4,'val')];
Hsize = get(0,'screensize');
pat = H.pat;
if isempty(dir([pat,filesep,'InfoofStatForShow']));
    mkdir([pat,filesep,'InfoofStatForShow']);
else
    rmdir([pat,filesep,'InfoofStatForShow'],'s');
    mkdir([pat,filesep,'InfoofStatForShow']);
end
if VALSHOW(1)
    HinfoT.fig = figure('pos',[Hsize(3)/2-300,Hsize(4)/2-150,600,300]);
    HinfoT.T1text = uicontrol('parent',HinfoT.fig,'units','norm','pos',[0.025 0.65 0.3 0.3],'style','text','string','Number of Group 1:');
    HinfoT.T2text = uicontrol('parent',HinfoT.fig,'units','norm','pos',[0.35 0.65 0.3 0.3],'style','text','string','Number of Group 2:');
    HinfoT.T1edit = uicontrol('parent',HinfoT.fig,'units','norm','pos',[0.025 0.3 0.3 0.3],'style','edit');
    HinfoT.T2edit = uicontrol('parent',HinfoT.fig,'units','norm','pos',[0.35 0.3 0.3 0.3],'style','edit');
    HinfoT.COVnumtext = uicontrol('parent',HinfoT.fig,'units','norm','pos',[0.675 0.65 0.3 0.3],'style','text','string','COV number:');
    HinfoT.COVnumedit = uicontrol('parent',HinfoT.fig,'units','norm','pos',[0.675 0.3 0.3 0.3],'style','edit','string','0');
    HinfoT.OK = uicontrol('parent',HinfoT.fig,'units','norm','pos',[0.4 0.1 0.2 0.1],'style','pushbutton','string','OK');
    set(HinfoT.OK,'callback',{@Tinfocol,HinfoT,H});
elseif VALSHOW(2)
    HinfoF.fig = figure('pos',[Hsize(3)/2-300,Hsize(4)/2-150,600,300]);
%     HinfoF.T1text = uicontrol('parent',HinfoF.fig,'units','norm','pos',[0.025 0.65 0.3 0.3],'style','text','string','Design of ANOVA:');
    HinfoF.T1text = uicontrol('parent',HinfoF.fig,'units','norm','pos',[0.025 0.65 0.3 0.3],'style','text','string','Group Number(only for 1-Way ANOVA):');
    HinfoF.T2text = uicontrol('parent',HinfoF.fig,'units','norm','pos',[0.35 0.65 0.3 0.3],'style','text','string','Number of Each Group:');
    HinfoF.T1edit = uicontrol('parent',HinfoF.fig,'units','norm','pos',[0.025 0.3 0.3 0.3],'style','edit');
    HinfoF.T2edit = uicontrol('parent',HinfoF.fig,'units','norm','pos',[0.35 0.3 0.3 0.3],'style','edit');
    HinfoF.COVnumtext = uicontrol('parent',HinfoF.fig,'units','norm','pos',[0.675 0.65 0.3 0.3],'style','text','string','COV number:');
    HinfoF.COVnumedit = uicontrol('parent',HinfoF.fig,'units','norm','pos',[0.675 0.3 0.3 0.3],'style','edit','string','0');
    HinfoF.OK = uicontrol('parent',HinfoF.fig,'units','norm','pos',[0.4 0.1 0.2 0.1],'style','pushbutton','string','OK');
    set(HinfoF.T1edit,'callback',{@Fedit,HinfoF,H});
    set(HinfoF.OK,'callback',{@Finfocol,HinfoF,H});
elseif VALSHOW(3)
    HinfoR.fig = figure('pos',[Hsize(3)/2-300,Hsize(4)/2-150,600,300]);
    HinfoR.T1text = uicontrol('parent',HinfoR.fig,'units','norm','pos',[0.1 0.65 0.3 0.3],'style','text','string','Particant Number:');
%     Hinfo.T2text = uicontrol('parent',Hinfo.fig,'units','norm','pos',[],'style','text','string','Number of Group 2:');
    HinfoR.T1edit = uicontrol('parent',HinfoR.fig,'units','norm','pos',[0.1 0.3 0.3 0.3],'style','edit');
%     Hinfo.T2edit = uicontrol('parent',Hinfo.fig,'units','norm','pos',[],'style','edit');
    HinfoR.COVnumtext = uicontrol('parent',HinfoR.fig,'units','norm','pos',[0.6 0.65 0.3 0.3],'style','text','string','COV number:');
    HinfoR.COVnumedit = uicontrol('parent',HinfoR.fig,'units','norm','pos',[0.6 0.3 0.3 0.3],'style','edit','string','0');
    HinfoR.OK = uicontrol('parent',HinfoR.fig,'units','norm','pos',[0.4 0.1 0.2 0.1],'style','pushbutton','string','OK');
    set(HinfoR.OK,'callback',{@Rinfocol,HinfoR,H});
elseif VALSHOW(4)
%     pat = H.pat;    
    dof1 = 10000;
    dof2 = [];
    styleofstat = 'P';
    save(fullfile(pat,'InfoofStatForShow','DOFinfo.mat'),'dof1','dof2','styleofstat');
end
close(H.fig)
end
function Fedit(varargin)
HinfoF = varargin{3};
H = varargin{4};
pat = H.pat;
info1 = get(HinfoF.T1edit,'string');
Groupnum = str2num(info1);
for i = 1:Groupnum
    Ansinfo = inputdlg(['GROUP',num2str(i),':'],['GROUP',num2str(i),':'],1);
    Groupnums(1,i) = str2num(Ansinfo{1});
end
% set(HinfoF.T2edit,'string',mat2str(Groupnums));
set(HinfoF.T2edit,'string',num2str(Groupnums));
end
function Tinfocol(varargin)
HinfoT = varargin{3};
H = varargin{4};
pat = H.pat;
if isempty(dir([pat,filesep,'InfoofStatForShow']));
    mkdir([pat,filesep,'InfoofStatForShow']);
else
    rmdir([pat,filesep,'InfoofStatForShow'],'s');
    mkdir([pat,filesep,'InfoofStatForShow']);
end
T1info = get(HinfoT.T1edit,'string');
T2info = get(HinfoT.T2edit,'string');
G1num = str2num(T1info);
G2num = str2num(T2info);
COVinfo = get(HinfoT.COVnumedit,'string');
COVnum = str2num(COVinfo);
dof1 = G1num+G2num-1-COVnum;
dof2 = [];
styleofstat = 'T';
save(fullfile(pat,'InfoofStatForShow','DOFinfo.mat'),'dof1','dof2','styleofstat');
close(HinfoT.fig);
end
function Finfocol(varargin)
HinfoF = varargin{3};
H = varargin{4};
pat = H.pat;
if isempty(dir([pat,filesep,'InfoofStatForShow']));
    mkdir([pat,filesep,'InfoofStatForShow']);
else
    rmdir([pat,filesep,'InfoofStatForShow'],'s');
    mkdir([pat,filesep,'InfoofStatForShow']);
end

T1info = get(HinfoF.T1edit,'string');
T2info = get(HinfoF.T2edit,'string');
G1num = str2num(T1info);
G2num = str2num(T2info);
COVinfo = get(HinfoF.COVnumedit,'string');
COVnum = str2num(COVinfo);
dof1 = G1num-1;
dof2 = sum(G2num)-G1num-COVnum;
styleofstat = 'F';
save(fullfile(pat,'InfoofStatForShow','DOFinfo.mat'),'dof1','dof2','styleofstat');
close(HinfoF.fig);
end
function Rinfocol(varargin)
HinfoR = varargin{3};
H = varargin{4};
pat = H.pat;
if isempty(dir([pat,filesep,'InfoofStatForShow']));
    mkdir([pat,filesep,'InfoofStatForShow']);
else
    rmdir([pat,filesep,'InfoofStatForShow'],'s');
    mkdir([pat,filesep,'InfoofStatForShow']);
end
T1info = get(HinfoR.T1edit,'string');
G1num = str2num(T1info);

COVinfo = get(HinfoR.COVnumedit,'string');
COVnum = str2num(COVinfo);
dof1 = G1num-2-COVnum;
dof2 = [];
styleofstat = 'R';
save(fullfile(pat,'InfoofStatForShow','DOFinfo.mat'),'dof1','dof2','styleofstat');
close(HinfoR.fig);
end