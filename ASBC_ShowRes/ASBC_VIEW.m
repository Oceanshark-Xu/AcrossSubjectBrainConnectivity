function ASBC_VIEW
Hsize = get(0,'screensize');
msize = min(Hsize(3:4))*0.8;
Hasview.fig = figure('pos',[Hsize(3)/2-msize/4,Hsize(4)/2-msize/8,msize/2,msize/4],'name','ASBC Viewer V1.0');
Hasview.FC = uicontrol('parent',Hasview.fig,'units','norm','pos',[0.025 0.6 0.3 0.3],'style','pushbutton','string','Map(R T F P)');
Hasview.GC = uicontrol('parent',Hasview.fig,'units','norm','pos',[0.025 0.1 0.3 0.3],'style','pushbutton','string','Map(GCA)');
Hasview.MatrixFC = uicontrol('parent',Hasview.fig,'units','norm','pos',[0.35 0.6 0.3 0.3],'style','pushbutton','string','Matrix(R T F P)');
Hasview.MatrixGC = uicontrol('parent',Hasview.fig,'units','norm','pos',[0.35 0.1 0.3 0.3],'style','pushbutton','string','Matrix(GCA)');
Hasview.WTA = uicontrol('parent',Hasview.fig,'units','norm','pos',[0.675 0.6 0.3 0.3],'style','pushbutton','string','Winner-Take-All');
Hasview.Return = uicontrol('parent',Hasview.fig,'units','norm','pos',[0.675 0.1 0.3 0.3],'style','pushbutton','string','Return');
set(Hasview.FC,'callback',{@FCmapshow,Hasview})
set(Hasview.GC,'callback',{@GCmapshow,Hasview})
set(Hasview.MatrixFC,'callback',{@FCmatrixshow,Hasview})
set(Hasview.MatrixGC,'callback',{@GCmatrixshow,Hasview})
set(Hasview.WTA,'callback',{@WTAshow,Hasview})
set(Hasview.Return,'callback',{@Return,Hasview})
end
function FCmapshow(varargin)
Hasview = varargin{3};
close(Hasview.fig);
AS_ShowFCGUI;
end
function GCmapshow(varargin)
Hasview = varargin{3};
% uiwait(msgbox('comming soon'));
close(Hasview.fig);
answ1 = questdlg('which type of GC?','which type of GC?','Residual','Coefficient','Residual');
if strcmp(answ1,'Residual')
    AS_ShowGCGUI_res;
else
    AS_ShowGCGUI_coef;
end
end
function FCmatrixshow(varargin)
Hasview = varargin{3};
close(Hasview.fig);
AS_ShowFCmatrixGUI;
end
function GCmatrixshow(varargin)
Hasview = varargin{3};
% uiwait(msgbox('comming soon'));
close(Hasview.fig);
AS_ShowGCmatrixGUI;
end
function WTAshow(varargin)
Hasview = varargin{3};
close(Hasview.fig);
AS_ShowWTAGUI;
end
function Return(varargin)
Hasview = varargin{3};
close(Hasview.fig);
ASBCmain;
end