function varargout = slice_overlay(action, varargin);
% Function to display + manage slice display 
% Slice display works on a global structure SO
% with fields 
%  - img - array of images to display
%        - img structs contain fields
%             vol - vol struct info (see spm_vol)
%                   can also be vol containing image as 3d matrix
%                   set with slice_overlay('AddBlobs'...) call
%             cmap - colormap for this image
%             nancol - color for NaN. If scalar, this is an index into
%                    the image cmap.  If 1x3 vector, it's a colour
%             prop - proportion of intensity for this cmap/img 
%                    if = Inf, gives split cmap effect where values of 
%                    this cmap override previous image cmap values
%             func - function to apply to image before scaling to cmap
%                    (and therefore before min/max thresholding. E.g. a func of
%                    'i1(i1==0)=NaN' would convert zeros to NaNs
%             range - 2x1 vector of values for image to distribute colormap across
%                    the first row of the colormap applies to the first
%                    value in 'range', and the last value to the second
%                    value in 'range'
%             outofrange - behavior for image values to the left and
%                    right of image limits in 'range'.  Left means
%                    colormap values < 1, i.e for image values <
%                    range(1), if (range(1)<range(2)), and image values >
%                    range(1) where (range(1)>range(2)). If missing,
%                    display min (for Left) and max (for Right) value from colormap. 
%                    Otherwise should be a 2 element cell array, where
%                    the first element is the colour value for image values
%                    left of 'range', and the second is for image values
%                    right of 'range'.  Scalar values for
%                    colour index the colormap, 3x1 vectors are colour
%                    values.  An empty array attracts default settings
%                    appropriate to the mode - i.e. transparent colour (where 
%                    SO.prop ~= Inf), or split colour.  Empty cells
%                    default to 0. 0 specifies that voxels with this
%                    colour do not influence the image (split =
%                    background, true = black)
%            hold  - resampling order for image (see spm_sample_vol) -
%                    default 1
%            background - value when resampling outside image - default
%                    NaN
%            
% - transform - either - 4x4 transformation to apply to image slice position,
%             relative to mm given by slicedef, before display
%               or     - text string, one of axial, coronal, sagittal
%                        These orientations assume the image is currently
%                        (after its mat file has been applied) axially
%                        oriented
% - slicedef - 2x3 array specifying dimensions for slice images in mm
%             where rows are x,and y of slice image, and cols are neg max dim,
%             slice separation and pos max dim
% - slices   - vector of slice positions in mm in z (of transformed image)
% - figure    - figure handle for slice display figure
% - refreshf  - flag - if set or empty, refresh axis info for figure
%             else assume this is OK
% - clf       - flag, non zero -> clear figure before display.  Redundant
%               if refreshf == 0
% - area      struct with fields
%                  position - bottom left, x size y size 1x4 vector of
%                      area in which to display slices
%                  units    - one of
%                    inches,centimeters,normalized,points,{pixels}
%                  halign - one of left,{center},right
%                  valign - one of top,{middle},bottom
% - xslices  - no of slices to display across figure (defaults to an optimum)
% - cbar      - if empty, missing, no colourbar.  If an array of integers, then 
%             indexes img array, and makes colourbar for each cmap for
%             that img.  Cbars specified in order of appearance L->R
% - labels - struct can be absent (-> default numerical labels)
%                  empty (SO.labels = []) (no labels) or contain fields 
%                  colour - colour for label text 
%                  size - font size in units normalized to slice axes 
%                  format - if = cell array of strings =
%                  labels for each slice in Z.  If is string, specifies
%                  sprintf format string for labelling in distance of the
%                  origin (Xmm=0, Ymm=0) of each slice from plane containing
%                  the AC, in mm, in the space of the transformed image
% - callback - callback string for button down on image panels.  E.g.
%              setting SO.callback to 'slice_overlay(''getpos'')' prints to
%              the matlab window the equivalent position in mm of the
%              position of a mouse click on one of the image slices
% - printstr - string for printing slice overlay figure window, e.g.
%              'print -dpsc -painters -noui' (the default)
% - printfile - name of file to print output to; default 'slices.ps'
%
% FORMAT slice_overlay
% Checks, fills SO struct (slice_overlay('checkso')), and 
% displays slice overlay (slice_overlay('display'))
%
% FORMAT slice_overlay('checkso')
% Checks SO structure and sets defaults
%
% FORMAT cmap = slice_overlay('getcmap',cmapname)
% Gets colormap named in cmapname string
%
% FORMAT [mx mn] = slice_overlay('volmaxmin', vol)
% Returns maximum and minimum finite values from vol struct 'vol'
%
% FORMAT slice_overlay('addspm',SPM,VOL,dispf)
% Adds SPM blobs as new img to SO struct, split effect, 'hot' colormap, 
% Structures SPM and VOL are generated by calls to SPM results
% if not passed, they are fetched from the workspace
% If dispf is not passed, or nonzero, displays resulting SO figure also
%
% FORMAT slice_overlay('addblobs', imgno, XYZ, vals, mat)
% adds SPM blobs to img no 'imgno', as specified in 
% XYZ  - 3xN voxel coordinates of N blob values
% vals - N blob intensity values
% mat  - 4x4 matrix specifying voxels -> mm
%
% FORMAT vol = slice_overlay('blobs2vol', XYZ, vals, mat)
% returns (pseudo) vol struct for 3d blob volume specified
% in matrices as above
%
% FORMAT slice_overlay('addmatrix', imgno, mat3d, mat)
% adds 3d matrix image vol to img imgno.  Optionally
% mat  - 4x4 matrix specifying voxels -> mm
%
% FORMAT vol = slice_overlay('matrix2vol', mat3d, mat)
% returns (pseudo) vol struct for 3d matrix 
% input matrices as above
%
% FORMAT mmpos = slice_overlay('getpos')
% returns equivalent position in mm of last click on current axes (gca)
% if the axes contain an image slice (empty otherwise)
%
% FORMAT vals = slice_overlay('pointvals', XYZmm, holdlist)
% returns IxN matrix with values of each image 1..I, at each
% point 1..N specified in 3xN mm coordinate matrix XYZmm
% If specified, 'holdlist' contains I values giving hold
% values for resampling for each image (see spm_sample_vol)
%
% FORMAT slice_overlay('display')
% Displays slice overlay from SO struct
% 
% FORMAT slice_overlay('print', filename, printstr) 
% Prints slice overlay figure, usually to file.  If filename is not
% passed/empty gets filename from SO.printfile.  If printstr is not
% passed/empty gets printstr from SO.printstr
% 
% V 0.8 2/8/00  
% More or less  beta - take care.  Please report problems to  
% Matthew Brett matthew@mrc-cbu.cam.ac.uk

global SO

if nargin < 1
  checkso;
  action = 'display';
else
  action = lower(action);
end

switch action
 case 'checkso'
  checkso;
 case 'getcmap'
  varargout = {getcmap(varargin{1})};
 case 'volmaxmin'
  [mx mn] = volmaxmin(varargin{1});
  varargout = {mx, mn};
 case 'addspm'
  if nargin < 2
    varargin(1) = {evalin('base', 'SPM', ['error(''Cannot find SPM' ...
		    ' struct'')'])};
  end
  if nargin < 3
    varargin(2) = {evalin('base', 'VOL', ['error(''Cannot find VOL' ...
		    ' struct'')'])};
  end
  if nargin < 4
    varargin{3} = 1;
  end
  newimg = length(SO.img)+1;
  SO.img(newimg).vol = blobs2vol(varargin{1}.XYZ,varargin{1}.Z, varargin{2}.M);
  SO.img(newimg).prop = Inf;
  SO.img(newimg).cmap = hot;
  SO.img(newimg).range = [0 max(varargin{1}.Z)];
  SO.cbar = [SO.cbar newimg];
  if varargin{3}
    checkso;
    slice_overlay('display');
  end
  
 case 'addblobs'
  addblobs(varargin{1},varargin{2},varargin{3},varargin{4});
 case 'blobs2vol'
  varargout = {blobs2vol(varargin{1},varargin{2},varargin{3})};
 case 'addmatrix'
  if nargin<3,varargin{2}='';end
  if nargin<4,varargin{3}='';end
  addmatrix(varargin{1},varargin{2},varargin{3});
 case 'matrix2vol'
  if nargin<3,varargin{2}=[];end
  varargout = {matrix2vol(varargin{1},varargin{2})};
 case 'getpos'
  varargout = {getpos};
 case 'pointvals'
  varargout = {pointvals(varargin{1})};
 case 'print'
  if nargin<2,varargin{1}='';end
  if nargin<3,varargin{2}='';end
  printfig(varargin{1}, varargin{2});
 case 'display'

% get coordinates for plane
X=1;Y=2;Z=3;
dims = SO.slicedef;
xmm = dims(X,1):dims(X,2):dims(X,3);
ymm = dims(Y,1):dims(Y,2):dims(Y,3);
zmm = SO.slices;
[y x] = meshgrid(ymm,xmm');
vdims = [length(xmm),length(ymm),length(zmm)];

% no of slices, and panels (an extra for colorbars)
nslices = vdims(Z);
minnpanels = nslices;
cbars = 0;
if is_there(SO,'cbar')
  cbars = length(SO.cbar);
  minnpanels = minnpanels+cbars;
end

% get figure data
% if written to, the axes may be specified already
figno = figure(SO.figure);

% (re)initialize axes and stuff

% check if the figure is set up correctly
if ~SO.refreshf
  axisd = flipud(findobj(SO.figure, 'Type','axes','Tag', 'slice overlay panel'));
  npanels = length(axisd);
  if npanels < vdims(Z)+cbars;
    SO.refreshf = 1;
  end
end
if SO.refreshf
  % clear figure, axis store
  if SO.clf, clf; end
  axisd = [];

  % prevent print inversion problems
  set(figno,'InvertHardCopy','off');
  
  % calculate area of display in pixels
  parea = SO.area.position;
  if ~strcmp(SO.area.units, 'pixels')
    ubu = get(SO.figure, 'units');
    set(SO.figure, 'units','pixels');
    tmp = get(SO.figure, 'Position');
    ascf = tmp(3:4);
    if ~strcmp(SO.area.units, 'normalized')
      set(SO.figure, 'units',SO.area.units);
      tmp = get(SO.figure, 'Position');
      ascf = ascf ./ tmp(3:4);
    end
    set(figno, 'Units', ubu);
    parea = parea .* repmat(ascf, 1, 2);
  end
  asz = parea(3:4);
  
  % by default, make most parsimonious fit to figure
  yxratio = length(ymm)*dims(Y,2)/(length(xmm)*dims(X,2));
  if ~is_there(SO, 'xslices')
    % iteration needed to optimize, surprisingly.  Thanks to Ian NS
    axlen(X,:)=asz(1):-1:1;
    axlen(Y,:)=yxratio*axlen(X,:);
    panels = floor(asz'*ones(1,size(axlen,2))./axlen);
    estnpanels = prod(panels);
    tmp = find(estnpanels >= minnpanels);
    if isempty(tmp)
      error('Whoops, cannot fit panels onto figure');
    end
    b = tmp(1); % best fitting scaling
    panels = panels(:,b);
    axlen = axlen(:, b);
  else
    % if xslices is specified, assume X is flush with X figure dimensions
    panels([X:Y],1) = [SO.xslices; 0];
    axlen([X:Y],1) = [asz(X)/panels(X); 0];
  end
  
  % Axis dimensions are in pixels.  This prevents aspect ratio rescaling
  panels(Y) = ceil(minnpanels/panels(X));
  axlen(Y) = axlen(X)*yxratio;
  
  % centre (etc) panels in display area as required
  divs = [Inf 2 1];the_ds = [0;0];
  the_ds(X) = divs(strcmp(SO.area.halign, {'left','center','right'}));
  the_ds(Y) = divs(strcmp(SO.area.valign, {'bottom','middle','top'}));
  startc = parea(1:2)' + (asz'-(axlen.*panels))./the_ds;
  
  % make axes for panels
  r=0;c=1;
  npanels = prod(panels);
  lastempty = npanels-cbars;
  for i = 1:npanels
    % panel userdata
    if i<=nslices
      u.type = 'slice';
      u.no   = zmm(i);
    elseif i > lastempty
      u.type = 'cbar';
      u.no   = i - lastempty;
    else
      u.type = 'empty';
      u.no   = i - nslices;
    end
    axpos = [r*axlen(X)+startc(X) (panels(Y)-c)*axlen(Y)+startc(Y) axlen'];
    axisd(i) = axes(...
	'Parent',figno,...
	'XTick',[],...
	'XTickLabel',[],...
	'YTick',[],...
	'YTickLabel',[],...
	'Box','on',...
	'XLim',[1 vdims(X)],...
	'YLim',[1 vdims(Y)],...
	'Units', 'pixels',...
	'Position',axpos,...
	'Tag','slice overlay panel',...
	'UserData',u);
    r = r+1;
    if r >= panels(X)
      r = 0;
      c = c+1;
    end
  end
end

% sort out labels
if is_there(SO,'labels')
  labels = SO.labels;
  if iscell(labels.format)
    if length(labels.format)~=vdims(Z)
      error(...
	  sprintf('Oh dear, expecting %d labels, but found %d',...
		  vdims(Z), length(labels.contents)));
    end
  else
    % format string for mm from AC labelling
    fstr = labels.format;
    labels.format = cell(vdims(Z),1);
    acpt = SO.transform * [0 0 0 1]';
    for i = 1:vdims(Z)
      labels.format(i) = {sprintf(fstr,zmm(i)-acpt(Z))};
    end
  end
end

% modify colormaps with any new colours
nimgs = length(SO.img);
lrn = zeros(nimgs,3);
cmaps = cell(nimgs);
for i = 1:nimgs
  cmaps(i)={SO.img(i).cmap};
  lrnv = {SO.img(i).outofrange{:}, SO.img(i).nancol};
  for j = 1:length(lrnv)
    if prod(size(lrnv{j}))==1
      lrn(i,j) = lrnv{j};
    else
      cmaps(i) = {[cmaps{i}; lrnv{j}(1:3)]};
      lrn(i,j) = size(cmaps{i},1);
    end
  end
end

% cycle through slices displaying images
nvox = prod(vdims(1:2));
pandims = [vdims([2 1]) 3]; % NB XY transpose for display

zimg = zeros(pandims);
for i = 1:nslices
  ixyzmm = [x(:)';y(:)';ones(1,nvox)*zmm(i);ones(1,nvox)];
  img = zimg;
  for j = 1:nimgs
    thisimg = SO.img(j);
    % to voxel space of image
    vixyz = inv(SO.transform*thisimg.vol.mat)*ixyzmm;
    % raw data 
    if is_there(thisimg.vol, 'imgdata')
      V = thisimg.vol.imgdata;
    else
      V = thisimg.vol;
    end
    i1 = spm_sample_vol(V,vixyz(X,:),vixyz(Y,:),vixyz(Z,:), ...
			 [thisimg.hold thisimg.background]);
    if is_there(thisimg, 'func')
      eval(thisimg.func);
    end
    % transpose to reverse X and Y for figure
    i1 = reshape(i1, vdims(1:2))';
    % rescale to colormap
    [csdata badvals]= scaletocmap(...
	i1,...
	thisimg.range(1),...
	thisimg.range(2),...
	cmaps{j},...
	lrn(j,:));
    % take indices from colormap to make true colour image
    iimg = reshape(cmaps{j}(csdata(:),:),pandims);
    tmp = repmat(logical(~badvals),[1 1 3]);
    if thisimg.prop ~= Inf % truecolor overlay
      img(tmp) = img(tmp) + iimg(tmp)*thisimg.prop;
    else % split colormap effect
      img(tmp) = iimg(tmp);
    end
  end
  % threshold out of range values
  img(img>1) = 1;
  
  image('Parent', axisd(i),...
	'ButtonDownFcn', SO.callback,...
	'CData',img);
  if is_there(SO,'labels')
    text('Parent',axisd(i),...
	 'Color', labels.colour,...
	 'FontUnits', 'normalized',...
	 'VerticalAlignment','bottom',...
	 'HorizontalAlignment','left',...
	 'Position', [1 1],...
	 'FontSize',labels.size,...
	 'ButtonDownFcn', SO.callback,...
	 'String', labels.format{i});
  end
end
for i = (nslices+1):npanels
   set(axisd(i),'Color',[0 0 0]);
end
% add colorbar(s) 
for i = 1:cbars
  axno = axisd(end-cbars+i);
  cbari = SO.img(SO.cbar(i));
  cml = size(cbari.cmap,1);
  p = get(axno, 'Position');; % position of last axis
  cw = p(3)*0.2;
  ch = p(4)*0.75;
  pc = p(3:4)/2;
  [axlims idxs] = sort(cbari.range);
  a=axes(...
      'Parent',figno,...
      'XTick',[],...
      'XTickLabel',[],...
      'Units', 'pixels',...
      'YLim', axlims,...   
      'FontUnits', 'normalized',...
      'FontSize', 0.075,...
      'YColor',[1 1 1],...
      'Tag', 'cbar',...
      'Box', 'off',...
      'Position',[p(1)+pc(1)-cw/2,p(2)+pc(2)-ch/2,cw,ch]...
      );
  ih = image('Parent', a,...
	'YData', axlims(idxs),...     
	'CData', reshape(cbari.cmap,[cml,1,3]));

end

 otherwise
  error(sprintf('Unrecognized action string %s', action));

% end switch action
end

return
