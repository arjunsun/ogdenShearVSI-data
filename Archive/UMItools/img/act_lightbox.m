function [im, h]= act_lightbox(root1, root2, wscale1, wscale2, rows)
%function [im, h]= act_lightbox(root1, root2, wscale1, wscale 2, rows)
%
%   (c) 2005 Luis Hernandez-Garcia
%   University of Michigan
%   report bugs to:  hernan@umich.edu
%
% This program displays the slices in a 3D data set OR a time series
% root :  either the name of a file in a time series (could be asingle file)
% 	or could also be a 3D matrix that you want to display in slices
%	the program checks to see if it's a string or a matrix...
% wscale: window scale factor for display
% rows: number of rows of slices in the lightbox
%

doColorBars=1;

if isstr(root1)
	[im,h] = read_img(root1);
	if isfield(h,'magic')
		h=nii2avw_hdr(h);
	end
	if h.tdim==1
		im1 = reshape(im,h.xdim, h.ydim, h.zdim);
	else
		fprintf('\n\nThis is a time series.  Cannot lightbox it');
		return
	end
else
	im1 = root1;
	h.xdim=size(im1,1);
	h.ydim=size(im1,2);
	h.zdim=size(im1,3);

end


if isstr(root2)
	[im2,h] = read_img(root2);
	if isfield(h,'magic')
		h=nii2avw_hdr(h);
	end
	if h.tdim==1
		im2 = reshape(im2,h.xdim, h.ydim, h.zdim);
	else
		fprintf('\n\nThis is a time series.  Cannot lightbox it');
		return
	end
else
	im2 = root2;
	h.xdim=size(im2,1);
	h.ydim=size(im2,2);
	h.zdim=size(im2,3);

end

if nargin==2
	rows=[];
	wscale1=[];

	wscale2=[];
end

if isempty(rows)
	rows = floor(sqrt(h.zdim));
end


% scale the two images
im1 = (im1 - wscale1(1))*255 / (wscale1(2) - wscale1(1));
im1(im1(:)<0) = 0; im1(im1(:)>255) = 250;

% neg. activations
nim2 = -im2;

im2 = (im2 - wscale2(1))*255 / (wscale2(2) - wscale2(1));
im2(im2<0) = 0; im2(im2>255) = 250;

nim2 = (nim2 - wscale2(1))*255 / (wscale2(2) - wscale2(1));
nim2(nim2<0) = 0; nim2(nim2>255) = 250;


% select voxels above threshold
inds = find(im2(:) > abs(wscale2(1)));
	
% make the overlay:
D = im1;
if ~isempty(inds)
	D(inds) = 257 + im2(inds);
end

% now map the neg. map into the blue range
inds = find(nim2 > abs(wscale2(1)));
if ~isempty(inds)
	D(inds) = 513 + nim2(inds);
end

M1=[];
cols=round(h.zdim/rows);

for r=1:rows
	Mrow = [];

	for c=1:cols
		sl = c + cols*(r-1);
		if sl<=h.zdim
			Mrow = [Mrow  D(:,:,sl)'];
		else
			Mrow = [Mrow  zeros(h.ydim, h.xdim)];
		end

	end
	M1 = [M1 ; Mrow];
end

image(M1)
make_colormap

if doColorBars
    ncbar_labels=3;
    Ncolors = 256;

    c1 = colorbar('SouthOutside');
    set(c1,'XTick',[linspace(1,Ncolors-1,ncbar_labels)  ],...
        'XTickLabel',round([linspace(wscale1(1),wscale1(2), ncbar_labels) ]), ...
        'XLim', [0 Ncolors-1],...
        'FontSize',12);

    c1 = colorbar('WestOutside');
    set(c1,'YTick',[linspace(Ncolors+1,2*Ncolors,ncbar_labels) ],...
        'YTickLabel',round([linspace(wscale2(1), wscale2(2), ncbar_labels)]), ...
        'YLim', [Ncolors+1 2*Ncolors], ...
        'FontSize',12);

    c1 = colorbar('EastOutside');
    set(c1,'YTick',[linspace(2*Ncolors+1,3*Ncolors,ncbar_labels) ],...
        'YTickLabel',round([linspace(-wscale2(1), -wscale2(2), ncbar_labels)]), ...
        'YLim', [2*Ncolors+1 3*Ncolors], ...
        'FontSize',12);
end

if isstr(root1)
	set(gcf,'Name',root1)
end


return

%%
function make_colormap
mygray = [0:256]' * [1 1 1];

myhot = [128.5:0.5:255]' * [1 0 0] ;
tmp =   [128.5:0.5:255]' * [0 1 0] ;
tmp(:,1) = 256;
myhot = [myhot; tmp];
tmp =   [128.5:0.5:255]' * [0 0 1];
tmp(:,1:2) = 256;
myhot =  [myhot;  tmp;];
myhot = myhot(1:3:end, :);

myblue = myhot;
myblue(:,1) = myhot(:,3);
myblue(:,3) = myhot(:,1);

mymap = [mygray; myhot; myblue]/256;
colormap(mymap);

Ncolors = size(myhot,1);
% if isstr(root)
%       set(gcf, 'Name', root)
% end

axis image
grid off
axis off
axis xy
return
