function [out, kfilter]=blur3d(varargin)
%[out, kfilter]=blur(in [,kfilter,options]);
%in is a structure with fields image, kspace, pars as returned by varianms
%out is a structure with the same fields, but, the image has been blurred in 3 directions
%an if optional kfilter matrix has the right size and is nonzero, then the
%input kfilter is used; otherwise kfilter is calculated. The idea is that
%if you want to blur many images, you can reuse the kfilter calculated the
%first time around.
%options:
%'vox', [sz_ro sz_pe sz_pe2]   three component size vector: 
%                              in mm or in voxels, depending on structure or matrix input
%'abs'                         blur absolute values, not complex numbers

in=varargin{1};

if ~isstruct(in);
    blurtype='voxels';  %optional voxel blurring units: voxels
    if isnumeric(in);
        Min=in;
    else
        display('abort: First input to blur3D.m must be a structure containing the data, or a 3darray.')
    end
else
    blurtype='mm'; %optional vovel blurring units: mm
    Min=in.image;
    out=in;
    display('blur3d got a structure as input, the blurring will be by mm');
end

%default blur coupling, applies if no vox is defined
sigmavec=0.7*[1 1 1];

if any(strcmp(varargin,'vox'));
    %virtual voxel size is given in mm
    ind=find(strcmp(varargin,'vox'));
    %errortrapping
    if nargin<ind+1 || ~isnumeric(varargin{ind+1});
        display('vox option requires a numeric vector input for the virtual voxel size.');
        return
    end
    
    %%%
    if strcmp(blurtype,'voxels'); %the input was a 3d array, not a varianms structure
        voxvec=varargin{ind+1};
        if numel(voxvec)==3;             %a voxel size has been handed down
            %sz is the voxel size handed down;
            sz=ones(1,3);
            sigmavec=voxvec./sz*0.8/1.7;
        end
    end
    
    %%%
    if strcmp(blurtype,'mm'); %the input was varianms structure, virtual voxel size in mm
        si=size(Min);
        voxvec=varargin{ind+1};
        if numel(voxvec)==3;             %a voxel size has been handed down
            %sz is the voxel size of the varianms data
            sz(1)=in.pars.lro*10/si(1);  
            sz(2)=in.pars.lpe*10/si(2);
            sz(3)=in.pars.lpe2*10/si(3);
            sigmavec=voxvec./sz*0.8/1.7;
        end
    end
end 

if any(strcmp(varargin,'abs'));
    Min=abs(Min);
end

Mout=Min;

simage=size(Min);
Min=reshape(Min, [simage(1:3) prod(simage(4:end))]);

%make the kspace filter
siMin=size(Min);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%pad the image cube in the slice direction, repeating the first
%and last slice twice; eliminate wrap around issues
siMin_padded=siMin;
siMin_padded(3)=siMin_padded(3)+4;
sfilter=zeros(siMin_padded(1:3));
[h w d]=size(sfilter);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fw=ceil(3*sigmavec);
d1=-fw(1):fw(1);
d2=-fw(2):fw(2);
d3=-fw(3):fw(3);

fmask=zeros([numel(d1),numel(d2),numel(d3)]);
for ii=1:length(d1);
    for jj=1:length(d2);
        for kk=1:length(d3);
            fmask(ii,jj,kk)=exp(-d1(ii)^2/(2*sigmavec(1)^2) -d2(jj)^2/(2*sigmavec(2)^2) -d3(kk)^2/(2*sigmavec(3)^2));
        end
    end
end
fmask=fmask/sum(fmask(:));

%fmask is a patch that will go into the middle of sfilter;
sfilter(h/2+1+d1,w/2+1+d2,d/2+1+d3)=fmask;
kfilter=conj(fftn(sfilter));


if ndims(Min)==4;
    nec=siMin(4);
else
    nec=1;
end

Min_padded=cat(3,Min(:,:,1,:),Min(:,:,1,:),Min,Min(:,:,end,:),Min(:,:,end,:));
Mout_padded=zeros(size(Min_padded));

for ne=1:nec;
    FFT_Min_padded=fftn(Min_padded(:,:,:,ne));
    Mout_padded(:,:,:,ne)=fftshift(ifftn(FFT_Min_padded.*kfilter));
    Mout(:,:,:,ne)=Mout_padded(:,:,3:(end-2),ne);
end


if isnumeric(varargin{1});
    out=Mout;
else
    
    out.image=Mout;
end

