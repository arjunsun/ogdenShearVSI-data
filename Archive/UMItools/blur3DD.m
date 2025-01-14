function [out, kfilter]=blur3DD(varargin)
%[out, kfilter]=blur(in [,kfilter,options]);
%in is a structure with fields image, kspace, pars as returned by varianms
%out is a structure with the same fields, but, the image has been blurred in 3 directions 
%an if optional kfilter matrix has the right size and is nonzero, then the
%input kfilter is used; otherwise kfilter is calculated. The idea is that 
%if you want to blur many images, you can reuse the kfilter calculated the
%first time around.

in=varargin{1};

if ~isstruct(in);
    display('abort: First input to blur3D.m must be a structure containing the data.')
    return
end
out=in;

if any(strcmp(varargin,'thk'))
    %a thickness has been specified, calculate the sigmaz required to reach
    ind=find(strcmp(varargin,'thk'));
    if nargin<ind+1 || ~isnumeric(varargin{ind+1});
        display('thk option requires a numeric input for the virtual slice thickness, in mm.');
        return
    end
    virtualslicethickness=varargin{ind+1};
    sigmaz=virtualslicethickness/in.pars.thk/2*0.8;
else
    sigmaz=0.8; %cross plane coupling, 0.8 means half the weight is on the pixel itself, the rest is on the side
end

ftest=exp(-(-10:10).^2/(2*sigmaz^2));
ftest=ftest/sum(ftest);
display(['z-weighting of parent slice =' num2str(round(ftest(11)*100)/100)]);


Min=in.image;
Mout=Min;
if any(strcmp(varargin,'abs'));
    Min=abs(Min);
end

if nargin>1 && isnumeric(varargin{2});
	kfilter=varargin{2};
	infilterflag=true; %default
	%1. check if the sizes of filter and inmatrix Min match
	if any(size(Min)~=size(kfilter));
		infilterflag=false;
		if isempty(kfilter);
			sigma=0.8;
		elseif isnumeric(kfilter);
			sigma=kfilter(1);
		end
	else
		%sizes do match, but also check if there are nonzero entries in the infilter
		if ~any(kfilter(:)~=0);
			infilterflag=false;
			sigma=0.8;
		end
	end

else
	infilterflag=false;
	sigma=0.8;
end


%if required, (re-)make the kspace filter
if ~infilterflag;
    siMin=size(Min);
    sfilter=zeros(siMin(1:3));
    [h w d]=size(sfilter);
	
    fw=ceil(6*sigma);
    x=[-fw:fw];
    y=x';
    z=x;
    fmask=zeros([numel(x),numel(y),numel(z)]);
    for jj=1:length(x);
        for ii=1:length(y);
            for kk=1:length(z);
                fmask(ii,jj,kk)=exp(-( (x(jj)^2+y(ii)^2) /(2*sigma^2) + z(kk)^2/(2*sigmaz^2)));
            end
        end
    end
    fmask=fmask/sum(fmask(:));
    
    %fmask is a patch that will go into the middle of sfilter;
    sfilter(h/2+1+y,w/2+1+x,d/2+1+z)=fmask;
    %kfilter=fftshift(fftn(fftshift(sfilter)));
    kfilter=ifftn(sfilter);
end

if ndims(Min)==4;
    nec=siMin(4);
else
    nec=1;
end

for ne=1:nec;
    FFT_Min=fftn(Min(:,:,:,ne));
    %Mout(:,:,:,ne)=fftshift(ifftn(fftshift(FFT_Min.*kfilter)));
    Mout(:,:,:,ne)=fftshift(ifftn(FFT_Min.*kfilter));
end

out.image=Mout;

