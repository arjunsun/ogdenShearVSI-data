function R2starout=R2star(varargin)
%out=R2star(in1[,in2,in3,..]))
% in1, in2 etc are structures returned by in=varianms(<filename.fid>)
% from dwgesege experiments, 3 or 4 echoes, GE(TE0) SE(TE) GE(TE00) GE(2*TE00)
display('T2star requires data to be hardwarephase.m corrected and softshimmed.m');
display('For noise reduction, consider blurring the data in 3D, with blur3D.m');

%identify which of the inputs are structures
datfileflag=false(nargin,1);
for jj=1:nargin;
    datfileflag(jj)=isstruct(varargin{jj});
end
numdatafiles=sum(+datfileflag);     %number of data files
di=find(datfileflag);               %datfileindex

[nro,npe,nsl,nec]=size(varargin{1}.image);

%stick all the data into a common data block
MM = zeros(nro,npe,nsl,nec,numdatafiles);
for nd=1:numdatafiles;
    findex=di(nd);
    MM(:,:,:,:,nd)=varargin{findex}.image;
    TR(nd)=varargin{findex}.pars.tr;
    TE(nd)=varargin{findex}.pars.te;
    TE2(nd)=varargin{findex}.pars.te2;
end

%check if all the data have the same TE
if std(TE)~=0; 
    display('Use data with same TE!');
    display(['TE=' num2str(TE)]);
    return;
end

%check if all the data have the same TE, and TE2
if std(TE2)~=0; 
    display('Use data with same TE2!');
    display(['TE2=' num2str(TE)]);
    return;
end

%sum echoes over all the data files, regardless of TR and Diffusion
aveE = squeeze(sum(MM,5));
decayE = abs(mean(aveE(:,:,:,3:end),4));
refE = abs(aveE(:,:,:,2));
echotime=mean(1:(nec-2))*TE2(1);

R2starout= -log(abs(decayE./aveE(:,:,:,2)))/echotime; 

%make a mask;
maxSig=max(abs(aveE(:)));
ave2nd=aveE(:,:,:,1)/maxSig;
R2starout(ave2nd(:)<0.2)=NaN;


