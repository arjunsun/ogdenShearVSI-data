function out=matrix_invrecfit3(indata,taxis,varargin)
%out=matrix_invrecfit(indata,taxis,['oddeven'])
%indata:        3D or 4D stack of image matrices (the same slice) acquired at different times (last dimension
%taxis:         times after which the images were acquired
%'oddeven':     option to average odd even slices, eliminates funky systematics
%'convergence': option to monitor the convergence of the amplitude plot
%
%out.T1, out.amplitude contain the fitted results 

options=varargin;

indata=squeeze(indata);

if ndims(indata) == 3;                          %one slice, multiple times
    out = matrix_invrecfit_3d(indata,taxis,options{:});
elseif ndims(indata) ==4;                       %multislice, multiple times
    [ntd nph nsl nti] = size(indata);
    out.amplitude=zeros(ntd,nph,nsl);
    out.err_amplitude=zeros(ntd,nph,nsl);
    out.T1=zeros(ntd,nph,nsl);
    out.err_T1=zeros(ntd,nph,nsl);
    out.err_chi2=zeros(ntd,nph,nsl);
    
    for ns=1:nsl;
        oout= matrix_invrecfit_3d(squeeze(indata(:,:,ns,:)),taxis,options{:});
        out.amplitude(:,:,ns)=oout.amplitude;
        out.T1(:,:,ns)=oout.T1;
        out.err_amplitude(:,:,ns)=oout.err_amplitude;
        out.err_T1(:,:,ns)=oout.err_T1;
        %out.err_chi2(:,:,ns)=oout.err_chi2;
    end
end
out.T1(abs(out.T1(:))>5)=0;

function out=matrix_invrecfit_3d(incube,taxis,varargin)
%out=matrix_invrecfit_3d(incube,taxis,weighting)
%fit to exponential recovery
%incube: nro X npe (x numslices) X numel(taxis) is a stack of images acquired at different times
%the times passed down in the taxis.
%out is structure with fields 'amplitude' and 'T1', result of a pixelwise
%
%the inversion recovery  ' y= a*(1-2*exp(bt)) '
%can't be done by straight linear regression, so this routine implements an
%iterative linear regression, where an initial amplitude matrix a is
%guessed, then updated iteratively

options=varargin;

[taxis,sortindex]=sort(taxis); %sort by ascending times
incube=incube(:,:,sortindex);  %sort by ascending times

[tmin,tminind]=min(taxis);

%Third dimension is the echotime dimension
phase_matrix=incube(:,:,tminind)./abs(incube(:,:,tminind)); %phase of the image with the shortest wait time
phased_incube=zeros(size(incube));

%1. rotate the phase, short times are positive, long times negative
for nt=1:numel(taxis);
    phased_incube(:,:,nt)= real(incube(:,:,nt)./phase_matrix);
end

%2. generate a starting guess for the amplitude, 
% option 1: interpolating backwards on the two shortest times
%a0= abs( phased_incube(:,:,1) - (phased_incube(:,:,2)-phased_incube(:,:,1))/(taxis(2)-taxis(1)) * taxis(1) );
% option 2: half the sum of the first and last image
a0= (phased_incube(:,:,1)-phased_incube(:,:,end))/2;

noiseestimate=2*estimate_noiselevel(a0);


TR=10;%default if TR is invoked without time



%3. iteratively fit and update the amplitude guess
dY= diff(phased_incube,1,3);
dtaxis=diff(taxis);

for jj=1:numel(dtaxis);
    dY(:,:,jj)=dY(:,:,jj)/dtaxis(jj);       %coarse derivative
    mtaxis(jj)=(taxis(jj)+taxis(jj+1))/2;
end
outY=matrix_expfit(dY,mtaxis);
outY.Tau(outY.Tau>5)=1e-5;
outY.Tau(outY.Tau<0)=5e-3;

initamp=blur(outY.amplitude/2.*outY.Tau);
initTau=blur(outY.Tau);

%update the TR correction
if any(strcmp(options,'TR'));
    ind=find(strcmp(options,'TR'));
    TR=options{ind+1};
    TRT1correction=initamp;
else
    TRT1correction=zeros(size(initamp));
end

correction_multiplier=0.9*ones(size(TRT1correction));


niter=30;
for jj=1:niter;
    for nt=1:numel(taxis);
        Y(:,:,nt)= phased_incube(:,:,nt)+initamp+TRT1correction;
    end
    
    outY=matrix_expfit(Y,taxis,options{:});
    initamp=outY.amplitude/2;
    
    if any(strcmp(options,'TR'));
        TRT1correction=TRT1correction.*correction_multiplier;
        correction_multiplier(outY.error.chi2.chi2>1)=1;
    end
end

% fit it once more, now taking the average of the expected monoexponential
% inversion recovery correction and the correction found by chi2
if any(strcmp(options,'TR'));
    %TRT1correction=(TRT1correction+initamp.*exp(-TR./outY.Tau))/2; 
    for nt=1:numel(taxis);
        Y(:,:,nt)= phased_incube(:,:,nt)+initamp+TRT1correction;
    end
    outY=matrix_expfit(Y,taxis,options{:});
end




out.amplitude=outY.amplitude/2;
out.T1=outY.Tau;
out.err_amplitude=outY.error.amplitude;
out.err_T1=outY.error.Tau;
out.err_chi2=outY.error.chi2.chi2;

