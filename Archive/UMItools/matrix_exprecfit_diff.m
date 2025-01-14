function [T1,R1]=matrix_exprecfit_diff(varargin)
%out=matrix_exprecfit(incube,taxis,maxamp[,errmaxamp])
%fit to exponential recovery
%incube: nro X npe (x numslices) X numel(taxis) is a stack of images acquired at different times
%the times passed down in the taxis.
%out is structure with fields 'amplitude' and 'tau', result of a pixelwise

incube=varargin{1};
taxis=varargin{2};


function [amplitude, tau, R] = matrix_fit(incube,taxis,noiselevel)
%this routine fits the log of the signal to a line
%log(y) = log(a) + bt
oneovers2 = (abs(incube)/noiselevel).^2; %scale noise to account for change of variables to log(y)
sincube=size(incube);
incube=log(abs(incube));

x=ones(sincube);
for jj=1:sincube(3);
    x(:,:,jj)=x(:,:,jj)*taxis(jj);
end
y=incube;
xy=x.*y;
x2=x.*x;

Sx2=sum(x2.*oneovers2,3);
Sy=sum(y.*oneovers2,3);
Sx=sum(x.*oneovers2,3);
Sxy=sum(xy.*oneovers2,3);
Soneovers2=sum(oneovers2,3);

Delta=Soneovers2.*Sx2-Sx.^2;
b= (Soneovers2.*Sxy - Sx.*Sy)./Delta;
a= (Sx2.*Sy-Sx.*Sxy)./Delta;

amplitude=exp(a);
tau=-1./b;
R=-b;

%errors in the fitted parameters:


