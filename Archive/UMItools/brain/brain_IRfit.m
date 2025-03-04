function [outfit,roicell]=brain_IRfit(varargin)
%spine multi echo fit
%[outfit,roicell,synthimage]=brain_mefit(structure_returned_by_varianms [, sliceindex-for_picking_ROI] [,roi] [,'oddeven']);

data=varargin{1};
imagedata=data.image;

xaxis=data.pars.xaxis;
yaxis=data.pars.yaxis;
inversiontimes=data.pars.ti;
if nargin>1 && isnumeric(varargin{2});
    sli=varargin{2};                %slice index
else
    sidata=size(imagedata);
    sli=round(sidata(3));   %central slice as default
end

%loop to check what has been passed down
for jj=1:nargin;
    structflag(jj)  =isstruct(varargin{jj});
    cellflag(jj)    =iscell(varargin{jj});
    realflag(jj)    =isreal((varargin{jj}));
    %optional other type determinants
end


if any(cellflag);
    dummy=find(cellflag);
    for jj=1:numel(dummy);
        testcell=varargin{dummy(jj)};
        if strcmp(testcell{1},'ROI limits');
            xlim=testcell{2};
            ylim=testcell{3};
            xwinvec=testcell{4};
            ywinvec=testcell{5};
            
            roicell=testcell;
        end
    end
else
    [digiout, xlim, ylim]=pickroi(imagedata(:,:,sli,1));
    xwinvec=xlim(1):xlim(2);
    ywinvec=ylim(1):ylim(2);
    roicell={'ROI limits'; xlim; ylim;xwinvec;ywinvec};
end

datatofit=squeeze(imagedata(ywinvec,xwinvec,:,:));  %default is fit all the data
sidata=size(datatofit);

%get the phase
[minti,minti_index]=min(inversiontimes);
[maxti,maxti_index]=max(inversiontimes);
phaserefdata=datatofit(:,:,:,maxti_index)-datatofit(:,:,:,minti_index);
for sl=1:sidata(3);
    [os,rotmatrix,aux]=gradshim(phaserefdata(:,:,sl));
    for ti=1:numel(inversiontimes);
        datatofit(:,:,sl,ti)=(datatofit(:,:,sl,ti).*rotmatrix);
    end
end

%outfit=matrix_invrecfit(datatofit,inversiontimes);
outfit=matrix_invrecfit2(datatofit,inversiontimes,'TR', data.pars.tr);




