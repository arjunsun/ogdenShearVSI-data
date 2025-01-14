function out = lig_Q_anatomyplot(varargin)
% out= lig_Q_anatomyplot(anatomy,Q) combines Q strain data and anatomical image to create a
% a hybrid color BW plot tinting the BW anatomy images

A=abs(varargin{1}); % anatomy
Q=varargin{2}; % Q

%(1) assume the 3D-matrices cover the same volume, but not necessarily the same grid
siQ=size(Q);
siA=size(A);

% set the size of the detailed matrices to be calculated by interpolaion
siX=max(siQ,siA);
ND=ndims(siX);

%generate normalized axes and grids for the data cubes
axQ1=(1:siQ(1))/siQ(1); axQ2=(1:siQ(2))/siQ(2); axQ3=(1:siQ(3))/siQ(3);
[QG1,QG2,QG3]=ndgrid(axQ1,axQ2,axQ3);

axA1=(1:siA(1))/siA(1); axA2=(1:siA(2))/siA(2); axA3=(1:siA(3))/siA(3);
[AG1,AG2,AG3]=ndgrid(axA1,axA2,axA3);

axis1=(1:siX(1))/siX(1);axis2=(1:siX(2))/siX(2);axis3=(1:siX(3))/siX(3);
[X1,X2,X3]=ndgrid(axis1,axis2,axis3);

iA=interpn(AG1,AG2,AG3,A,X1,X2,X3);
iQ=interpn(QG1,QG2,QG3,Q,X1,X2,X3);

% iA and iQ are the interpolated anatomy and Q matrix
iA=2*iA/max(iA(:));
[sigA,noiseA]=estimate_noiselevel(iA);

iQ=iQ/max(iQ(:));
iQ(isnan(iQ(:)))=0;
[sigQ,noiseQ]=estimate_noiselevel(iQ);

%create a a red mask for the Q data
siQ=size(iQ);
red=cat(3,ones(siQ(1:2)), 0.1*ones(siQ(1:2)), zeros(siQ(1:2)));

%%
close all;
for sl=1:64;
    figure; colormap gray;
    hA=imagesc(iA(:,:,sl),[0 1]);
    hold on;
    hi=imshow(red);
    opacity_of_red=iQ(:,:,sl).*(iA(:,:,sl))*3;
    set(hi,'alphadata',opacity_of_red);
end




end

