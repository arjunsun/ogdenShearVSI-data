function mems=spine_mems(varargin)
%out=spine_mems(intimes)
%intimes is either a 4 digit time string, or a cell array of 4-digit time strings

inputwords={'plot','nav','vox','pix'}; %optional input words 
parseinput(inputwords,varargin);

%default values
planarpixels=384;   %zerofill to get a minimum pixel number of 192

%% input checking
jj=1;
loopend=false;
for kk=1:nargin;
    if (~any(strcmp(inputwords,varargin{kk})) && ~loopend);
        timelist{jj}=varargin{kk};
        jj=jj+1;
    else
        loopend=true;
    end
end

%input checking that all the listed times exist
d=dir;
names={d(:).name}; %cell array of file and directory names
for jj=1:numel(timelist);
    existflag=false;
    for kk=1:numel(names);
        if ~isempty(strfind(names{kk},timelist{jj}));
            existflag=true;
        end
    end
    if ~existflag;
        display(['No directory for time ' timelist{jj}]);
    end
end

%%
for jj=1:numel(timelist);
    data(jj)=varianms(timelist{jj}, num2str(planarpixels),'filter'); %doubles resolution in the 2nd PE dimension
end

%average the data;
mems=data(1);
fsems.image=zeros(size(data(1).image));
for jj=1:numel(data);
    mems.image=mems.image+abs(data(jj).image);
end

