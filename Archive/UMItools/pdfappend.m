function out=pdfappend(handlearray, filename,varargin)
% out=pdfappend(currentfigurehandlearray, filename,['size',[sizeX sizeY])
% sizes in inches
% appends the figure specified by the input handle to the pdf file
% specified in name. if the file does not exist, creates it


numfigs=numel(handlearray);

if any(strcmp(varargin,'size'));
    ind=find(strcmp(varargin,'size'));
    posvec=varargin{ind+1};
    if numel(posvec==2);
        %print at the specified size in inches, center it
        xmargin=(11-posvec(1))/2;
        ymargin=(8.5-posvec(2))/2;
        posvec=[xmargin ymargin posvec(1) posvec(2)];
    else
        posvec=[0.25 0.25 10.5 8];
        display('size must have two elements, the x size and the y size on a landscape page, in inches.');
    end
else
    posvec=[0.25 0.25 10.5 8];
end
    

for jj=1:numfigs;
    figurehandle=handlearray(jj);
    figure(figurehandle);
    pause(0.2)
    set(figurehandle,'units','inches','PaperOrientation','landscape','Paperunits','inches'); 
    set(figurehandle,'paperpositionmode','auto',...
        'Renderer','painters','paperposition', posvec),
    
    %print the figure to a tempfile
    stringofnumbers= ['000000000' num2str(jj)];
    stringofnumbers=stringofnumbers(end-3:end);
    tempfilename{jj}=['temp_pdf_' stringofnumbers '.pdf'];
    
    
    print(figurehandle,'-dpdf', tempfilename{jj});
    
    close(figurehandle);
end

%concatenate the pdf slides

%commandline:
execute_string=['pdfjam --landscape --outfile ' filename ' temp_pdf*'];
out=system(execute_string);
system('rm temp_pdf_*');
