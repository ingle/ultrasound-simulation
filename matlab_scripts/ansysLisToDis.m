function [DX, DY] = interpsave2(in_xfile, in_yfile, in_cfile)

%given 3 .lis files generated in ANSYS.  This MATLAB code will create the
%.dis file the C++ program needs to generate a post-compression phantom
%from a pre-compression phantom.  It also creates plots of the displacement
%fields

%the 3 .lis files are the x displacement, y displacement, and nodal
%coordinates.  They are text files.

% close all; 
% clear all;



%%%%%%%%%first read in all the nodal positions%%%%%%%%
h = fopen(in_cfile);
nodeNo = 1;

while(1)
    tline = fgetl(h);
    if ~ischar(tline), break, end
    
    %%check to see if line is empty
    if(~isempty(tline))
       if( numel(tline) >= 8)
        %%node number is eighth character in a line
        if(~isempty(str2num(tline(8))))
        
            x(nodeNo) = str2double(tline(14:26));
            y(nodeNo) = str2double(tline(34:46));
            
            nodeNo = nodeNo + 1;
        
        end
       end
    end
end

fclose(h);

%%%read in the x displacements
h = fopen(in_xfile);
nodeNo = 1;

while(1)
    tline = fgetl(h);
    if ~ischar(tline), break, end
    
    %%check to see if line is empty
    if(~isempty(tline))
       if(numel(tline) >=8)
        %%node number is eighth character in a line
        if(~isempty(str2num(tline(8))))
        
            Ux(nodeNo) = str2double(tline(10:21));
         
            nodeNo = nodeNo + 1;
        end
       end
    end
end

fclose(h);

%%%%%read in the y displacements
h = fopen(in_yfile);
nodeNo = 1;

while(1)
    tline = fgetl(h);
    if ~ischar(tline), break, end
    
    %%check to see if line is empty
    if(~isempty(tline))
       
        %%node number is eighth character in a line
        if(numel(tline) >= 8)
        if(~isempty(str2num(tline(8))))
        
            Uy(nodeNo) = str2double(tline(10:21));
            
            nodeNo = nodeNo + 1;
        end
        end
    end
end

fclose(h);






%creating a line with regularly sampled points: make sure it matches ansys geometry
xi=linspace(min(x(:)),max(x(:)),400); 
yi=linspace(min(y(:)),max(y(:)),400);

%creating a grid of x positions and y positions
[XI,YI] = meshgrid(xi,yi);


%the nodal displacements are rearranged and interpolated
%to be indexed by x and y positions in a regular grid
xDataInterp = TriScatteredInterp(x',y',Ux', 'linear');
yDataInterp = TriScatteredInterp(x',y',Uy', 'linear');

DX = xDataInterp(XI,YI);
DY = yDataInterp(XI,YI);

DY = flipud(DY);
DX = flipud(DX);
