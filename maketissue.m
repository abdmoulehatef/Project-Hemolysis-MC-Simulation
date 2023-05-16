% maketissue_18apr17.m
% maketissue.m
%   Creates a cube of optical property pointers,T(y,x,z), saved in
%       myname_T.bin = a tissue structure file
%   which specifies a complex tissue for use by mcxyz.c.
%
%   Also prepares a listing of the optical properties at chosen wavelength
%   for use by mcxyz.c, [mua, mus, g], for each tissue type specified
%   in myname_T.bin. This listing is saved in
%       myname_H.mci = the input file for use by mcxyz.c.
%
%   Will generate a figure illustrating the tissue with its various
%   tissue types and the beam being launched.
%
%   Uses
%       makeTissueList.m
%
%   To use, 
%       1. Prepare makeTissueList.m so that it contains the tissue
%   types desired.
%       2. Specify the USER CHOICES.
%       2. Run this program, maketissue.m.
%
%   Note: mcxyz.c can use optical properties in cm^-1 or mm^-1 or m^-1,
%       if the bin size (binsize) is specified in cm or mm or m,
%       respectively.
%
%  Steven L. Jacques. updated Aug 21, 2014.
%       

clear
format compact
clc
home

%%% USER CHOICES %%%%%%%% <-------- You must set these parameters ------
SAVEON      = 1;        % 1 = save myname_T.bin, myname_H.mci 
                        % 0 = don't save. Just check the program.

myname         = 'skinvessel'; % name for files: myname_T.bin, myname_H.mci  
Nphotons       = 1000000;         % number of photons in the simulation 
nm             = 700;   	   % desired wavelength of simulation
Nbins          = 600;    	   % # of bins in each dimension
binsize        = 0.005; 	   % size of each bin, [cm] 
          
% Set Monte Carlo launch flags
mcflag      = 0;     	% launch: 0 = uniform beam, 1 = Gaussian, 2 = isotropic pt. 
                        % 3 = rectangular beam (use xfocus,yfocus for x,y halfwidths)
launchflag  = 0;        % 0 = let mcxyz.c calculate launch trajectory
                        % 1 = manually set launch vector.
boundaryflag = 1;       % 0 = no boundaries, 1 = escape at boundaries
                        % 2 = escape at surface only. No x, y, bottom z
                        % boundaries
                        
% only used if mcflag == 0 or 1 or 3 (not 2=isotropic pt.)
radius      = 1;             % 1/e radius of cylinder model
theta       = atan(0.02/1);  % angle between edge of beam and center of cylinder
waist       = 0.02;  	       % 1/e radius of beam at focus
width       = 0.04;            % width of source at cylinder surface, look at zy

% Sets position of source
xs          = -0.2;     	                        % x of source
ys          = 0;                                    % y of source
zs          = -sqrt(radius*radius-xs*xs)+1.5;    	% z of source
alpha       = acos(abs(xs)/radius);                 % Einfallswinkel zwischen beam und x Achse


% Set position of focus, so mcxyz can calculate launch trajectory
xfocus      = 0;                   % set x,position of focus
yfocus      = 0;                       % set y,position of focus
zfocus      = 1.5;    	               % set z,position of focus (=inf for collimated beam)

% only used if launchflag == 1 (manually set launch trajectory):
ux0         = 0.7;      % trajectory projected onto x axis
uy0         = 0.4;      % trajectory projected onto y axis
uz0         = sqrt(1 - ux0^2 - uy0^2); % such that ux^2 + uy^2 + uz^2 = 1
%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%% 
% Prepare Monte Carlo 
%%%

% Create tissue properties
tissue = makeTissueList(nm); % also --> global tissue(1:Nt).s
Nt = length(tissue);
for i=1:Nt
    muav(i)  = tissue(i).mua;
    musv(i)  = tissue(i).mus;
    gv(i)    = tissue(i).g;
end

% Specify Monte Carlo parameters    
Nx = Nbins;
Ny = Nbins/2;
Nz = Nbins;
dx = binsize;
dy = binsize;
dz = binsize;
x  = ([1:Nx]'-Nx/2)*dx;
y  = ([1:Ny]'-Ny/2)*dy;
z  = [1:Nz]'*dz;
zmin = min(z);
zmax = max(z);
xmin = min(x);
xmax = max(x);
ymin = min(y);
ymax = max(y);

if isinf(zfocus), zfocus = 1e12; end

%%%%%%
% CREATE TISSUE STRUCTURE T(y,x,z)
%   Create T(y,x,z) by specifying a tissue type (an integer)
%   for each voxel in T.
%
%   Note: one need not use every tissue type in the tissue list.
%   The tissue list is a library of possible tissue types.

T = double(zeros(Ny,Nx,Nz)); 

T = T + 1;      % fill background with air

zsurf = 0.0100;  % position of air/skin surface

for iz= 1: Nz % for every depth z(iz)

    % air
    if iz<=round(zsurf/dz)
        T(:,:,iz) = 2; 
    end

    % Living epidermis 
%     xc      = 0;            % [cm], center of epidermis
%    
%     zc      = Nz*0.5*dz;     	% [cm], center of epidermis
%     epidermisradius  = 2.5;      	% epidermis radius [cm]
%     for ix=1:Nx
%             xd = x(ix) - xc;	
%             
%             zd = z(iz) - zc;   	              
%             r  = sqrt(xd^2 + zd^2);	% r from epidermis center
%             if (r<=epidermisradius)     	% if r is within epidermis
%                 T(:,ix,iz) = 10; %  Living epidermis
%             end
%     end
%     
    % Papillary dermis
%     xc      = 0;            % [cm]
%    
%     zc      = Nz*3/4*dz;     	% [cm]
%     p_dermisradius  = 2.485;      	
%     for ix=1:Nx
%             xd = x(ix) - xc;	
%             
%             zd = z(iz) - zc;   	              
%             r  = sqrt(xd^2 + zd^2);	
%             if (r<=p_dermisradius)     	
%                 T(:,ix,iz) = 9; 
%             end
%     end
%    

    % upper blood net @ xc, zc, radius, oriented along y axis

    
     % Reticular dermis
%     xc      = 0;            % [cm]
%    
%     zc      = Nz*3/4*dz;     	% [cm]
%     dermisradius  = 2.46;      	
%     for ix=1:Nx
%             xd = x(ix) - xc;	
%             
%             zd = z(iz) - zc;   	              
%             r  = sqrt(xd^2 + zd^2);	
%             if (r<=dermisradius)     	
%                 T(:,ix,iz) = 7; 
%             end
%     end
%     
%      % Deep blood net
%     xc      = 0;            % [cm]
%    
%     zc      = Nz*3/4*dz;     	% [cm]
%     bloodradius  = 2.365;      	
%     for ix=1:Nx
%             xd = x(ix) - xc;	
%             
%             zd = z(iz) - zc;   	              
%             r  = sqrt(xd^2 + zd^2);	
%             if (r<=bloodradius)     	
%                 T(:,ix,iz) = 6; 
%             end
%     end
%     
%      % Fettgewebe
%     xc      = 0;            % [cm]
%    
%     zc      = Nz*3/4*dz;     	% [cm]
%     fatradius  = 2.35;      	
%     for ix=1:Nx
%             xd = x(ix) - xc;	
%             
%             zd = z(iz) - zc;   	              
%             r  = sqrt(xd^2 + zd^2);	
%             if (r<=fatradius)     	
%                 T(:,ix,iz) = 5; 
%             end
%     end
%     
%     % Muskel 
%     xc      = 0;            % [cm]
%     zc      = Nz*3/4*dz;     	% [cm]
%     muscle_radius  = 2.1 ;      	
%     for ix=1:Nx
%             xd = x(ix) - xc;	
%             zd = z(iz) - zc;   	       
%             r  = sqrt(xd^2 + zd^2);	
%             if (r<=muscle_radius)     	
%                 T(:,ix,iz) = 4; 
%             end
% 
%     end %ix
%     
    % Silicon @ xc, zc, radius, oriented along y axis
    xc      = 0;            % [cm], center of Silicon
    zc      = Nz*0.5*dz;     	% [cm], center of Silicon
    siliconradius  = 1;      	% fat radius [cm]
    for ix=1:Nx
            xd = x(ix) - xc;	% vessel, x distance from vessel center
            zd = z(iz) - zc;   	% vessel, z distance from vessel center                
            r  = sqrt(xd^2 + zd^2);	% r from vessel center
            if (r<=siliconradius)     	% if r is within vessel
                T(:,ix,iz) = 3; % Silicon
            end

    end %ix
    
    xc      = 0;            % [cm], center of blood vessel
    zc      = Nz*0.5*dz;     	% [cm], center of blood vessel
    upperblood  = 0.85 ;      	% blood vessel radius [cm]
    for ix=1:Nx
            xd = x(ix) - xc;	% vessel, x distance from vessel center
            zd = z(iz) - zc;   	% vessel, z distance from vessel center                
            r  = sqrt(xd^2 + zd^2);	% r from vessel center
            if (r<=upperblood)     	% if r is within vessel
                T(:,ix,iz) = 8; % Upper blood net
            end

    end %ix
    
end % iz


%%
if SAVEON
    tic
    % convert T to linear array of integer values, v(i)i = 0;
    v = uint8(reshape(T,Ny*Nx*Nz,1));

    %% WRITE FILES
    % Write myname_H.mci file
    %   which contains the Monte Carlo simulation parameters
    %   and specifies the tissue optical properties for each tissue type.
    commandwindow
    disp(sprintf('--------create %s --------',myname))
    filename = sprintf('%s_H.mci',myname);
    fid = fopen(filename,'w');
        % run parameters
        fprintf(fid,'%d\n',Nphotons);
        fprintf(fid,'%d\n'   ,Nx);
        fprintf(fid,'%d\n'   ,Ny);
        fprintf(fid,'%d\n'   ,Nz);
        fprintf(fid,'%0.4f\n',dx);
        fprintf(fid,'%0.4f\n',dy);
        fprintf(fid,'%0.4f\n',dz);
        % launch parameters
        fprintf(fid,'%d\n'   ,mcflag);
        fprintf(fid,'%d\n'   ,launchflag);
        fprintf(fid,'%d\n'   ,boundaryflag);
        fprintf(fid,'%0.4f\n',xs);
        fprintf(fid,'%0.4f\n',ys);
        fprintf(fid,'%0.4f\n',zs);
        fprintf(fid,'%0.4f\n',xfocus);
        fprintf(fid,'%0.4f\n',yfocus);
        fprintf(fid,'%0.4f\n',zfocus);
        fprintf(fid,'%0.4f\n',ux0); % if manually setting ux,uy,uz
        fprintf(fid,'%0.4f\n',uy0);
        fprintf(fid,'%0.4f\n',uz0);
        fprintf(fid,'%0.4f\n',radius);
        fprintf(fid,'%0.4f\n',theta);
        fprintf(fid,'%0.4f\n',alpha);
        fprintf(fid,'%0.4f\n',waist);
        fprintf(fid,'%0.4f\n',width);
        % tissue optical properties
        fprintf(fid,'%d\n',Nt);
        for i=1:Nt
            fprintf(fid,'%0.4f\n',muav(i));
            fprintf(fid,'%0.4f\n',musv(i));
            fprintf(fid,'%0.4f\n',gv(i));
        end
    fclose(fid);

    %% write myname_T.bin file
    filename = sprintf('%s_T.bin',myname);
    disp(['create ' filename])
    fid = fopen(filename,'wb');
    fwrite(fid,v,'uint8');
    fclose(fid);

    toc
end % SAVEON


%% Look at structure of Tzx at iy=Ny/2
Txzy = shiftdim(T,1);   % Tyxz --> Txzy
Txz  = Txzy(:,:,Ny/2)'; % Tzy

%%
figure(1); clf
sz = 12;  fz = 10; 
imagesc(x,z,Txz,[1 Nt])
hold on
set(gca,'fontsize',sz)
xlabel('x [cm]')
ylabel('z [cm]')

colorbar
cmap = makecmap(Nt);
colormap(cmap)
set(colorbar,'fontsize',1)
% label colorbar
zdiff = zmax-zmin;
%%%

for i=1:Nt
    yy = (Nt-i)/(Nt-1)*Nz*dz;
    text(max(x)*1.33,yy, tissue(i).name,'fontsize',fz)
end
text(xmax,zmin - zdiff*0.06, 'Tissue types','fontsize',fz)
axis equal image
axis([xmin xmax zmin zmax])

%%% draw launch
N = 20; % # of beam rays drawn
switch mcflag
    case 0 % uniform
        
        for i=0:N
           plot([-radius*cos(alpha-theta+2*theta*i/N),-radius*cos(alpha-theta+2*theta*i/N)+radius*sin(pi/2-alpha-theta+2*theta*i/N)],...
           [zs, zc],'y-')
        end
        
    case 1 % Gaussian
        for i=0:N
            plot([-radius*cos(pi/2-theta/2+theta*i/N), xfocus],[-radius*sin(pi/2-theta/2+theta*i/N)+2.5, zfocus],'y-')
        end

    case 2 % iso-point
        for i=1:N
            th = (i-1)/19*2*pi;
            xx = Nx/2*cos(th) + xs;
            zz = Nx/2*sin(th) + zs;
            plot([xs xx],[zs zz],'y-')
        end
        
    case 3 % rectangle
        zz = max(z);
        for i=0:N
            xx = -radius*cos(pi/2-theta/2)+2*radius*cos(pi/2-theta/2)*i/N;
            plot([xx xx],[-radius*sin(pi/2-theta/2+theta*i/N)+2.5 zz],'y-')
        end
end

%% Look at structure of Tzy at ix=Nx/2

Tzy = Txzy(Nx/2,:,:);
TzyNew = squeeze(Tzy);

%%
figure(2); clf
sz = 12;  fz = 10; 
imagesc(y,z,TzyNew,[1 Nt])
hold on
set(gca,'fontsize',sz)
xlabel('y [cm]')
ylabel('z [cm]')

colorbar
cmap = makecmap(Nt);
colormap(cmap)
set(colorbar,'fontsize',1)
% label colorbar
zdiff = zmax-zmin;
%%%

for i=1:Nt
    xx = (Nt-i)/(Nt-1)*Nz*dz;
    text(max(y)*1.4,xx, tissue(i).name,'fontsize',fz)
end

text(ymax,zmin - zdiff*0.06, 'Tissue types','fontsize',fz)
axis equal image
axis([ymin ymax zmin zmax])

%%% draw launch
N = 20; % # of beam rays drawn
switch mcflag
    case 0 % uniform
        for i=0:N
            plot((-width + ys + 2*width*i/N)*[1 1],[0.5 0.5+radius],'y-')
        end

    case 1 % Gaussian
        for i=0:N
            plot([(-width + 2*width*i/N) xfocus],[zs zfocus],'y-')
        end

    case 2 % iso-point
        for i=1:N
            th = (i-1)/19*2*pi;
            xx = Nx/2*cos(th) + xs;
            zz = Nx/2*sin(th) + zs;
            plot([xs xx],[zs zz],'y-')
        end
        
    case 3 % rectangle
        zz = max(z);
        for i=1:N
            xx = -width + 2*width*i/20;
            plot([xx xx],[zs zz],'y-')
        end
end


disp('done')

