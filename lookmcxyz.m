% lookmcxyz.m
%   Looks at myname_F.bin, created by mcxyz.c 
%   where myname is the name of the run: myname_T.bin, myname_H.mci
%   Makes figures:
%       myname_tissue.jpg   = tissue structure (shows tissue types)
%       myname_Fzx.jpg      = fluence rate vs z,x
%       myname_Fzy.jpg      = fluence rate vs z,y
%   Uses:
%       myname_H.mci    = input file from maketissue.m
%       myname_T.bin    = tissue input file from maketissue.m
%       myname_F.bin    = fluence rate output from Monte Carlo
%       reportH_mci.m   = lists input parameters in myname_H.mci
%       makecmap.m      = makes colormap for tissue types
%       makec2f.m       = makes colormap for fluence rate
%
%   This example sets myname = 'skinvessel'.
%
% 7/feb/2017, add boundaryflag (see A(10)).
% 1/june/2017 , no major changes, just clean up display outputs.
% Steven L Jacques

home; clear
format compact
commandwindow


SAVEPICSON = 1;
if SAVEPICSON
    sz = 10; fz = 10; fz2 = 10; % to use savepic.m
else
    sz = 12; fz = 10; fz2 = 10; % for screen display
end

%%%% USER CHOICES <---------- you must specify -----
myname = 'skinvessel'; 
nm = 940;
%%%%


disp(sprintf('------ mcxyz %s -------',myname))

% Load header file
filename = sprintf('%s_H.mci',myname);
disp(['loading ' filename])
fid = fopen(filename, 'r');
A = fscanf(fid,'%f',[1 Inf])';
fclose(fid);

%% parameters
Nphotons = A(1);
Nx = A(2);
Ny = A(3);
Nz = A(4);
dx = A(5);
dy = A(6);
dz = A(7);
mcflag = A(8);
launchflag = A(9);
boundaryflag = A(10);
xs = A(11);
ys = A(12);
zs = A(13);
xfocus = A(14);
yfocus = A(15);
zfocus = A(16);
ux0 = A(17);
uy0 = A(18);
uz0 = A(19);
radius = A(20);
theta = A(21);
alpha = A(22);
waist = A(23);
width = A(24);
Nt = A(25);
j = 25;
for i=1:Nt
    j=j+1;
    muav(i,1) = A(j);
    j=j+1;
    musv(i,1) = A(j);
    j=j+1;
    gv(i,1) = A(j);
end

reportHmci(myname)

%% Load Fluence rate F(y,x,z) 
filename = sprintf('%s_F.bin',myname);
disp(['loading ' filename])
tic
    fid = fopen(filename, 'rb');
    [Data count] = fread(fid, Ny*Nx*Nz, 'float');
    fclose(fid);
toc
F = reshape(Data,Ny,Nx,Nz); % F(y,x,z)


%%
% Load tissue structure in voxels, T(y,x,z) 
filename = sprintf('%s_T.bin',myname);
disp(['loading ' filename])
tic
    fid = fopen(filename, 'rb');
    [Data count] = fread(fid, Ny*Nx*Nz, 'uint8');
    fclose(fid);
toc
T = reshape(Data,Ny,Nx,Nz); % T(y,x,z)

clear Data

%%
x = ([1:Nx]-Nx/2-1/2)*dx;
y = ([1:Ny]-Ny/2-1/2)*dy;
z = ([1:Nz]-1/2)*dz;
ux = [2:Nx-1];
uy = [2:Ny-1];
uz = [2:Nz-1];
zmin = min(z);
zmax = max(z);
zdiff = zmax-zmin;
xmin = min(x);
xmax = max(x);
xdiff = xmax-xmin;

% %% Look at structure, Tzx
% Tzx = reshape(T(Ny/2,:,:),Nx,Nz)';
% tissue = makeTissueList(nm);
% Nt = length(tissue);
% 
% % figure(1);clf
% % imagesc(x(ux),z(uz),Tzx(uz,ux),[1 Nt])
% % hold on
% % cmap = makecmap(Nt);
% % colormap(cmap)
% % colorbar
% % set(gca,'fontsize',sz)
% % set(colorbar,'fontsize',1)
% % xlabel('x [cm]')
% % ylabel('z [cm]')
% % title('Tissue','fontweight','normal','fontsize',fz2)
% 
% figure(1); clf
% sz = 12;  fz = 10; 
% imagesc(x(ux),z(uz),Tzx(uz,ux),[1 Nt])
% hold on
% %set(gca,'fontsize',sz)
% %xlabel('x [cm]')
% %ylabel('z [cm]')
% 
% colorbar
% cmap = makecmap(Nt);
% colormap(cmap)
% %set(colorbar,'fontsize',1)
% % label colorbar
% zdiff = zmax-zmin;
% 
% for i=1:Nt
% %     yy = zmin + (Nt-i)/(Nt-1)*zdiff;
% %     text(xmax*1.33,yy, sprintf('%s',tissue(i).name),'fontsize',fz2)
%        yy = (Nt-i)/(Nt-1)*Nz*dz;
%    % text(max(x)*1.23,yy, tissue(i).name,'fontsize',fz)
% end
% %text(xmax,zmin - zdiff*0.06, 'Tissue types','fontsize',fz)
% axis equal image
% %axis([xmin xmax zmin zmax])
%% Look at structure, Tzx
Tzx = reshape(T(Ny/2,:,:),Nx,Nz)';
tissue = makeTissueList(nm);
Nt = length(tissue);

figure(1);clf
imagesc(x(ux),z(uz),Tzx(uz,ux),[1 Nt])
hold on
cmap = makecmap(Nt);
colormap(cmap)
colorbar
set(gca,'fontsize',sz)
set(colorbar,'fontsize',1)
xlabel('x [cm]')
ylabel('z [cm]')
title('Tissue','fontweight','normal','fontsize',fz2)
for i=1:Nt
    yy = zmin + (Nt-i)/(Nt-1)*zdiff;
    text(xmax*1.2,yy, sprintf('%s',tissue(i).name),'fontsize',fz2)
end
% draw launch
N = 20; % # of beam rays drawn
switch mcflag
    case 0 % uniform  YL
        for i=0:N
         plot([-radius*cos(alpha-theta+2*theta*i/N),-radius*cos(alpha-theta+2*theta*i/N)+radius*sin(pi/2-alpha-theta+2*theta*i/N)],...
           [zs, 1.5],'y-')
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
            plot([xs xx],[zs zz],'r-')
        end
        
    case 3 % rectangle
        zz = max(z);
        for i=0:N
            xx = -radius*cos(pi/2-theta/2)+2*radius*cos(pi/2-theta/2)*i/N;
            plot([xx xx],[-radius*sin(pi/2-theta/2+theta*i/N)+2.5 zz],'y-')
        end
end

hold on;

% YL 
x_det = 0;
z_det = 0.5;
Center = [0,1.5];
Q2 = [x_det, z_det];   % position of detector
Q1 = [xs,zs];          % position of source
Nphotons_det = 175;
            

N_blood = 0;             % number of detected photons in epidermis
N_Silicon = 0;            % number of detected photons in papillary dermis
N_upper = 0;           % number of detected photons in upper blood net
N_Rder = 0;            % number of detected photons in reticular dermis
N_deep = 0;            % number of detected photons in deep blood net
N_fat = 0;             % number of detected photons in fat
N_muscle = 0;          % number of detected photons in muscle
N_bone = 0;            % number of detected photons in bone
 
circle(0,0.5,0.02,3)  % draw detector
hold on;

%  draw trace of detected photons
position = load('detected.txt');
row = load('row.txt');

   
      % first detected photon
         X = position(1:row(1,1),1);
         Y = position(1:row(1,1),2);
         Z = position(1:row(1,1),3);
   
       plot(X,Z,'b');
       hold on 
       
       P = [X, Z];   % position of lauched photon
       OP = zeros(Nphotons_det,1);  % optical path 
       
       count = 1;
       sl = zeros(length(P)-1,1);                 % pathlength in each step
       for is = 2:length(P)
           sl(is-1,1)=sqrt((X(is)-X(is-1)).^2 + (Y(is)-Y(is-1)).^2 + (Z(is)-Z(is-1)).^2);
       end
       OP(1,1) = sum(sl(:));
       %fprintf('Total pathlength of #%d photon is: %f [cm]\n',count,OP(1,1));
       
       depth = zeros(length(P),1); % initial depth
       d = zeros(length(P),1);     % initial distance between the deepst position and center
       depth_max = zeros(Nphotons_det,1);
       d_min = zeros(Nphotons_det,1);
       for ip = 1:length(P)
            depth(ip,1) = abs(det([Q2-Q1;[P(ip,1),P(ip,2)]-Q1]))/norm(Q2-Q1);  % maximal Eindringtiefe
            d(ip,1) = sqrt((P(ip,1)-0)^2 + (P(ip,2)-3.75)^2);
       end
       
       depth_max(1,1) = max(depth(:));
       d_min(1,1) = min(d(:));
   
%        if    d_min(1,1) <= 0.65  % in bone
%            
%           N_blood = N_blood + 1;
%           N_Silicon = N_Silicon + 1;
%           N_upper = N_upper + 1;
%           N_Rder = N_Rder + 1;
%           N_deep = N_deep + 1;
%           N_fat = N_fat + 1;
%           N_muscle = N_muscle + 1;
%           N_bone = N_bone + 1;
% 
%        elseif   d_min(1,1) <= 2.1  % in muscle
%           
%           N_blood = N_blood + 1;
%           N_Silicon = N_Silicon + 1;
%           N_upper = N_upper + 1;
%           N_Rder = N_Rder + 1;
%           N_deep = N_deep + 1;
%           N_fat = N_fat + 1;
%           N_muscle = N_muscle + 1;
%        
%        elseif d_min(1,1) <= 2.35  % in fat
%            
%           N_blood = N_blood + 1;
%           N_Silicon = N_Silicon + 1;
%           N_upper = N_upper + 1;
%           N_Rder = N_Rder + 1;
%           N_deep = N_deep + 1;
%           N_fat = N_fat + 1;
%       
%        elseif  d_min(1,1) <= 2.365  % in deep blood net
%           
%           N_blood = N_blood + 1;
%           N_Silicon = N_Silicon + 1;
%           N_upper = N_upper + 1;
%           N_Rder = N_Rder + 1;
%           N_deep = N_deep + 1;
%          
%        elseif  d_min(1,1) <= 2.46  % in reticular dermis
%           
%           N_blood = N_blood + 1;
%           N_Silicon = N_Silicon + 1;
%           N_upper = N_upper + 1;
%           N_Rder = N_Rder + 1;
% 
%        elseif   d_min(1,1) <= 2.47  % in upper blood net 
%           
%           N_blood = N_blood + 1;
%           N_Silicon = N_Silicon + 1;
%           N_upper = N_upper + 1;
 
       if  d_min(1,1) > 0.95  % in papillary dermis

          N_blood = N_blood + 1;
          N_Silicon = N_Silicon + 1;
           
       elseif  d_min(1,1) >= 1  % in living epidermis

          N_blood = N_blood + 1;
  
       end
       
  
   % from 2nd detected photon to the last
   for i = 1:(size(row,1)-1)
       
       count = count + 1;
       
       X = position((row(i,1)+1):row(i+1,1),1);
       Y = position((row(i,1)+1):row(i+1,1),2);
       Z = position((row(i,1)+1):row(i+1,1),3);
       plot(X,Z,'b');
       hold on
       
       P = [X, Z];   % position of lauched photon
       depth = zeros(length(P),1);  % initial depth
       d = zeros(length(P),1);     % initial distance between the deepst position and center
       sl = zeros(length(P)-1,1); % pathlength in each step                
       for is = 2:length(P)
           sl(is-1,1)=sqrt((X(is)-X(is-1)).^2 + (Y(is)-Y(is-1)).^2 + (Z(is)-Z(is-1)).^2);
       end
       OP(1+i,1) = sum(sl(:));
       %fprintf('Total pathlength of #%d photon is: %f [cm]\n',count, OP(1+i,1));
  
       for ip = 1:length(P)
            depth(ip,1) = abs(det([Q2-Q1;[P(ip,1),P(ip,2)]-Q1]))/norm(Q2-Q1);  % maximal Eindringtiefe
            d(ip,1) = sqrt((P(ip,1)-0)^2 + (P(ip,2)-3.75)^2);
       end
       
       depth_max(i+1,1) = max(depth(:));
       d_min(i+1,1) = min(d(:));
   
%        if     d_min(i+1,1) < 0.65  % in bone
%            
%           N_blood = N_blood + 1;
%           N_Silicon = N_Silicon + 1;
%           N_upper = N_upper + 1;
%           N_Rder = N_Rder + 1;
%           N_deep = N_deep + 1;
%           N_fat = N_fat + 1;
%           N_muscle = N_muscle + 1;
%           N_bone = N_bone + 1;
%        
%        elseif   d_min(i+1,1) < 2.1  % in muscle
%           
%           N_blood = N_blood + 1;
%           N_Silicon = N_Silicon + 1;
%           N_upper = N_upper + 1;
%           N_Rder = N_Rder + 1;
%           N_deep = N_deep + 1;
%           N_fat = N_fat + 1;
%           N_muscle = N_muscle + 1;
%           
%        elseif d_min(i+1,1) < 2.35  % in fat
%           
%           N_blood = N_blood + 1;
%           N_Silicon = N_Silicon + 1;
%           N_upper = N_upper + 1;
%           N_Rder = N_Rder + 1;
%           N_deep = N_deep + 1;
%           N_fat = N_fat + 1;
%       
%        elseif  d_min(i+1,1) < 2.365  % in deep blood net 
%           
%           N_blood = N_blood + 1;
%           N_Silicon = N_Silicon + 1;
%           N_upper = N_upper + 1;
%           N_Rder = N_Rder + 1;
%           N_deep = N_deep + 1;
%  
%        elseif  d_min(i+1,1) < 2.46 % in reticular dermis
%            
%           N_blood = N_blood + 1;
%           N_Silicon = N_Silicon + 1;
%           N_upper = N_upper + 1;
%           N_Rder = N_Rder + 1;
% 
%        elseif   d_min(i+1,1) < 2.47  % in upper blood net
% 
%           N_blood = N_blood + 1;
%           N_Silicon = N_Silicon + 1;
%           N_upper = N_upper + 1;

       if  d_min(i+1,1) > 0.95 % in papillary dermis
          
          N_blood = N_blood + 1;
          N_Silicon = N_Silicon + 1;

       elseif  d_min(i+1,1) > 1  % in living epidermis
          
           N_blood = N_blood + 1;

       end
   end
   
   MOP = (1/count)*sum(OP(:));       % mean optical path
   MD  = (1/count)*sum(depth_max(:));    % mean penetration depth
   
   fprintf('------- number of detected photons through each tissue------\n');
   fprintf('Total detected photons: %d \n',N_blood);
   fprintf('%d detected photons through blood\n',N_blood);
   fprintf('%d detected photons through silicon\n',N_Silicon);
%    fprintf('%d detected photons through upper blood net\n',N_upper);
%    fprintf('%d detected photons through reticular dermis\n',N_Rder);
%    fprintf('%d detected photons through deep blood net\n',N_deep);
%    fprintf('%d detected photons through fat\n',N_fat);
%    fprintf('%d detected photon through muscle\n',N_muscle);
%    fprintf('%d detected photon through bone\n\n',N_bone);
   fprintf('Mean optical path: %f [cm]\n',MOP);
  % fprintf('Optical path: %f [cm]\n',sum(OP(:)));
   fprintf('Mean penetration depth: %f [cm]\n',MD);
   fprintf('Max penetration depth: %f [cm]\n',max(depth_max));
   
axis equal image

if SAVEPICSON
    name = sprintf('%s_tissue.jpg',myname);
    savepic(1,[4 3],name)
end


%% Look at Fluence Fzx @ launch point
% 
% Fzx = reshape(F(Ny/2,:,:),Nx,Nz)'; % in z,x plane through source
% 
% figure(2);clf
% imagesc(x,z,log10(Fzx),[.5 2.8])
% hold on;
% text(max(x)*1.2,min(z)-0.04*max(z),'log_{10}( \phi )','fontsize',fz)
% colorbar
% set(gca,'fontsize',sz)
% xlabel('x [cm]')
% ylabel('z [cm]')
% title('Fluence \phi [W/cm^2/W.delivered] ','fontweight','normal','fontsize',fz)
% colormap(makec2f)
% axis equal image
% axis([min(x) max(x) min(z) max(z)])
% text(min(x)-0.2*max(x),min(z)-0.08*max(z),sprintf('number of photons = %d',Nphotons),...
%     'fontsize',fz2)
% 
% if SAVEPICSON
%     name = sprintf('%s_Fzx.jpg',myname);
%     savepic(2,[4 3],name)
% end

%% Look at Fluence Fxy @ launch point   YL

%  Nz_s = round(zs/dz);
  Nz_det = round(0.5/dz);
  Nz_source = round(zs /dz);
%  Fxy = reshape(F(:,:,Nz_s),Ny,Nx);% in x,y plane through source
% figure(3);clf

xdet_c = 0;  % center of detector in x
ydet_c = 0;  % center of detector in y

xdet = (xdet_c-0.2):dx:(xdet_c+0.2);  % [cm] length of detector in x axis
ydet = (ydet_c-0.2):dy:(ydet_c+0.2);  % [cm] length of detector in y axis
Nxdet = round((xdet_c-0.2)/dx):1:round((xdet_c+0.2)/dx);  % voxel of detector in x
Nydet = round((ydet_c-0.2)/dy):1:round((ydet_c+0.2)/dy);  % voxel of detector in y


% fx = ([((xdet_c-0.2)/dx):((xdet_c+0.2)/dx)]-xdet_c/dx-1/2)*dx*2;
% fy = ([((ydet_c-0.2)/dy):((ydet_c+0.2)/dy)]-ydet_c/dy-1/2)*dy;

% Fdet = zeros(Ny,Nx,Nz);
% for iy = 1:length(Nydet)
%     for ix = 1:length(Nxdet)
%         zdet= -sqrt(radius*radius-(xdet(ix)-2.5)*(xdet(ix)-2.5))+3.75;
%         Nzdet = round(zdet/dz);
%         Fdet(Ny/2-round(length(Nydet)/2)+iy,Nx/2-round(length(Nxdet)/2)+ix,Nzdet)= F(Nydet(iy),Nxdet(ix),Nzdet); 
%         
%     end
% end

Fdet = F((Ny/2-0.2/dy):(Ny/2+0.2/dy),(Nx/2-0.2/dx):(Nx/2+0.2/dx),Nz_det);
value = sum(Fdet(:));
fprintf('Total Fluence rate at detector = %f [W/cm^2/W.delivered]\n',value);

% imagesc(x,y,log10(Fdet(:,:,Nz_det)),[.5 2.8])
% 
% 
% hold on
       
% text(1.25*max(x),1.1*min(y),'log_{10}( \phi )','fontsize',fz)
% colorbar
% set(gca,'fontsize',sz)
% xlabel('x [cm]')
% ylabel('y [cm]')
% title('Fluence in detector \phi ','fontweight','normal','fontsize',fz)
% colormap(makec2f)
% axis equal image
% 
% text(min(x),1.3*min(y),sprintf('number of photons = %d',Nphotons), 'fontsize',fz2)
% text(min(x),1.48*min(y),sprintf('Total Fluence rate = %f [W/cm^2/W.delivered]',value),...
%     'fontsize',fz2)
% 
% if SAVEPICSON
%     name = sprintf('%s_Fxy_det.jpg',myname);
%     savepic(3,[4 3],name)
% end




% %% look Fzy
%  Nx_new = round((radius+xs)/dx);
%  Fzy = reshape(F(:,Nx_new,:),Ny,Nz)'; % in z,y plane through source
% 
% figure(4);clf
% imagesc(y,z,log10(Fzy),[.5 2.8])
% hold on
% text(max(x)*0.7,min(z)-0.04*max(z),'log_{10}( \phi )','fontsize',fz)
% colorbar
% set(gca,'fontsize',sz)
% xlabel('y [cm]')
% ylabel('z [cm]')
% title('Fluence \phi [W/cm^2/W.delivered] ','fontweight','normal','fontsize',fz)
% colormap(makec2f)
% axis equal image
% text(min(x)-0.2*max(x),min(z)-0.07*max(z),sprintf('number of photons = %d',Nphotons),...
%     'fontsize',fz2)
% 
% if SAVEPICSON
%     name = sprintf('%s_Fzy.jpg',myname);
%     savepic(4,[4 3],name)
% end

%% look Azx
Fzx = reshape(F(Ny/2,:,:),Nx,Nz)'; % in z,x plane through source
mua = muav(reshape(T(Ny/2,:,:),Nx,Nz)');
Azx = Fzx.*mua;
figure(2);clf

NX = round((radius+position(1:row(1,1),1))/dx);
NZ = round(position(1:row(1,1),3)/dz);
A = zeros(Nx,Nz);

for in = 1:length(NX)
    A(NX(in),NZ(in)) =  Azx(NX(in),NZ(in));
    
end

for i = 1:(size(row,1)-1)

       NX = round((radius+position((row(i,1)+1):row(i+1,1),1))/dx);
       NZ = round(position((row(i,1)+1):row(i+1,1),3)/dz);
     
       for in = 1:length(NX)
           A(NX(in),NZ(in)) =  Azx(NX(in),NZ(in));  
            
       end
     
end
%imagesc(x,z,log10(A'));
imagesc(x,z,log10(Azx))
hold on
circle(0,0.5,0.02,1.5)  % draw detector
hold on;
circle(xs,zs,0.02,1.5)  % draw source
hold on;

 
text(max(x)*1.2,min(z)-0.04*max(z),'log_{10}( A )','fontsize',fz)
colorbar
set(gca,'fontsize',sz)
xlabel('x [cm]')
ylabel('z [cm]')
title('Deposition A [W/cm^3/W.delivered] ','fontweight','normal','fontsize',fz)
colormap(makec2f)
axis equal image
%axis([min(x) max(x) min(z) max(z)])
text(min(x)-0.2*max(x),min(z)-0.07*max(z),sprintf('number of photons = %d',Nphotons),...
    'fontsize',fz2)

if SAVEPICSON
    name = sprintf('%s_Azx.jpg',myname);
    savepic(2,[4 3],name)
end

% %% look Azy
% Nx_new = round((radius+xs)/dx);
% Fzy = reshape(F(:,Nx_new,:),Ny,Nz)'; % look at source point
% mua = muav(reshape(T(:,Nx_new,:),Ny,Nz)');
% 
% Azy = Fzy.*mua;
% 
% figure(6);clf
% imagesc(y,z,log10(Azy))
% hold on
% text(max(x)*0.7,min(z)-0.04*max(z),'log_{10}( A )','fontsize',fz)
% colorbar
% set(gca,'fontsize',sz)
% xlabel('y [cm]')
% ylabel('z [cm]')
% title('Deposition A [W/cm^3/W.delivered] ','fontweight','normal','fontsize',fz)
% colormap(makec2f)
% axis equal image
% %axis([min(x) max(x) min(z) max(z)])
% text(min(x)-0.2*max(x),min(z)-0.07*max(z),sprintf('number of photons = %d',Nphotons),...
%     'fontsize',fz2)
% 
% if SAVEPICSON
%     name = sprintf('%s_Azy.jpg',myname);
%     savepic(6,[4 3],name)
% end


  

drawnow

disp('done')


