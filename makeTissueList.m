function tissue = makeTissueList(nm)
%function tissueProps = makeTissueList(nm)
%   Returns the tissue optical properties at the wavelength nm:
%       tissueProps = [mua; mus; g]';
%       global tissuenames(i).s
%   Uses 
%       SpectralLIB.mat

%% Load spectral library
load spectralLIB.mat
%   muadeoxy      701x1              5608  double              
%   muamel        701x1              5608  double              
%   muaoxy        701x1              5608  double              
%   muawater      701x1              5608  double              
%   musp          701x1              5608  double              
%   nmLIB         701x1              5608  double              
MU(:,1) = interp1(nmLIB,muaoxy,nm);
MU(:,2) = interp1(nmLIB,muadeoxy,nm);
MU(:,3) = interp1(nmLIB,muawater,nm);
MU(:,4) = interp1(nmLIB,muafat,nm); % getauscht mit 5 IB
MU(:,5) = interp1(nmLIB,muamel,nm);
MU(:,6) = interp1(nmLIB,muasilicon,nm);
MU(:,7) = interp1(nmLIB,mussilicon,nm);
LOADED = 1;



%% Create tissueList
H = 0.45;       % Hematocrit Vrbc/Vblood
SaO2 = 0.99;    % Arterial oxygen saturation
SvO2 = 0.72;    % Venous oxygen saturation
illum = 0;      % 0 for collimated illumination; 1 for diffuse

j=1;
tissue(j).name  = 'air';
tissue(j).mua   = 0.01; % Negligible absorption yet still tracks, 
tissue(j).mus   = 1.0;    % but take steps in air
tissue(j).g     = 1.0;    % and don't scatter.

j=2;
tissue(j).name  = 'water';
tissue(j).mua   = MU(3);
tissue(j).mus   = 10;   % Take steps in water,
tissue(j).g     = 1.0;  % but don't scatter.


j = 10;
tissue(j).name = 'Living epidermis';
thickness(j)    = 0.015;    % Thickness of tissue in cm 
mus0(j)         = 300;      % Scattering coefficients at the reference wavelength 577 nm 
dvessels(j)     = 0;        % Mean diameter of vessels [cm]
Vblood(j)       = 0;        % Blood volume fraction 
VbloodA(j)      =  Vblood(j)*0.7;           % Arterial blood volume fraction   
VbloodV(j)      =  Vblood(j)-VbloodA(j);    % Venous blood volume fraction
Vmel(j)         = 0.1;      % Melanin volume fraction (1-10%)
Vwater(j)       = 0.2;      % Water volume fraction 
Vfat(j)         = 0.151;    % Fat volume fraction
n(j)            = 1.44;     % Refractive index 

j=9;
tissue(j).name  = 'Papillary dermis';
thickness(j)    = 0.015;      % Thickness of tissue in cm from 2002 Meglinski
mus0(j)         = 120;        % Scattering coefficients at the reference wavelength 577 nm
dvessels(j)     = 0.0006;        % Mean diameter of vessels [cm]
Vblood(j)       = 0.04;
VbloodA(j)      =  Vblood(j)*0.7;           % Arterial blood volume fraction   
VbloodV(j)      =  Vblood(j)-VbloodA(j);    % Venous blood volume fraction
Vmel(j)         = 0;        % Melanin volume fraction
Vwater(j)       = 0.5;      % Water volume fraction from 2019 Chatterjee, Kyriacou
Vfat(j)         = 0.1733;        % Fat volume fraction
n(j)            = 1.39;     % Refractive index from Maxim Integrated    

j=8;
tissue(j).name  = 'Blood';
thickness(j)    = 0.01;       % Thickness of tissue in cm from 2002 Meglinski
mus0(j)         = 0;        % Scattering coefficients at the reference wavelength 577 nm
dvessels(j)     = 0.004;        % Mean diameter of vessels  [cm]
Vblood(j)       = 1;        % Blood volume fraction from 2019 Chatterjee, Kyriacou
VbloodA(j)      =  Vblood(j)*0.7;           % Arterial blood volume fraction   
VbloodV(j)      =  Vblood(j)-VbloodA(j);    % Venous blood volume fraction
Vmel(j)         = 0;        % Melanin volume fraction
Vwater(j)       = 0.6;      % Water volume fraction from 2019 Chatterjee, Kyriacou
Vfat(j)         = 0;        % Fat volume fraction
n(j)            = 0;        % Refractive index  

j=7;
tissue(j).name  = 'Reticular dermis';
thickness(j)    = 0.095;      % Thickness of tissue in cm from 2002 Meglinski
mus0(j)         = 120;        % Scattering coefficients at the reference wavelength 577 nm
dvessels(j)     = 0.0015;        % Mean diameter of vessels [cm]
Vblood(j)       = 0.04;        % Blood volume fraction from 2019 Chatterjee, Kyriacou
VbloodA(j)      =  Vblood(j)*0.7;           % Arterial blood volume fraction   
VbloodV(j)      =  Vblood(j)-VbloodA(j);    % Venous blood volume fraction
Vmel(j)         = 0;        % Melanin volume fraction
Vwater(j)       = 0.7;      % Water volume fraction from 2019 Chatterjee, Kyriacou
Vfat(j)         = 0.1733;        % Fat volume fraction
n(j)            = 1.41;     % Refractive index

j=6;
tissue(j).name = 'Deep blood net';
thickness(j)    = 0.015;      % Thickness of tissue in cm from 2002 Meglinski
mus0(j)         = 0;        % Scattering coefficients at the reference wavelength 577 nm
dvessels(j)     = 0.009;        % Mean diameter of vessels [cm]
Vblood(j)       = 0.1;        % Blood volume fraction from 2019 Chatterjee, Kyriacou
VbloodA(j)      =  Vblood(j)*0.7;           % Arterial blood volume fraction   
VbloodV(j)      =  Vblood(j)-VbloodA(j);    % Venous blood volume fraction
Vmel(j)         = 0;        % Melanin volume fraction
Vwater(j)       = 0.7;        % Water volume fraction from 2019 Chatterjee, Kyriacou
Vfat(j)         = 0;        % Fat volume fraction
n(j)            = 1.4;        % Refractive index

j=5;
tissue(j).name  = 'Fat';
thickness(j)    = 0.25;      % Thickness of tissue in cm from 2002 Meglinski
mus0(j)         = 130;        % Scattering coefficients at the reference wavelength 577 nm
dvessels(j)     = 0.0075;        % Mean diameter of vessels [cm]
Vblood(j)       = 0.05;        % Blood volume fraction from 2019 Chatterjee, Kyriacou
VbloodA(j)      =  Vblood(j)*0.7;           % Arterial blood volume fraction   
VbloodV(j)      =  Vblood(j)-VbloodA(j);    % Venous blood volume fraction
Vmel(j)         = 0;        % Melanin volume fraction
Vwater(j)       = 0.05;        % Water volume fraction from 2019 Chatterjee, Kyriacou
Vfat(j)         = 0.8;        % Fat volume fraction
n(j)            = 1.44;        % Refractive index

j=4;
tissue(j).name  = 'Muscle';     % 2020 Chatterjee, Kiriacou
% nm = 940                     
% tissue(j).mua = 1.61;        
% tissue(j).mus = 79.8;
% tissue(j).g   = 0.9112;          
% nm = 810
tissue(j).mua = 1.11;
tissue(j).mus = 81.1;
tissue(j).g   = 0.9088;
% nm = 770
% tissue(j).mua = 1.06;
% tissue(j).mus = 83.8;
% tissue(j).g   = 0.9013;
% nm = 660
% tissue(j).mua = 1.28;
% tissue(j).mus = 85.6;
% tissue(j).g   = 0.8813;
% nm = 530
% tissue(j).mua = 4.15;
% tissue(j).mus = 85.2;
% tissue(j).g   = 0.7949;


j=3;
tissue(j).name  = 'Silicon';
            
tissue(j).mua = MU(:,6);      
tissue(j).mus = MU(:,7);  
tissue(j).g   = 0.93;      

% j=11;
% tissue(j).name  = 'standard tissue';
% tissue(j).mua   = 1;
% tissue(j).mus   = 100;
% tissue(j).g     = 0.90;

%% %% Thicknesses in cm
% for j = 3:1:9
%     tissue(j).dum = thickness(j);
% end

%% Absorption coefficients mua
% from Althschuler: https://iopscience.iop.org/article/10.1088/0022-3727/38/15/027/pdf

% Background tissue absorption coefficient
muabase = 78.4*10^7*nm^(-3.255); % [cm-1]  % Literatur Chatterjee  YL

% V(j) = Volume fraction
% H = Hematocrit
% mua = muaeumel*Vmel + muabil*Vblood +
%   + H*(SaO2*muaoxy+(1-SaO2)*muadeoxy)*VbloodA + 
%   + H*(SvO2*muaoxy+(1-SvO2)*muadeoxy)*VbloodV + (1-H)*muwater*Vblood
%   +  muawater*Vwater + muafat*Vfat +
%   + (muabase+muacarot)*(1-Vmel-(VbloodA+VbloodV)-Vwater-Vfat)
for j=5:1:10
    % Assumption: Vbil = Vcarot = 0
    tissue(j).mua = MU(:,5)*Vmel(j) + ...
                    H*(SaO2*MU(:,1) + (1-SaO2)*MU(:,2))*VbloodA(j) + ...
                    H*(SvO2*MU(:,1) + (1-SvO2)*MU(:,2))*VbloodV(j) + ...
                    (1-H)*MU(:,3)*(VbloodA(j)+VbloodV(j))+ ...
                     MU(:,3)*Vwater(j) + MU(:,4)*Vfat(j) + ...
                    (muabase)*(1-Vmel(j)-(VbloodA(j)+VbloodV(j))-Vwater(j)-Vfat(j));
end

 %% Scaterring coefficient mus
% from Althschuler: https://iopscience.iop.org/article/10.1088/0022-3727/38/15/027/pdf

musblood0 = 4407.2; % cm^(-1)
musblood = musblood0*H*(1-H)*(1.4-H)*(685/nm);  % H = hematocrit defined at the beginning
% disp("musblood = ")
% disp(musblood)
% Correction factor. For very thin vessels = 1; for very thick vessels << 1
if (illum == 0)
    a = 1.007;
    b = 1.228;
elseif (illum == 1)
    a = 1.482;
    b = 1.151;
end
for j=5:1:10  
    % Correction factor
    Ccorr(j) = 1/(1+a*(0.5*musblood*dvessels(j))^b);
%     disp("Ccorr")
%     disp(Ccorr(j))
    % mus0 = scattering coefficients at reference wavelength 577nm defined in
    % each layer
    musTissue(j) = mus0(j)*(577/nm); % For nm<950nm
    tissue(j).mus = Vblood(j)*Ccorr(j)*musblood+(1-Vblood(j))*musTissue(j);
end

%% %% Anisotropy g
% from Althschuler: https://iopscience.iop.org/article/10.1088/0022-3727/38/15/027/pdf

% Anisotropy factor of blood sacterring assumed constant over visible and
% NIR ranges
gb = 0.995;
if (nm>=1125)
    gT = 0.9;
else
    gT = 0.7645+0.2355*(1-exp((-(nm-500)/729.1)));
%     disp("gT = ")
%     disp(gT)
end
for j=5:1:10  
    tissue(j).g = (Vblood(j)*Ccorr(j)*musblood*gb+(1-Vblood(j))*musTissue(j)*gT)/tissue(j).mus;
%     disp("tissue.g")
%     disp(tissue(j).g)
end

%% Print coefficients in command window

disp(sprintf('--- tissueList [1/cm] ---- \tmua   \tmus  \tg  \tmusp'))
for i=1:length(tissue)
    disp(sprintf('%d\t%15s\t%0.4f\t%0.1f\t%0.3f\t%0.1f',...
        i,tissue(i).name,   tissue(i).mua,   tissue(i).mus,   tissue(i).g,...
        tissue(i).mus*(1-tissue(i).g)))
end
disp(' ')

