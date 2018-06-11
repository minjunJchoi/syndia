%% global parameters

global e
global me
global eps
global c
global mc2
global verbose 
global ds
global Rmin
global Rmax
global zmin
global zmax
global Nz
global Nf
global zstart
global zend
global fstart
global fend
global pstart
global pend
global pint
global Bfactor

e = 1.602*10^-19;
me = 9.109*10^-31;
eps = 8.854*10^-12;
c = 299792458;
mc2 = me*c^2; 

Bfactor = 1; % 1 = include Br and Bz from EFIT

verbose = 0; % plot 
ds = 0.001; % resolution [m]
Rmin = 1.2; % [m]
Rmax = 2.4; % [m]
zmin = -0.5; % [m]
zmax = 0.5; % [m]

Nz = 10; % number of vertical rays for a single channel
Nf = 10; % number of frequency rays for a single channel
zstart = -10; % [mm] first ray at the mini lens
zend = 10; % [mm] last ray at the mini lens
fstart = -0.3; % [GHz] frequency bandwidth
fend = 0.3; % [GHz] frequency bandwidth

pstart = -0.14; % start point = cold resonance + pstart
pend = 0.03; % end point = cold resonance + pend
pint = 0.1; % interpolation step. 0.1 = 10%
%% physical parameters

% calculat path
calpath = 0;

% ECEI shot
fdir = 'Z:\KSTAR\ecei_data\013250';

% output filename
outputf = fullfile(fdir,'eceipath.mat');

% EFIT
eq.shotname = '13250';
eq.edir = 'Z:\KSTAR\keydiag_data\13250_EFIT';
eq.tCur = 2300;

% electron density
eq.ne0 = 4.0; % [10^19 m^-3]
eq.np1 = 0.95; 
eq.np2 = 5; 
eq.ne_eqn = '(0.05 + (ne0/2)*(1 + tanh(np2^2*(np1^2 - (r/(bb)).^2))).*exp(-(r/(bb)).^2/2))'; % [10^-19 m^-3]

% electron temperature
eq.Te0 = 3.6; % [keV]
eq.Tp1 = 0.95; %
eq.Tp2 = 5; %
eq.Te_eqn  = '(0.05 + (Te0/2)*(1 + tanh(Tp2^2*(Tp1^2 - (r/(bb)).^2))).*exp(-(r/(bb)).^2/2))'; % [keV]

% choose device and channels : vacuum propagation 
dn = 1;
fidx = 5; % low 1--high 8 frerquency channel numbering
zidx = 17; % low 1--high 24 vertical channel numbering

Lcz = 9; % gaussian optics coupling e^2-fallding distance [mm] 
Bcf = 0.3; % gaussian optics coupling e^2-fallding distance [mm] 

% load optics settings
opt = load_hdf5_attributes(fdir, eq.shotname);

%% make propfile

[R2d, z2d, F_Br, F_Bz, F_Bt, F_B, F_ne, F_Te] = make_profile(eq);

%% YOU CAN ADD PERTURBATION HERE FOR ACCURATE CALCULATION %%
%%%% YOU CAN ADD PERTURBATION HERE FOR ACCURATE CALCULATION %%%%

%% calculate ece path for selected channels

Rp = cell(1,length(fidx));
zp = cell(1,length(fidx));
theta = cell(1,length(fidx));
fs = cell(1,length(fidx)); % sub ray in frequency
dz = cell(1,length(fidx)); % sub ray in z @ mini lens
if verbose
    figure('color',[1 1 1]) %%%% PLOT %%%%
end
if calpath 
    for k = 1:length(fidx)
        [fs{k}, dz{k}, Rp{k}, zp{k}, theta{k}] = ece_path(R2d, z2d, F_Br, F_Bz, F_Bt, F_B, opt, dn, fidx, zidx);
    end
    save(outputf,'fidx','zidx','fs','dz','Rp','zp','theta');
else
    load(outputf,'fidx','zidx','fs','dz','Rp','zp','theta');
end

%% YOU CAN ADD PERTURBATION HERE FOR FAST CALCULATION %%
%%%% YOU CAN ADD PERTURBATION HERE FOR FAST CALCULATION %%%%

%% calculate ece intensity for selected channels

tic
Imeas = zeros(size(fidx));
Rch = zeros(size(fidx));
zch = zeros(size(fidx));
Is = cell(1,length(fidx));
S = cell(1,length(fidx));
for k = 1:length(fidx)
    Is{k} = zeros(Nf,Nz);
    S{k} = zeros(Nf,Nz);
    for i = 1:length(fs{k})
        for j = 1:length(dz{k})
            if verbose
                plot(Rp{k}{i,j}, zp{k}{i,j}, 'k'); hold all; %%%% PLOT %%%%
            end
            
            % re-define functions for speed
            F_nes = scatteredInterpolant(Rp{k}{i,j}', zp{k}{i,j}', F_ne(Rp{k}{i,j},zp{k}{i,j})'); % [m^-3]
            F_Tes = scatteredInterpolant(Rp{k}{i,j}', zp{k}{i,j}', F_Te(Rp{k}{i,j},zp{k}{i,j})'); % [J]
            F_Bs = scatteredInterpolant(Rp{k}{i,j}', zp{k}{i,j}', F_B(Rp{k}{i,j},zp{k}{i,j})'); % [T]
            
            % calculate ECE intensity
            [Rm, zm, ~, ~, ~, ~, Iece] = ece_intensity(Rp{k}{i,j}, zp{k}{i,j}, theta{k}{i,j}, 0, fs{k}(i)*2*pi()*10^9, opt.harmonic(dn), F_Bs, F_Tes, F_nes);
            
%             % flat response function average
%             Imeas(k) = Imeas(k) + Iece/(Nf*Nz); 

            % Gaussian response function
            dS = exp(-2*(dz{k}(j)/Lcz).^4) * exp(-2*((fs{k}(i)-mean(fs{k}))/Bcf).^4);
            S{k}(i,j) = dS;
            Imeas(k) = Imeas(k) + Iece * dS;
            
            % channel position
            Rch(k) = Rch(k) + Rm/(Nf*Nz); 
            zch(k) = zch(k) + zm/(Nf*Nz); 
            
            if verbose
                plot(Rm, zm, 'o') %%%% PLOT %%%%
            end
            
            Is{k}(i,j) = Iece;
        end
    end
    Imeas(k) = Imeas(k)/sum(sum(S{k}));
end
toc

Imeas(1) % measured ECE intensity
(mean(fs{1})*2*pi()*10^9/(2*pi*c)).^2*F_Te(Rch(1), zch(1)) % Ibb

if verbose
    plot(Rch, zch, 'sq') %%%% PLOT %%%%
end


%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ece intensity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         if get(h.check_cpc_rel, 'value') == 1 
%             % re-define functions for speed
%             F_nes = scatteredInterpolant(Rp', zp', F_ne(Rp,zp)'); % [m^-3]
%             F_Tes = scatteredInterpolant(Rp', zp', F_Te(Rp,zp)'); % [J]
%             F_Bs = scatteredInterpolant(Rp', zp', F_B(Rp,zp)'); % [T]
% 
%             idx2 = find(Rp >= P(1) - 0.14, 1, 'last'); %%%%%%%%%%%%%%%%%%%%% Rp range 
%             idx1 = find(Rp >= P(1) + 0.03, 1, 'last'); 
%             nidx = idx2:-0.1:idx1; %%%%%%%%%%%%%%%%%%%%% Rp resolution 
%             Rpi = interp1(idx2:-1:idx1, Rp(idx2:-1:idx1), nidx);
%             zpi = interp1(idx2:-1:idx1, zp(idx2:-1:idx1), nidx);
%             thetai = interp1(idx2:-1:idx1, theta(idx2:-1:idx1), nidx);
% %             if verbose
% %                 plot(Rpi, zpi, 'LineWidth',1, 'color','r'); drawnow
% %             end
%             
%             [Rc, zc, ~, ~, ~, ba] = ece_intensity(Rpi, zpi, thetai, P(1), omega, hm, F_Bs, F_Tes, F_nes);
%         end            
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ece intensity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%