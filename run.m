%%
% This code calculates 2D T_rad from the given ne, Te, and EFIT
% M.J. Choi (mjchoi@nfri.re.kr)
% CC BY-NC-SA
%% global parameters
%% hahaha

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
ds = 0.001; % profile resolution [m]
Rmin = 1.2; % [m]
Rmax = 2.4; % [m]
zmin = -0.5; % [m]
zmax = 0.5; % [m]

Nz = 3; % number of vertical rays for a single channel
Nf = 4; % number of frequency rays for a single channel
zstart = -10; % [mm] first ray at the mini lens -14
zend = 10; % [mm] last ray at the mini lens 14
fstart = -0.2; % [GHz] frequency bandwidth -0.3
fend = 0.2; % [GHz] frequency bandwidth +0.3

pstart = -0.14; % [m] start point = cold resonance + pstart
pend = 0.03; % [m] end point = cold resonance + pend
pint = 0.1; % ECE integration path inter step. 0.1 = 10%
%% physical parameters

% calculate path
calpath = 1;

% ECEI shot
fdir = 'C:\data\ecei_data\013728';

% analysis time
eq.tCur = 3.4; % [s]

% EFIT
eq.shotname = '13728';
eq.edir = 'C:\Gdrive\summary\Choi_avalanche_mode_coupling\ECEI_cpc';

% profile data
eq.pdir = 'C:\Gdrive\summary\Choi_avalanche_mode_coupling\ECEI_cpc';


% 1 density functional form (psi2d - g.simag coordinate) or 2 data (nomralized psi coordinate) % [10^-19 m^-3]
eq.ne_opt = 2;

% electron density (functional form; ne_opt = 1)
eq.ne0 = 4.0; % [10^19 m^-3]
eq.np1 = 0.95;
eq.np2 = 5;
eq.ne_eqn = '(0.05 + (ne0/2)*(1 + tanh(np2^2*(np1^2 - (r/(bb)).^2))).*exp(-(r/(bb)).^2/2))';

% electron density (data; ne_opt = 2)
eq.nefname = fullfile(eq.pdir, sprintf('ne_%d_smooth.mat',eq.tCur*1000));

% 1 temperature functional form (psi2d - g.simag coordinate) or 2 data (nomralized psi coordinate) % [keV]
eq.Te_opt = 2;

% electron temperature (functional form; Te_opt = 1)
eq.Te0 = 3.6; % [keV]
eq.Tp1 = 0.95; %
eq.Tp2 = 5; %
eq.Te_eqn  = '(0.05 + (Te0/2)*(1 + tanh(Tp2^2*(Tp1^2 - (r/(bb)).^2))).*exp(-(r/(bb)).^2/2))';

% electron temperature (data; Te_opt = 2)
eq.Tefname = fullfile(eq.pdir, sprintf('Te_%d_smooth.mat',eq.tCur*1000));



% choose device and channels : vacuum propagation
dn = 1;
fidx = (1:8); % low 1--high 8 frerquency channel numbering
zidx = (1:24); % low 1--high 24 vertical channel numbering

Lcz = 9; % gaussian optics coupling e^2-fallding distance [mm]
Bcf = 0.3; % gaussian optics coupling e^2-fallding distance [GHz]


% output filename
dev = 'LHG';
outputf = fullfile(fdir, sprintf('%cD_ecei_modeling_%05d.mat', dev(dn), eq.tCur*1000));


% load optics settings
opt = load_hdf5_attributes(fdir, eq.shotname);

%% make propfile

[R2d, z2d, F_Br, F_Bz, F_Bt, F_B, F_ne, F_Te] = make_profile(eq);

%% YOU CAN ADD PERTURBATION HERE FOR ACCURATE CALCULATION %%
%%%% YOU CAN ADD PERTURBATION HERE FOR ACCURATE CALCULATION %%%%


%% calculate ece path for selected channels

Rp = cell(length(zidx),length(fidx));
zp = cell(length(zidx),length(fidx));
theta = cell(length(zidx),length(fidx));
fs = cell(length(zidx),length(fidx)); % sub ECE ray frequency [GHz]
dz = cell(length(zidx),length(fidx)); % sub ECE ray in z @ mini lens
if verbose
    figure('color',[1 1 1]) %%%% PLOT %%%%
end
if calpath
    for f = 1:length(fidx)
        for z = 1:length(zidx)
            tic

            [fs{z,f}, dz{z,f}, Rp{z,f}, zp{z,f}, theta{z,f}] = ece_path(R2d, z2d, F_Br, F_Bz, F_Bt, F_B, opt, dn, fidx(f), zidx(z));

            fprintf('##################### z %d - f %d PATH done ################# \n', z, f)

            toc
        end
    end
    save(outputf,'fidx','zidx','fs','dz','Rp','zp','theta');
else
    load(outputf,'fidx','zidx','fs','dz','Rp','zp','theta');
end

%% YOU CAN ADD PERTURBATION HERE FOR FAST CALCULATION %%
%%%% YOU CAN ADD PERTURBATION HERE FOR FAST CALCULATION %%%%

%% calculate ece intensity for selected channels

Imeas = zeros(length(zidx),length(fidx));
Trad  = zeros(length(zidx),length(fidx));
Te    = zeros(length(zidx),length(fidx));
Rch   = zeros(length(zidx),length(fidx));
zch   = zeros(length(zidx),length(fidx));
Is     = cell(length(zidx),length(fidx));
S      = cell(length(zidx),length(fidx));
for f = 1:length(fidx)
    for z = 1:length(zidx)
        tic

        Is{z,f} = zeros(Nf,Nz);
        S{z,f} = zeros(Nf,Nz);
        for j = 1:length(fs{z,f})
            for i = 1:length(dz{z,f})
                if verbose
                    plot(Rp{z,f}{i,j}, zp{z,f}{i,j}, 'k'); hold all; %%%% PLOT %%%%
                end

                % re-define functions for speed
                F_nes = scatteredInterpolant(Rp{z,f}{i,j}', zp{z,f}{i,j}', F_ne(Rp{z,f}{i,j},zp{z,f}{i,j})'); % [m^-3]
                F_Tes = scatteredInterpolant(Rp{z,f}{i,j}', zp{z,f}{i,j}', F_Te(Rp{z,f}{i,j},zp{z,f}{i,j})'); % [J]
                F_Bs =  scatteredInterpolant(Rp{z,f}{i,j}', zp{z,f}{i,j}',  F_B(Rp{z,f}{i,j},zp{z,f}{i,j})'); % [T]

                % calculate ECE intensity
                [Rm, zm, ~, ~, ~, ~, Iece] = ece_intensity(Rp{z,f}{i,j}, zp{z,f}{i,j}, theta{z,f}{i,j}, 0, fs{z,f}(j)*2*pi()*10^9, opt.harmonic(dn), F_Bs, F_Tes, F_nes);

    %             % flat response function average
    %             Imeas(f) = Imeas(f) + Iece/(Nf*Nz);

                % Gaussian response function [Rathgeber PPCF2013 (15--17)]
                dS = exp(-2*(dz{z,f}(i)/Lcz).^4) * exp(-2*((fs{z,f}(j)-mean(fs{z,f}))/Bcf).^4);
                S{z,f}(i,j) = dS;
                Imeas(z,f) = Imeas(z,f) + Iece * dS;

                % channel position
                Rch(z,f) = Rch(z,f) + Rm/(Nf*Nz);
                zch(z,f) = zch(z,f) + zm/(Nf*Nz);

                if verbose
                    plot(Rm, zm, 'o') %%%% PLOT %%%%
                end

                Is{z,f}(i,j) = Iece;
            end
        end
        Imeas(z,f) = Imeas(z,f)/sum(sum(S{z,f})); % spatial average  [Rathgeber PPCF2013 (15--17)]
        Trad(z,f) = Imeas(z,f) / ((mean(fs{z,f})*2*pi()*10^9/(2*pi*c)).^2) / (1000*1.602*10^-19); % [keV]

        Te(z,f)   = F_Te(Rch(z,f), zch(z,f)) / (1000*1.602*10^-19); % [keV]
        fprintf('##################### z %d - f %d Imeas done ################# \n', z, f)

        toc
    end
end

% Imeas(1) % measured ECE intensity over frequency and vertical area
% (mean(fs{1})*2*pi()*10^9/(2*pi*c)).^2 * F_Te(Rch(1), zch(1)) % Ibb

% fprintf('first selected channel T_rad [keV] : %g \n', Trad(1,1) );
% fprintf('first selected channel T_e [keV]   : %g \n', Te(1,1) );

save(outputf,'Imeas', 'Trad', 'Te', 'Is', 'S', 'Rch', 'zch', '-append');

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
