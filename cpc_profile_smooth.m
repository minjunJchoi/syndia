
tCur = 4200;
Nint = 100;

% read EFIT
edir = 'Z:\KSTAR\keydiag_data\17245_EFIT02_HSKim\';
[g, t_efit, gname] = efit_readg(edir, tCur, 1); 

% profile
pdir = 'Z:\KSTAR\keydiag_data\17245\L-mode_4200';
%% density
% % ne data (R -> normalized_psi_coordinate)
% RawNeFname = fullfile(pdir, sprintf('ne_%d.mat',tCur));
% load(RawNeFname, 'R_ne', 'ne'); % LOAD R_ne [m] and ne [1e19 m^-3]
% [rr, zz] = meshgrid(g.r, g.z);
% psin_ne = interp2(rr, zz, g.psinorm, R_ne, zeros(size(R_ne)));   
% 
% figure; plot(psin_ne, ne,'o'); hold on;
% 
% % need smoothing
% pp = polyfit(psin_ne, ne, 7);
% 
% psin_axis = ((linspace(g.simag,g.sibry,100) - g.simag)/(g.sibry - g.simag)).'; 
% nefit = polyval(pp, psin_axis);
% nefit(nefit < 0) = 0.01;
% 
% plot(psin_axis, nefit, 'x')
% psin_ne = psin_axis;
% ne = nefit;
% NewNeFname = fullfile(pdir, sprintf('ne_%d_smooth.mat',tCur));
% save(NewNeFname, 'psin_ne', 'ne');
% fprintf('save new data %s\n', NewNeFname);
% 
% % No smoothing; just check
% % psin_axis = ((linspace(g.simag,g.sibry,100) - g.simag)/(g.sibry - g.simag)).'; 
% % nefit = interp1(psin_ne, ne, psin_axis, 'spline').*(psin_axis <= 1) + 0.01*(psin_axis > 1);
% 
% plot(psin_axis, nefit);
% 
% [rr, zz] = meshgrid(g.r, g.z);
% ds = 0.001; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1 mm resolution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [R2d, z2d] = meshgrid((1.2:ds:2.4),(-0.5:ds:0.5)); % [m] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% psi2d = interp2(rr, zz, g.psinorm, R2d, z2d);    
% ne2d = interp1(psin_ne, ne, psi2d, 'spline').*(psi2d <= 1) + 0.01*(psi2d > 1);
% 
% figure; 
% imagesc(R2d(1,:), z2d(:,1), ne2d);




% nefunc = @(r) (10.02*(r/bb).^5-34.1*(r/bb).^4+41.74*(r/bb).^3-17.97*(r/bb).^2-3.278*(r/bb)+4.623).*(r/bb<1)+0.1*(r/bb>=1);
% ne_fit = nefunc(psin_ne*bb);
% plot(psin_ne, ne_fit);

%% Temperature 

% Te data (R -> normalized_psi_coordinate)
Tefname = fullfile(pdir, sprintf('Te_%d.mat',tCur));
load(Tefname, 'R_Te', 'Te'); % LOAD R_Te [m] and Te [keV]
[rr, zz] = meshgrid(g.r, g.z);
psin_Te = interp2(rr, zz, g.psinorm, R_Te, zeros(size(R_Te)));

figure; plot(psin_Te, Te,'o'); hold on;

% need smoothing
pp = polyfit(psin_Te, Te, 7);

psin_axis = ((linspace(g.simag,g.sibry,100) - g.simag)/(g.sibry - g.simag)).'; 
Tefit = polyval(pp, psin_axis);
Tefit(Tefit < 0) = 0.01;

plot(psin_axis, Tefit, 'x')
psin_Te = psin_axis;
Te = Tefit;
NewTeFname = fullfile(pdir, sprintf('Te_%d_smooth.mat',tCur));
save(NewTeFname, 'psin_Te', 'Te');
fprintf('save new data %s\n', NewTeFname);

% No smoothing; just check
% psin_axis = ((linspace(g.simag,g.sibry,100) - g.simag)/(g.sibry - g.simag)).'; 
% Tefit = interp1(psin_Te, Te, psin_axis, 'spline').*(psin_axis <= 1) + 0.01*(psin_axis > 1);

plot(psin_axis, Tefit);

[rr, zz] = meshgrid(g.r, g.z);
ds = 0.001; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1 mm resolution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[R2d, z2d] = meshgrid((1.2:ds:2.4),(-0.5:ds:0.5)); % [m] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psi2d = interp2(rr, zz, g.psinorm, R2d, z2d);    
Te2d = interp1(psin_Te, Te, psi2d, 'spline').*(psi2d <= 1) + 0.01*(psi2d > 1);

figure; 
imagesc(R2d(1,:), z2d(:,1), Te2d);


% Tefname = fullfile(edir, sprintf('Te_%d.mat',t_efit));
% load(Tefname, 'psin_Te', 'Te'); % LOAD R_Te and Te [keV]
% 
% figure; plot(psin_Te, Te,'o'); hold on;
% 
% Tefit = interp1(psin_Te, Te, psin_axis, 'spline').*(psin_axis <= 1) + 0.01*(psin_axis > 1);
% 
% plot(psin_axis, Tefit);
% 
% Te2d = interp1(psin_Te, Te, psi2d, 'spline').*(psi2d <= 1) + 0.01*(psi2d > 1);
% figure; 
% imagesc(R2d(1,:), z2d(:,1), Te2d);
% 
% % R coordinate -> normalized psi coordinate
% psirz = g.psirz';
% [rr, zz] = meshgrid(g.r, g.z);
% psin_Te = (interp2(rr, zz, psirz, R_Te, zeros(size(R_Te))) - g.simag)/bb;    
% 
% figure; plot(psin_Te, Te, 'o'); hold on;
% 
% Tefunc = @(r) (-34.8*(r/bb).^6+160.8*(r/bb).^5-287.1*(r/bb).^4+243.5*(r/bb).^3-91.54*(r/bb).^2+4.181*(r/bb)+5.054).*(r/bb<1)+0.1*(r/bb>=1);
% 
% Te_fit = Tefunc(psin_Te*bb);
% plot(psin_Te, Te_fit);


%%(10.02*(r/bb).^5-34.1*(r/bb).^4+41.74*(r/bb).^3-17.97*(r/bb).^2-3.278*(r/bb)+4.623).*(r/bb<=1)+0.1*(r/bb>1);

%%(-34.8*(r/bb).^6+160.8*(r/bb).^5-287.1*(r/bb).^4+243.5*(r/bb).^3-91.54*(r/bb).^2+4.181*(r/bb)+5.054).*(r/bb<=1)+0.1*(r/bb>1);