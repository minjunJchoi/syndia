function [R2d, z2d, F_Br, F_Bz, F_Bt, F_B, F_ne, F_Te] = make_profile(eq)

global verbose 
global ds
global Rmin
global Rmax
global zmin
global zmax
global Bfactor

%% use EFIT only

% function make_profile(fname, in, r, R2d, z2d, verbose, dosave)

% %%%%%%%%%%%%%%%%% shaping factors
% if in.eq_opt == 0 %  MODEL
%     R0 = ones(size(r))*in.Rc;
%     z0 = ones(size(r))*in.zc;
%     delta = ones(size(r))*in.delta;
%     kappa = ones(size(r))*in.kappa;
%     rb = in.rb;
% elseif in.eq_opt == 1 % MODEL + EFIT
    % g = efit_readg(edir, shotname, tCur*1000, verbose);    
    % [rr, zz] = meshgrid(g.r, g.z);
    % [edge.rb, edge.delta, edge.kappa, edge.zc, edge.Rc] = get_flux_shape(rr*100, zz*100, g.psirz', g.sibry); % [cm]

    % r = (0.001:0.001:edge.rb);
    % R0 = linspace(in.Rc, edge.Rc, length(r));    
    % z0 = linspace(in.zc, edge.zc, length(r));
    
    % ridx = find(r > in.rin*0.01, 1, 'first');
    % delta = zeros(size(r));
    % delta(1:ridx) = linspace(min(0,in.delta), in.delta, ridx);
    % delta(ridx:end) = linspace(in.delta, edge.delta, length(delta)-ridx+1);
    % kappa = zeros(size(r));
    % kappa(1:ridx) = linspace(min(1,in.kappa), in.kappa, ridx);
    % kappa(ridx:end) = linspace(in.kappa, edge.kappa, length(kappa)-ridx+1);    

    % % radial boundary    
    % rb = edge.rb;    
% elseif in.eq_opt == 2 % EFIT
    % g = efit_readg(edir, shotname, tCur*1000, verbose);    
    % % [rr, zz] = meshgrid(g.r, g.z);
    % % get_flux_shape(rr*100, zz*100, g.psirz', g.sibry); % [cm]

    % % radial boundary
    % rb = g.sibry - g.simag;
% end

% %%%%%%%%%%%%%% 2-D density, temperature, magnetic field equilibrium
% %% density
% nefunc = eval(['@(r) ',in.prof_eqn]);
% 
% if in.eq_opt < 2 % make 2D density profile from delta kappa
%     R = [];
%     z = [];
%     ne = [];
%     for i = 1:length(r)
%         thi = linspace(0, 2*pi, fix(2*pi*r(i)/in.ds))';
%         Ri = R0(i) + r(i)*cos(thi + delta(i)*sin(thi)); % horizontal, starting from (R0 + r, z0)
%         zi = z0(i) + kappa(i)*r(i)*sin(thi); % vertical        
%         nei = ones(size(Ri))*nefunc(r(i)); % 10e19 m^-3
% 
%         R = [R; Ri];  
%         z = [z; zi];  
%         ne = [ne; nei];
% 
%         if 0
%             plot(Ri,zi); hold on;
%             input(sprintf('%d : ',i),'s')
%         end 
%     end
% 
%     F_ne = scatteredInterpolant(R, z, ne); % 10e19 m^-3
%     ne2d = F_ne(R2d,z2d);
%     ne2d(ne2d < 0) = 0.01;        
%     ne2d(isnan(ne2d)) = 0.01;
% 
% % elseif in.eq_opt == 2 % make 2D density profile from EFIT psi contour
%     % psi2d = interp2(rr, zz, g.psirz', R2d, z2d);    
%     % ne2d = nefunc(psi2d - g.simag); % 10e19 m^-3
%     % ne2d(ne2d < 0) = 0.01;
%     % ne2d(isnan(ne2d)) = 0.01;
% end
% 
% %% temperature
% Te2d = ne2d*in.Trn; % [keV]
% 
% %% magnetic field
% % toroidal field 
% Bt2d = repmat(in.Bt*1.8./R2d(1,:), [size(R2d,1) 1]);
% % poloidal field
% Bp2d = zeros(size(Bt2d));
% 
% % % poloidal field from EFIT
% % if eq.cpc_eqopt == 2
% %     [dpsidr, dpsidz] = gradient(g.psirz');
% %     dpsidr = dpsidr/(g.r(2)-g.r(1));
% %     dpsidz = dpsidz/(g.z(2)-g.z(1));
% %     dpsidr = interp2(rr, zz, dpsidr, R2d, z2d);
% %     dpsidz = interp2(rr, zz, dpsidz, R2d, z2d);
% %     Bp2d = sqrt(dpsidr.^2 + dpsidz.^2)./R2d;
% % else
% %     Bp2d = zeros(size(Bt2d));
% % end
% 
% % total field
% B2d = sqrt(Bt2d.^2 + Bp2d.^2);

%% input parameters for equilibrium

% % EFIT
% eq.edir = 'Z:\KSTAR\keydiag_data\13250_EFIT';
% eq.tCur = 2300;
% 
% % electron density
% eq.ne0 = 4.0; % [10^19 m^-3]
% eq.np1 = 0.95; 
% eq.np2 = 5; 
% eq.ne_eqn = '(0.05 + (ne0/2)*(1 + tanh(np2^2*(np1^2 - (r/(bb)).^2))).*exp(-(r/(bb)).^2/2))'; % [10^-19 m^-3]
% 
% % electron temperature
% eq.Te0 = 3.6; % [keV]
% eq.Tp1 = 0.95; %
% eq.Tp2 = 5; %
% eq.Te_eqn  = '(0.05 + (Te0/2)*(1 + tanh(Tp2^2*(Tp1^2 - (r/(bb)).^2))).*exp(-(r/(bb)).^2/2))'; % [keV]

%% coordinates and spatial resolution for interpolation

[R2d, z2d] = meshgrid((Rmin:ds:Rmax),(zmin:ds:zmax));         % [m] 

%% read EFIT

[g, ~, gname] = efit_readg(eq.edir, eq.tCur*1000, 0);     
psirz = g.psirz';
[rr, zz] = meshgrid(g.r, g.z);

% boundary
bb = g.sibry - g.simag;

%% read geqdsk file from the selected EFIT run time and convert it to topfile for TORBEAM

gname(strfind(gname,'\')) = '/';                                                   
system(sprintf('cp %s C:/cygwin64/home/mjchoi/eqdsk2topfile/EFIT.geqdsk', gname)); 
fprintf('     Use the EFIT file %s for TORBEAM\n', gname);
oldFolder = cd('C:\cygwin64\home\mjchoi\eqdsk2topfile');
system('readeqdsk<EFIT.geqdsk');
cd(oldFolder);
% send topfile to the iKSTAR server for magnetic field equilibrium to run TORBEAM
system('scp -P 2201 /home/mjchoi/eqdsk2topfile/topfile mjchoi@172.17.250.11:~/torbeam_ifortran/run_torbeam');

%% 2D B profile [T]

B0 = abs(g.bcentr); % [T]  

[dpsidr,dpsidz] = gradient(psirz, g.r(2) - g.r(1), g.z(2) - g.z(1));
F_Br = @(x,y) -interp2(rr, zz, dpsidz, x, y)./x*Bfactor;
F_Bz = @(x,y) interp2(rr, zz, dpsidr, x, y)./x*Bfactor;
F_Bt = @(x,y) B0*1.8./x;
F_B = @(x,y) sqrt(F_Br(x,y).^2 + F_Bz(x,y).^2 + F_Bt(x,y).^2);

%% 2D electron density profile [1e19]

ne0 = eq.ne0; % [10^19 m^-3] 
np1 = eq.np1; %
np2 = eq.np2; %
nefunc = eval(['@(r)',eq.ne_eqn]);

% scatteredInterpolant function
psi2d = interp2(rr, zz, psirz, R2d, z2d);    
ne2d = nefunc(psi2d - g.simag);
F_ne = scatteredInterpolant(reshape(R2d, numel(R2d), 1), reshape(z2d, numel(z2d), 1), reshape(ne2d*10^19, numel(ne2d), 1)); % [m^-3]

%% 1D electron density data for TORBEAM [1e19]

psi_axis = linspace(g.simag,g.sibry,80) - g.simag; 
ne1d = nefunc(psi_axis);

% save
fileID = fopen(fullfile('C:\cygwin64\home\mjchoi\eqdsk2topfile', 'ne.dat'),'w');
fprintf(fileID,'%d\n',80);
fprintf(fileID,'%f %f\n',[sqrt(psi_axis(end:-1:1)/bb); ne1d(end:-1:1)]);
fclose(fileID);

% send ne.dat to the iKSTAR server to run TORBEAM
system('scp -P 2201 /home/mjchoi/eqdsk2topfile/ne.dat mjchoi@172.17.250.11:~/torbeam_ifortran/run_torbeam');

%% 2D electron temperature profile [keV]

Te0 = eq.Te0; % [keV]  
Tp1 = eq.Tp1; %
Tp2 = eq.Tp2; %
Tefunc = eval(['@(r)',eq.Te_eqn]); % [keV]

psi2d = interp2(rr, zz, psirz, R2d, z2d);    
Te2d = Tefunc(psi2d - g.simag);
F_Te = scatteredInterpolant(reshape(R2d, numel(R2d), 1), reshape(z2d, numel(z2d), 1), reshape(Te2d*1000*1.602*10^-19, numel(Te2d), 1)); % [J]

%% 1D electron temperature data for TORBEAM [keV]

Te1d = Tefunc(psi_axis);

% save
fileID = fopen(fullfile('C:\cygwin64\home\mjchoi\eqdsk2topfile', 'Te.dat'),'w');
fprintf(fileID,'%d\n',80);
fprintf(fileID,'%f %f\n',[sqrt(psi_axis/bb); Te1d]);
fclose(fileID);

% send Te.dat to the iKSTAR server to run TORBEAM
system('scp -P 2201 /home/mjchoi/eqdsk2topfile/Te.dat mjchoi@172.17.250.11:~/torbeam_ifortran/run_torbeam');


if verbose
    figure('color',[1 1 1]);
    subplot(2,1,1); contour(R2d(1,:),z2d(:,1),ne2d,50); title('ne'); axis equal;
    subplot(2,1,2); plot(R2d(1,:),F_ne(R2d(1,:),zeros(size(R2d(1,:))))); xlabel('R'); ylabel('ne [1e19]');
    figure('color',[1 1 1]);
    subplot(2,1,1); contour(R2d(1,:),z2d(:,1),Te2d,50); title('Te'); axis equal;
    subplot(2,1,2); plot(R2d(1,:),F_Te(R2d(1,:),zeros(size(R2d(1,:))))); xlabel('R'); ylabel('Te [J]');
end

% %%%%%%%%%%%%%% save
% 
% if verbose
%     midx = round(size(ne2d,1)/2);
% 
%     figure('color',[1 1 1]);
%     subplot(2,1,1); contour(R2d(1,:),z2d(:,1),ne2d,50); title('ne'); axis equal;
%     subplot(2,1,2); plot(R2d(1,:),ne2d(midx,:)); xlabel('R');
%     figure('color',[1 1 1]);
%     subplot(2,1,1); contour(R2d(1,:),z2d(:,1),Te2d,50); title('Te'); axis equal;
%     subplot(2,1,2); plot(R2d(1,:),Te2d(midx,:)); xlabel('R');
% end
% 
% save(fname,'R0','z0','kappa','delta','ne2d','Te2d','B2d','-append');

% end

end