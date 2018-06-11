% FUNCTION: g = efit_readg(fdir, shot, t, verbose)
%
% read EFIT G file
% fdir: directory containing the EFIT files
% shot: shot number
% t: time in msec
% g: EFIT G structure
% t_efit: the time point at which the EFIT is calculated (in msec). 
% 
% 2010/12/10, Gunsu S. Yun 
%

%% Referece:
% Variables
% CASE 	Identification character string
% NW	Number of horizontal R grid points
% NH	Number of vertical Z grid points
% 
% Namelist OUT1:
% 
%     BCENTR 	Vacuum toroidal magnetic field in Tesla at RCENTR
%     CURRENT 	Plasma current in Ampere
%     FPOL	Poloidal current function in m-T, F = RBT on flux grid
%     FFPRIM	FF'(psi) in (mT)2 / (Weber/rad) on uniform flux grid
%     LIMITR	Number of limiter points
%     NBBBS	Number of boundary points
%     NQPSI
%     PPRIME	P'(psi) in (nt/m2) / (Weber/rad) on uniform flux grid
%     PRES	Plasma pressure in nt / m2 on uniform flux grid
%     PSIZR	Poloidal flux in Weber/rad on the rectangular grid points
%     QPSI	q values on uniform flux grid from axis to boundary
%     RBBBS	R of boundary points in meter
%     RCENTR	R in meter of vacuum toroidal magnetic field BCENTR
%     RDIM	Horizontal dimension in meter of computational box
%     RLEFT	Minimum R in meter of rectangular computational box
%     RLIM	R of surrounding limiter contour in meter
%     RMAXIS	R of magnetic axis in meter
%     SIMAG	poloidal flux at magnetic axis in Weber/rad
%     SIBRY	poloidal flux at the plasma boundary in Weber/rad
%     ZBBBS	Z of boundary points in meter
%     ZDIM	Vertical dimension in meter of computational box
%     ZLIM	Z of surrounding limiter contour in meter
%     ZMAXIS	Z of magnetic axis in meter
%     ZMID	Z of center of computational box in meter 


function [g, t_efit, fname] = efit_readg(fdir, t, verbose)

% fdir = sprintf('%s/%s_EFIT',fdir,shot);
flist = ls(fdir);
tlist = zeros(size(flist,1),1);
for i=3:size(flist,1)
    st = sscanf(flist(i,:),'g%d.%d');
    if length(st) > 1
        tlist(i) = st(2);
    end
end

fname = strtrim(flist(find(abs(tlist - t)==min(abs(tlist - t)), 1, 'first'),:));
t_efit = tlist(find(abs(tlist - t)==min(abs(tlist - t)), 1, 'first'));
fname = fullfile(fdir,fname)

if exist(sprintf('%s.mat',fname), 'file') 
    load(sprintf('%s.mat',fname));
%     disp(fname)
else
    if exist(fname, 'file') 
        fid = fopen(fname, 'rt');

        %% basic EFIT parameters
        %line01:   EFITD    01/23/2002    #  4362  1701ms  3    65  65
        s = fgetl(fid);
        a = sscanf(s, '%*s %*s # %d %dms %d %d %d'); % shotid, time, mw, mh
        if length(a) == 2
            a = sscanf(s, '%*s %*s # %d %d %d %d %d'); % shotid, time, mw, mh
        end
        g.shot = a(1);
        g.t = a(2);
        g.mw = a(4);
        g.mh = a(5);

        %line02: g.xdim,g.zdim,g.rzero,g.rgrid1,g.zmid
        s = fgetl(fid);
        a = sscanf(s, '%f %f %f %f %f');
        g.xdim = a(1);
        g.zdim = a(2);
        g.r0 = a(3);
        g.rgrid1 = a(4);
        g.zmid = a(5);

        %line03: g.rmaxis,g.zmaxis,g.simag,g.sibry,g.bcentr
        s = fgetl(fid);
        a = sscanf(s, '%f %f %f %f %f');
        g.rmaxis = a(1);
        g.zmaxis = a(2);
        g.simag = a(3);
        g.sibry = a(4);
        g.bcentr = a(5);

        %line04: g.cpasma,g.simag,0.,g.rmaxis,0.
        s = fgetl(fid);
        a = sscanf(s, '%f %f %f %f %f');
        g.cpasma = a(1);
        g.simag = a(2);
        g.rmaxis = a(4);

        %line05: g.zmaxis,0.,g.sibry,0.,0.
        s = fgetl(fid);
        a = sscanf(s, '%f %f %f %f %f');
        g.zmaxis = a(1);
        g.sibry = a(3);

        %% Profiles: Fpol, Pressure, ff', P'
        %profile01: g.fpol: (FPOL) Poloidal current function in m-T, F = R*BT on flux grid
        cnt = 0;
        y = zeros(1, g.mw);
        while cnt < (g.mw)
            s = fgetl(fid);
            [a, c] = sscanf(s, '%f %f %f %f %f');
            y((1:c)+cnt) = a;
            cnt = cnt + c;
        end
        g.fpol = y;

    %profile02: g.pres: (PRES) Plasma pressure in nt / m2 on uniform flux grid
        cnt = 0;
        y = zeros(1, g.mw);
        while cnt < (g.mw)
            s = fgetl(fid);
            [a, c] = sscanf(s, '%f %f %f %f %f');
            y((1:c)+cnt) = a;
            cnt = cnt + c;
        end
        g.pres = y;

    %profile03: g.ffp: (FFPRIM)	FF'(psi) in (mT)2 / (Weber/rad) on uniform flux grid
        cnt = 0;
        y = zeros(1, g.mw);
        while cnt < (g.mw)
            s = fgetl(fid);
            [a, c] = sscanf(s, '%f %f %f %f %f');
            y((1:c)+cnt) = a;
            cnt = cnt + c;
        end
        g.ffp = y;

    %profile04: g.pp: (PPRIME) P'(psi) in (nt/m2) / (Weber/rad) on uniform flux grid
        cnt = 0;
        y = zeros(1, g.mw);
        while cnt < (g.mw)
            s = fgetl(fid);
            [a, c] = sscanf(s, '%f %f %f %f %f');
            y((1:c)+cnt) = a;
            cnt = cnt + c;
        end
        g.pp = y;

    %% Profile: poloidal flux on rectangular grid points
    %profile05: g.psirz: 
        cnt = 0;
        y = zeros(1, g.mw*g.mh);
        while cnt < (g.mw*g.mh)
            s = fgetl(fid);
            [a, c] = sscanf(s, '%f %f %f %f %f');
            y((1:c)+cnt) = a;
            cnt = cnt + c;
        end
        g.psirz = reshape(y, g.mw, g.mh);

    %% Profile: safety factor
    %profile06: safety factor on uniform flux grid  
    % g.qpsi
        cnt = 0;
        y = zeros(1, g.mw);
        while cnt < (g.mw)
            s = fgetl(fid);
            [a, c] = sscanf(s, '%f %f %f %f %f');
            y((1:c)+cnt) = a;
            cnt = cnt + c;
        end
        g.qpsi = y;

    %% boundary and limiter
        % g.nbdry,g.limitr   
        s = fgetl(fid);
        a = sscanf(s, '%d %d');
        g.nbdry = a(1);
        g.limitr = a(2);

    % boundary:
        cnt = 0;
        y = zeros(1, g.nbdry*2);
        while cnt < (g.nbdry*2)
            s = fgetl(fid);
            [a, c] = sscanf(s, '%f %f %f %f %f');
            y((1:c)+cnt) = a;
            cnt = cnt + c;
        end
        g.bdry = reshape(y, 2, g.nbdry);

    % limiter
        cnt = 0;
        y = zeros(1, g.limitr*2);
        while cnt < (g.limitr*2)
            s = fgetl(fid);
            [a, c] = sscanf(s, '%f %f %f %f %f');
            y((1:c)+cnt) = a;
            cnt = cnt + c;
        end
        g.lim = reshape(y, 2, g.limitr);

    % kvtor,rvtor,nmass   
        s = fgetl(fid);
        a = sscanf(s, '%d %f %d');
        g.kvtor = a(1);
        g.rvtor = a(2);
        g.nmass = a(3);

        fclose(fid);

        dR = g.xdim/(g.mw-1);
        dZ = g.zdim/(g.mh-1);
        g.r = g.rgrid1 + (0:g.mw-1)*dR; % R grid coordinates
        g.z = g.zmid - 0.5*g.zdim + (0:g.mh-1)*dZ; % R grid coordinates

        g.lim = g.lim;
        g.bdry = g.bdry;
        g.r = g.r;
        g.z = g.z;

        save(sprintf('%s.mat',fname), 'g');
    else
        fprintf('%s file does not exist\n',fname)
    end
end

%%
if exist('verbose', 'var') && verbose == 1
    figure('position', [100 100 400 800], 'color', [1 1 1]); 
    plot(g.lim(1,:), g.lim(2,:), '-b'); hold on;
    plot(g.bdry(1,:), g.bdry(2,:), '-r', 'LineWidth',2);

    contour(g.r, g.z, g.psirz', 135, ':k');

    axis equal;
    xlim([min(g.lim(1,:)), max(g.lim(1,:))]);
    ylim([min(g.lim(2,:)), max(g.lim(2,:))]);
    set(gca, 'fontsize', 15)
    xlabel('R [m]', 'fontsize', 15)
    ylabel('z [m]', 'fontsize', 15)
end

if exist('verbose', 'var') && verbose == 2
    figure('color',[1 1 1]); 
    subplot(2,2,[1,3]);
    plot(g.lim(1,:), g.lim(2,:), '-b'); hold on;
    plot(g.bdry(1,:), g.bdry(2,:), '-r', 'LineWidth',2);

    vcdelta = (g.sibry-g.simag)/(g.mw-1)*2;
    vc = g.simag:vcdelta:g.sibry;
    % [C, h]= contour(g.r, g.z, g.psirz', vc, ':k');
    [C, h]= contour(g.r, g.z, g.psirz', 50);
%     clabel(C,h);

    axis equal;
    xlim([min(g.lim(1,:)), max(g.lim(1,:))]);
    ylim([min(g.lim(2,:)), max(g.lim(2,:))]);
    fprintf('%s\n',fname);
    set(gca, 'FontSize', 14);

    
    subplot(2,2,2);
%     psi = g.psirz((g.mh+1)/2, :);
    psi = g.psirz(:,(g.mh+1)/2);
    bp = zeros(1, g.mw);
    for i=1:(g.mw-1)
        bp(i) = (psi(i+1) - psi(i))/(g.r(i+1)-g.r(i))*2/(g.r(i+1)+g.r(i));
    end
    plot(g.r, bp);
    xlim([min(g.lim(1,:)), max(g.lim(1,:))]); 
    xlabel('R (m)', 'FontSize', 14);
    ylabel('Bpol (T)', 'FontSize', 14);
    set(gca, 'FontSize', 14);
    
    subplot(2,2,4);
    bz = abs(g.bcentr)*g.rmaxis./g.r;
    alpha = bp.^2 ./ (bz.^2 + bp.^2);
    plot(g.r, g.r.*(sqrt(1+alpha)-1)*100);
    xlim([min(g.lim(1,:)), max(g.lim(1,:))]); 
    xlabel('R (m)', 'FontSize', 14);
    ylabel('Delta R (cm)', 'FontSize', 14);
    set(gca, 'FontSize', 14);
    
end



