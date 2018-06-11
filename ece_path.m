function [fs, dz, Rpi, zpi, thetai] = ece_path(R2d, z2d, F_Br, F_Bz, F_Bt, F_B, opt, dn, fidx, zidx)

%% global parameters

global e
global me
global Nz
global Nf
global zstart
global zend
global fstart
global fend
global pstart
global pend
global pint
global verbose

%% characteristic frequencies 

hm = opt.harmonic(dn); % O1 or X2

wce = @(x,y) e*F_B(x,y)/me; % [rad/s]

%% determine zs, as, fs of a selected single channel (zidx, fidx)

Rinit = 2.39; % [m]
dz = linspace(zstart,zend,Nz); % [mm] effective vertical range at mini lens 
zs = zeros(size(dz)); % [m]
as = zeros(size(dz)); % [rad]
for i = 1:length(dz)
    [zR, aR] = vac_path_now(dn, opt.LensFocus(dn), opt.LensZoom(dn), Rinit, dz(i)); % find vertical beam position and angle at initial point R = 2.32 m [m] [rad]
    zs(i) = zR(zidx);
    as(i) = aR(zidx);
end

fs = linspace((fidx-1)*0.9 + 2.6 + opt.LO(dn) + fstart, (fidx-1)*0.9 + 2.6 + opt.LO(dn) + fend, Nf); % [GHz]  

%% calculate

Rpi = cell(length(fs), length(zs));
zpi = cell(length(fs), length(zs));
thetai = cell(length(fs), length(zs));

% inside single channel 
for fn = 1:length(fs) 
    % find EC resonant position (cold)
    [fC] = contourc(R2d(1,:), z2d(:,1), wce(R2d,z2d)/(2*pi()*1e9)*hm, [fs(fn) fs(fn)]); 
    fR = fC(1,2:end); 
    fz = fC(2,2:end);
    
    if verbose
        plot(fR, fz, 'k'); hold all; % plot
    end
    
    for vn = 1:length(zs) % low:1 ~ high:24
        % TORBEAM  % write inbeam.dat for TORBEAM with beam parameters 
        write_inbeam(hm, fs(fn)*10^9, as(vn)/pi*180, zs(vn)*100); % [GHz, angle, cm]
        
        % run TORBEAM
        system('ssh -p 2201 mjchoi@172.17.250.11 /home/users/mjchoi/torbeam_ifortran/run.sh');

        % get outfile for beam path and width
        system('scp -P 2201 mjchoi@172.17.250.11:~/torbeam_ifortran/run_torbeam/t1_LIB.dat ./');
        beampath = dlmread(fullfile('./', 't1_LIB.dat'));
        cpath = beampath(:,1:2)/100; % for central position [m]
        Rp = cpath(:,1)'; % beam path starts from (calculation start point(xxb) - 1) cm 
        zp = cpath(:,2)';
        theta = zeros(size(Rp));            
        for i=2:length(Rp) % calculate angle between emission path and B-field
            Rvec = [-(Rp(i)-Rp(i-1)), -(zp(i)-zp(i-1)), 0]; % opposite direction for emission path                
            Bvec = [F_Br(Rp(i),zp(i)), F_Bz(Rp(i),zp(i)), F_Bt(Rp(i),zp(i))];
            theta(i) = acos((Bvec*Rvec')/(sqrt(sum(Bvec.^2))*sqrt(sum(Rvec.^2)))); % [rad]                
        end
        
%         plot(Rp, zp, 'LineWidth',1, 'color','k'); drawnow % plot
        
        % find a cross with cold resonance
        P = InterX([Rp;zp],[fR;fz]); % [m]
        
        if verbose
            plot(P(1), P(2), 'x'); hold all; %%%% PLOT %%%%
        end

        idx2 = find(Rp >= P(1) + pstart, 1, 'last'); 
        idx1 = find(Rp >= P(1) + pend, 1, 'last'); 
                
        nidx = idx2:-pint:idx1; 
        Rpi{fn,vn} = interp1(idx2:-1:idx1, Rp(idx2:-1:idx1), nidx);
        zpi{fn,vn} = interp1(idx2:-1:idx1, zp(idx2:-1:idx1), nidx);
        thetai{fn,vn} = interp1(idx2:-1:idx1, theta(idx2:-1:idx1), nidx);        

        fprintf('# %d - %d sub_ray done\n', fn, vn)     
    end
end 


end