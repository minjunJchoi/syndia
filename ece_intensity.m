function [Rm, zm, s, tau, jms, theta_max, Iece] = ece_intensity(Rp, zp, theta, Rc, omega, m, F_B, F_Te, F_ne)
% all mks units except Te
% Rp : R coordinates on beam path [m]
% zp : z coordinate on beam path [m]
% theta : angle between field and beam path [rad]
% omega : emission frequency [rad/s]
% m : harmonic number
% B : B field function in R-z coordinate [T]
% F_Te : 2d TriScatteredInterp Te function in R-z coordinates [J]
% F_ne : 2d TriScatteredInterp ne function in R-z coordinates [m^-3]
% @(x,y) : (R,z) coordinate

% it will calculate emission profile of frequency omega along s
% Rm : maximum emission position [m]
% zm : maximum emission position [m]
% Iece : ece intensity [W/m^2 Hz] 
% % all mks units, Temperature in [J]
global e
global me
global eps
global c
global mc2
global zeta
global mu
global eta

% 
N = length(Rp);

% characteristic frequencies
wce = @(x,y) e*F_B(x,y)/me; % [rad/s]
wpe = @(x,y) sqrt(F_ne(x,y)*e^2/(eps*me)); % [rad/s]

% characteristic functions
zeta = @(x,y) mc2./F_Te(x,y)/2;
mu = @(x,y,w) (w./(m.*wce(x,y))).^2;
eta = @(x,y,w,theta) 1 + mu(x,y,w).*(cos(theta)).^2;

% refractive index 
X = @(x,y) (wpe(x,y)./(m*wce(x,y))).^2;  
if m == 1 % O-mode
%     Nperp2C = @(x,y,w) 1 - (wpe(x,y)./w).^2; % B(3.1.22)    
    N1OCsq = @(x,y,theta) 1 - X(x,y).*(1-X(x,y))./(1 - X(x,y) - sin(theta).^2/(2*m^2) + sqrt((sin(theta).^2/(2*m^2)).^2 - (1 - X(x,y)).*cos(theta).^2/m^2));
    N1OCre = @(x,y,theta) real(sqrt(N1OCsq(x,y,theta)));    
elseif m == 2 % X-mode
%     Nperp2C = @(x,y,w) 1 - (wpe(x,y)./w).^2.*(w.^2 - wpe(x,y).^2)./(w.^2 - wce(x,y).^2 - wpe(x,y).^2); % B(3.1.12)
    N2XCsq = @(x,y,theta) 1 - X(x,y).*(1-X(x,y))./(1 - X(x,y) - sin(theta).^2/(2*m^2) - sqrt((sin(theta).^2/(2*m^2)).^2 - (1 - X(x,y)).*cos(theta).^2/m^2)); % H(5.2.48)
    N2XCre = @(x,y,theta) real(sqrt(N2XCsq(x,y,theta)));    
end
    
% absorption coefficient
if m == 1 % O-mode
    amO = @(x,y,theta) pi/(2*c).*wpe(x,y).^2.*N1OCre(x,y,theta).*(1 + 2*cos(theta)^2)^2*sin(theta)^4/(1 + cos(theta)^2)^3.*F_Te(x,y)/mc2; % H(5.2.52) %
    
    amRz = @(x,y,theta) amO(x,y,theta); 
elseif m == 2 % X-mode
    am = @(x,y) e^2*F_ne(x,y)/(4*c*me*eps) .* m^(2*m-1)/factorial(m-1) .* (F_Te(x,y)/(2*mc2)).^(m-1); % H(5.2.39)

    a2sq = @(x,y,theta) (1 + (1 - X(x,y)).*N2XCre(x,y,theta).^2*cos(theta)^2./(1 - X(x,y) - N2XCre(x,y,theta).^2*sin(theta)^2).^2 ...
        *m^2.*(1 - (m^2 -1)/m^2./X(x,y).*(1 - N2XCre(x,y,theta).^2)).^2).^2*sin(theta)^2;
    b2sq = @(x,y,theta) (1 + (1 - X(x,y))./(1 - X(x,y) - N2XCre(x,y,theta).^2*sin(theta)^2) ...
        *m^2.*(1 - (m^2 -1)/m^2./X(x,y).*(1 - N2XCre(x,y,theta).^2)).^2).^2*cos(theta)^2;
    eta2X = @(x,y,theta) N2XCre(x,y,theta).^(2*m-3)*(m - 1)^2.*(1 - (m+1)/m./X(x,y).*(1 - N2XCre(x,y,theta).^2)).^2 ...
        ./((1 + cos(theta)^2).*(a2sq(x,y,theta) + b2sq(x,y,theta)).^(1/2)); % H(5.2.47)    
    
    amX = @(x,y,theta) am(x,y) .* eta2X(x,y,theta); % H(5.2.54)

    amRz = @(x,y,theta) amX(x,y,theta); 
end

% shape Maxwellian (relativistic + Doppler)
ba1 = @(x,y,w,theta) (mu(x,y,w)*cos(theta) - sqrt(1 - mu(x,y,w)*sin(theta)^2))./(1 + mu(x,y,w)*cos(theta)^2);
ba2 = @(x,y,w,theta) (mu(x,y,w)*cos(theta) + sqrt(1 - mu(x,y,w)*sin(theta)^2))./(1 + mu(x,y,w)*cos(theta)^2);

shape = @(x,y,w,theta) 2*pi*zeta(x,y).^(7/2).*w./(sqrt(pi).*(m*wce(x,y)).^2) ...
    .*integral(@(ba) (1 - ba.^2 - (1 - ba*cos(theta)).^2.*mu(x,y,w)).^2 ...
    .*(1 - ba.*cos(theta)).*exp(-zeta(x,y).*(1 - (1 - ba.*cos(theta)).^2.*mu(x,y,w))), ...
    ba1(x,y,w,theta), ba2(x,y,w,theta));

% absorption coefficent along beam path s (from cold position Rc + 3 cm to Rc - 12 cm)
ams = zeros(1,N);
for i=1:N
%     if F_Te(Rp(i), zp(i))/(1.602*10^-19) <= 100; % less than 100 eV %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% consider all range
%         ams(i) = 0;
%         fprintf('Te is less than 100 eV, so neglected \n');
%     else    
%         ams(i) = real(amRz(Rp(i), zp(i), theta(i)))*real(shape(Rp(i), zp(i), omega, theta(i)));
%     end
    if mu(Rp(i),zp(i),omega)*sin(theta(i))^2 < 1 % the integral blows up with the imaginary beta
        ams(i) = real(amRz(Rp(i), zp(i), theta(i)))*real(shape(Rp(i), zp(i), omega, theta(i)));
    else
        ams(i) = 0;
    end
end

% %% interpolation for a better accuracy in intensity //
% Rp0 = Rp; zp0 = zp;
% Rp = Rp0(1):0.00001:Rp0(end); %%%%%%%%%%%%%%%%%%%%% Rp resolution 
% zp = interp1(Rp0, zp0, Rp);
% theta = interp1(Rp0, theta, Rp);
% ams = interp1(Rp0, ams, Rp);
% N = length(Rp);

% emission profile along beam path s [Rp, zp]
Ibb = @(x,y,w) (w/(2*pi*c)).^2.*F_Te(x,y); % H(5.2.37)

% define path from the inside
s = zeros(1,N);
ds = zeros(1,N);
for i=2:N
    ds(i) = sqrt((Rp(i) - Rp(i-1))^2 + (zp(i) - zp(i-1))^2);
    s(i) = s(i-1) + ds(i);
end

% calculate optical depth and transfer intensity
tau = trapz(s,ams) - cumtrapz(s, ams);
Is = zeros(1,N);
for i=1:N
    if i > 2
        Is(i) = Is(i-1) + (ams(i-1)*Ibb(Rp(i-1),zp(i-1),omega) - ams(i-1)*Is(i-1))*ds(i);
%         Is(i) = Is(i-1) + (ams(i-1)*Ibb(Rp(i-1),zp(i-1),omega) - ams(i-1)*Is(i-1))/N2XCsq(Rp(i-1),zp(i-1),theta(i-1))*ds(i); %%%%%%%%%%%%%%%%%% Nsq effect on intensity
    end 
end
% Is = Is.*N2XCsq(Rp,zp,theta); %%%%%%%%%%%%%%%%%% Nsq effect on intensity

jms = ams.*Ibb(Rp,zp,omega).*exp(-tau); % emissivity after reabsorption. B(2.2.13), B(2.2.15)

% return maximum emissivity position
midx = jms == max(jms);
Rm = Rp(midx);
zm = zp(midx);
theta_max = theta(midx);
if length(Rm) ~= 1
    Rm = 0;
    zm = 0;
    theta_max = pi/2;
end

% total intensity at the outside 
Iece = sum(jms.*ds);

% % plot
% figure;
% plot(s, jms)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%% save image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h = figure; 
% subplot(3,1,1)
% Cidx = find(Rp <= Rc, 1, 'last'); 
% plot(s-s(Cidx),jms./exp(-tau)); hold on; % pure
% plot(s-s(Cidx),jms); % after reabsorption
% % plot(s-s(Cidx),ams) % absorption coefficient
% % xlim([-0.03, 0.12])
% % title(sprintf('jms, theta max min = %g %g', max(theta)/pi*180, min(theta)/pi*180))
% subplot(3,1,2)
% plot(s - s(Cidx), Ibb(Rp,zp,omega)/(1e-11)); hold all % IBB from Te [1e-11 W/msq Hz] 
% % plot(s - s(Cidx), Ibb(Rp,zp,omega).*N2XCsq(Rp,zp,theta)/(1e-11)); hold all % IBB from Te [1e-11 W/msq Hz] 
% plot(s - s(Cidx), Is/(1e-11)); hold all % Is 
% plot(s(Cidx) - s(Cidx), Ibb(Rp(Cidx),zp(Cidx),omega)/1e-11,'o') % cold position IBB 
% plot(s(midx) - s(Cidx), sum(jms.*ds)/1e-11, 'x') % measured intensity
% 
% subplot(3,1,3)
% % plot(Rp, mu(Rp,zp,omega))
% plot(Rp, zp)
% % plot(s-s(Cidx),theta);
% % title('theta')
% % plot(s - s(Cidx), Rp); hold all
% % plot(0, Rc, 'o')
% 
% % export as an image
% f = getframe(h);
% [im, ~] = frame2im(f); 
% filename = sprintf('C:/data/ecei_data/016150/jms_%g_%g.png', Rp(1), zp(1));
% imwrite(im, filename, 'png')
% close(h) 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%% save image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% function val = shape(x, y, theta, w, m, F_B, F_Te, F_ne)
% % Maxwellian velocity distribution
% global e
% global me
% global c
% global mc2
% global zeta
% global mu
% global eta
% global alp1
% global alp2
% 
% % characteristic frequencies
% wce = @(x,y) e*F_B(x,y)/me; % [rad/s]
% 
% % shape condition
% Cth = abs(cos(theta)); 
% vth = sqrt(F_Te(x, y)/me);
% 
% 
% % % if Cth/(vth/c) > 1 && abs(theta/pi*180-90) >= 3 % H(5.2.19)
% % if abs(theta/pi*180-90) >= 3 % H(5.2.19)    
% %     % shape 2 Appendix equation II in Rathgeber PPCF2013
% %     val = 2*pi*zeta(x,y).^(7/2).*w./(2*sqrt(pi).*(m*wce(x,y)).^2.*(cos(theta)^5)) ...
% %         .*integral(@(alp) (sin(theta)^4 - 4.*sqrt(alp)*sin(theta)^2 + 2.*alp ...
% %         .*(2 + eta(x,y,w,theta).*sin(theta)^2) - 4.*eta(x,y,w,theta).*alp.^(3/2) + eta(x,y,w,theta).^2.*alp.^2).*exp(-zeta(x,y).*(1 - mu(x,y,w).*alp)), ...
% %         alp2(x,y,w,theta), alp1(x,y,w,theta)); % 2*pi factor for rad/s -> 1/s; 
% %     fprintf('shape 2 angle = %g, Te = %g, ne = %g, Cth/(vth/c) = %g \n', theta/pi*180, F_Te(x,y), F_ne(x,y), Cth/(vth/c))
% % else
% %     % shape 1 H(5.2.30)
% %     val = sqrt(pi)*w./(m*wce(x,y)).^2 .* (2*mc2./F_Te(x,y)).^(m+3/2) * factorial(m)/factorial(2*m+1)...
% %         .*(1-(w./(m*wce(x,y))).^2).^(m+1/2) .* exp(-mc2./(2*F_Te(x,y)).*(1-(w./(m*wce(x,y))).^2));    
% %     fprintf('shape 1 angle = %g, Te = %g, ne = %g, Cth/(vth/c) = %g \n', theta/pi*180, F_Te(x,y), F_ne(x,y), Cth/(vth/c))
% % end
% 
% if abs(theta/pi*180-90) < 3 % H(5.2.19)
% % if Cth < beta % H(5.2.19)
%     % shape 1 H(5.2.30)
%     val = sqrt(pi)*w./(m*wce(x,y)).^2 .* (2*mc2./F_Te(x,y)).^(m+3/2) * factorial(m)/factorial(2*m+1)...
%         .*(1-(w./(m*wce(x,y))).^2).^(m+1/2) .* exp(-mc2./(2*F_Te(x,y)).*(1-(w./(m*wce(x,y))).^2));
% %     % check shape 1 
% %     domega = 0.01*10^9*2*pi;
% %     oaxis = ((omega-10*10^9*2*pi):domega:(omega+3*10^9*2*pi));
% %     ts = real(shape(mean(Rp)*ones(size(oaxis)), mean(zp)*ones(size(oaxis)), oaxis)); %% Rp and zp should be resonant position of omega
% %     % figure; plot(oaxis, ts); hold on; plot(omega,0,'o')
% %     fprintf('angle = %g // shape 1 integration = %g \n', theta/pi*180, sum(ts*domega)/(2*pi));
% %     fprintf('shape 1 angle = %g \n', theta/pi*180)
% elseif abs(theta/pi*180-90) >= 3
% % elseif Cth > beta
%     % shape 2 Appendix equation II in Rathgeber PPCF2013
% 
%     val = 2*pi*zeta(x,y).^(7/2).*w./(sqrt(pi).*(m*wce(x,y)).^2) ...
%         .*integral(@(ba) (1 - ba.^2 - (1 - ba*cos(theta)).^2.*mu(x,y,w)).^2 ...
%         .*(1 - ba.*cos(theta)).*exp(-zeta(x,y).*(1 - (1 - ba.*cos(theta)).^2.*mu(x,y,w))), ...
%         ba1(x,y,w,theta), ba2(x,y,w,theta));
%     
%     % test 1
% %     val = 2*pi*zeta(x,y).^(7/2).*w./(2*sqrt(pi).*(m*wce(x,y)).^2.*(cos(theta)^5)) ...
% %         .*integral(@(alp) (sin(theta)^4 - 4.*sqrt(alp)*sin(theta)^2 + 2.*alp ...
% %         .*(2 + eta(x,y,w,theta).*sin(theta)^2) - 4.*eta(x,y,w,theta).*alp.^(3/2) + eta(x,y,w,theta).^2.*alp.^2).*exp(-zeta(x,y).*(1 - mu(x,y,w).*alp)), ...
% %         alp2(x,y,w,theta), alp1(x,y,w,theta)); % 2*pi factor for rad/s -> 1/s; 
%     
% %     F = @(x) exp(-x.^2).*integral(@(y) exp(y.^2), 0, x); % integral form of dawson function
% % %     F = @(x) 1/2*sqrt(pi).*exp(-x.^2).*erfi(x); % alternative form of dawson function need symbolic math toolbox
% % %     k = (0:100);
% % %     F2 = @(x) sum((-1).^k.*factorial(k).*2.^(2*k).*x.^(2*k+1)./factorial(2*k+1)); % approximated form McCabe Math.Comput 1974. diverges with x
% %     rshape = @(x,y,w,alp) zeta(x,y).^(3/2)./(sqrt(pi)*w.*mu(x,y,w)*(cos(theta)^5))...
% %         .*exp(-zeta(x,y).*(1-mu(x,y,w).*alp))...
% %         .*(-(sin(theta))^2./sqrt(alp) + eta(x,y,w,theta).^2./(zeta(x,y).*mu(x,y,w))...
% %         - 1./(sqrt(zeta(x,y).*mu(x,y,w))).*(3*eta(x,y,w,theta) - 2*zeta(x,y).*mu(x,y,w)*(sin(theta)^2)).*F(sqrt(zeta(x,y).*mu(x,y,w).*alp)));
% %     shape = @(x,y,w) 2*pi*(rshape(x,y,w,alp1(x,y,w,theta)) - rshape(x,y,w,alp2(x,y,w,theta)));
%     
% %     % check shape 2
% %     domega = 0.01*10^9*2*pi;
% %     oaxis = ((omega-10*10^9*2*pi):domega:(omega+3*10^9*2*pi));
% %     ts = zeros(size(oaxis));
% %     for i = 1:length(oaxis)
% %         ts(i) = shape(Rp(Cidx),zp(Cidx),oaxis(i)); %% Rp and zp should be resonant position of omega
% %     end
% %     figure; plot(oaxis, real(ts)); hold on; plot(omega,0,'o')    
% %     fprintf('angle = %g // shape 2 integration = %g \n', theta/pi*180, sum(ts*domega)/(2*pi));
% %     fprintf('shape 2 angle = %g \n', theta/pi*180)
% end
