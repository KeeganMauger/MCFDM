% Initialization
clear all
close all
figure(2)
set(0,'DefaultFigureWindowStyle','docked')
set(0,'defaultaxesfontsize',10)
set(0,'defaultaxesfontname','Times New Roman')
set(0,'DefaultLineLineWidth', 0.5);



global C

C.q_0 = 1.60217653e-19;             % electron charge
C.hb = 1.054571596e-34;             % Dirac constant
C.h = C.hb * 2 * pi;                % Planck constant
C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.eps_0 = 8.854187817e-12;          % vacuum permittivity
C.mu_0 = 1.2566370614e-6;           % vacuum permeability
C.c = 299792458;                    % speed of light
C.g = 9.80665;                      % metres (32.1740 ft) per s²
C.m_n = 0.26 * C.m_0;               % effective electron mass
C.am = 1.66053892e-27;              % atomic mass unit
C.T = 300;
C.vth = sqrt(2*C.kb * C.T / C.m_n);


temp = C.T;
subplot(2,1,1);
rectangle('Position',[0 0 200e-9 100e-9])
hold on

%--------------------------------------------------------------------------
% Initializing Positions
%--------------------------------------------------------------------------


N = 10000;        % Number of electrons
i = 0;
j = 0;

for i=1:N
    px(i) = 0 + (200e-9 - 0).*rand(1,1);
    py(i) = 0 + (100e-9 - 0).*rand(1,1);
end

%--------------------------------------------------------------------------
% Thermal Velocity and Direction
%--------------------------------------------------------------------------

vth = C.vth;

for j=1:N
    vx(j) = (vth/sqrt(2))*randn();                           
    vy(j) = (vth/sqrt(2))*randn();
    vth_calc(j) = sqrt(vx(j)^2 + vy(j)^2);
end


t = 0;
T(1) = 0;
dt = 1e-14;     % time step

for l=1:N           %Scattering time step
    ndt(l) = dt;
end
P_scat = 0;
Tmn = 0.2e-12;

px_prev = 0;
py_prev = 0;
T_prev = 0;

sampleidx = randi(N,10,1);
for t=2:200
    for k=1:N
        
        P_scat(k) = 1 - exp(-(dt/Tmn));
        if P_scat(k) > rand()
            vx(k) = (vth/sqrt(2))*randn();
            vy(k) = (vth/sqrt(2))*randn();
        else
            ndt(k) = ndt(k) + dt;
        end
        
        px_prev(k) = px(k);
        px(k) = px(k) + vx(k)*dt;
        py_prev(k) = py(k);
        py(k) = py(k) + vy(k)*dt;
        
        if py(k) >= 100e-9 || py(k) <= 0
            vy(k) = -vy(k);
            if py(k) >= 100e-9
                py(k) = 100e-9;
            end
            if py(k) <= 0
                py(k) = 0;
            end
        end
        if px(k) >= 200e-9
            px(k) = 0;
            px_prev(k) = px(k);
        elseif px(k) <= 0
            px(k) = 200e-9;
            px_prev(k) = px(k);
        else
            px(k) = px(k);
        end
        
        v(k) = sqrt(vx(k)^2 + vy(k)^2);
        v2(k) = v(k).*v(k);
        
    end
    for h=1:length(sampleidx)
        subplot(2,1,1);
        plot([px_prev(sampleidx(h)) px(sampleidx(h))],[py_prev(sampleidx(h)) py(sampleidx(h))],'SeriesIndex',h)
        hold on 
    end
    
    KE = 0.5 * C.m_n * mean(v2);
    T_prev = T;
    T = KE / C.kb;
    subplot(2,1,2);
    plot([t-1 t], [T_prev T],'r')
    hold on
    
    pause(0.01)
end

T_T_C = mean(ndt);
MFP = mean(v)*mean(ndt);

subplot(2,1,1);
title('Maxwell-Boltzmann Velocity Distributions of Electrons with Scattering')
xlabel('Region Width (m)')
ylabel('Region Height (m)')
subplot(2,1,2);
title('Mean Temperature over Time')
xlabel('Timesteps (1e-14 s per step)') 
ylabel('Temperature (K)')
disp('Section 2')
fprintf('\nThe given mean free time between collisions is %e seconds.',Tmn)
fprintf('\nThe mean calculated time between collisions is %e seconds.',T_T_C)
fprintf('\nThe mean calculated thermal velocity is %e meters per second.',mean(vth_calc))
fprintf('\nThe thermal velocity is %e meters per second.',vth)
fprintf('\nThe mean free path is %e meters.\n\n\n',MFP)
pause(0.1)
saveas(gcf,'Figure2')

figure(3)
histogram(vth_calc,10)
title('Maxwell-Boltzmann Distribution of Electron Velocities')
ylabel('Number of Electrons in Velocity Bins')
xlabel('Velocity of Electrons (meters per second)')
saveas(gcf,'Figure3')
