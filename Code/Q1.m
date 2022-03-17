% Testing

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

figure(1)
temp = C.T;
subplot(2,1,1);
rectangle('Position',[0 0 200e-9 100e-9])
hold on

%--------------------------------------------------------------------------
% Initializing Positions
%--------------------------------------------------------------------------


N = 30000;        % Number of electrons
i = 0;
j = 0;

for i=1:N
    px(i) = 0 + (200e-9 - 0).*rand(1,1);
    py(i) = 0 + (100e-9 - 0).*rand(1,1);
end

%--------------------------------------------------------------------------
% Voltage Applied Across x-Dimension to Find Electric Field
%--------------------------------------------------------------------------

V0x = 0.1;
V0y = 0;
L = 200e-9;
W = 100e-9;
E0x = V0x / L;
E0y = V0y / W;
fMesh = 1;
nx = fMesh*200;
ny = fMesh*100;
% G = sparse(nx,ny);
% F = sparse(1,nx*ny);

La = linspace(0,L,nx);
Wa = linspace(0,W,ny);

Emapx = zeros(ny,nx);
Emapy = zeros(ny,nx);
for i = 1:width(La)
    for j = 1:width(Wa)
        Emapx(j,i) = E0x;
        Emapy(j,i) = E0y;
    end
end
%surf(La,Wa,Emapx)

Fex = abs(C.q_0*E0x);
aex = Fex/C.m_n;
Fey = abs(C.q_0*E0y);
aey = Fey/C.m_n;


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
vx_total = 0;
vy_total = 0;

sampleidx = randi(N,10,1);
Ix = 0;
Jx = 0;

for t=2:300
    vx_total = 0;
    vy_total = 0;
    for k=1:N
        
        P_scat(k) = 1 - exp(-(dt/Tmn));
        if P_scat(k) > rand()
            vx(k) = (vth/sqrt(2))*randn();
            vy(k) = (vth/sqrt(2))*randn();
        else
            ndt(k) = ndt(k) + dt;
        end
        
        px_prev(k) = px(k);
        px(k) = px(k) + vx(k)*dt + aex*dt^2;        % Adding acceleration
        vx(k) = vx(k) + aex*dt;
        py_prev(k) = py(k);
        py(k) = py(k) + vy(k)*dt + aey*dt^2;
        vy(k) = vy(k) + aey*dt;
        
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
        
        vx_total = vx_total + vx(k);        % Drift velocity x
        vy_total = vy_total + vy(k);        % Drift velocity y
    end
    vx_total_alt = sum(vx);
    for h=1:length(sampleidx)
        subplot(2,1,1);
        plot([px_prev(sampleidx(h)) px(sampleidx(h))],[py_prev(sampleidx(h)) py(sampleidx(h))],'SeriesIndex',h)
        hold on 
    end
    
    
    vx_drift = 1/N * vx_total;
    vy_drift = 1/N * vy_total;
    Jx = C.q_0 * 1e17 * vx_drift;
    Ix_prev = Ix;
    Ix = Jx * (W*L);
    
    
    subplot(2,1,2);
    plot([t-1 t], [Ix_prev Ix], 'r')
    hold on
    
    KE = 0.5 * C.m_n * mean(v2);
    T_prev = T;
    T = KE / C.kb;
%     subplot(2,1,2);
%     plot([t-1 t], [T_prev T],'r')
%     hold on
    
    pause(0.01)
end

T_T_C = mean(ndt);
MFP = mean(v)*mean(ndt);

subplot(2,1,1);
title('Maxwell-Boltzmann Velocity Distributions of Electrons with Scattering')
xlabel('Region Width (m)')
ylabel('Region Height (m)')
subplot(2,1,2);
title('Current over Time')
xlabel('Timesteps (1e-14 s per step)') 
ylabel('Current (A)')
% disp('Section 2')
% fprintf('\nThe given mean free time between collisions is %e seconds.',Tmn)
% fprintf('\nThe mean calculated time between collisions is %e seconds.',T_T_C)
% fprintf('\nThe mean calculated thermal velocity is %e meters per second.',mean(vth_calc))
% fprintf('\nThe thermal velocity is %e meters per second.',vth)
% fprintf('\nThe mean free path is %e meters.\n\n\n',MFP)
pause(0.1)
saveas(gcf,'Figure1')




E_map = [reshape(px,[N,1]),reshape(py,[N,1])];
figure(2)
hist3(E_map,'CDataMode','auto','FaceColor','interp')
%view(2)
title('Electron Density in Modeled Region')
xlabel('Region x-Axis (m)')
ylabel('Region y-Axis (m)')
zlabel('Electron Density')
saveas(gcf,'Figure2')

Nbins = 21;
d = 1;
u = 1;
vtm = 0;
vbm = 0;
vbm2 = 0;
T_map = zeros(Nbins);
[X,Xe] = discretize(px,Nbins);
[Y,Ye] = discretize(py,Nbins);

for e=1:Nbins
    for f=1:Nbins
        for g=1:N
            if X(g) == e && Y(g) == f
                vtm(d) = v(g);
                vtm2(d) = vtm(d).*vtm(d);
                d = d+1;
            end
        end
    vbm2 = mean(vtm2);
    d = 1;
    T_map(e,f) = (0.5*C.m_n*vbm2)/C.kb;
    end
end

[V,W] = meshgrid(0:1e-8:2e-7,0:0.5e-8:1e-7);
figure(3)
surf(V,W,T_map,'FaceColor','interp')
%view(2)
title('Temperature Map of Modeled Region')
xlabel('Region x-Axis (m)')
ylabel('Region y-Axis (m)')
zlabel('Temperature (K)')
saveas(gcf,'Figure3')


        
        


