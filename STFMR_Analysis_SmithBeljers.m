%% STFMR model
clearvars
% close all


% Here, an in-plane magnetization thin film of sample_length (y axis) 
% and sample_width (x axis) is placed in an external magnetic field which
% is apllied in the x-y plane, its orientation is determined by phi_Hext
% where theta_Hext is contant at the x-y plane and slightly tilted. 
% the axis are determined by the sample itself, where the (0,0,0) is the
% middle of the sample. Current is transported in the y axis, 
% This makes h_Oe, S at a constant x orientation. 
% The only vector that changes orientation is H_ext.
%% physical constants
q=1.60217646e-19;                        %[C]
miu0=4*pi*1e-7;                          %[H/m]
hbar=6.62606885e-34/2/pi;                %[J*sec]
gamma_divide_by_2pi=28.024;              %[GHz/T]    TODO - is this with a minus sign?
gamma=gamma_divide_by_2pi*2*pi;          %Giga*[(rad/s)/T]
rho_py = 35e-8;                          %[ohm*m] for 7.5nm Py

%% Do Plots
plot_Mx_My_Mz = 0;
plot_Lorentzians = 1;
plot_linewidth_vs_current = 1;
plot_f_vs_h0 = 0;
plot_dH_vs_f = 0;

%% Frequency loop parameters:
w_loop = 0;

w_vec = 0;
if w_loop
    plot_f_vs_h0 = 1;
    plot_linewidth_vs_current = 0;
    w_vec = 2 * pi * linspace(0.1, 1, 9); % [rad/s]
end

%% DC Current Loop
Ic_dc_vec=(1e-3)*linspace(-30,30,7); % [A]
% Ic_dc_vec = 0e-3;
% 
%% sample parameters
Theta_sh = 0.09;        % conversion efficiency of charge to spin current
Ms_Oe = 661;            % [Oe]
Ms_Am = Ms_Oe/4/pi*1e3; % [A/m]
Ms_T = Ms_Am*miu0;      % [T]
alpha = 0.019;          % damping factor
Ku = 0;                 % Anisotropy const [J / m^3]

t_FM = 7.5e-9;          % [m]
t_NM = 5e-9;           % [m]
sample_width = 25e-6;   % [m]
sample_length = 55e-6;  % [m]

% Experiment parameters
f = 6;        %[GHz]
w = f*2*pi;

%% AC current
ZL= 64.7;     % load resistance [ohm]
Z0=50;        %[Ohm*m] Transmission line impedance.

Z_Py = rho_py * sample_length / (sample_width * t_FM); % Permalloy resistance.
I_ac_power=20; %[dBm]
gma_reflection=abs((ZL-Z0)./(ZL+Z0));  %reflectivity for voltage coefficient  at z=0 of the transmission line (Pozar book pp 59).
Pmw=10.^(I_ac_power./10);    % [mW]
V_rms = sqrt(Pmw*1e-3*Z0);   % [V]
V_peak = V_rms*sqrt(2);      %[V] Incident amplitude (Vp) from the source.
ratio_of_currents = 1 - ZL/Z_Py;  % the relative amount of current going through the OHE layer.
I_ac=V_peak.*(1-gma_reflection)/Z0 * ratio_of_currents;    % [A] Current through load.
Jc_ac=I_ac/(sample_width*t_NM);   %[A/m^2]   Current density throguh OHE layer.

%% Fields
Hext_vec=linspace(40e-3,200e-3,2^12); %[T] If simulation shows wierd resonances at the wrong h0, limit the field to where you expect the resonance field.

z=t_NM+0.5*t_FM;
h_Oe_AC = -Jc_ac/(4*pi).* (F(-sample_width/2,z)-F(sample_width/2,z)-...
        F(-sample_width/2,z-t_NM)+F(sample_width/2,z-t_NM))*miu0; % from amir matlab code, [T] units

% Current is in y direction so h_Oe_AC in x direction:
h_Oe_AC_x = h_Oe_AC;
h_Oe_AC_y = 0;
h_Oe_AC_z = 0;

Hstt_AC=(hbar*Theta_sh*Jc_ac)/(2*q*Ms_Am*t_FM);  %[T] spin Hall effect field parameter 

% Demagnetizing Field
shape_type=3;% sphere=1, rod=2, disk=3
[Ny, Nx, Nz] = shape_anisotropy(sample_length,sample_width,t_FM,shape_type); %Nx+Ny+Nz=1 [SI] units
% Nx = 0.3;
% Ny=0;
% Nz = 0.6;


%% System Orientation - This is setup like the sample 
% External field
theta_Hext=pi/2+0*(pi/180); %[rad] if =pi/2 DC Hext is in x-y plane.
phi_Hext=(pi/180)*45; %[rad] if =0, DC Hext is along X axis, in this regime the simulation breaks.

amr_factor = pi/4;

% S Orientation
Sx=1;
Sy=0;
Sz=0; %|[Sx,Sy,Sz]| have to be 1

%% DC Current
Ic_effective = Ic_dc_vec * ratio_of_currents;
Jc_dc_vec=Ic_effective./(sample_width*t_NM);   %[A/m^2] 
%% Loops
res_field_vec_currents=zeros(1,length(Ic_dc_vec)); % [T] this is H0 - resonance field
res_field_vec_freqs = zeros(1,length(w_vec));

res_width_vec_currents=zeros(1,length(Ic_dc_vec)); % [T] resonance linewidth
res_width_vec_freqs=zeros(1,length(w_vec)); % [T] resonance linewidth

% Frequency Loop
if w_loop
    Ic_dc_vec = 0;
    Jc_dc_vec = 0;
end
for l=1:length(w_vec)
    if w_loop
        w = w_vec(l);
        disp("------------- " + l + "/" + length(w_vec) + " -------------")
    end

    % DC current Loop
    for n=1:length(Ic_dc_vec)    
        Hstt_DC=-(hbar*Theta_sh*Jc_dc_vec(n))/(q*Ms_Am*t_FM);  %[T] spin Hall effect field parameter 
        
        % Oersted field from DC current
        H_Oe_DC = -Jc_dc_vec(n)/(4*pi).* (F(-sample_width/2,z)-F(sample_width/2,z)-...
            F(-sample_width/2,z-t_NM)+F(sample_width/2,z-t_NM))*miu0; % from amir matlab code, [T] units
        H_Oe_cartesian_vec=[H_Oe_DC,0,0];
    
        % Magnetisation calculation
        delta_phi_vec=zeros(1,length(Hext_vec));
        delta_R_vec=zeros(1,length(Hext_vec));
        phase_M_vec=zeros(1,length(Hext_vec));

        delta_Mx_vec = zeros(1,length(Hext_vec));
        delta_My_vec = zeros(1,length(Hext_vec));
        delta_Mz_vec = zeros(1,length(Hext_vec));
    
        % Field sweep: For each Hext data point, calculates total H = H_ext +
        % H_Oe, only DC components - to calculate steady state of magnetization
        % M0 = Mx0, My0, Mx0.
        test_phi0 = zeros(1,length(Hext_vec));
        test_theta0 = zeros(1,length(Hext_vec));
        test_Mx = zeros(1,length(Hext_vec));
        test_My = zeros(1,length(Hext_vec));
        test_Mz = zeros(1,length(Hext_vec));
        

        %% Magnetic Field Loop
        for k=1:length(Hext_vec)
            Hext=Hext_vec(k);
            Hext_cartesian_vec=Hext*[sin(theta_Hext)*cos(phi_Hext), sin(theta_Hext)*sin(phi_Hext), cos(theta_Hext)];
            H_cartesian=H_Oe_cartesian_vec+Hext_cartesian_vec;

            Hx = H_cartesian(1);
            Hy = H_cartesian(2);
            Hz = H_cartesian(3);


            phi_H = atan(Hy/Hx);
            theta_H = atan(sqrt(Hx^2 + Hy^2)/Hz);
            H0 = sqrt(Hx^2 + Hy^2 + Hz^2); % [T]

            H0_Am = H0/miu0; % [A/m]
            
            if k == 1
                initial_guess = [pi/2+1*pi/180, pi/2-1*pi/180]; % initial guess for optimization: [pi/2, pi/4] is for m in x-y plane at 45 degrees.
                % start at y.
            else
                initial_guess = [theta0, phi0]; % Update guess from last iter.
            end
            % The solver is always going 'downhill' for example for Hext in phi=0 and initial
            % guess is phi=89 degrees, it will converge to phi=0, but if
            % initial guess is phi=91 degrees it will converge to phi=180 -
            % the unstable soluton.
            [theta0, phi0] = find_steady_state(H0_Am, theta_Hext, phi_Hext, Ms_Am, Nx, Ny, Nz, Sx, Sy, Sz, Hstt_DC, Ku, initial_guess);
            % disp("Theta0: " + theta0 * 180/pi + " deg, Phi0: " + phi0 * 180/pi + " deg")
            Mx0 = sin(theta0)*cos(phi0);
            My0 = sin(theta0)*sin(phi0);
            Mz0 = cos(theta0);

            h_Oe_theta = h_Oe_AC * cos(theta0) * cos(phi0);
            h_Oe_phi = -h_Oe_AC * sin(phi0);
            [delta_theta, delta_phi] = find_perturbations(Ku, miu0, gamma, w, Ms_Am*4*pi, H0_Am, alpha, theta0, theta_H, phi0, phi_H, Nx, Ny, Nz, Sx, Sy, Sz, Hstt_DC, Hstt_AC, h_Oe_theta, h_Oe_phi);


            % ---------------------------- TEST --------------------------
            % E_Zeeman_theta = -miu0 * Ms_Am * H0_Am * (cos(theta0) * sin(theta_H) * cos(phi0 - phi_H) - sin(theta0) * cos(theta_H)) % [J / m^3]
            % E_Demag_theta = 0.5*miu0 * Ms_Am^2 * sin(2*theta0) * (Nx * cos(phi0)^2 + Ny * sin(phi0)^2 - Nz) % [J / m^3]
            % E_SHE_theta = Ms_Am * Hstt_DC * (Sx * sin(phi0) - Sy * cos(phi0)) % [J / m^3]
            % E_Crystalline_theta = Ku * sin(2 * theta0) % [J / m^3]


            % E_Zeeman_phi = miu0 * Ms_Am * H0_Am * sin(theta0) * sin(theta_H) * sin(phi0 - phi_H) % [J / m^3]
            % E_Demag_phi = 0.5*miu0 * Ms_Am^2 * sin(theta0)^2 * sin(2*phi0) * (-Nx + Ny) % [J / m^3]
            % E_SHE_phi = Ms_Am * sin(theta0) * Hstt_DC * (Sx * cos(theta0) * cos(phi0) - Sy * cos(theta0) * sin(phi0) + Sz * sin(theta0)) % [J / m^3]

            %-------------------------------------------------------------
    
            delta_phi_vec(k) = delta_phi;
            delta_R_vec(k) = AMR_factor(amr_factor)*delta_phi;
    
            test_phi0(k) = phi0;
            test_theta0(k) = theta0;
            test_Mx(k) = Mx0;
            test_My(k) = My0;
            test_Mz(k) = Mz0;
    
        end
    
        V_STFMR_vec=0.5*delta_R_vec*I_ac;
        
        % The initial guess of M along y axis makes the first few iterations
        % invalid.
        num_sample_to_reject = 35; % This is a factor of H resolution, high resolution demands rejecting more points.
        valid_H = Hext_vec(num_sample_to_reject:end);
        valid_V = V_STFMR_vec(num_sample_to_reject:end);
 
        param=fit_func_STFMR_ranen(valid_H,valid_V);
        res_width_vec_currents(n) = param(1); % [T]
        res_width_vec_freqs(l) = param(1);

        res_field_vec_currents(n) =  param(2); % [T] this is H0 - resonance field
        res_field_vec_freqs(l) = param(2);
        
        disp("------------- \phi_{H_{ext}} = " + phi_Hext*180/pi + ", Ic = "+Ic_dc_vec(n)*1e3 + " [mA], f = " + w/2/pi + " [GHz] -------------")
        disp("\H_{res} = " + param(2)*1e4 + " [Oe], \DeltaH = " + param(1)*1e4 + " [Oe], S,A,C = " + param(3) + ", " + param(4) + ", " + param(5))
        
        if plot_Mx_My_Mz
            figure()
            plot(Hext_vec, test_phi0, 'LineWidth', 3, 'DisplayName', "phi0")
            hold on
            plot(Hext_vec, test_theta0, 'LineWidth', 3, 'DisplayName', "theta0")
            xline(Ms_T, 'DisplayName', "Ms")
            plot(Hext_vec, test_Mx, 'LineWidth', 3, 'DisplayName', "Mx0")
            plot(Hext_vec, test_My, 'LineWidth', 3, 'DisplayName', "My0")
            plot(Hext_vec, test_Mz, 'LineWidth', 3, 'DisplayName', "Mz0")
            legend
            xlabel('H_{ext} [T]')
            hold off
        end
    
        if plot_Lorentzians
            figure()
            ft = param(3).*(((0.5*param(1)).^2)./((valid_H-param(2)).^2+(0.5*param(1)).^2))+param(4).*(((0.5*param(1))^2)./((valid_H-param(2)).^2+(0.5*param(1)).^2)).*((valid_H-param(2))./(0.5*param(1)))+param(5);
            plot(valid_H * 1e4, valid_V, 'o', LineWidth=4)
            hold on
            plot(valid_H * 1e4, ft)
            title("Ic = "+Ic_dc_vec(n)*1e3 + "[mA] f = " + w/2/pi + " [GHz] \phi_{H_{ext}} = " + phi_Hext*180/pi)
            legend("simul", "fit")
            xlabel('H_{ext} [Oe]')
            ylabel('V_{mix} [V]')
            hold off
        end
    
    
        %deletes all variables except
        %clearvars -except Hext_vec V_STFMR_vec delta_theta_vec delta_phi_vec phase_M_vec width_param res_param
    end
end


%% fit to I VS resonance width - linear fit


%plots
if plot_linewidth_vs_current
    ft = @(b,x) b(1)*x+b(2);
    start_points=[0, 0];
    res_width_vec_Oe=res_width_vec_currents*1e4; %[T]-->[Oe]
    param = nlinfit(Jc_dc_vec(:),res_width_vec_Oe(:),ft, start_points);
    slope=param(1);
    figure()
    plot(Jc_dc_vec,res_width_vec_Oe, "o")
    hold on
    plot(Jc_dc_vec, ft(param, Jc_dc_vec))
    hold off
    title("slope = " + slope)
    legend("Simul", "Linear fit")
    xlabel('J_{DC} [A/m^2]')
    ylabel('resonance width [Oe]')
end

if plot_f_vs_h0
    figure()
    plot(res_field_vec_freqs, w_vec/2/pi, 'o', LineWidth=4)
    hold on
    title("f vs Hres at \phi_{H_{ext}} = " + phi_Hext*180/pi)
    xlabel("H_{res} [T]")
    ylabel("f [GHz]")
    hold off
end

if plot_dH_vs_f
    figure()
    plot(w_vec/2/pi, res_width_vec_freqs , 'o', LineWidth=4)
    hold on
    title("\DeltaH vs f at \phi_{H_{ext}} = " + phi_Hext*180/pi)
    ylabel("\DeltaH [T]")
    xlabel("f [GHz]")
    hold off
end

% clearvars -except slope