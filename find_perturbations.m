function [delta_theta, delta_phi] = find_perturbations(Ku, miu0, gamma, w, Ms_Am, H0_Am, alpha, theta0, theta_H, phi0, phi_H, Nx, Ny, Nz, Sx, Sy, Sz, H_STT_DC, H_STT_AC, h_Oe_theta, h_Oe_phi)
% Gamma in rad/(s*T)
% H_STT_DC, H_STT_AC, h_Oe [T]
% H0_Am, Ms_Am [A/m]

% Second derivatives of energies:
E_zeeman_theta_theta = miu0 * Ms_Am * H0_Am * (sin(theta0) * sin(theta_H) * cos(phi0-phi_H) + cos(theta0)*cos(theta_H)); % [J / m^3]
E_zeeman_theta_phi = miu0 * Ms_Am * H0_Am * cos(theta0) * sin(theta_H) * sin(phi0-phi_H); % [J / m^3]
E_zeeman_phi_theta = E_zeeman_theta_phi; % [J / m^3]
E_zeeman_phi_phi = miu0 * Ms_Am * H0_Am * sin(theta0) * sin(theta_H) * cos(phi0-phi_H); % [J / m^3]

E_demag_theta_theta = miu0 * Ms_Am^2 * cos(2*theta0) * (Nx*cos(phi0)^2 + Ny*sin(phi0)^2 - Nz); % [J / m^3]
E_demag_theta_phi = 0.5*miu0 * Ms_Am^2 * sin(2*theta0) * sin(2*phi0) * (-Nx+Ny); % [J / m^3]
E_demag_phi_theta = E_demag_theta_phi; % [J / m^3]
E_demag_phi_phi = miu0 * Ms_Am^2 * sin(theta0)^2 * cos(2*phi0) * (-Nx+Ny); % [J / m^3]

E_Anis_theta_theta = Ku*2*cos(2*theta0); % [J / m^3]
% All other E_Anis derivatives are 0.

E_theta_theta = E_zeeman_theta_theta + E_demag_theta_theta + E_Anis_theta_theta; % [J / m^3]
E_theta_phi = E_zeeman_theta_phi + E_demag_theta_phi; % [J / m^3]
E_phi_theta = E_zeeman_phi_theta + E_demag_phi_theta; % [J / m^3]
E_phi_phi = E_zeeman_phi_phi + E_demag_phi_phi; % [J / m^3]

epsilon_theta_theta = E_theta_theta; % [J / m^3]
epsilon_theta_phi = E_theta_phi + Ms_Am * H_STT_DC * (Sx*cos(phi0) - Sy*sin(phi0)); % [J / m^3]
epsilon_phi_theta = E_phi_theta - Ms_Am * sin(theta0) * H_STT_DC * (Sx*sin(theta0)*cos(phi0) + Sy*sin(theta0)*sin(phi0) - Sz*cos(theta0)); % [J / m^3]
epsilon_phi_phi = E_phi_phi + Ms_Am * sin(theta0) * H_STT_DC * (-Sx*sin(phi0)*cos(theta0) + Sy*cos(theta0)*cos(phi0)); % [J / m^3]

% A matrix determinant equals w^2 - iwdw - w0^2:
w0 = gamma / (sin(theta0)*(alpha^2+1)*Ms_Am) * sqrt(abs(epsilon_phi_phi*epsilon_theta_theta - epsilon_phi_theta*epsilon_theta_phi));
dw = gamma / (sin(theta0)*(alpha^2+1)*Ms_Am) * (epsilon_phi_theta - epsilon_theta_phi + alpha * (1 / sin(theta0) * epsilon_phi_phi + sin(theta0)*epsilon_theta_theta));


% RF Fields [T]:
H_rf_theta = h_Oe_theta - H_STT_AC * (Sx*sin(phi0) - Sy*cos(phi0));
H_rf_phi = h_Oe_phi - H_STT_AC * (Sx*cos(theta0)*cos(phi0) - Sy*cos(theta0)*sin(phi0) + Sz*sin(theta0));

% Solution:
denum = 1 / ((w^2-w0^2)^2 + w^2*dw^2); % [rad/s ^-4]
% disp("denum = " + denum)
% disp("w0 = " + w0 + " and w = " + w)
% disp("dw = " + dw)


% delta_theta = gamma / ((w^2 - w0^2) +w^2*dw^2) * ( gamma / Ms * (w0^2-w^2)*(h_rf_phi*e_theta_phi+ 1/sin(theta) *h_rf_theta*e_phi_phi) - w^2 *dw * (h_rf_phi -alpha*h_rf_theta)*sin(theta) ) * Jc_ac; 
delta_theta = denum * ...
            (gamma^2 / Ms_Am * (w0^2-w^2) * ...
            (epsilon_theta_phi * H_rf_phi + epsilon_phi_phi * H_rf_theta / sin(theta0)) ...
            - w^2*dw*gamma*(H_rf_phi - alpha*H_rf_theta)*sin(theta0));

delta_phi = denum * ...
            (gamma^2 / Ms_Am * (w^2-w0^2) * ...
            (epsilon_theta_theta * H_rf_phi + epsilon_phi_theta * H_rf_theta / sin(theta0)) ...
            - w^2*dw*gamma*(H_rf_theta + alpha*H_rf_phi));


end