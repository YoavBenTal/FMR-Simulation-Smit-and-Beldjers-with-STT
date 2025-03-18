function [theta0, phi0] = find_steady_state(H0_, theta_H_, phi_H_, Ms_Am_, Nx_, Ny_, Nz_, Sx_, Sy_, Sz_, H_STT_dc_, Ku_, initial_guess_)
    % This function finds the steady state angles of the magnetization theta0,
    % phi from the energy equations in SI units. for E_Zeeman, E_Crystalline,
    % E_demag.
    % Ms, H0 in [A/m]
    % H_STT_dc_ [T]

    % Define constants
    miu0 = 4 * pi * 1e-7; % Permeability of free space [H/m]

    % Function to solve E_theta = 0 and E_phi = 0
    function F = energy_equations(vars)
        theta = vars(1);
        phi = vars(2);

        % Energy derivatives
        E_Zeeman_theta = -miu0 * Ms_Am_ * H0_ * (cos(theta) * sin(theta_H_) * cos(phi - phi_H_) - sin(theta) * cos(theta_H_)); % [J / m^3]
        E_Demag_theta = 0.5*miu0 * Ms_Am_^2 * sin(2*theta) * (Nx_ * cos(phi)^2 + Ny_ * sin(phi)^2 - Nz_); % [J / m^3]
        E_Crystalline_theta = Ku_ * sin(2 * theta); % [J / m^3]
        E_SHE_theta = Ms_Am_ * H_STT_dc_ * (Sx_ * sin(phi) - Sy_ * cos(phi)); % [J / m^3]
        E_theta = E_Zeeman_theta + E_Demag_theta + E_Crystalline_theta + E_SHE_theta;

        E_Zeeman_phi = miu0 * Ms_Am_ * H0_ * sin(theta) * sin(theta_H_) * sin(phi - phi_H_); % [J / m^3]
        E_Demag_phi = 0.5*miu0 * Ms_Am_^2 * sin(theta)^2 * sin(2*phi) * (-Nx_ + Ny_); % [J / m^3]
        E_SHE_phi = Ms_Am_ * sin(theta) * H_STT_dc_ * (Sx_ * cos(theta) * cos(phi) - Sy_ * cos(theta) * sin(phi) + Sz_ * sin(theta)); % [J / m^3]
        E_phi = E_Zeeman_phi + E_Demag_phi + E_SHE_phi;

        % Output the equations to be solved
        F = [E_theta; E_phi];
    end

    % Set optimization options
    options = optimoptions('lsqnonlin', ...
    'Algorithm', 'trust-region-reflective', ...
    'TolFun', 1e-6, ...
    'TolX', 1e-6, ...
    'Display', 'off');

    % Bounds
    lb = [0, 0];       % Lower bounds
    ub = [pi, 2*pi]; % Upper bounds

    solution = lsqnonlin(@energy_equations, initial_guess_, lb, ub, options);
    % Extract results
    theta0 = solution(1);
    phi0 = solution(2);
end
