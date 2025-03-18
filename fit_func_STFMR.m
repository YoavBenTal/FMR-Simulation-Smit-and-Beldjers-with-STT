function [param] = fit_func_STFMR(field, volt)
    % Width parameter
    [max_val, max_index] = max(volt);
    [min_val, min_index] = min(volt);
    width_param = abs(field(max_index) - field(min_index)) / 2;
    
    % Resonance parameter
    res_param = field(round((max_index + min_index) / 2));
    
    % Fit function
    ft = @(b, x) b(3) .* (((0.5 * b(1)) .^ 2) ./ ((x - b(2)) .^ 2 + (0.5 * b(1)) .^ 2)) + ...
                  b(4) .* (((0.5 * b(1)) ^ 2) ./ ((x - b(2)) .^ 2 + (0.5 * b(1)) .^ 2)) .* ((x - b(2)) ./ (0.5 * b(1))) + b(5);
    
    % Initial guess
    start_points = [width_param, res_param, 0, 1, 0];

    % Lower and upper bounds
    lower_bounds = [1e-6, min(field), -Inf, -Inf, -Inf]; % Ensure dH > 0
    upper_bounds = [Inf, max(field), Inf, Inf, Inf];

    % Fit using lsqcurvefit
    options = optimset('Display', 'off'); % Suppress output
    param = lsqcurvefit(ft, start_points, field(:), volt(:), lower_bounds, upper_bounds, options);

    % Plot results (optional)
    % H = linspace(min(field), max(field), 1000);
    % v = param(3) .* (((0.5 * param(1)) .^ 2) ./ ((H - param(2)) .^ 2 + (0.5 * param(1)) .^ 2)) + ...
    %     param(4) .* (((0.5 * param(1)) ^ 2) ./ ((H - param(2)) .^ 2 + (0.5 * param(1)) .^ 2)) .* ((H - param(2)) ./ (0.5 * param(1))) + param(5);
    % figure; hold on;
    % plot(field, volt, 'o', 'DisplayName', 'Data');
    % plot(H, v, 'r', 'DisplayName', 'Fit');
    % legend; hold off;
end
