% Define common parameters and variables
% (assuming these are defined before the CVX model)

% Objective function
cvx_begin
    variable SgenR(Nbus, Ngen)
    variable SgenC(Nbus, Ngen)
    variable v(Nbus)
    variable theta(Nbus)
    
    % Generation cost
    minimize(sum(sum(C(:, :, 2) .* SgenR.^2 + C(:, :, 1) .* SgenR + C(:, :, 0))))
    
    % Power flow balance constraints
    subject to
    for b = 1:Nbus
        % Real power flow balance
        SDR(b) + sum(YffR(b, a, i) * v(b)^2 + YftC(b, a, i) * v(b) * v(a) * sin(theta(b) - theta(a)) + YftR(b, a, i) * v(b) * v(a) * cos(theta(b) - theta(a))) + sum(YtfC(a, b, i) * v(b) * v(a) * sin(theta(b) - theta(a)) + YtfR(a, b, i) * v(b) * v(a) * cos(theta(b) - theta(a)) + YttR(a, b, i) * v(b)^2) == -shR(b) * v(b)^2 + sum(SgenR(b, :))
        
        % Imaginary power flow balance
        SDC(b) + sum(-YffC(b, a, i) * v(b)^2 - YftC(b, a, i) * v(b) * v(a) * cos(theta(b) - theta(a)) + YftR(b, a, i) * v(b) * v(a) * sin(theta(b) - theta(a))) + sum(-YtfC(a, b, i) * v(b) * v(a) * cos(theta(b) - theta(a)) + YtfR(a, b, i) * v(b) * v(a) * sin(theta(b) - theta(a)) - YttC(a, b, i) * v(b)^2) == shC(b) * v(b)^2 + sum(SgenC(b, :))
    end
    
    % Power bounds on lines
    for b = 1:Nbus
        for a = 1:Nbus
            for i = 1:Nlines
                if SU(b, a, i) > 0 && SU(b, a, i) < Inf
                    % Power bound 1
                    YffC(b, a, i)^2 * v(b)^4 + 2 * YffC(b, a, i) * YftC(b, a, i) * v(a) * v(b)^3 * cos(theta(b) - theta(a)) - 2 * YffC(b, a, i) * YftR(b, a, i) * v(a) * v(b)^3 * sin(theta(b) - theta(a)) + YffR(b, a, i)^2 * v(b)^4 + 2 * YffR(b, a, i) * YftC(b, a, i) * v(a) * v(b)^3 * sin(theta(b) - theta(a)) + 2 * YffR(b, a, i) * YftR(b, a, i) * v(a) * v(b)^3 * cos(theta(b) - theta(a)) + YftC(b, a, i)^2 * v(a)^2 * v(b)^2 + YftR(b, a, i)^2 * v(a)^2 * v(b)^2 <= SU(b, a, i)^2;
                    
                    % Power bound 2
                    YtfC(a, b, i)^2 * v(a)^2 * v(b)^2 + 2 * YtfC(a, b, i) * YttC(a, b, i) * v(a) * v(b)^3 * cos(theta(b) - theta(a)) + 2 * YtfC(a, b, i) * YttR(a, b, i) * v(a) * v(b)^3 * sin(theta(b) - theta(a)) + YtfR(a, b, i)^2 * v(a)^2 * v(b)^2 - 2 * YtfR(a, b, i) * YttC(a, b, i) * v(a) * v(b)^3 * sin(theta(b) - theta(a)) + 2 * YtfR(a, b, i) * YttR(a, b, i) * v(a) * v(b)^3 * cos(theta(b) - theta(a)) + YttC(a, b, i)^2 * v(b)^4 + YttR(a, b, i)^2 * v(b)^4 <= SU(b, a, i)^2;
                end
            end
        end
    end
    
    % Bounds on phase difference
    for b = 1:Nbus
        for a = 1:Nbus
            for i = 1:Nlines
                pdLB(b, a, i) <= theta(b) - theta(a);
                theta(b) - theta(a) <= pdUB(b, a, i);
            end
        endApologies, but I made an error in the previous translation. Here's the corrected version:

```matlab
% Define common parameters and variables
% (assuming these are defined before the CVX model)

% Objective function
cvx_begin
    variable SgenR(Nbus, Ngen)
    variable SgenC(Nbus, Ngen)
    variable v(Nbus)
    variable theta(Nbus)
    
    % Generation cost
    minimize(sum(sum(C(:, :, 2) .* SgenR.^2 + C(:, :, 1) .* SgenR + C(:, :, 0))))
    
    % Power
