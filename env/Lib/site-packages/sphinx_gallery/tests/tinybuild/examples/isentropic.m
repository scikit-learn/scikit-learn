function isentropic(g)
    %% isentropic, adiabatic flow example
    %
    % In this example, the area ratio vs. Mach number curve is computed for a
    % hydrogen/nitrogen gas mixture.

    if nargin == 1
        gas = g;
    else
        gas = Solution('gri30.yaml', 'gri30');
    end

    %% Set the stagnation state
    gas.TPX = {1200.0, 10.0 * OneAtm, 'H2:1,N2:0.1'};
    gas.basis = 'mass';
end
