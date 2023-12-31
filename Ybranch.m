function [Yff,Yft,Ytf,Ytt] = Ybranch(mpc)
%Returns Yff, Yft,Ytf,Ytt  for every branch in the network
%As defined in: Bienstock OPF math programminf formulations
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

branch = mpc.branch;
bus = mpc.bus;
nb = size(bus, 1);          %% number of buses
nl = size(branch, 1);



stat = branch(:, BR_STATUS); % ones at in-service branches
Ys = stat ./ (branch(:, BR_R) + 1j * branch(:, BR_X));        % series admittance
Bc = stat .* branch(:, BR_B);                                % line charging susceptance
tap = ones(nl, 1);                                           % default tap ratio = 1
i = find(branch(:, TAP));                                    % indices of non-zero tap ratios
tap(i) = branch(i, TAP);                                     % assign non-zero tap ratios
tap = tap .* exp(1j*pi/180 * branch(:, SHIFT));              % add phase shifters
Ytt = Ys + 1j*Bc/2;
Yff = Ytt ./ (tap .* conj(tap));
Yft = - Ys ./ conj(tap);
Ytf = - Ys ./ tap;

end

