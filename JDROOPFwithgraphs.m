clc; clear all; close all;
warning off

disp(['Running Data-based distributionally robust stochastic optimal power flow (OPF)'])
disp(['Stochastic Optimal power flow problem with wind farms'])

define_constants;
mpc = loadcase("case118wind3.m");

tic 

%% Input Data Section
% System MVA base
baseMVA = mpc.baseMVA;
% Number of Nodes
bus = mpc.bus;
branch = mpc.branch;
gen = mpc.gen;
gencost = mpc.gencost;

n_bus = size(bus,1); % number of buses
n_branch = size(branch,1); % number of lines
n_gen = size(gen,1); % number of traditional generators
refBus = mpc.bus(mpc.bus(:,2) == 3,1); %reference bus index


Ns = 6;% set Number of scenarios for distributionally robust optimization (30)
Nxi = 3; % set Number of stochastic variables (number of renewable generators)
beta = 0.01; % CVar level

rho_Matrix = [10,100,300,500,700,1000,1300,1500,2000];
epsilon_Matrix = [0.01,0.02,0.05,0.1,0.2]; % wasserstein distance
n_rho = length(rho_Matrix);
n_Var = length(epsilon_Matrix);
formt = 'TestJabr2.mat' %saved file name
%% define named indices into bus, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus; 
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%% load wind energy forecast errors data
% Before import the wind power forecast errors, Please shift the forecast 
% errors data to make sure zero mean.
% load three vector for each individual wind farm

%from CESI wind farm
file_name = "C:\Users\Gabor\Documents\tesi_magistrale\data\previsioni\eolico\error_scenarios.csv";
error_scenarios = readmatrix(file_name);
error_scenarios = error_scenarios(1:Ns);
e_mean = mean(error_scenarios);
e_std = std(error_scenarios);
G_error_1 = (error_scenarios - e_mean)/e_std * 200;
G_error_2 = G_error_1;
G_error_3 = G_error_1*3/2;
S_error = transpose([G_error_1, G_error_2, G_error_3]);

% generator located matrix
Cg = zeros(n_bus,n_gen);
for i = 1:n_gen
    Cg(gen(i,1),i) = 1;
end

% Got the load demand at each bus
Pd = bus(:,3);

% Got nominal wind power injection at each bus
wind = mpc.wind;
Nwind = size(wind,1);
Pwind = zeros(n_bus,1);
wingIdx = size(mpc.wind,1);
n_rho = length(rho_Matrix);
n_epsilon = length(epsilon_Matrix);
for i = 1:Nwind
    Pwind(wind(i,1)) = wind(i,2);
end

%% Costructing matrices

[Yff, Yft, Ytf, Ytt] = Ybranch(mpc);
Gff = real(Yff);
Gft = real(Yft);
Gtf = real(Ytf);
Gtt = real(Ytt);
Bff = imag(Yff);
Bft = imag(Yft);
Btf = imag(Ytf);
Btt = imag(Ytt);


%% parameter check before start
disp(['program and system setup done !']);
disp(['Please check the following settings before start optimization']);
disp(['number of buses:    ', num2str(n_bus)]);
disp(['number of lines:   ', num2str(n_branch)]);
disp(['number of generators:    ', num2str(n_gen)]);
disp(['number of scenarios for distributionally robust stochastic OPF:  ', num2str(Ns)]);
disp(['number of renewable generatos:   ', num2str(Nxi)]);
disp(['range of weight factors:', 'start from ', num2str(min(rho_Matrix)),', up to ', num2str(max(rho_Matrix))]);
disp(['wind injection locations: ', '  bus 1,  ','    bus 9,  ','    bus 26']);
disp(['wind nominal injection:   ', '  500 (MW), ','  500 (MW), ','  800 (MW)']);
disp(['..................................................................................'])

%disp(['press any keys to continue .......']);
%pause

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  data-based distributionally robust stochastic OPF   %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Distribution Robust Optimization
% initialize data saving matrices
Cost_Data = zeros(n_rho,n_epsilon);   
CVaR_Data = zeros(n_rho,n_epsilon); % 

%PGn at each generator
PG = zeros(n_gen,n_rho,n_epsilon);
QG = zeros(n_gen,n_rho,n_epsilon);

%un,cn,sn
U = zeros(n_bus,n_rho,n_epsilon);
S = zeros(n_branch,n_rho,n_epsilon);
C = zeros(n_branch,n_rho,n_epsilon);

%Pn,Qn, Active and reactive nominal power flow
Pn_data = zeros(n_branch,2,n_rho,n_epsilon);
Qn_data = zeros(n_branch,2,n_rho,n_epsilon);

%save reserve policies
PG_reserve = zeros(n_gen,Nxi,n_rho,n_epsilon);
C_reserve = zeros(n_gen,Nxi,n_rho,n_epsilon);
U_reserve = zeros(n_gen,Nxi,n_rho,n_epsilon);
P_reserve = zeros(n_branch,2,Nxi,n_rho,n_epsilon);

for id_eps = 1:length(epsilon_Matrix)
    Varepsilon = epsilon_Matrix(id_eps)
    for  id_rho = 1:length(rho_Matrix)
        rho = rho_Matrix(id_rho);
        disp(['..................................................................................'])
        disp(['Wasserstein distance:    ', num2str(Varepsilon)]);
        disp(['Current weight factor:   ', num2str(rho)]);
        
        [PGn, QGn, PGe, QGe, cn, sn, un, ce, se, ue, Pn, Qn, Pe, Qe, Fop, Frisk, F] = JDROOPF3("case118wind3.m","winderrors.mat", Ns, Nxi, beta, rho, Varepsilon)
        

        %% Data save and plot
        % save total operational cost and sum-up CVaR values of all line
        % constraints violation
        Cost_Data(id_rho,id_eps) = Fop;   
        CVaR_Data(id_rho,id_eps) = Frisk/rho; % 

        %PGn at each generator
        PG(1:n_gen,id_rho,id_eps) = PGn;
        QG(1:n_gen,id_rho,id_eps) = QGn;
        
       
        %un,cn,sn
        U(1:n_bus,id_rho,id_eps) = un;
        S(1:n_branch,id_rho,id_eps) = sn;
        C(1:n_branch,id_rho,id_eps) = cn;

        %Pn,Qn, Active and reactive nominal power flow
        Pn_data(1:n_branch,1:2,id_rho,id_eps) = Pn;
        Qn_data(1:n_branch,1:2,id_rho,id_eps) = Qn;

        %save reserve policies
        PG_reserve(1:n_gen, 1:Nxi, id_rho, id_eps) = PGe;
        C_reserve(1:n_branch, 1:Nxi, id_rho, id_eps) = ce;
        U_reserve(1:n_bus, 1:Nxi, id_rho, id_eps) = ue;
        P_reserve(:, :, :, id_rho, id_eps) = Pe;
        
%         variable PGn(n_gen) %nominal real power production
%         variable QGn(n_gen) %nominal reactive power production 
%         variable PGe(n_gen, Nxi) %real power reserve policy
%         variable QGe(n_gen, Nxi) %reactive power reserve policy
% 
%         variable cn(n_branch)
%         variable sn(n_branch)
%         variable un(n_bus)
%         variable ce(n_branch,Nxi)
%         variable se(n_branch,Nxi)
%         variable ue(n_bus,Nxi)


    % save optimal status for each weight factor \rho


    % average actual line flow (mean of forecast errors is zero)
    
    
    %     cvx_status
    %     cvx_optval
    %     Cost_Data_0
    %     CVaR_Data_0
        

    end

end

toc
elapsedTime = toc;
fprintf('Elapsed time: %f seconds\n', elapsedTime);

%% PLOTS
date = datestr(now, 'yyyy-mm-dd_HH-MM-SS')

%plot nominal power generators of 6 generators chosen at random
% Selecting 6 random generator indices
sixgens = randi([1, n_gen], 1, 6);
figure;

for i = 1:6
    subplot(3,2,i)
    
    % Plotting the power generation values of 3 scenarios (or cases) for the selected generator
    for e = 1:n_Var
        plot(rho_Matrix, PG(sixgens(i),:,e),'LineWidth',1.5, "DisplayName", ['Epsilon ' num2str(epsilon_Matrix(e))] )
        hold on
    end
    
    % Plotting the nominal power generation limit of the selected generator
    nominal_limit = mpc.gen(sixgens(i), PMAX); % Fetching the nominal power generation limit
    plot(rho_Matrix, nominal_limit .* ones(size(rho_Matrix)), 'k--', 'LineWidth', 2, 'color', [0.5 0.5 0.5], "DisplayName", "PGmax (MW)")
    nominal_limit = mpc.gen(sixgens(i), PMIN); % Fetching the nominal power generation limit
    plot(rho_Matrix, nominal_limit .* ones(size(rho_Matrix)), 'k--', 'LineWidth', 2, 'color', [0.5 0.5 0.5], "DisplayName", "PGmin (MW)")
    
    
    xlabel('weight factor \rho')
    ylabel('Power (MW)')  
    title(['Power generation at bus ', num2str(sixgens(i))])  
    
end
set(gca,'FontSize',13);
grid on
legend('Location','best')
hold off

filename = strcat(date, formt,'PG_Plot.png'); 
saveas(gcf, filename);
%% plot nominal voltage of 6 nodes chosen at random
% Selecting 6 random bus indices
sixbus = randi([1, n_bus], 1, 6);
figure;

for i = 1:6
    
    subplot(3,2,i)
    
    for e = 1:n_Var
        % Plotting the voltage values of 3 scenarios (or cases) for the selected bus
        plot(rho_Matrix, U(sixbus(i),:,e)./U(refBus,:,e),'LineWidth',1.5, 'DisplayName', ['Epsilon ' num2str(epsilon_Matrix(e))])
        hold on
    end
    
    % Plotting the nominal voltage limits of the selected bus
    nominal_VMlimit = mpc.bus(sixbus(i), VMAX)^2; % Fetching the max nominal voltage
    plot(rho_Matrix, nominal_VMlimit .* ones(size(rho_Matrix)), 'k--', 'LineWidth', 2, 'DisplayName', "V Max (p.u.)")
    nominal_Vmlimit = mpc.bus(sixbus(i), VMIN)^2; % Fetching the min nominal voltage
    plot(rho_Matrix, nominal_Vmlimit .* ones(size(rho_Matrix)), 'r--', 'LineWidth', 2, 'DisplayName', "V Min (p.u.)")  % Changed color to red for distinction
    
     
    xlabel('weight factor \rho')
    ylabel('Voltage^2 p.u.')  
    title(['Voltage at bus ', num2str(sixbus(i))])
    
end
set(gca,'FontSize',13);
legend('Location', 'best'); % Display legend to identify each row
grid on
hold off

filename = strcat(date, formt,'U_Plot.png');
saveas(gcf, filename);
%% Plot Cost
figure; 
plotStyles = {'o-', '*-', 's-'};  % Different styles for each series

for j=1:3
    plot(rho_Matrix, Cost_Data(:,j), plotStyles{j},'LineWidth',1.5)
    hold on
end

legend([char(949),'  = 0.02'], [char(949), '  = 0.04'], [char(949), '  = 0.06'])
xlabel('weight factor \rho')
ylabel('Nominal Cost')  
title('Cost')
set(gca,'FontSize',13);
grid on
filename = strcat(date, formt,'Cost_Plot.png');
saveas(gcf, filename ); % Saves the current figure (gcf = get current figure) to a file named 'Cost_Plot.png'

%% Plot CVaR
figure;
for j=1:3
    plot(rho_Matrix, CVaR_Data(:,j)./rho_Matrix, plotStyles{j},'LineWidth',1.5)
    hold on
end

legend([char(949),'  = 0.02'], [char(949), '  = 0.04'], [char(949), '  = 0.06'])
xlabel('weight factor \rho')
ylabel('CVar')  
title('CVar')
set(gca,'FontSize',13);
grid on
filename = strcat(date, formt,'CVaR_Plot.png');
saveas(gcf, filename); % Saves the current figure to a file named 'CVaR_Plot.png'


%% Plot Line Flow Constraint (Constraint not enforced directly)
sixlines = randi([1, n_branch], 1, 6);
Sn_data = Pn_data.^2 + Qn_data.^2;
for i = 1:6
    subplot(3,2,i)
    
    for e = 1:n_Var
        % Squeezing data to remove singleton dimensions, conversion to p.u.
        data = squeeze(Sn_data(sixlines(i),1,:,e))/(un(refBus)^2);
        % Plotting the power flow values of 3 scenarios (or cases) for the selected line
        plot(rho_Matrix, data,'DisplayName', ['Epsilon ' num2str(epsilon_Matrix(e))],'LineWidth',1.5);
        hold on
    end
    
    % Plotting the nominal power flow limit of the selected line
    nominal_PMlimit = mpc.Fmax(sixlines(i)); % Fetching the max power flow
    plot(rho_Matrix, nominal_PMlimit .* ones(size(rho_Matrix)), 'k--', 'LineWidth', 2, 'DisplayName', 'Max Power Flow (MW)'); % Using PMlimit, since PLlimit isn't defined in the code you provided
    
    %legend( 'Vmax', 'Vmin')  % Updated legend
    xlabel('weight factor \rho')
    ylabel('Power flow through line (MW).')  
    title(['Power flow through line ', num2str(sixlines(i))]);
    
end
set(gca,'FontSize',13);
legend('Location', 'best'); % Display legend to identify each row
grid on
hold off
filename = strcat(date, formt,'Powerflow_plot.png');
saveas(gcf, filename);

%%OutOfSample Tests

%% Check for out of sample constraints violations
% load out of sample samples
errors_name = "winderrors.mat"
data = load(errors_name);
S_error = data.S_error;
outErr = S_error(Ns+1:size(S_error,1),:); %gets out of sample errors
n_samples = length(outErr);

osPG = zeros(n_rho, n_Var, n_samples, n_gen); %Power generation at each bus at each scenario using the found policy
osU = zeros(n_rho, n_Var, n_samples, n_bus); %for voltage magnitude contraints
osP = zeros(n_rho, n_Var, n_samples, 2, n_branch); %for power flow contraints (non imposed stochasticly on model, see how it goes)

P_reserve = permute(P_reserve, [1,3,2,4,5]); %reshape matrix to fit multiplication by \xi
for r = 1:n_rho
    for e = 1:n_epsilon
        for j = 1:n_samples
            %Calculate value of PG and U for the out of samples scenario j
            %at each bus.
            osPG(r,e,j,:) = PG(:,r,e) + PG_reserve(:,:,r,e)*transpose(outErr(j,:));
            osU(r,e,j,:) = U(:,r,e) + U_reserve(:,:,r,e)*transpose(outErr(j,:));
            osP(r,e,j,1,:) = Pn_data(:,1,r,e) + P_reserve(:,:,1,r,e)*transpose(outErr(j,:));
            osP(r,e,j,2,:) = Pn_data(:,2,r,e) + P_reserve(:,:,2,r,e)*transpose(outErr(j,:));
        end
    end
end

filename = strcat([date, formt]);
save(filename);

%Compute violations of constraints for various out of sample scencarios
%Violations on power generation
NPGmax_viol = zeros(n_rho, n_Var, n_gen); %number of violations
tot_PGmax_viol = zeros(n_rho, n_Var, n_gen); %sum of magnitute of violations
NPGmin_viol = zeros(n_rho, n_Var, n_gen);
tot_PGmin_viol = zeros(n_rho, n_Var, n_gen);
%Violations on voltage magnitude
n_UMax_viol = zeros(n_rho, n_Var, n_bus);
tot_UMax_viol = zeros(n_rho, n_Var, n_bus);
n_UPos_viol = zeros(n_rho, n_Var, n_bus);
tot_UPos_viol = zeros(n_rho, n_Var, n_bus);
%Power loss violations
P_loss = osP(:,:,:,1,:) + osP(:,:,:,2,:);
n_Ploss_viol = P_loss < 0;
tot_Ploss_viol = P_loss(n_Ploss_viol);


for r = 1:n_rho
    for e = 1:n_Var
        %constraints violations on gen contraints
        for g = 1:n_gen
            for j = 1:n_samples
                if osPG(r,e,j,g) > mpc.gen(g,PMAX)
                    NPGmax_viol(r,e,g) = NPGmax_viol(r,e,g) +1;
                    tot_PGmax_viol(r,e,g) = tot_PGmax_viol(r,e,g) + osPG(r,e,j,g) - mpc.gen(g,PMAX);
                end
                if osPG(r,e,j,g) < mpc.gen(g,PMIN)
                    NPGmin_viol(r,e,g) = NPGmin_viol(r,e,g) +1;
                    tot_PGmin_viol(r,e,g) = tot_PGmin_viol(r,e,g) - osPG(r,e,j,g) + mpc.gen(g,PMIN);
                end
            end
        end
        
        %contraint violtations on bus constraints
        for b = 1:n_bus
            for j = 1:n_samples
                if osU(r,e,j,b) < 0
                    n_UPos_viol(r,e,b) = n_UPos_viol(r,e,b) +1;
                    tot_UPos_viol(r,e,b) = tot_UPos_viol(r,e,b) - osU(r,e,j,b);
                end
                if osU(r,e,j,b) > mpc.bus(b,VMAX)^2*un(refBus)
                    n_UMax_viol(r,e,b) = n_UMax_viol(r,e,b) + 1;
                    tot_UMax_viol(r,e,b) = tot_UMax_viol(r,e,b) + osU(r,e,j,b) - mpc.bus(b,VMAX)^2*un(refBus)^2;
                end
                    
            end
        end
        
        
    end
end


%% Plot Error Violations

avNPGmax_viol = sum(NPGmax_viol,3)/n_samples;
avtot_PGmax_viol = sum(tot_PGmax_viol,3)/n_samples;
avNPGmin_viol = sum(NPGmin_viol,3)/n_samples;
avtot_PGmin_viol = sum(tot_PGmin_viol,3)/n_samples;
avTn_UMax_viol = sum(n_UMax_viol,3)/n_samples;
avtot_UMax_viol = sum(tot_UMax_viol,3)/n_samples;
avn_UPos_viol = sum(n_UPos_viol,3)/n_samples;
avtot_UPos_viol = sum(tot_UPos_viol,3)/n_samples;
avn_Ploss_viol = sum(sum(n_Ploss_viol,5),3)/n_samples;
atot_Ploss_viol = -sum(sum(tot_Ploss_viol,5),3)/n_samples;
% P_loss = osP(:,:,:,1,:) + osP(:,:,:,2,:)
% n_Ploss_viol = P_loss < 0;
% tot_Ploss_viol = P_loss(n_Ploss_viol);
% Define matrices
% I'm assuming these matrices have already been defined and computed in your code.

% Create a cell array of matrices for easy looping
matrices = {avNPGmax_viol, avtot_PGmax_viol, avNPGmin_viol, avtot_PGmin_viol, ...
            avTn_UMax_viol, -avtot_UMax_viol, avn_UPos_viol, avtot_UPos_viol};

% Names for each matrix (for the title of the graph)
matrixNames = {'TNPGmax\_violations', 'TPGmaxtot\_violation', 'TNPGmin\_violations', 'TPGmintot\_violation', ...
               'Tn\_UMaxViolation', 'TTot\_UMaxViolation', 'Tn\_UPosViolation', 'TTot\_UPosViolation'};

ylabelNames = {'n_MaxMagnitutePGViolations', 'Total_MaxMagnitutePGViolations', 'TNPGmin\_violations', 'TPGmintot\_violation', ...
               'Tn\_UMaxViolation', 'TTot\_UMaxViolation', 'Tn\_UPosViolation', 'TTot\_UPosViolation'};           
           
% Loop through each matrix
for m = 1:length(matrices)
    figure; % Creates a new figure for each matrix
    currentMatrix = matrices{m};
    
    % Loop through each row in the matrix and plot it as a separate line
    for e = 1:size(currentMatrix, 2)
        plot(rho_Matrix,currentMatrix(:, e), 'DisplayName', ['Epsilon ' num2str(epsilon_Matrix(e))],'LineWidth',1.5);
        hold on;
    end
    
    title(['Average:',matrixNames{m}]);
    xlabel('rho');
    ylabel(ylabelNames{m});
    legend('Location', 'best'); % Display legend to identify each row
    grid on;
    hold off;
end

%% Save
filename = strcat([date, 'OutOfSample',formt]);
save(filename);






