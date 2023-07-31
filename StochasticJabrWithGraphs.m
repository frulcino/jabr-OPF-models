%%TODOS
%todo variables varyiing respect to rho 8 valori di rho


clc; clear all; close all;
warning off

disp(['Running Data-based distributionally robust stochastic optimal power flow (OPF)'])
disp(['Stochastic Optimal power flow problem with wind farms'])
disp(['Latest update: 28/07/2023'])

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

Nbus = size(bus,1); % number of buses
Nbranch = size(branch,1); % number of lines
Ngen = size(gen,1); % number of traditional generators
refBus = mpc.bus(mpc.bus(:,2) == 3,1); %reference bus index

%% define named indices into bus, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus; 
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%% load wind energy forecast errors data
% Before import the wind power forecast errors, Please shift the forecast 
% errors data to make sure zero mean.
load('wind.mat')


% load three vector for each individual wind farm
%PROVA
%v = (1:52)';
%v_mean = mean(v);
%v_std = std(v);
%G_error_1 = (v - v_mean)/v_std * 200;
%G_error_2 = (v - v_mean)/v_std * 200;
%G_error_3 = (v - v_mean)/v_std * 300;

%from CESI wind farm
file_name = "C:\Users\Gabor\Documents\tesi_magistrale\data\previsioni\eolico\error_scenarios.csv";
error_scenarios = readmatrix(file_name);
G_error_1 = error_scenarios;
G_error_2 = error_scenarios;
G_error_3 = error_scenarios;
S_error = transpose([G_error_1, G_error_2, G_error_3]);

% generator located matrix
Cg = zeros(Nbus,Ngen);
for i = 1:Ngen
    Cg(gen(i,1),i) = 1;
end

% Got the load demand at each bus
Pd = bus(:,3);

% Got nominal wind power injection at each bus
wind = mpc.wind;
Nwind = size(wind,1);
Pwind = zeros(Nbus,1);
wingIdx = size(mpc.wind,1)

for i = 1:Nwind
    Pwind(wind(i,1)) = wind(i,2)
end


%% Get the number of loads
loadlocation = find(bus(:,3)>0);
load = mpc.bus(find(bus(:,3)>0),[1,3]);
Nload = size(load,1);

Ns = 3;% set Number of scenarios for distributionally robust optimization (30)
Nxi = 3; % set Number of stochastic variables (number of renewable generators)
beta = 0.01; % CVar level
rho_Matrix = [10,300,900,3000];

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
disp(['number of buses:    ', num2str(Nbus)]);
disp(['number of lines:   ', num2str(Nbranch)]);
disp(['number of generators:    ', num2str(Ngen)]);
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
epsilonMatrix = [0.02,0.04,0.06]; % wasserstein distance

for id_eps = 1:3
    Varepsilon = epsilonMatrix(id_eps)
    for  id_rho = 1:length(rho_Matrix)
        rho = rho_Matrix(id_rho);

        disp(['..................................................................................'])
        disp(['Wasserstein distance:    ', num2str(Varepsilon)]);
        %disp(['Current weight factor:   ', num2str(rho)]);
        %disp(["round number:",i,"out of", size])
        cvx_begin
        cvx_solver GUROBI 

        %policy variables
        disp(['defining policy variables'])
        variable PGn(Ngen) %nominal real power production
        variable QGn(Ngen) %nominal reactive power production 
        variable PGe(Ngen, Nxi) %real power reserve policy
        variable QGe(Ngen, Nxi) %reactive power reserve policy

        variable cn(Nbranch)
        variable sn(Nbranch)
        variable un(Nbus)
        variable ce(Nbranch,Nxi)
        variable se(Nbranch,Nxi)
        variable ue(Nbus,Nxi)

        %dummy variables
        variable Pn(Nbranch,2) %matrix having for each row P_ij and ji with ij a branc
        variable Qn(Nbranch,2) %Q_ij
        variable Pe(Nbranch,2,Nxi) %matrix having for each row P_ij and ji with ij a branc
        variable Qe(Nbranch,2,Nxi)

        %cvar auxiliary variables (sigma)
        variable oPL(Nbranch)
        variable oVM(Nbus)
        variable oVm(Nbus)
        variable oPGM(Ngen)
        variable oPGm(Ngen)


        disp(['defining finite reduction variables...'])
        %power loss contraint finite reduction variable
        variable lambdaPL(Nbranch)
        variable sPL(Nbranch,Ns)
        %variable gammaPL(Nbranch,Ns,2)
        %variable APL(Nbranch,Ns)
        %variable bPL(Nbranch)


        %voltage mangitude constraint finite reduction variable
        variable lambdaVm(Nbus) %for minimum voltage
        variable sVm(Nbus,Ns)
        %variable gammaVm(Nbus,Ns,2)
        %variable AVm(Nbus,Ns)
        %variable bVm(Nbus)

        variable lambdaVM(Nbus)
        variable sVM(Nbus,Ns)
        %variable gammaVM(Nbus,Ns,2)
        %variable AVM(Nbus,Ns)
        %variable bVM(Nbus) 

        %PG magnitude constraing finite reduction variable
        variable lambdaPGm(Ngen) %for minimum voltage
        variable sPGm(Ngen,Ns)
        %variable gammaPGm(Ngen,Ns,2)
        %variable APGm(Ngen,Ns)
        %variable bPGm(Ngen) 

        variable lambdaPGM(Ngen)
        variable sPGM(Ngen,Ns)
        %variable gammaPGM(Ngen,Ns,2)
        %variable APGM(Ngen,Ns)
        %variable bPGM(Ngen) 

        Fop = 0;
        % operation costs (avarage over dataset)
        for i = 1:Ngen
            Fop = Fop + sum( gencost(i,5)*(PGe(i,1).*G_error_1(:,1) + PGe(i,2).*G_error_2(:,1) + PGe(i,3).*G_error_3(:,1) + PGn(i)).^2 +(gencost(i,6))*(PGe(i,1).*G_error_1(:,1) + PGe(i,2).*G_error_2(:,1) + PGe(i,3).*G_error_3(:,1) + PGn(i)) + gencost(i,7));
        end
        Fop = Fop./Ns;
        % cost of constraints violation (finite reduction)
        Frisk = (sum(sum(sPL)) + sum(sum(sVM + sVm)) + sum(sum(sPGM+sPGm)))/Ns + Varepsilon*(sum(sum(lambdaPL)) + sum(sum(lambdaVM+lambdaVm)) + sum(sum(lambdaPGM+lambdaPGm)));
        F=Fop +  Frisk;
        minimize F
        subject to
        disp(['Adding line constraints...'])

        for b = 1:Nbranch

            f = mpc.branch(b,1);%from
                t = mpc.branch(b,2); %to %[Gff(b),Gft(b),Gtf(b),Gtt(b)]
                %s ha segni diversi paper per contare che s(i,j) = -s(j,i)
                %invece che imporlo come constraint

                %nominal power flow constraint
                Pn(b,1) == Gff(b)*un(f) + Gft(b)*cn(b) + Bft(b)*sn(b); % P_ft b=ft
                Qn(b,1) == -Bff(b)*un(f) - Bft(b)*cn(b) + Gft(b)*sn(b);
                Pn(b,2) == Gtt(b)*un(t) + Gtf(b)*cn(b) - Btf(b)*sn(b); % P_tf b=ft
                Qn(b,2) == -Btt(b)*un(t) - Btf(b)*cn(b) - Gtf(b)*sn(b);

                %policy power flow constraint
                Pn(b,1,:) == Gff(b)*ue(f,:) + Gft(b)*ce(b,:) + Bft(b)*se(b,:); % P_ft b=ft
                Qn(b,1,:) == -Bff(b)*ue(f,:) - Bft(b)*ce(b,:) + Gft(b)*se(b,:);
                Pn(b,2,:) == Gtt(b)*ue(t,:) + Gtf(b)*ce(b,:) - Btf(b)*se(b,:); % P_tf b=ft
                Qn(b,2,:) == -Btt(b)*ue(t,:) - Btf(b)*ce(b,:) - Gtf(b)*se(b,:);

                %power loss nominal constraints
                Pn(b,1,:) + Pn(b,2,:) >= 0; %todo: matrix constraint

                %power loss finite reduction constraints (H=0 d=0)
                %APL(b) == (-((Gft(b)+Gtf(b))*ce(b,:)+(Bft(b)-Btf(b))*se(b,:)-Gff(b)*ue(f,:)-Gtt(b)*ue(t,:))); 
                %bPL(b) == (-((Gft(b)+Gtf(b))*cn(b,:)+(Bft(b)-Btf(b))*sn(b,:)-Gff(b)*un(f,:)-Gtt(b)*un(t,:)));
                for i = 1:Ns %todo: substitute with matrix multiplication constraint
                    rho*((-((Gft(b)+Gtf(b))*cn(b,:)+(Bft(b)-Btf(b))*sn(b,:)-Gff(b)*un(f,:)-Gtt(b)*un(t,:)))+(1-beta)*oPL(b)+(-((Gft(b)+Gtf(b))*ce(b,:)+(Bft(b)-Btf(b))*se(b,:)-Gff(b)*ue(f,:)-Gtt(b)*ue(t,:)))*S_error(:,i)) <= sPL(b,i);
                    rho*(-oPL(b)*beta) <= sPL(b,i);
                    norm(-rho*(-((Gft(b)+Gtf(b))*ce(b,:)+(Bft(b)-Btf(b))*se(b,:)-Gff(b)*ue(f,:)-Gtt(b)*ue(t,:))),Inf) <= lambdaPL(b)
                end

        end

        disp(['Adding bus constraints...'])
        for b = 1:Nbus

            %policy constraints
                N = Nmat(b,mpc);

                %power flow constraint at bus b
                isgen = 0; %is there a gen at bus b?
                for g = 1:Ngen
                    if b == mpc.gen(g,1);
                        isgen=1;
                        gIdx = g; %bus index of g
                        break
                    end
                end

                     %add se è un generatore eolico:
                iswind = 0; %is there a gen at bus b?
                for g = 1:Nwind
                    if b == mpc.wind(g,1);
                        iswind =1;
                        wingIdx = g;
                        break
                    end
                end

                if isgen == 1
                    %nominal injection constraints
                   sum(sum(N .* Pn)) == PGn(gIdx) - mpc.bus(b,PD); %    segni?
                   sum(sum(N .* Qn)) == QGn(gIdx) - mpc.bus(b,QD);
                   %nominal power generation constraints
                   PGn(gIdx) <= mpc.gen(gIdx,PMAX)
                   PGn(gIdx) >= mpc.gen(gIdx,PMIN)
                   if mpc.gen(g,QMIN) ~= Inf
                       QGn(gIdx) <= mpc.gen(gIdx,QMAX)
                   end
                   if mpc.gen(g,QMIN) ~= -Inf
                       QGn(gIdx) >= mpc.gen(gIdx,QMIN)
                   end

                   for s = 1:Nxi
                    %policy constraints
                       sum(sum(N .* Pe(:,:,s))) == PGe(gIdx,s); %    segni?
                       sum(sum(N .* Qe(:,:,s))) == QGe(gIdx,s);
                   end

                   % PG max magnitude finite reduction constraints
                   %APGM(gIdx) == PGe(gIdx) ;
                   %bPGM(gIdx) == (PGn(gIdx) - mpc.gen(gIdx,PMAX));
                   for i = 1:Ns %todo: substitute with matrix multiplication constraint
                       rho*((PGn(gIdx) - mpc.gen(gIdx,PMAX))+(1-beta)*oPGM(gIdx)+PGe(gIdx)*S_error(:,i)) <= sPGM(gIdx,i);
                       rho*(-oPGM(gIdx)*beta) <= sPGM(gIdx,i);
                       norm(-rho*PGe(gIdx),Inf) <= lambdaPGM(gIdx);
                   end

                   % PG min magnitude finite reduction constraints
                   %APGm(gIdx) == (-PGe(gIdx));
                   %bPGm(gIdx) == (-PGn(gIdx) + mpc.gen(gIdx,PMIN));
                   for i = 1:Ns %todo: substitute with matrix multiplication constraint
                       rho*((-PGn(gIdx) + mpc.gen(gIdx,PMIN))+(1-beta)*oPGm(gIdx)+(-PGe(gIdx))*S_error(:,i)) <= sPGm(gIdx,i);
                       rho*(-oPGM(gIdx)*beta) <= sPGm(gIdx,i);
                       norm(-rho*(-PGe(gIdx)),Inf) <= lambdaPGm(gIdx);
                   end


                elseif iswind == 1
                   %nominal constraint with wind power generator
                   sum(sum(N .* Pn)) == mpc.wind(wingIdx,2) - mpc.bus(b,PD); 
                   sum(sum(N .* Qn)) == - mpc.bus(b,QD);
                   for s = 1:Nxi
                   %policy constraints with wind power generator P^Ge_g is a
                   %e_g
                       if s == wingIdx
                           sum(sum(N .* Pe(:,:,s))) == 1; %    segni?
                           sum(sum(N .* Qe(:,:,s))) == 0;
                       else
                           sum(sum(N .* Pe(:,:,s))) == 0; %    segni?
                           sum(sum(N .* Qe(:,:,s))) == 0;
                       end

                   end

                else
                    sum(sum(N .* Pn)) == - mpc.bus(b,PD); %sto considerando injected power or outjected power?
                    sum(sum(N .* Qn)) == - mpc.bus(b,QD);
                    for s = 1:Nxi
                    %policy constraints (with no generator at bus b)
                       sum(sum(N .* Pe(:,:,s))) == 0; %    segni?
                       sum(sum(N .* Qe(:,:,s))) == 0;
                   end
                end

            %V magnitude nominal constraints
            un(b) >= 0;
            un(b) >= mpc.bus(b,VMIN)^2*un(refBus);
            un(b) <= mpc.bus(b,VMAX)^2*un(refBus);

            %V max magnitude finite reduction constraints
            %AVM(b) == (ue(b) - mpc.bus(b,VMAX)^2*ue(refBus));
            %bVM(b) == (un(b) - mpc.bus(b,VMAX)^2*un(refBus));
                for i = 1:Ns %todo: substitute with matrix multiplication constraint
                    rho*((un(b) - mpc.bus(b,VMAX)^2*un(refBus))+(1-beta)*oVM(b)+(ue(b) - mpc.bus(b,VMAX)^2*ue(refBus))*S_error(:,i)) <= sVM(b,i);
                    rho*(-oVM(b)*beta) <= sVM(b,i);
                    norm(-rho*(ue(b) - mpc.bus(b,VMAX)^2*ue(refBus)),Inf) <= lambdaVM(b);
                end

            %V min magnitude finite reduction constraints
            %AVm(b) == (- ue(b) + mpc.bus(b,VMIN)^2*ue(refBus));
            %bVm(b) == (- un(b) + mpc.bus(b,VMIN)^2*un(refBus));
                for i = 1:Ns %todo: substitute with matrix multiplication constraint
                    rho * ((- un(b) + mpc.bus(b,VMIN)^2*un(refBus)) + (1-beta) * oVm(b)+(- ue(b) + mpc.bus(b,VMIN)^2*ue(refBus))*S_error(:,i)) <= sVm(b,i);
                    rho*(-oVm(b)*beta) <= sVm(b,i);
                    norm(-rho*(- ue(b) + mpc.bus(b,VMIN)^2*ue(refBus)),Inf) <= lambdaVm(b);
                end

        end
        disp(['Finished adding constraints, solving model...'])
        cvx_end 

        disp(['Finished solving model, saving results... ', num2str(rho)]);

         %% Data save and plot
        % save total operational cost and sum-up CVaR values of all line
        % constraints violation
        Cost_Data(id_rho,id_eps) = Fop;   
        CVaR_Data(id_rho,id_eps) = Frisk/rho; % 

        %PGn at each generator
        PG(1:Ngen,id_rho,id_eps) = PGn;
        QG(1:Ngen,id_rho,id_eps) = QGn;

        %un,cn,sn
        U(1:Nbus,id_rho,id_eps) = un;
        S(1:Nbranch,id_rho,id_eps) = sn;
        C(1:Nbranch,id_rho,id_eps) = cn;

        % save reserve policies


        % save schedule power output


        % save CVaR auxillary variables


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

%plot nominal power generators of 6 generators chosen at random
% Selecting 6 random generator indices
sixgens = randi([1, Ngen], 1, 6);
figure;

for i = 1:6
    subplot(3,2,i)
    
    % Plotting the power generation values of 3 scenarios (or cases) for the selected generator
    plot(rho_Matrix, PG(sixgens(i),:,1),'o-','LineWidth',1.5)
    hold on
    plot(rho_Matrix, PG(sixgens(i),:,2),'*-','LineWidth',1.5)
    hold on
    plot(rho_Matrix, PG(sixgens(i),:,3),'s-','LineWidth',2)
    
    % Plotting the nominal power generation limit of the selected generator
    nominal_limit = mpc.gen(sixgens(i), PMAX); % Fetching the nominal power generation limit
    plot(rho_Matrix, nominal_limit .* ones(size(rho_Matrix)), 'k--', 'LineWidth', 2, 'color', [0.5 0.5 0.5])
    
    legend([char(949),'  = 0.00'], [char(949), '  = 0.02'], [char(949), '  = 0.04'], 'nominal power generation limit')
    xlabel('weight factor \rho')
    ylabel('Power (MW)')  
    title(['Power generation at bus ', num2str(sixgens(i))])  
    set(gca,'FontSize',13);
    grid on
end
saveas(gcf, 'PG_Plot.png');
%plot nominal voltage of 6 nodes chosen at random
% Selecting 6 random bus indices
sixbus = randi([1, Nbus], 1, 6);
figure;

for i = 1:6
    subplot(3,2,i)
    
    % Plotting the voltage values of 3 scenarios (or cases) for the selected bus
    plot(rho_Matrix, U(sixbus(i),:,1)./U(refBus,:,1),'o-','LineWidth',1.5)
    hold on
    plot(rho_Matrix, U(sixbus(i),:,2)./U(refBus,:,2),'*-','LineWidth',1.5)
    plot(rho_Matrix, U(sixbus(i),:,3)./U(refBus,:,3),'s-','LineWidth',2)
    
    % Plotting the nominal voltage limits of the selected bus
    nominal_VMlimit = mpc.bus(sixbus(i), VMAX)^2; % Fetching the max nominal voltage
    plot(rho_Matrix, nominal_VMlimit .* ones(size(rho_Matrix)), 'k--', 'LineWidth', 2)
    nominal_Vmlimit = mpc.bus(sixbus(i), VMIN)^2; % Fetching the min nominal voltage
    plot(rho_Matrix, nominal_Vmlimit .* ones(size(rho_Matrix)), 'r--', 'LineWidth', 2)  % Changed color to red for distinction
    
    legend([char(949),'  = 0.02'], [char(949), '  = 0.04'], [char(949), '  = 0.06'], 'Vmax', 'Vmin')  % Updated legend
    xlabel('weight factor \rho')
    ylabel('Voltage^2 p.u.')  
    title(['Voltage at bus ', num2str(sixbus(i))])
    set(gca,'FontSize',13);
    grid on
end
saveas(gcf, 'U_Plot.png');
% Plot Cost
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
saveas(gcf, 'Cost_Plot.png'); % Saves the current figure (gcf = get current figure) to a file named 'Cost_Plot.png'

% Plot CVaR
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
saveas(gcf, 'CVaR_Plot.png'); % Saves the current figure to a file named 'CVaR_Plot.png'



