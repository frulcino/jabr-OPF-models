function [PGn, QGn, PGe, QGe, cn, sn, un, ce, se, ue, Pn, Qn, Pe, Qe, Fop, Frisk, F] = JDROOPF(casename,errors_name,Ns, Nxi, beta, rho, Varepsilon)
% Jabr distributionally Robust OPtimal Power Flows
%Input:
% casename: str containing the casename for the network
%    -also contains the additions field "wind" containing two columns the
%    first one with the windpower location (bus index) and the second one
%    with the windpower nominal value
% errors_name: str containing the file name containing the samples of the
% error of prediction of the wind power production
%   -each row has the prediction error for each wind turbine
% Ns: number of samples of error of prediction to be used
% Nxi: number of renewable generators
% beta: CVar level
% rho: risk aversion coeficient
% Varepsilon: Wasserstein ball radious
% [casename, errors_name, Ns, Nxi, beta, rho, Varepsilon] = deal("case118wind3.m","winderrors.mat", 4, 3, 0.01, 20, 0.02)
define_constants;
mpc = loadcase(casename);

tic 

%% Input Data Section
% Number of Nodes
bus = mpc.bus;
branch = mpc.branch;
gen = mpc.gen;
gencost = mpc.gencost;
n_bus = size(bus,1); % number of buses
n_branch = size(branch,1); % number of lines
n_gen = size(gen,1); % number of traditional generators
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
%from CESI wind farm
data = load(errors_name); %extracts variable S_error
S_error = transpose(data.S_error);
S_error = S_error(:,1:Ns);
% Got the load demand at each bus
Pd = bus(:,3);

% Got nominal wind power injection at each bus
wind = mpc.wind;
Nwind = size(wind,1);
Pwind = zeros(n_bus,1);


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
%disp(['range of weight factors:', 'start from ', num2str(min(rho_Matrix)),', up to ', num2str(max(rho_Matrix))]);
disp(['wind injection locations: ', '  bus 1,  ','    bus 9,  ','    bus 26']);
disp(['wind nominal injection:   ', '  500 (MW), ','  500 (MW), ','  800 (MW)']);
disp(['..................................................................................'])

%disp(['press any keys to continue .......']);
%pause

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  data-based distributionally robust stochastic OPF   %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Distribution Robust Optimization (wasserstein distance = 0.00)

 
    disp(['..................................................................................'])
    disp(['Wasserstein distance:    ', num2str(Varepsilon)]);
    %disp(['Current weight factor:   ', num2str(rho)]);
    %disp(["round number:",i,"out of", size])
    
    %% Set up Gurobi parameters
    %gurobi_params = struct();
    %gorobi_params.BarHomogeneous = 1;
    %% Build Model
    cvx_begin
    cvx_solver GUROBI
    
    %policy variables
    disp(['defining policy variables'])
    variable PGn(n_gen) %nominal real power production
    variable QGn(n_gen) %nominal reactive power production 
    variable PGe(n_gen, Nxi) %real power reserve policy
    variable QGe(n_gen, Nxi) %reactive power reserve policy
    
    variable cn(n_branch)
    variable sn(n_branch)
    variable un(n_bus)
    variable ce(n_branch,Nxi)
    variable se(n_branch,Nxi)
    variable ue(n_bus,Nxi)
    
    %dummy variables
    variable Pn(n_branch,2) %matrix having for each row P_ij and ji with ij a branc
    variable Qn(n_branch,2) %Q_ij
    variable Pe(n_branch,2,Nxi) %matrix having for each row P_ij and ji with ij a branc
    variable Qe(n_branch,2,Nxi)
    
    %cvar auxiliary variables (sigma)
    variable oPL(n_branch)
    variable oVM(n_bus)
    variable oVm(n_bus)
    variable oPGM(n_gen)
    variable oPGm(n_gen)
   
    
    disp(['defining finite reduction variables...'])
    %power loss contraint finite reduction variable
    variable lambdaPL(n_branch)
    variable sPL(n_branch,Ns)

    %voltage mangitude constraint finite reduction variable
    variable lambdaVm(n_bus) %for minimum voltage
    variable sVm(n_bus,Ns)
    
    variable lambdaVM(n_bus)
    variable sVM(n_bus,Ns)
    
    %PG magnitude constraing finite reduction variable
    variable lambdaPGm(n_gen) %for minimum voltage
    variable sPGm(n_gen,Ns)
    variable lambdaPGM(n_gen)
    variable sPGM(n_gen,Ns)
    
    Fop = 0;
    % operation costs (avarage over dataset)
    for i = 1:n_gen
        Fop = Fop + sum( gencost(i,5)*( PGn(i) + sum(transpose(S_error)*diag(PGe(i,:)),2) ).^2 + gencost(i,6)*( PGn(i) +  sum(transpose(S_error)*diag(PGe(i,:)),2)) + gencost(i,7));
    end
    Fop = Fop./Ns;
    % cost of constraints violation (finite reduction)
    Frisk = (sum(sum(sPL)) + sum(sum(sVM + sVm)) + sum(sum(sPGM+sPGm)))/Ns + Varepsilon*(sum(sum(lambdaPL)) + sum(sum(lambdaVM+lambdaVm)) + sum(sum(lambdaPGM+lambdaPGm)));
    F = Fop +  Frisk;
    minimize F
    subject to
    disp(['Adding line constraints...'])
    %%nominal line constraints:
    Pn(:,1).^2 + Qn(:,1).^2 <= mpc.Fmax*(un(refBus))*20;
    %Pn(:,2).^2 + Qn(:,2).^2 <= mpc.Fmax*(un(refBus))*1000;
    
   %power loss nominal constraints
   Pn(:,1,:) + Pn(:,2,:) >= 0;  %78043
   %nominal power generation constraints
   PGn(:) <= mpc.gen(:,PMAX)
   PGn(:) >= mpc.gen(:,PMIN)
   %V magnitude nominal constraints
   un(:) >= 0;
   un(:) >= mpc.bus(:,VMIN).^2*un(refBus);
   un(:) <= mpc.bus(:,VMAX).^2*un(refBus);
   %un(refBus) == 1
    
for b = 1:n_branch

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
    Pe(b,1) == Gff(b)*ue(f) + Gft(b)*ce(b) + Bft(b)*se(b); % P_ft b=ft
    Qe(b,1) == -Bff(b)*ue(f) - Bft(b)*ce(b) + Gft(b)*se(b);
    Pe(b,2) == Gtt(b)*ue(t) + Gtf(b)*ce(b) - Btf(b)*se(b); % P_tf b=ft
    Qe(b,2) == -Btt(b)*ue(t) - Btf(b)*ce(b) - Gtf(b)*se(b);



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
for b = 1:n_bus

    %policy constraints
        N = Nmat(b,mpc);

        %power flow constraint at bus b
        isgen = 0; %checks if there is a gen at bus b
        for g = 1:n_gen
            if b == mpc.gen(g,1);
                isgen=1;
                gIdx = g; %bus index of g
                break
            end
        end

             
        iswind = 0; %is there a gen at bus b
        for g = 1:Nwind
            if b == mpc.wind(g,1);
                iswind =1;
                windIdx = g;
                break
            end
        end

        if isgen == 1
            %nominal injection constraints
           sum(sum(N .* Pn)) == PGn(gIdx) - mpc.bus(b,PD); %    segni?
           sum(sum(N .* Qn)) == QGn(gIdx) - mpc.bus(b,QD);
           %nominal power generation constraints
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
           sum(sum(N .* Pn)) == mpc.wind(windIdx,2) - mpc.bus(b,PD); 
           sum(sum(N .* Qn)) == - mpc.bus(b,QD);
           for s = 1:Nxi
           %policy constraints with wind power generator P^Ge_g is a
           %e_g
               if s == windIdx
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



    %V max magnitude finite reduction constraints
    %AVM(b) == (ue(b) - mpc.bus(b,VMAX)^2*ue(refBus));
    %bVM(b) == (un(b) - mpc.bus(b,VMAX)^2*un(refBus));
        for i = 1:Ns %todo: substitute with matrix multiplication constraint
            rho*((un(b) - mpc.bus(b,VMAX)^2*un(refBus))+(1-beta)*oVM(b)+(ue(b) - mpc.bus(b,VMAX)^2*ue(refBus))*S_error(:,i)) <= sVM(b,i);
            rho*(-oVM(b)*beta) <= sVM(b,i);
            norm(-rho*(ue(b) - mpc.bus(b,VMAX)^2*ue(refBus)),Inf) <= lambdaVM(b);
        end

    %V min magnitude finite reduction constraints
        for i = 1:Ns %todo: substitute with matrix multiplication constraint
            rho * ((- un(b) + mpc.bus(b,VMIN)^2*un(refBus)) + (1-beta) * oVm(b)+(- ue(b) + mpc.bus(b,VMIN)^2*ue(refBus))*S_error(:,i)) <= sVm(b,i);
            rho*(-oVm(b)*beta) <= sVm(b,i);
            norm(-rho*(- ue(b) + mpc.bus(b,VMIN)^2*ue(refBus)),Inf) <= lambdaVm(b);
        end

end
disp(['Finished adding constraints, solving model...'])
cvx_end 


elapsedTime = toc;
fprintf('Elapsed time: %f seconds\n', elapsedTime);

end


