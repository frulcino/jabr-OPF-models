%function [PGn, QGn, PGe, QGe, cn, sn, un, ce, se, ue, Pn, Qn, Pe, Qe, Fop, Frisk, F] = StochasticJabrOPF(casename,errors_name,Ns, Nxi, beta, rho, Varepsilon)
clc; clear all; close all;
warning off
[casename, errors_name, Ns, Nxi, beta, rho, Varepsilon] = deal("case118wind3.m","winderrors.mat", 4, 3, 0.01, 20, 0.02)

%StochasticJabrOPF 
% casename: string name of matpower case. Must also contain
%   -mpc.Fmax contains the power flow limit in p.u. for each branch in the
%             network, in the same order as mpc.branch
%   -mpc.wind has two columns, each row correspons to a windgenerator
%             at bus indicated in the first column and generating a nominal
%             value located in the second column
% errors_name: string name of .mat containing error samples for
%                       each renewable generator.
%   - has one column for every eolic generator
%   - every row is a sample (has at least Ns rows)
%   - has variable calles S_error
define_constants;
mpc = loadcase(casename);
tic 

%% Input Data Section
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
%from CESI wind farm
data = load(errors_name);
error_scenarios = data.S_error
Err = error_scenarios(1:Ns,:); %has as Nxi columns (wind turbunes) and Ns rows (scenarios)
%S_error = transpose(error_scenarios);

% Got the load demand at each bus
Pd = bus(:,3);

% Got nominal wind power injection at each bus
wind = mpc.wind;
Nwind = size(wind,1);
Pwind = zeros(Nbus,1);


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
disp(['number of buses:    ', num2str(Nbus)]);
disp(['number of lines:   ', num2str(Nbranch)]);
disp(['number of generators:    ', num2str(Ngen)]);
disp(['number of scenarios for distributionally robust stochastic OPF:  ', num2str(Ns)]);
disp(['number of renewable generatos:   ', num2str(Nxi)]);
disp(['wind injection locations: ', '  bus 1,  ','    bus 9,  ','    bus 26']); %todo change
disp(['wind nominal injection:   ', '  500 (MW), ','  500 (MW), ','  800 (MW)']); %todo change
disp(['Wasserstein distance:    ', num2str(Varepsilon)]);

%% Set up Gurobi parameters
%gurobi_params = struct();
%gorobi_params.BarHomogeneous = 1;
%% Build Model
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

%voltage mangitude constraint finite reduction variable
variable lambdaVm(Nbus) %for minimum voltage
variable sVm(Nbus,Ns)

variable lambdaVM(Nbus)
variable sVM(Nbus,Ns)

%PG magnitude constraing finite reduction variable
variable lambdaPGm(Ngen) %for minimum voltage
variable sPGm(Ngen,Ns)
variable lambdaPGM(Ngen)
variable sPGM(Ngen,Ns)

Fop = 0;
% operation costs (avarage over dataset)
for i = 1:Ngen
    %{
    for j = 1:Ns
        Fop = Fop + sum( gencost(i,5)*(PGn(i) + PGe(i,:).*Err(j,:)).^2 + gencost(i,6)*(PGn(i) + PGe(i,:).*Err(j,:)) + gencost(i,7));
    end
    %}
    % sum(PGn(i) + Err*diag(PGe(i,:)),2) the second term in the sum corresponds to
    % the matrix having as i-throw the correction to PGn for the i-th scenario 
    Fop = Fop + sum( gencost(i,5)*(PGn(i) + sum(Err*diag(PGe(i,:)),2) ).^2 + gencost(i,6)*(PGn(i) + sum(Err*diag(PGe(i,:)),2))  + gencost(i,7)); %Err*diag(PGe(i,:)) is a component wise multiplication of the colums of Err
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
Pn(:,1,:) + Pn(:,2,:) >= 0; 
%nominal power generation constraints
PGn(:) <= mpc.gen(:,PMAX)
PGn(:) >= mpc.gen(:,PMIN)
%V magnitude nominal constraints
un(:) >= 0;
un(:) >= mpc.bus(:,VMIN).^2*un(refBus);
un(:) <= mpc.bus(:,VMAX).^2*un(refBus);
%un(refBus) == 1

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
    Pe(b,1) == Gff(b)*ue(f) + Gft(b)*ce(b) + Bft(b)*se(b); % P_ft b=ft
    Qe(b,1) == -Bff(b)*ue(f) - Bft(b)*ce(b) + Gft(b)*se(b);
    Pe(b,2) == Gtt(b)*ue(t) + Gtf(b)*ce(b) - Btf(b)*se(b); % P_tf b=ft
    Qe(b,2) == -Btt(b)*ue(t) - Btf(b)*ce(b) - Gtf(b)*se(b);



    %power loss finite reduction constraints (H=0 d=0)
    %APL(b) == (-((Gft(b)+Gtf(b))*ce(b,:)+(Bft(b)-Btf(b))*se(b,:)-Gff(b)*ue(f,:)-Gtt(b)*ue(t,:))); 
    %bPL(b) == (-((Gft(b)+Gtf(b))*cn(b,:)+(Bft(b)-Btf(b))*sn(b,:)-Gff(b)*un(f,:)-Gtt(b)*un(t,:)));
    for i = 1:Ns %todo: substitute with matrix multiplication constraint
        rho*((-((Gft(b)+Gtf(b))*cn(b,:)+(Bft(b)-Btf(b))*sn(b,:)-Gff(b)*un(f,:)-Gtt(b)*un(t,:)))+(1-beta)*oPL(b)+(-((Gft(b)+Gtf(b))*ce(b,:)+(Bft(b)-Btf(b))*se(b,:)-Gff(b)*ue(f,:)-Gtt(b)*ue(t,:)))*transpose(Err(i,:))) <= sPL(b,i);
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
            windIdx = g;
            break
        end
    end

    if isgen == 1
       %nominal injection constraints
       sum(sum(N .* Pn)) == PGn(gIdx) - mpc.bus(b,PD); 
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
           rho*((PGn(gIdx) - mpc.gen(gIdx,PMAX))+(1-beta)*oPGM(gIdx)+PGe(gIdx)*transpose(Err(i,:))) <= sPGM(gIdx,i);
           rho*(-oPGM(gIdx)*beta) <= sPGM(gIdx,i);
           norm(-rho*PGe(gIdx),Inf) <= lambdaPGM(gIdx);
       end

       % PG min magnitude finite reduction constraints
       %APGm(gIdx) == (-PGe(gIdx));
       %bPGm(gIdx) == (-PGn(gIdx) + mpc.gen(gIdx,PMIN));
       for i = 1:Ns %todo: substitute with matrix multiplication constraint
           rho*((-PGn(gIdx) + mpc.gen(gIdx,PMIN))+(1-beta)*oPGm(gIdx)+(-PGe(gIdx))*transpose(Err(i,:))) <= sPGm(gIdx,i);
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
            rho*((un(b) - mpc.bus(b,VMAX)^2*un(refBus))+(1-beta)*oVM(b)+(ue(b) - mpc.bus(b,VMAX)^2*ue(refBus))*transpose(Err(i,:))) <= sVM(b,i);
            rho*(-oVM(b)*beta) <= sVM(b,i);
            norm(-rho*(ue(b) - mpc.bus(b,VMAX)^2*ue(refBus)),Inf) <= lambdaVM(b);
        end

    %V min magnitude finite reduction constraints
        for i = 1:Ns %todo: substitute with matrix multiplication constraint
            rho * ((- un(b) + mpc.bus(b,VMIN)^2*un(refBus)) + (1-beta) * oVm(b)+(- ue(b) + mpc.bus(b,VMIN)^2*ue(refBus))*transpose(Err(i,:))) <= sVm(b,i);
            rho*(-oVm(b)*beta) <= sVm(b,i);
            norm(-rho*(- ue(b) + mpc.bus(b,VMIN)^2*ue(refBus)),Inf) <= lambdaVm(b);
        end

end
disp(['Finished adding constraints, solving model...'])
cvx_end 


elapsedTime = toc;
fprintf('Elapsed time: %f seconds\n', elapsedTime);

%end

