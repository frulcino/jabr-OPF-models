%TODO quello che ho fatto con P_ij ha senso?
%TODO sto usando la corretta admittance matrix? Devo fare P_ij - P_ji? Cosa
%è delta(k)?
%TODO: recover original variables from c and 
%TODO:come considerare reference bus?
%Ab?
%TODO: non riesce a passare corrente?

define_constants;
mpc = loadcase("case9.m");
%mpc = ext2int(mpc);
[YBUS, YF, YT] = makeYbus(mpc);
G = real(YBUS);
B = complex(YBUS);
refBus = mpc.bus(mpc.bus(:,2) == 3,1); %reference


size_gen = size(mpc.gen);
size_bus = size(mpc.bus);
size_branch = size(mpc.branch);
Ngen = size_gen(1);
Nbus = size_bus(1);
Nbranch = size_branch(1);


%[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
%[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;


%%

cvx_begin
cvx_solver MOSEK
    variable P_G(Ngen)
    variable Q_G(Ngen)
    variable u(Nbus) %c_ii
    variable c(Nbranch) %c_ij
    variable s(Nbranch) %s_ij
    variable P(Nbranch,2) %matrix having for each row P_ij and ji with ij a branc
    variable Q(Nbranch,2) %Q_ij
    
    %total cost
    % total cost
   
    F = sum(mpc.gencost(:, 7)) + sum(mpc.gencost(:, 6) .* P_G) + sum(mpc.gencost(:, 5) .* P_G.^2);
    
    minimize F
    subject to
    
        for b = 1:Nbranch %constraints over branch b
            f = mpc.branch(b,1);%from
            t = mpc.branch(b,2); %to
            [G(f,f),G(f,t),G(t,f),G(t,t)]
            P(b,1) == G(f,f)*u(f) + G(f,t)*c(b) + B(f,t)*s(b) %P_ft b=ft
            Q(b,1) == -B(f,f)*u(f) - B(f,t)*c(b) + G(f,t)*s(b)
            P(b,2) == G(t,t)*u(t) + G(t,f)*c(b) - B(t,f)*s(b)  %P_tf b=ft 
            Q(b,2) == -B(t,t)*u(t) - B(t,f)*c(b) - G(t,f)*s(b)
            norm([2*c(b),2*s(b),u(f)-u(t)]) <= u(f) + u(t); %jabr inequality in convex form

            
            %max power flow constraint (non in case description)?
            
        end
    
        for b = 1:Nbus %constraints over buses
            u(b) >= 0;
            u(b) >= mpc.bus(b,VMIN)^2*u(refBus);
            u(b) <= mpc.bus(b,VMAX)^2*u(refBus); %^2?
            
            N = Nmat(b,mpc);
            
            %power flow constraint at bus b
            isgen = 0; %is there a gen at bus b?
            for g = 1:Ngen
                if b == mpc.gen(g,1);
                    isgen=1;
                    ngen = g;
                    break
                end
            end
            
            %constraints su Q
            if isgen == 1
               sum(sum(N .* P)) == P_G(g) - mpc.bus(b,PD); %    segni?
               sum(sum(N .* Q)) == Q_G(g) - mpc.bus(b,QD)
               P_G(g) >= mpc.gen(g,PMIN);
               P_G(g) <= mpc.gen(g,PMAX);
            else
                sum(sum(N .* P)) == - mpc.bus(b,PD); %che sia qui il problrms
                %erchè non c'è un magico pg?
                %sto considerando injected power or outjected power?
                sum(sum(N .* Q)) == - mpc.bus(b,QD);
            end
            
        end
            
            %{
            %calculates power injection at bus b
            sumP = 0;
            sumQ = 0;
            neighbours = k(b,mpc);
            for t = neighbours
                si = idside(b,t,mpc) %branch index of (b,t)
                ssbt = ss(b,t,mpc) %sign of branch (b,t)
                sumP = sumP + G(b,b)*u(b) + G(b,t)*c(si) + ssbt*B(b,t)*s(si);
                sumQ = sumQ + -B(b,b)*u(b) - B(b,t)*c(si) + ssbt*B(b,t)*s(si);
            end
               
            
            
            %power flow constraint at bus b
            isgen = 0; %is there a gen at bus b?
            for g = 1:Ngen
                if b == mpc.gen(g,1);
                    isgen=1;
                    ngen = g;
                    break
                end
            end
            
            %constraints su Q
            if isgen == 1
               sumP == P_G(g) - mpc.bus(b,PD); %    segni?
            else
                sumP == - mpc.bus(b,PD);
            end
            
                
            
        end

        for g = 1:Ngen
            P_G(g) >= mpc.gen(g,PMIN);
            P_G(g) <= mpc.gen(g,PMAX);
        end
        
            
        %}

cvx_end


%% helper functions


