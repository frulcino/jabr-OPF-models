
define_constants;
mpc = loadcase("case118.m");
mpc = ext2int(mpc);
[Yff, Yft, Ytf, Ytt] = Ybranch(mpc);
Gff = real(Yff);
Gft = real(Yft);
Gtf = real(Ytf);
Gtt = real(Ytt);
Bff = imag(Yff);
Bft = imag(Yft);
Btf = imag(Ytf);
Btt = imag(Ytt);
refBus = mpc.bus(mpc.bus(:,2) == 3,1); %reference Bus

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
cvx_solver GUROBi
    disp("defining variables...")
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
    disp("adding branch constraints...")
        for b = 1:Nbranch %constraints over branch b
            f = mpc.branch(b,1);%from
            t = mpc.branch(b,2); %to %[Gff(b),Gft(b),Gtf(b),Gtt(b)]
            %s ha segni diversi paper per contare che s(i,j) = -s(j,i)
            %invece che imporlo come constraint
            
            %questi contraint sono coerenti con la posizione di
            %trasformatori?
            P(b,1) == Gff(b)*u(f) + Gft(b)*c(b) + Bft(b)*s(b); % P_ft b=ft
            Q(b,1) == -Bff(b)*u(f) - Bft(b)*c(b) + Gft(b)*s(b);
            P(b,2) == Gtt(b)*u(t) + Gtf(b)*c(b) - Btf(b)*s(b); % P_tf b=ft
            Q(b,2) == -Btt(b)*u(t) - Btf(b)*c(b) - Gtf(b)*s(b);
            
            %jabr inequality Bienstock Linear Form:
            %norm([2*c(b),2*s(b),u(f)-u(t)]) <= u(f) + u(t); becomes:
            P(b,1) + P(b,2) >= 0;
            
            
            %voltage angle difference constraints (if angles are smaller than
            %pi/4
            %tan(mpc.branc(b,ANGMIN)/(2*pi))*c(b) <= s(b)
            %s(b) <= tan(mpc.branc(b,ANGMAC)/(2*pi))*c(b)
            c(b) >= 0 %voltage angle difference constraints
            
            %not sure if necessary constraints:
            %c(b) <= mpc.bus(f,VMAX)*mpc.bus(t,VMAX) %magnitude constraint, ridondante?
            %s(b) <= mpc.bus(f,VMAX)*mpc.bus(t,VMAX)
            %s(b) >= -mpc.bus(f,VMAX)*mpc.bus(t,VMAX)
            
           
            
        end
    
       disp("adding bus constraints...")
       
        for b = 1:Nbus 
            %voltage constraints (with conversion form p.u. constraints)
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
               %disp("QMIN:" + num2str(mpc.gen(g,QMIN)));
               if mpc.gen(g,QMIN) ~= -Inf
                    Q_G(g) >= mpc.gen(g,QMIN);
               end
               %disp("QMAX:" + num2str(mpc.gen(g,QMAX)));
               if mpc.gen(g,QMAX) ~= Inf
                    Q_G(g) <= mpc.gen(g,QMAX);
               end
               
               
            else
                sum(sum(N .* P)) == - mpc.bus(b,PD); %sto considerando injected power or outjected power?
                sum(sum(N .* Q)) == - mpc.bus(b,QD);
            end
            
        end
            

cvx_end


%% helper functions


