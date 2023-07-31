# jabr-OPF-models
Implemented in matlab using cvx the Jabr inequality model (jabr2.m) as definied in "Mathematical Programming formulations for the Alternating Current Optimal Power Flow problem" by Bienstock. Also implemented the linearization introduced by Bienstock in the Gurobi seminar.

# Scripts descrioption
# Mains scripts
Jabr2.m: Implementation of the Jabr inequality model  as definied in "Mathematical Programming formulations for the Alternating Current Optimal Power Flow problem" by Bienstock.

Jabr2LinBiestock.m the Jabr inequality model using the power loss inequality introduced by Bienstock and Munos in "On linear relaxation of OPF problems"

StochasticJabr.m Stochastic OPF model based on Jabr OPF, following "Data-driven distributionally robust optimization using the Wassestein metric: performance guarantees and tractable reformulation" by Esfahani and Kuhn and "Data_Based Distributionally Robust Stochastic Optimal Power Flow" by Yi Guo et al.
