%create error scenarios file

% file_name = "C:\Users\Gabor\Documents\tesi_magistrale\data\previsioni\eolico\error_scenarios.csv";
% error_scenarios = readmatrix(file_name);
% error_scenarios = error_scenarios(:);
% e_mean = mean(error_scenarios);
% e_std = std(error_scenarios);
% %oss questi errori fanno un po pena perchè sono lin dipendenti
% G_error_1 = (error_scenarios - e_mean)/e_std * 200;
% G_error_2 = G_error_1;
% G_error_3 = G_error_1*3/2;
% S_error = [G_error_1, G_error_2, G_error_3];
% save("winderrors.mat", "S_error");

clc; clear all; close all;%
%warning off
[PGn, QGn, PGe, QGe, cn, sn, un, ce, se, ue, Pn, Qn, Pe, Qe, Fop, Frisk, F] = JDROOPF2("case118wind3.m","winderrors.mat",4, 3, 0.01, 20, 0.02)