% Toolbox DNA_MELT
% Copyright (C) 2005 Tom Michoel
% Version 2006/10/12
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%%
%% This set of programs is free software; you can redistribute it and/or 
%% modify it under the terms of the GNU General Public License as
%% published by the Free Software Foundation; either version 2 of the
%% License, or (at your option) any later version.
%%
%% This set of programs is distributed in the hope that it will be
%% useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%% General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License along
%% with this program; if not, write to the Free Software Foundation, Inc.,
%% 675 Mass Ave, Cambridge, MA 02139, USA.
%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%%
%% You are welcome to use the code for your research under the terms of the
%% license. However, please acknowledge its use with the following
%% citation:
%%
%%   Tom Michoel and Yves van de Peer, "A helicoidal transfer matrix model
%%      for inhomogeneous DNA melting", Phys. Rev. E 73 011908 (2006),
%%      arXiv:q-bio/0507036
%%
%% Author contact information:
%%
%%  Tom Michoel
%%
%%  Bioinformatics & Evolutionary Genomics
%%  VIB/Ghent University
%%  Technologiepark 927
%%  B-9000 Gent, Belgium
%%
%%  Email: tom.michoel@psb.ugent.be
%%  URL:   <http://www.psb.ugent.be/~tomic/>
%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
% Before using this toolbox, add DNA_melt directory to path and set the
% model parameters:
%
% >> addpath DIR;
% >> DNA_melt_preferences;
%
% To compute the opening probabilities for sequence "seq" at temperature
% "temp", torque "Gamma", and boundary conditions "bc":
%
% >> prob = melt_prob(temp,Gamma,seq,bc);
%
% "seq" can be a fasta-file or sequence of {1,2,3,4}; for a fasta-file,
% the Bioinformatics toolbox needs to be installed.
%
% If the sequence is about 10^6 base pairs or more, use
%
% >> prob = melt_prob_long(temp,Gamma,seq, N0, d);
%
% "N0" is the window length and "2d" the overlap.
%
% If you want to analyze different sequences at the same temperature and
% torque, first run
%
% >> T = transfer_matrix(temp, Gamma);
%
% followed by (for every sequence)
%
% >> prob = melt_prob(temp,Gamma,seq,T);
%
% To compute the data for plotting a melting curve and melt map, use
%
% >> [temp,gamma,diff,tempr] = melt_curve(tempI, tempF, step, Gamma, seq);
%
% To compute more general expectation values (such as base pair opening
% fluctuations, degrees of untwisting, ...) for testfunctions
% F(r)*G(theta):
%
% >> prob = exp_val_fg(temp,Gamma,seq,F,G,bc);
%
% If more than one sequence is analyzed at the same temperature and
% torque:
% 
% >> T = transfer_matrix(temp, Gamma);
% >> Tx = transfer_matrix_fg(temp,Gamma,F,G);
% >> prob = exp_val_fg(temp,Gamma,seq,F,G,T,Tx,bc);
%
% If the degrees of untwisting are stored in "prob_theta", run
%
% >> [dev_theta, link_diff, hel_dens, prob_phi] = angles(prob_theta,seq);
%
% to obtain the deviations of theta from its zero temperature average,
% the linking difference, the helical density, and the expectation of phi.
%
% To compute the melting probability at fixed helical density sigma:
%
% >> [prob,u,un] = melt_prob_lk(temp,sigma,seq,bc,Gamma,omega);
%
% "prob" contains the probability and "u" and "un" the integrands for the
% denominator, resp. nominators in eq. (10) of the paper. In the input,
% Gamma is a vector containing the torque values for which to compute the
% free energies (the maximum is found by differentiating a cubic splines
% interpolation), and omega is a vector containing the imaginary part of
% the complex integration variables (see eq. (9-10)).
%
% Computation of the melting probability at fixed helical density relies
% on free energy computations for fixed torque. The following functions
% are provided 
%
% >> fe = free_energy(temp,Gamma,seq,bc);
%
% computes the free energy of the system
%
% >> [fe,fen] = free_energy_n(temp,Gamma,seq,bc);
%
% computes the free energy fe as well as an array of position dependent
% free energies fen such that prob(n) = exp(-beta(fen(n)-fe)) is the melting
% probabiliy.
%
% >> [fe,fen] = free_energy_n_fg(temp,Gamma,seq,F,G,bc);
%
% same computation for general test function F*G, i.e., exp_val_fg(n) =
% exp(-beta(fen(n)-fe))
%
% Like for melt_prob and exp_val_fg, all 3 functions can be called with
% optional transfer matrix arguments T, Tx.
%
% To compute the free energy in a temperature interval use
%
% >> [temp,fe] = free_energy_curve(tempI, tempF, step, Gamma, seq);
%
% To have a score of how well the opening probability identifies a
% biological property such as coding/non-coding, write the 
% property as a vector "cod" of 0 and 1, e.g., coding=1, non-coding=0,
% then run  
%
% >> score = overlap(prob,cod);
%
% Specifically for the coding property, the vector "cod" can be extracted
% from a GenBank file if the Bioinformatics toolbox is installed.
%
% >> [cod,A] = extract_coding_gb(GenbankFileName);
%
%
% All functions have a short help section with more info:
%
% >> help functionname;
%
% For an elaborate usage example, see the file "DNAmeltPaperFigs.m" which
% generates all the figures in our paper referenced above.
%
