function [temp,gamma,diff,tempr] = melt_curve(tempI,tempF,step,Gamma,seq,bc)
% MELT_CURVE - Compute a melting curve 
%   [temp,TT,diff,tempr] = melt_curve(tempI, tempF, step, Gamma, seq)
%   computes a melting curve for sequence "seq" between temperatures
%   "tempI" and "tempF" with step size "step" and torque "Gamma" and
%   boundary conditions "bc".
%   It returns:
%     * temp  : temperature values
%     * gamma : fraction of open base pairs, piecewise polynomial
%     * diff  : derivative of gamma w.r.t. temp, piecewise polynomial
%     * tempr : for each sequence position, the temperature where the
%               melting probability exceeds 50%

sz = size(seq);
% sequence length
N = sz(2);

nstep = (tempF - tempI)/step;
TT = zeros(nstep+1,1);
tempr = zeros(N,1);
done = zeros(N,1);

for k = 0:nstep
  temp = tempI + (k/nstep)*(tempF-tempI);
  prob = melt_prob(temp,Gamma,seq,bc);
  TT(k+1) = sum(prob)/N;
  tempr(find(prob>=0.5 & done==0)) = temp;
  done(find(prob>=0.5 & done==0)) = 1;
end
tempr = [tempF;tempr];
tempr(end) = tempF;
temp = tempI + (0:nstep)*(tempF-tempI)/nstep;
diff = (TT(2:nstep+1)-TT(1:nstep))*nstep/(tempF-tempI);
diff = [diff; diff(end)];
gamma = spline(temp, TT);
%diff2 = mmppder(gamma);

