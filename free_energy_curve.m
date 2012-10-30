function [temp,fe] = free_energy_curve(tempI, tempF, step, Gamma, seq)
% FREE_ENERGY_CURVE - Compute free energy in temperature interval
%  [temp,fe] = free_energy(tempI, tempF, step, Gamma, seq) computes the
%  free energy for all temperatures between "tempI" and "tempF" with step
%  size "step" at torque "Gamma" for the sequence "seq", and returns a
%  cubic splines interpolation "fe" of the computed values (i.e., a
%  piecewise polynomial).
%  The result can be plotted by >>plot(temp,ppval(fe,temp));
  
  % if seq is not a numeric array, it should be a fasta file which is
  % then converted to a numeric array
  if ~isnumeric(seq)
    [h seq] = fastaread(seq);
    seq = nt2int(seq,'ACGTOnly', true);
    if length(find(seq==0)) ~= 0
      warning('%d nucleotides of unknown type are deleted\n', ...
              length(find(seq==0)));
      seq = seq(find(seq~=0));
    end
  end
  N = length(seq);
  
  nstep = (tempF - tempI)/step;
  fev = zeros(nstep+1,1);
  for k = 0:nstep
    temp = tempI + (k/nstep)*(tempF-tempI)
    fev(k+1) = free_energy(temp,Gamma,seq,bc);
  end

  temp = tempI + (0:nstep)*(tempF-tempI)/nstep;
  fe = spline(temp,fe);
