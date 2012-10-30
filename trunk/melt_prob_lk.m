function [prob,u,un] = melt_prob_lk(temp,sigma,seq,bc,Gamma,omega)
% MELT_PROB_LK - computes the melting probability
% melt_prob_lk(temp,sigma,seq,bc) computes the melting probability at
%     temperature "temp" and superhelical density "sigma" for the
%     sequence "seq" with boundary conditions bc.
%     Gamma is a vector of torque values over which to maximize the free
%     energy
%     omega is a vector of values to compute the complex integrals
% seq is either a 1-dim array of {1,2,3,4} or the name of a fasta file.
% Output: "prob" contains the melting probability, "in" the inequivalence
% of ensembles measure (prob_lk / prob_tq)  

  
 
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
  % get model and integration parameters
  p = getpref('DNA_melt');
  beta = 1./(p.kB.*temp);
  % compute linking number
  Lk0 = sum(p.theta0(seq(1:N-1)+4*(seq(2:N)-1)));
  alpha = (1+sigma)*Lk0/N;

  % find maximum of free energy
  f=zeros(size(Gamma));
  fn=zeros(size(Gamma));
  for k=1:length(Gamma)
    f(k) = free_energy(temp,Gamma(k),seq,bc) + alpha*Gamma(k);
  end
  % cubic splines interpolation of f
  ppf = spline(Gamma,f);
  % derivative
  ppfd = mmppder(ppf);
  % maximum
  Gamma0 = fzero(@(Gamma) ppval(ppfd,Gamma), -0.05);
  
  % integrate along same line for all n
  z = Gamma0 + i*omega;
  
  [f0, f0n] = free_energy_n(temp,Gamma0,seq,bc);
  %[f0, f0n] = free_energy_n_theta(temp,Gamma0,seq);
  f0 = f0 + alpha*Gamma0;
  f0n = f0n + alpha*Gamma0;

  ur=zeros(size(z));
  ui=zeros(size(z));
  urn=zeros(length(z),N);
  uin=zeros(length(z),N);
  for k=1:length(ur)
    [fe, fen] = free_energy_n(temp,z(k),seq,bc);
    %[fe, fen] = free_energy_n_theta(temp,z(k),seq);
    ur(k) = real(fe) + alpha*Gamma0 - f0;
    ui(k) = imag(fe) + alpha*omega(k);
    
    urn(k,:) = real(fen) + alpha*Gamma0 - f0n;
    uin(k,:) = imag(fen) + alpha*omega(k);
  end
  % denominator integral
  pr = spline(omega,exp(-beta*N*ur).*cos(beta*N*ui));
  pri = mmppint(pr,0);
  int = ppval(pri,omega(end));
  % nominator integrals
  intn = zeros(N,1);
  for n=1:N
    prn = spline(omega,exp(-beta*N*urn(:,n)).*cos(beta*N*uin(:,n)));
    prni = mmppint(prn,0);
    intn(n) = ppval(prni,omega(end));
  end
  prob = exp(-beta*N*(f0n-f0)).*intn/int;

  u = exp(-beta*N*ur).*exp(i*beta*N*ui);
  un = exp(-beta*N*urn).*exp(i*beta*N*uin);
