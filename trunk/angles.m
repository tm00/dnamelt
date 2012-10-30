function [dev_theta, link_diff, hel_dens] = angles(prob_theta,seq)
% ANGLES - compute various angular quantities given the expectation of theta
%   [dev_theta, link_diff, hel_dens, prob_phi] = angles(prob_theta,seq) 
%   takes the expectation value of the angular differences theta_n =
%   phi_(n+1)-phi_n (in radians) for the sequence seq and computes:
%     * dev_theta: deviation of theta from (step dependent) theta0
%     * link_diff: the linking difference
%     * hel_dens: the helical density
%     * prob_phi: the corresponding expectation of phi, in degrees
%   seq is either a 1-dim array of {1,2,3,4} or the name of a fasta file.
  
  % if seq is not a numeric array, it should be a fasta file which is
  % then converted to a numeric array
  if ~isnumeric(seq)
    [h seq] = fastaread(seq);
    seq = nt2int(seq,'ACGTOnly', true);
    if length(find(seq==0)) ~= 0
      seq = seq(find(seq~=0));
      sprintf('%d nucleotides of unknown type were deleted\n',...
              length(find(seq==0)));
    end
  end
  N = length(seq);
  
  % get model and integration parameters
  p = getpref('DNA_melt');
  
  % deviation of theta
  dev_theta = zeros(N,1);
  dev_theta(1:N-1) = (prob_theta(1:N-1) - p.theta0(seq(1:N-1)+4*(seq(2:N)-1))) ...
      ./p.theta0(seq(1:N-1)+4*(seq(2:N)-1));
  dev_theta(N) = dev_theta(N-1);
  
  
  % linking difference
  link_diff = sum((prob_theta(1:N-1) ... 
                   -p.theta0(seq(1:N-1)+4*(seq(2:N)-1))))/(2*pi); 
  
  % helical density
  hel_dens = 2*pi*link_diff/sum(p.theta0(seq(1:N-1)+4*(seq(2:N)-1)));
  
  % expectation of phi
  prob_phi = zeros(N,1);
  for l=2:N
    prob_phi(l) = rem(prob_phi(l-1)+prob_theta(l-1),2*pi);
  end
  prob_phi = (180/pi)*prob_phi;
