function [fe,fen] = free_energy_n(varargin)
% FREE_ENERGY_N - compute free energies
% free_energy_n(temp,Gamma,seq,bc) computes the free energy as well as the
%  `expectation value free energies' at temperature "temp" and 
%  torque "Gamma" for the sequence "seq" and boundary conditions "bc".
% free_energy_z(temp,z,seq,T,bc) computes the free energies using the
%   4x4 cell array of transfer matrices "T".
    
  temp = varargin{1};
  Gamma = varargin{2};
  seq = varargin{3};
  switch nargin
   case 4
     T = transfer_matrix(temp,Gamma);
     bc = varargin{4};
   case 5
    T = varargin{4};
    bc = varargin{5};
   otherwise
    error('Wrong number of input arguments');
  end
  
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
  
  % transfer matrix for n=N-1
  TR = T{seq(N-1),seq(N)} .* ...
       sparse(diag(p.legw.*p.xi.* Tmorse(p.xi,p.D(seq(N-1)),p.a(seq(N-1)),...
                                         p.r0,beta)));
  % <xi> expectation value function
  %df = diag(p.xi);
  % melting probability function
  df = diag(heaviside(p.xi-12));
  
  % boundary vector
  if strcmp(bc, 'free')
    v = ones(p.ML,1)/sqrt(p.ML);
  elseif strcmp(bc, 'closed')
    v = zeros(p.ML,1);  v(find(p.xi<=13))=1;  v = v/norm(v);
  elseif strcmp(bc, 'open')
    v = zeros(p.ML,1);  v(find(p.xi>13))=1;  v = v/norm(v);
  else
    error([bc ' boundary conditions not supported, choose free, open or closed.'])
  end
  
  VL = zeros(N-1,p.ML);
  VR = zeros(p.ML,N-1);
  VL(1,:) = v';
  VR(:,N-1) = TR*v/norm(TR*v);
  fe_vec = zeros(N-1,1);
  for l=2:N-1
    VL(l,:) = VL(l-1,:)*T{seq(l-1),seq(l)};
    fe_vec(l-1) = (VL(l,:)*v)/(VL(l-1,:)*v);
    VL(l,:) = VL(l,:)/norm(VL(l,:));
    VR(:,N-l) = T{seq(N-l),seq(N-l+1)}*VR(:,N-l+1); 
    VR(:,N-l) = VR(:,N-l)/norm(VR(:,N-l));
  end
  fe_vec(N-1) = (VL(N-1,:)*TR*v)/(VL(N-1,:)*v);
  fe = -sum(log(fe_vec))/(beta*N);
  
  fen = zeros(N,1);
  for n=1:N-1
    fen(n) = -(beta*N)^(-1)*log((VL(n,:)*df*VR(:,n))/(VL(n,:)*VR(:,n))) + ...
             fe;
  end
  fen(N) = -(beta*N)^(-1)*log((VL(N-1,:)*TR*df*v)/(VL(N-1,:)*TR*v)) + fe;

  
function T = Tmorse(x,D,a,R,beta)
  T = exp(-beta.*D.*(exp(-a.*(x-R))-1).^2);

function T = Tstack(x,y,K,alpha,R,beta)
  T = exp(-beta.*K.*(x-y).^2.*exp(-alpha.*(x+y-2*R)));
  
function T = Ttwist(x,y,z,E,h,l0,beta)
  l = sqrt(h^2+x.^2+y.^2-2.*x.*y.*z);
  T = exp(-beta.*E.*(l-l0).^2);
