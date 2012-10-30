function fe = free_energy(varargin)
% FREE_ENERGY - Compute the free energy 
%  fe = free_energy(temp, Gamma, seq, bc) computes the free energy for
%  sequence "seq" at temperature "temp", torque "Gamma" and boundary conditions "bc".
%  fe = free_energy(temp, Gamma, seq, T, bc) computes the free energy using
%  the 4x4 cell array "T" of transfer matrices
  
  
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

  fe_vec = zeros(N-1,1);
  V = zeros(N-1,p.ML);
  V(1,:) = v';
  for l=2:N-1
    V(l,:) = V(l-1,:)*T{seq(l-1),seq(l)};
    fe_vec(l-1) = (V(l,:)*v)/(V(l-1,:)*v);
    V(l,:) = V(l,:)/norm(V(l,:));
  end
  fe_vec(N-1) = (V(N-1,:)*TR*v)/(V(N-1,:)*v);
  fe = -sum(log(fe_vec))/(beta*N);

function T = Tmorse(x,D,a,R,beta)
  T = exp(-beta.*D.*(exp(-a.*(x-R))-1).^2);

function T = Tstack(x,y,K,alpha,R,beta)
  T = exp(-beta.*K.*(x-y).^2.*exp(-alpha.*(x+y-2*R)));
  
function T = Ttwist(x,y,z,E,h,l0,beta)
  l = sqrt(h^2+x.^2+y.^2-2.*x.*y.*z);
  T = exp(-beta.*E.*(l-l0).^2);
