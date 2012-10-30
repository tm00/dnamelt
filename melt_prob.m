function prob = melt_prob(varargin)
% MELT_PROB computes the melting probability
% melt_prob(temp,Gamma,seq,bc) computes the melting probability at
%     temperature "temp" and torque "Gamma" for the sequence "seq" with
%     boundary conditions "bc" .
% melt_prob(temp,Gamma,seq,T,bc) computes the melting probability using the
%     4x4 cell array of transfer matrices "T" and with
%     boundary conditions "bc".
% seq is either a 1-dim array of {1,2,3,4} or the name of a fasta file.

  
  temp = varargin{1};
  Gamma = varargin{2};
  seq = varargin{3};
  switch nargin
   case 4
    T = transfer_matrix(temp, Gamma);
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
  % store left and right matrix products
  VL = zeros(N-1,p.ML);
  VR = zeros(p.ML,N-1);
  VL(1,:) = v';
  VR(:,N-1) = TR*v/norm(TR*v);
  for l=2:N-1
    VL(l,:) = VL(l-1,:)*T{seq(l-1),seq(l)};
    VL(l,:) = VL(l,:)/norm(VL(l,:));
    VR(:,N-l) = T{seq(N-l),seq(N-l+1)}*VR(:,N-l+1); 
    VR(:,N-l) = VR(:,N-l)/norm(VR(:,N-l));
  end
  % melting probability function
  df = diag(heaviside(p.xi-12));
  % melting probabilities
  prob = zeros(N,1);
  for l=1:N-1
    prob(l) = (VL(l,:)*df*VR(:,l))/(VL(l,:)*VR(:,l));
  end
  prob(N) = (VL(N-1,:)*TR*df*v)/(VL(N-1,:)*TR*v);
  
function T = Tmorse(x,D,a,R,beta)
  T = exp(-beta.*D.*(exp(-a.*(x-R))-1).^2);

function T = Tstack(x,y,K,alpha,R,beta)
  T = exp(-beta.*K.*(x-y).^2.*exp(-alpha.*(x+y-2*R)));
  
function T = Ttwist(x,y,z,E,h,l0,beta)
  l = sqrt(h^2+x.^2+y.^2-2.*x.*y.*z);
  T = exp(-beta.*E.*(l-l0).^2);
  
function h = heaviside(x)
  h = ones(size(x));
  h(find(x<0)) = 0;