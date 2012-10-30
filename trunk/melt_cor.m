function [prob,cor] = melt_cor(n,temp,Gamma,seq,bc)
% MELT_COR computes the melting correlation between one site and all the others
% melt_cor(n,temp,Gamma,seq,bc) computes the melting correlation of site
%     "n" with all the others at temperature "temp" and torque "Gamma"
%     for the sequence "seq" with boundary conditions "bc" .
% correlation is defined as <r_n r_m>-<r_n><r_m>
% seq is either a 1-dim array of {1,2,3,4} or the name of a fasta file.

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

  % transfer matrices
  T1 = transfer_matrix(temp,Gamma);
  T2 = transfer_matrix_fg(temp,Gamma,'x','1');
  
  % expectation values <r_n>
  %prob = exp_val_fg(temp,Gamma,seq,'x','1',bc);
  prob=melt_prob(temp,Gamma,seq,bc);
  
  % next compute <r_m>' (see paper for explanation)
  % transfer matrix for n=N-1
  TR = T1{seq(N-1),seq(N)} .* ...
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

%  df = diag(p.xi);
  df = diag(heaviside(p.xi-12));
  id = speye(N);
  % store left and right matrix products
  VL = zeros(N-1,p.ML);
  VR = zeros(p.ML,N-1);
  VL(1,:) = v';
  VR(:,N-1) = TR*v/norm(TR*v);
  for l=2:N-1
    VL(l,:) = VL(l-1,:)*(id(l-1,n)*(T2{seq(l-1),seq(l)}-T1{seq(l-1),seq(l)}) ...
                         + T1{seq(l-1),seq(l)});
    VL(l,:) = VL(l,:)/norm(VL(l,:));
    VR(:,N-l) = (id(N-l,n)*(T2{seq(N-l),seq(N-l+1)}-T1{seq(N-l),seq(N-l+1)}) ...
                 + T1{seq(N-l),seq(N-l+1)})*VR(:,N-l+1); 
    VR(:,N-l) = VR(:,N-l)/norm(VR(:,N-l));
  end
  probpr = zeros(N,1);
  for l=1:N-1
    probpr(l) = (VL(l,:)*df*VR(:,l))/(VL(l,:)*VR(:,l));
  end
  probpr(N) = (VL(N-1,:)*TR*df*v)/(VL(N-1,:)*TR*v);
  
  % correlation
  cor = prob(n)*(probpr - prob);
  
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