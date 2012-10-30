function T = transfer_matrix_fg(temp, Gamma, F, G)
% TRANSFER_MATRIX_FG - Create transfer matrix for every possible dinucleotide step
% transfer_matrix_fg(temp, Gamma, F, G) creates a 4x4 cell array of transfer 
%    matrices at temperature "temp" and torque "Gamma" for the test
%    function "F(r)*G(cos theta)" (F,G are given as strings).

  if isstr(F)
    f = inline(F);
  else
    f = F;
  end
  
  if isstr(G)
    g = inline(G);
  else
    g = G;
  end
 % get model and integration parameters
  p = getpref('DNA_melt');
  beta = 1./(p.kB.*temp);
  cutoff = 10^(-6); % values below cutoff are set to 0 to obtain sparse matrices
  
  % define some indices
  n1 = 1:p.ML;
  n2 = 1:p.ML^2;
  [i1,i2] = ind2sub([p.ML p.ML],n2); 
  n3 = 1:p.ML^2*p.MC;
  [j1,j2,j3] = ind2sub([p.ML p.ML p.MC],n3);
  n4 = 1:p.ML*p.MC;
  [k1,k2] = ind2sub([p.MC p.MC],n4);
  
  % create transfer matrices
  T = cell(4,4);
  for s=1:4
    S1diag = p.legw(n1).*p.xi(n1).* feval(f,p.xi(n1)) ...
             .* Tmorse(p.xi(n1),p.D(s),p.a(s),p.r0,beta);
    S1 = sparse(n1,n1,S1diag,p.ML,p.ML);
    for t=1:4
      S2 = Tstack(p.xi(i1),p.xi(i2),p.K(s,t),p.alpha(s,t),p.r0,beta);
      S2 = reshape(S2,[p.ML,p.ML]);
      
      S3 = p.chebw .* feval(g,p.chebz(j3)) ...
           .*Ttwist(p.xi(j1),p.xi(j2),p.chebz(j3),p.E(s,t),p.h,...
                    p.l0(s,t),beta).* exp(beta.*Gamma.*acos(p.chebz(j3)));
      S3 = reshape(S3, [p.ML p.ML p.MC]);
      S3 = sum(S3,3);
      tt = (S1*S2).*S3;
      tt(find(tt<cutoff)) = 0;
      T{s,t} = sparse(tt);
    end
  end
  
function T = Tmorse(x,D,a,R,beta)
  T = exp(-beta.*D.*(exp(-a.*(x-R))-1).^2);
  
function T = Tstack(x,y,K,alpha,R,beta)
  T = exp(-beta.*K.*(x-y).^2.*exp(-alpha.*(x+y-2*R)));
  
function T = Ttwist(x,y,z,E,h,l0,beta)
  l = sqrt(h^2+x.^2+y.^2-2.*x.*y.*z);
  T = exp(-beta.*E.*(l-l0).^2);