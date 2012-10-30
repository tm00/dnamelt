function ppi = mmppint(pp,c)
% MMPPINT - Cubic spline integral interpolation
%   Copied from "Mastering Matlab, 6th Ed." p. 292   
  if prod(size(c))~=1
    error('C must be a scalar.')
  end
  [br,co,npy,nco]=unmkpp(pp);
  sf=nco:-1:1;
  ico=[co./sf(ones(npy,1),:) zeros(npy,1)];
  nco=nco+1;
  ico(1,nco)=c;
  for k=2:npy
    ico(k,nco)=polyval(ico(k-1,:),br(k)-br(k-1));
  end
  ppi=mkpp(br,ico);
