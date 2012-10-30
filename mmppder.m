function ppd = mmppder(pp)
% MMPPDER - Cubic Spline Derivative Interpolation
%   Copied from "Mastering Matlab, 6th Ed." p. 294
  [br,co,npy,nco] = unmkpp(pp);
  sf = nco-1:-1:1;
  dco = sf(ones(npy,1),:).*co(:,1:nco-1);
  ppd = mkpp(br,dco);
