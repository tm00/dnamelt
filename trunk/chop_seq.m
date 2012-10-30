function A = chop_seq(seq, N0, d)
% CHOP_SEQ - Chop up a sequence
%   A = chop_seq(seq, N0, d) chops up the sequence seq in windows of
%   length N0 and overlap 2d; A is an (N0+2d) x n array, n the number of
%   windows; the last window will have larger overlap if
%   rem(length(seq),N0) is not 0;
  
  N = length(seq); 
  if rem(N,N0) == 0
    Nwin = N/N0;
  else
    Nwin = floor(N/N0)+1;
  end
  A = zeros(N0+2*d,Nwin);
  A(:,1) = seq(1:N0+2*d);
  for k = 2:Nwin-1
    A(:,k) = seq((k-1)*N0-d+1:k*N0+d);
  end
  A(:,Nwin) = seq(end-N0-2*d+1:end);