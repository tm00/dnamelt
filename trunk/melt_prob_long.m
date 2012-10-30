function prob = melt_prob_long(varargin)
% MELT_PROB_LONG computes the melting probability for long sequences 
% melt_prob_long(temp,Gamma,seq,N0,d,bc) computes the melting probability at
%     temperature "temp" and torque "Gamma for the sequence "seq" using
%     windows of length "N0" and overlap "2d" and boundary conditions "bc".
% melt_prob(temp,Gamma,seq,T,N0,d,bc) computes the melting probability using
%     the 4x4 cell array of transfer matrices "T" and boundary conditions "bc".
% seq is either a 1-dim array of {1,2,3,4} or the name of a fasta file.

  temp = varargin{1};
  Gamma = varargin{2};
  seq = varargin{3};
  switch nargin
   case 6
    T = transfer_matrix(temp, Gamma);
    N0 = varargin{4};
    d = varargin{5};
    bc = varargin{6};
   case 6
    T = varargin{4};
    N0 = varargin{5};
    d = varargin{6};
    bc = varargin{7};
   otherwise
    error('Wrong number of input arguments');
  end
    
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
  
  prob = zeros(N,1);
  
  % first window
  tmpprob = melt_prob(temp,Gamma,seq(1:N0+d),T,bc);
  prob(1:N0) = tmpprob(1:N0);
  
  % second to next to last window
  for k=2:floor((N-d)/N0)-1
    tmpprob = melt_prob(temp,Gamma,seq((k-1)*N0-d+1:k*N0+d),T,bc);
    prob((k-1)*N0+1:k*N0)=tmpprob(d+1:N0+d);
  end
  
  % last window
  tmpprob = melt_prob(temp,Gamma,seq((floor((N-d)/N0)-1)*N0-d+1:end),T,bc);
  prob((floor((N-d)/N0)-1)*N0+1:end) = tmpprob(d+1:end);
