function c = overlap(prob,cod)
% OVERLAP - Compute the overlap between opening probability and sequence property
%   c = overlap(prob,cod) computes the overlap between the opening
%   probabilities "prob" and a sequence property "cod".
%   Example: property is coding/non-coding; cod is a vector of 0/1's,
%   0=non-coding, 1=coding.
%   The score function is the standard normalized vector product, i.e.,
%   lies between 0 and 1 and is 1 if and only if prob=cod.
  
  if length(prob) == length(cod)
    c = (1-prob)'*cod/(norm(1-prob)*norm(cod));
  else
    error('Input arguments should be vectors of the same length.')
  end
