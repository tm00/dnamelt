function [ind, A] = extract_coding_gb(GenbankFile)
% EXTRACT_CODING_GB - extract coding regions out of Genbank file
%   [ind,A]=extract_coding_gb(GenbankFile) extracts the coding regions
%   out of the file "Genbankfile" and puts them in the array A. A(:,1)
%   contains the start positions, A(:,2) the end positions; the array
%   "ind" has the same length as the sequence with 1 in coding regions, 0
%   otherwise. This function uses "genbankread" from the  bioinformatics 
%   toolbox  
  
  A = [];
  L = size(A);
  % read the file
  data = genbankread(GenbankFile);
      
  % the field data.CDS.location unfortunately only returns the first CDS
  % line; use data.CDS.text instead
  for k = 1:length(data.CDS)
    % write data.CDS(k).text on one line
    sz = size(data.CDS(k).text);
    str = reshape(data.CDS(k).text', [1 sz(1)*sz(2)]);
    % the format is "CDS  -coding sites- /locustag="; hence read to first
    % "/"; this returns a 1x1 cell containing the char array
    str = textscan(str,'%s',1,'delimiter','/');
    str = char(str{1});
    % get rid of all whitespace
    str = strrep(str,' ','');
    % get rid of junk in front of coding positions
    str2 = textscan(str,'%[^123456789]');
    str2 = str2{1};
    str = strrep(str,str2,'');
    % get rid of junk brackets at the end of coding positions
    str = strrep(str,')','');
    % get rid of > signs
    str = strrep(str,'>','');
    % now str has the form A(1,1)..A(1,2),A(2,1)..A(2,2),etc
    % first split at the ','
    str = textscan(char(str),'%s','delimiter',',');
    str = str{1};
    for l = 1:length(str)
      str2 = textscan(char(str{l}), '%d %d', 'delimiter','..');
      A(L(1)+l,1) = str2{1};
      A(L(1)+l,2) = str2{2};
    end
    % current row length of A;
    L = size(A);
  end
  ind = zeros(length(data.Sequence),1);
  for l = 1: L(1)
    ind(A(l,1):A(l,2)) = 1;
  end
