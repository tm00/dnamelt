function dna_melt_preferences()
% DNA_MELT_PREFERENCES - Set preferences for DNA_melt programs
%   set model and integration parameters at will, then run this function.   

% MODEL PARAMETERS
kB = 8.617385*10^(-5); % Boltzmann constant [eV/K]
Ds = 0.18; % Morse potential depth C-G [eV]
as = 6.9; % Morse potential inv. width C-G [1/Angstrom]
Dw = 0.12; % Morse potential depth A-T [eV]
aw = 4.2; % Morse potential inv. width A-T [1/Angstrom]
D =[Dw Ds Ds Dw]';
a=[aw as as aw]';
r0 = 10.0; % closed bp radius
sigmaslide = [0.28 0.71 0.82 0.48;
              1.23 1.17 1.02 0.82;
              0.69 0.86 1.17 0.71;
              1.09 0.69 1.23 0.28]; % slide variance
K = 0.1*sigmaslide.^(-1); % anharmonic stacking strength [eV]
%K = (sum(sum(K))/16)*ones(4);
%K = 0.65*ones(4);

alpha = repmat(0.5, [4 4]); % anharmonic stacking inv. width [1/Angstrom]
sigmatwist = [3.3 3.8 4.8 2.8;
              9.5 3.7 5.3 4.8;
              3.8 4.0 3.7 3.8;
              6.7 3.8 9.5 3.3]; % twist variance
E =0.4*sigmatwist.^(-1); % twist energy strength [eV]
%E = (sum(sum(E))/16)*ones(4);
%E = 0.04*ones(4);

h = 3.4; % vertical distance between bps [Angstrom]
theta0 = (pi/180).*[35.9 32.9 34.8 32.4;
                    37.4 31.9 35.1 34.8;
                    37.8 37.4 31.9 32.9;
                    30.6 37.8 37.4 35.9]; % average twist angles
%theta0 = (sum(sum(theta0))/16)*ones(4);
%theta0 = 0.60707*ones(4);


l0 = sqrt(h^2*ones(4) + 4*r0^2*(sin(0.5*theta0)).^2); 
% ground state length between nucleotides on same strand [Angstrom]

% NUMERICAL INTEGRATION PARAMETERS
ML = 36; % size of Gauss-Legendre grid
% if you change ML, replace legz by "legz=legpolzeros(ML)"
legz = [-0.997830462484086;
        -0.988586478902212;
        -0.972027691049698;
        -0.948272984399508;
        -0.917497774515659;
        -0.879929800890397;
        -0.835847166992475;
        -0.785576230132207;
        -0.729489171593557;
        -0.668001236585521;
        -0.601567658135981;
        -0.530680285926245;
        -0.45586394443342;
        -0.377672547119689;
        -0.296684995344028;
        -0.213500892316866;
        -0.128736103809385;
        -0.0430181984737086;
        0.0430181984737086;
        0.128736103809385;
        0.213500892316866;
        0.296684995344028;
        0.377672547119689;
        0.45586394443342;
        0.530680285926245;
        0.601567658135981;
        0.668001236585521;
        0.729489171593557;
        0.785576230132207;
        0.835847166992475;
        0.879929800890397;
        0.917497774515659;
        0.948272984399508;
        0.972027691049698;
        0.988586478902212;
        0.997830462484086]; % roots of Legendre polynomial P_ML


% $$$ ML = 70;
% $$$ legz = [ -0.999418285973576
% $$$          -0.99693625196168
% $$$          -0.99247605521169
% $$$         -0.986045558070399
% $$$         -0.977657405957592
% $$$         -0.967328223664986
% $$$         -0.955078509114293
% $$$         -0.940932579003815
% $$$         -0.924918516897934
% $$$         -0.907068116260923
% $$$         -0.887416816863348
% $$$         -0.866003634213859
% $$$          -0.84287108199898
% $$$         -0.818065087625441
% $$$         -0.791634901007893
% $$$           -0.7636329967719
% $$$         -0.734114970060943
% $$$         -0.703139426151529
% $$$         -0.670767864094077
% $$$         -0.637064554609778
% $$$         -0.602096412485356
% $$$         -0.565932863718808
% $$$         -0.528645707679711
% $$$         -0.490308974557637
% $$$         -0.450998778381648
% $$$         -0.410793165902631
% $$$         -0.369771961638462
% $$$         -0.328016609389643
% $$$         -0.285610010540038
% $$$         -0.242636359463741
% $$$         -0.199180976364858
% $$$          -0.15533013788207
% $$$         -0.111170905794299
% $$$        -0.0667909541675513
% $$$        -0.0222783952861403
% $$$         0.0222783952861403
% $$$         0.0667909541675513
% $$$          0.111170905794299
% $$$           0.15533013788207
% $$$          0.199180976364858
% $$$          0.242636359463741
% $$$          0.285610010540038
% $$$          0.328016609389643
% $$$          0.369771961638462
% $$$          0.410793165902631
% $$$          0.450998778381648
% $$$          0.490308974557637
% $$$          0.528645707679711
% $$$          0.565932863718808
% $$$          0.602096412485356
% $$$          0.637064554609778
% $$$          0.670767864094077
% $$$          0.703139426151529
% $$$          0.734114970060943
% $$$            0.7636329967719
% $$$          0.791634901007893
% $$$          0.818065087625441
% $$$           0.84287108199898
% $$$          0.866003634213859
% $$$          0.887416816863348
% $$$          0.907068116260923
% $$$          0.924918516897934
% $$$          0.940932579003815
% $$$          0.955078509114293
% $$$          0.967328223664986
% $$$          0.977657405957592
% $$$          0.986045558070399
% $$$           0.99247605521169
% $$$           0.99693625196168
% $$$          0.999418285973576];

legw = 2.*(1-legz.^2)./((ML+1).*legz.*legpol(legz,ML) ...
                        - (ML+1).*legpol(legz,ML+1)).^2; % weights
inta = 9.3;
intb = 40;  % integration interval [inta,intb]
xi = (intb-inta)/2.*legz + (intb+inta)/2;   % legendre nodes in [inta,intb]

MC = 24; % size of Gauss-Chebyshev grid
mc=1:MC;
chebz = cos(pi.*(mc-0.5)./MC)'; % roots of Chebyshev polynomial T_MC
chebw = pi/MC; % weights 


% DO NOT CHANGE BELOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
setpref('DNA_melt',{'kB','D','a','r0','K','alpha','E','h','theta0','l0',...
                    'ML','legw','xi','MC','chebz','chebw'},...
                   {kB,D,a,r0,K,alpha,E,h,theta0,l0,ML,legw,xi,MC,chebz,chebw});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function z = legpolzeros(n)
% LEGPOLZEROS - find the zeros of the Legendre polynomial of order n
%   
  global ord
  ZR = zeros(ceil(n/2)+1,n);
  leg = inline('legpol(x)');
  ZR(2,1)=1;
  for l = 2:n
    ZR(ceil(l/2)+1,l)=1;
    ord = l;
    for k = 1+mod(l,2):ceil(l/2)
      ZR(k,l) = fzero(leg, [ZR(k-mod(l,2),l-1) ZR(k+1-mod(l,2),l-1)]);
    end
  end
   
  if mod(n,2) == 0
    z = sort([-ZR(1:ceil(n/2),n).' ZR(1:ceil(n/2),n).'].');
  else
    z = sort([-ZR(2:ceil(n/2),n).' 0  ZR(2:ceil(n/2),n).'].');
  end
  
function l = legpol(x,n)
% LEGPOL - legendre polynomial of order n evaluated at x
%   x is either real value in [-1,1] or a (column) vector of such values
    
  global ord
  
  if nargin == 1
    tmp = legendre(ord,x);
    s = size(x);
    if s(1,1) == 1
      l = tmp(1);
    else
      l = tmp(1,:).';
    end
  else
    tmp = legendre(n,x);
    s = size(x);
    if s(1,1) == 1
      l = tmp(1);
    else
      l = tmp(1,:).';
    end
  end
  