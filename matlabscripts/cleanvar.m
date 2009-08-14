

% 1:  log10 PofU (rho0,U)
% 2:  log10 HofU (rho0,U)
% 3:  log10 UofP (rho0,P)
% 4:  log10 PofCHI (rho0,CHI)
% 5:  log10 PofS (rho0,S)               % S is specific entropy
% 6:  log10 UofS (rho0,S)               % S is specific entropy
% 7:  log10 extraofU (rho0,U)
% 8:  log10 SofU (rho0,U)               % S is entropy density
% 9:  lin   dPofUdrho0cS (rho0,U) | S   % S is specific entropy
% 20: lin   dPofUdrho0out (rho0,U)
% 21: lin   dPofUduout (rho0,U)
% 22: lin   cs2/csq (rho0,U)
% 23: lin   cs2cgs (rho0,U)
% 24: lin   dSofUdrho0out (rho0,U)  % S is entropy density
% 25: lin   dSofUduout (rho0,U)     % S is entropy density
% 26: lin   dPofCHIdrho0out (rho0,CHI)
% 27: lin   dPofCHIdchiout (rho0,CHI)
% 28: log10 extraofUout (rho0,U)
% 29: log10 Tk(rho0,U)
% 30: log10 Tk(rho0,P)
% 31: log10 Tk(rho0,CHI)
% 200: log10(UofUdiff)
% 201: log10(PofPdiff)
% 202: log10(CHIofCHIdiff)
% 203: log10(SofSdiff)
% 204: log10(SSofSSdiff)
% 65 : UofSdiff


% input is full 4-D quantities
function funout = cleanvar(whichvar, funin, X, Y)

  % check if nan and fix if so
  myisnan=isnan(funin);
  % check if log10() is not finite and fix if so and if should be able to have finite log10() value
  myisnotfinite=~isfinite(log10(funin));

  if whichvar==1
    funin(myisnan) = idealPofU(X(myisnan),Y(myisnan));
    funin(myisnotfinite) = idealPofU(X(myisnotfinite),Y(myisnotfinite));
  end
  if whichvar==2
    funin(myisnan) = idealHofU(X(myisnan),Y(myisnan));
    funin(myisnotfinite) = idealHofU(X(myisnotfinite),Y(myisnotfinite));
  end
  if whichvar==3
    funin(myisnan) = idealUofP(X(myisnan),Y(myisnan));
    funin(myisnotfinite) = idealUofP(X(myisnotfinite),Y(myisnotfinite));
  end
  if whichvar==65
    funin(myisnan) = idealUofSdiff(X(myisnan),Y(myisnan));
    funin(myisnotfinite) = idealUofSdiff(X(myisnotfinite),Y(myisnotfinite));
  end
  if whichvar==4
    funin(myisnan) = idealPofCHI(X(myisnan),Y(myisnan));
    funin(myisnotfinite) = idealPofCHI(X(myisnotfinite),Y(myisnotfinite));
  end
  if whichvar==400
    funin(myisnan) = idealSSofCHI(X(myisnan),Y(myisnan));
    funin(myisnotfinite) = idealSSofCHI(X(myisnotfinite),Y(myisnotfinite));
  end
  if whichvar==5
    funin(myisnan) = idealPofS(X(myisnan),Y(myisnan));
    funin(myisnotfinite) = idealPofS(X(myisnotfinite),Y(myisnotfinite));
  end
  if whichvar==6
    funin(myisnan) = idealUofS(X(myisnan),Y(myisnan));
    funin(myisnotfinite) = idealUofS(X(myisnotfinite),Y(myisnotfinite));
  end
  if whichvar==7
    funin(myisnan) = idealextraofU(X(myisnan),Y(myisnan));
    funin(myisnotfinite) = idealextraofU(X(myisnotfinite),Y(myisnotfinite));
  end
  if whichvar==8
    funin(myisnan) = idealSdenofU(X(myisnan),Y(myisnan));
    funin(myisnotfinite) = idealSdenofU(X(myisnotfinite),Y(myisnotfinite));
  end
  if whichvar==9
    funin(myisnan) = idealdPofUdrho0cS(X(myisnan),Y(myisnan));
    % should be positive because eventually comes into sound speed with all other positive terms
    funin(myisnotfinite) = idealdPofUdrho0cS(X(myisnotfinite),Y(myisnotfinite));
  end

  if whichvar==29
    funin(myisnan) = idealtkofU(X(myisnan),Y(myisnan));
    funin(myisnotfinite) = idealtkofU(X(myisnotfinite),Y(myisnotfinite));
  end

  if whichvar==30
    funin(myisnan) = idealtkofP(X(myisnan),Y(myisnan));
    funin(myisnotfinite) = idealtkofP(X(myisnotfinite),Y(myisnotfinite));
  end

  if whichvar==31
    funin(myisnan) = idealtkofCHI(X(myisnan),Y(myisnan));
    funin(myisnotfinite) = idealtkofCHI(X(myisnotfinite),Y(myisnotfinite));
  end

  if whichvar==40
    funin(myisnan) = idealcs2ofU(X(myisnan),Y(myisnan));
    funin(myisnotfinite) = idealcs2ofU(X(myisnotfinite),Y(myisnotfinite));
  end

  if whichvar==200
    funin(myisnan) = idealUofU(X(myisnan),Y(myisnan));
    funin(myisnotfinite) = idealUofU(X(myisnotfinite),Y(myisnotfinite));
  end

  if whichvar==201
    funin(myisnan) = idealUofU(X(myisnan),Y(myisnan));
    funin(myisnotfinite) = idealUofU(X(myisnotfinite),Y(myisnotfinite));
  end

  if whichvar==202
    funin(myisnan) = idealUofU(X(myisnan),Y(myisnan));
    funin(myisnotfinite) = idealUofU(X(myisnotfinite),Y(myisnotfinite));
  end

  if whichvar==203
    funin(myisnan) = idealSofS(X(myisnan),Y(myisnan));
    funin(myisnotfinite) = idealSofS(X(myisnotfinite),Y(myisnotfinite));
  end

  if whichvar==204
    funin(myisnan) = idealSSofSS(X(myisnan),Y(myisnan));
    funin(myisnotfinite) = idealSSofSS(X(myisnotfinite),Y(myisnotfinite));
  end

  % assume rest are derivative quantities (i.e. derived from the above)


  % assume we fixed funin and set to funout
  funout = funin;

end




% whether c^2 is with density or not for the below functions
% cleanvar(1,PofU): rhob
% cleanvar(2,PofU): rhobcsq
% cleanvar(3,PofU): rhob
% cleanvar(4,PofU): rhob
% cleanvar(5,PofU): rhobcsq
% cleanvar(6,PofU): rhobcsq
% cleanvar(7,PofU): rhob
% cleanvar(8,PofU): rhobcsq
% cleanvar(29,tkofU): rhob
% cleanvar(30,tkofP): rhob
% cleanvar(31,tkofCHI): rhob



function out = idealPofU(rho0, U)

GAMMA=(4.0/3.0);

%out = (GAMMA-1.0).*U;
%		out = OUTBOUNDSVALUE;

% Mignone EOS
c = 2.99792458E10;
rhobcsq=rho0.*(c.*c);
out = U.*(2.0.*rhobcsq+U)./(3.0.*(rhobcsq+U));

end

% ratios of energies assumed to be dimensionless
function out = idealHofU(rhobcsq, U)

GAMMA=(4.0/3.0);
%out = 1.0 + GAMMA.*U./(rhobcsq);
%		out = OUTBOUNDSVALUE;

% Mignone EOS
c = 2.99792458E10;
rho0=rhobcsq./(c.*c);
PofU = idealPofU(rho0, U);

% specific enthalpy
out = 1.0 + (U+PofU)./rhobcsq;

end


function out = idealUofP(rho0, P)

GAMMA=(4.0/3.0);

%out = P./(GAMMA-1.0);
%		out = OUTBOUNDSVALUE;

% Mignone EOS
c = 2.99792458E10;
rhobcsq=rho0.*(c.*c);
out = 1.5.*(P + 3.0.*P.*P./(2.0.*rhobcsq+sqrt(9.0.*P.*P+4.0.*rhobcsq.*rhobcsq)));

end


function out = idealUofSdiff(rho0, Sdiff)

OUTBOUNDSVALUE=1E-20;
out = OUTBOUNDSVALUE;

end

  
function out = idealPofCHI(rho0, CHI)

GAMMA=(4.0/3.0);

%out = ((GAMMA-1.0)./GAMMA).*CHI;
%		out = OUTBOUNDSVALUE;

% Mignone EOS
c = 2.99792458E10;
rhobcsq=rho0.*(c.*c);
wmrho0=CHI;
Q=wmrho0./rhobcsq;
delta=9.0./25.0.*wmrho0.*(2.0+Q);
delta2=delta./rhobcsq;

out=(5.0./8.0).*(wmrho0 - delta./(1.0+sqrt(1.0+delta2)));




end



function out = idealSSofCHI(rho0, CHI)

OUTBOUNDSVALUE=1E-20;
out = OUTBOUNDSVALUE;
  
end



% Sden 
function out = idealPofSden(rhobcsq, Sden)

GAMMA=(4.0/3.0);
C = 0.00039324313707644323;
%out = (rho0.^GAMMA) .* exp(Sden./(kb.*rho0./mN));

out = (GAMMA-1.0).*((Sden./C).^(GAMMA));
%		out = OUTBOUNDSVALUE;

% not setup for Mignone EOS yet


end



% Sden
function out = idealUofSden(rhobcsq, Sden)

GAMMA=(4.0/3.0);

%out = (Sden./C).^(4.0/3.0);


out		= idealPofSden(rhobcsq,Sden)./(GAMMA-1.0);



%		out = OUTBOUNDSVALUE;
% not setup for Mignone EOS yet


end


% here S is specific entropy
function out = idealPofS(rhobcsq, S)

% consistent with how defined specific entropy above (sspec)
Sden = S.*rhobcsq;

out = idealPofSden(rhobcsq, Sden);


%		out = OUTBOUNDSVALUE;

% not setup for Mignone EOS yet

end



% here S is specific entropy
function out = idealUofS(rhobcsq, Ss)

GAMMA=(4.0/3.0);

% idealPofS takes in specific entropy
out = idealPofS(rhobcsq,Ss)./(GAMMA-1.0);

%		out = OUTBOUNDSVALUE;
% not setup for Mignone EOS yet


end





function out = idealextraofU(rho0, U)

GAMMA=(4.0/3.0);

% put in something small so can use log-log interpolation
%out = 1E-20;
OUTBOUNDSVALUE=1E-20;
		out = OUTBOUNDSVALUE;

end


function out = idealtkofU(rho0, U)

% put in something small so know that out of Kaz-EOS for function(rho0,u)
%out = 1E-20;
OUTBOUNDSVALUE=1E-20;
		out = OUTBOUNDSVALUE;

end

function out = idealtkofP(rho0, P)

% put in something small so know that out of Kaz-EOS for function(rho0,u)
%out = 1E-20;
OUTBOUNDSVALUE=1E-20;
		out = OUTBOUNDSVALUE;

end

function out = idealtkofCHI(rho0, CHI)

% put in something small so know that out of Kaz-EOS for function(rho0,u)
%out = 1E-20;
OUTBOUNDSVALUE=1E-20;
		out = OUTBOUNDSVALUE;

end


function out = idealtkofS(rho0, S)

% put in something small so know that out of Kaz-EOS for function(rho0,u)
%out = 1E-20;
OUTBOUNDSVALUE=1E-20;
		out = OUTBOUNDSVALUE;

end


function out = idealUofU(rho0, U)

  out = U;

end

function out = idealSofS(rho0, S)

  out = S;

end

function out = idealSSofSS(rho0, SS)

  out = SS;

end


% S = entropy density
function out = idealSdenofU(rhobcsq, U)

GAMMA=(4.0/3.0);
C = 0.00039324313707644323;

% from Kaz's formula for photon entropy and see kaz_entropy_conversion.nb
% for how converted to cgs
out = C.*U.^(3.0/4.0);
%		out = OUTBOUNDSVALUE;
% not setup for Mignone EOS yet

end



% S = specific entropy
function out = idealSofU(rhobcsq, U)

GAMMA=(4.0/3.0);
P = (GAMMA-1.0)*U;
Sden = idealSdenofU(rhobcsq, U);

% true specific entropy?
%out = Sden./(rhobcsq + U + P);
% makes no sense if no baryons
out = Sden./rhobcsq;
%		out = OUTBOUNDSVALUE;
% not setup for Mignone EOS yet

end



% ratios from derivatives assumed to be dimensionless
% holding constant: S = specific entropy
function out = idealdPofUdrho0cS(rhobcsq,U)

GAMMA=(4.0/3.0);
hspec = idealHofU(rhobcsq,U);
P = (GAMMA-1.0).*U;

% dP/drho|S = h c_s^2 = h GAMMA P/(\rho_0+u+P) = h GAMMA(GAMMA-1)u/(\rho_0+u+P)
out = hspec.*GAMMA.*P./(rhobcsq+U+P);
%		out = OUTBOUNDSVALUE;
% not setup for Mignone EOS yet


end



function out = idealcs2ofU(rhobcsq, U)

GAMMA=(4.0/3.0);
P = (GAMMA-1.0).*U;
out = GAMMA.*P./(rhobcsq+U+P);
% not setup for Mignone EOS yet

end







function out = pwfutot(rho0,u)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PWF99 u_tot
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
% baryon mass (m_n+m_p)/2
mb = 1.67377E-24;
  
R = kb/mb;
K=1.24e15;
arad=7.56641E-15;
out=3/2.*rhob*R.*tk+11/4*arad.*tk.^4+3*K.*rhob.^(4/3)/2^(1/3);

%for ii=1:48 fprintf('%21.15g %21.15g %21.15g %21.15g\n',log10(utot(28,ii,1,1)/(c*c)),log10(utotnew(28,ii,1,1)/(c*c)),log10(out(28,ii,1,1)/(c*c)),log10(tk(28,ii,1,1))); end


end

