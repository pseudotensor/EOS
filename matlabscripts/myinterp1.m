function varargout = myinterp1(whichvar, varargin)


%extrap1(lsspecptotx, lsspecptoty, lsspecgrid','linear');

% do extrapolation
% as designed, derived quantities may obtain awkward values

%if 0


%if whichvar==8

% probably should include extra since often reach out of Kaz EOS
%if whichvar==5 || whichvar==6 || whichvar==7

% NORMAL: (when was not outputting Tk for each output type so didn't know if valid inversion)
%if whichvar==5 || whichvar==6


% DEBUG:
%if whichvar==5 || whichvar==6 || whichvar==1 || whichvar==2 || whichvar==3 || whichvar==4 || whichvar==7 || whichvar==8 || ...
%        whichvar == 20 || whichvar == 21 || whichvar==22 || whichvar==23 || whichvar==24 || whichvar==25 || whichvar==26 || whichvar==27 || whichvar==28

% DEBUG:
% if 1 || whichvar==5 || whichvar==6


% NORMAL: (when using Tk to indicate whether proper inversion)
% seems like extrapolation shouldn't be too bad as long as near boundaries of valid EOS
% extrapolation seems possibly best for iterative inversion too so can temporarily step out of validity -- but no guarantee that doesn't keep going away from valid EOS
if (whichvar~=29 && whichvar~=30 && whichvar~=31 && whichvar~=32)


  % DO EXTRAPOLATION
  varargout{:}=extrap1(varargin{:});


else

		%whichvar
    
% don't do extrapolation
% returns NaN's so can extrapolate PER quantity rather than deriving things
% from extrapolated quantities
    varargout{:}=interp1(varargin{:});
end
% try just setting NaN's to 0
%varargout{1}(isnan(varargout{1}(:))) = 0;

end
