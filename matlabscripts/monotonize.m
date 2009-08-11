%http://physics.gac.edu/~huber/matlab/mtlbasic.htm


function out = monotonize(inputlist)
% check if list is not sorted
% If not sorted, then enforce the function to be monotonically *increasing*
% Fill-in non-monotone regions with linear interpolation over smallest domain needed for monotonicity to be secured
%
  
  if(issorted(inputlist))
    out = inputlist;
    % then already sorted
    %		fprintf('got here0\n');
    return;
  end
  
  % we know it's not sorted
  originallist=inputlist;

  
  %lowest = min(inputlist);
  sizelisttemp=size(inputlist(:));
  sizelist=sizelisttemp(1);

  
  
  %default
  nextii=sizelist;
  iihigh=sizelist;
  highest=inputlist(iihigh);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % enforce that from lowest to highest that monotonic.  Fill-in with linear interpolation if not monotonically increasing
  for ii=sizelist-1:-1:1
   
    if(inputlist(ii+1)>inputlist(ii))
      % use as new highest so linear interpolation is most local
      iihigh=ii;
      highest=inputlist(ii);
    else
      % Then fix
      iilow=-1;
      % First find next monotonic point
      for jj=ii-1:-1:1
        if(inputlist(jj)<highest)
          lowest=inputlist(jj);
          iilow=jj;
          break; % stop then
        end
      end
      if iilow==-1
        % then no more lower points, so just extrapolate from last 2 points that were monotonic
        jj=ii;
          if highest>0.0 && inputlist(ii+2)>0.0 && inputlist(ii+1)>0.0
            inputlist(jj) = 10.^(log10(inputlist(ii+2)) + ( jj - (ii+2) )*(log10(inputlist(ii+1))-log10(inputlist(ii+2)))/((ii+1) - (ii+2)));
          else
            inputlist(jj) = inputlist(ii+2) + ( jj - (ii+2) )*(inputlist(ii+1)-inputlist(ii+2))/((ii+1) - (ii+2));
          end
      else
        % fix all between lowest and highest that are part of non-monotonic values
        jj=ii;
          % Second use current lowest and current highest to do linear interpolation
          if highest>0.0 && lowest>0.0 && highest>0.0
            inputlist(ii) = 10.^(log10(lowest) + (jj-iilow)*(log10(highest)-log10(lowest))/(iihigh-iilow));
          else
            inputlist(ii) = lowest + (jj-iilow)*(highest-lowest)/(iihigh-iilow);
          end
      end  %end if iilow==-1
    end % end if checking mono
    
  end
  
  
  %%%%%%%%%%%%%%%%%%%%%
  % Ensure really monotonic
  if(issorted(inputlist))
    % output modified inputlist
    out = inputlist;
    return;
  else
    fprintf('Failed to monotonize\n');
    
    fprintf('%g ',originallist);
    fprintf('\n');
    fprintf('%g ',inputlist);
    return;
  end
  
		

end






% test function
function out = testfunc()

inputlist = 1:.1:10;
sizelisttemp=size(inputlist(:));
sizelist=sizelisttemp(1);
inputlist(50)=1;
inputlist(51)=1.5;
inputlist(53)=1.2;
inputlist(sizelist)=1;

inputlist(50)
inputlist=monotonize(inputlist);
inputlist(50)


end



