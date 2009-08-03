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
  
  originallist=inputlist;

  
  % we know it's not sorted
  lowest = inputlist(1);
  iilow=1;
  %lowest = min(inputlist);
  iihigh=iilow; % default
  sizelisttemp=size(inputlist(:));
  sizelist=sizelisttemp(1);

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % find overall lowest value
  lowest = inputlist(1);
  iilow=1;
  for ii=2:sizelist
    if(inputlist(ii)<inputlist(ii-1))
      lowest = inputlist(ii);
      iilow=ii;
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % find 2nd lowest value beyond first non-monotonic region
  lowest2nd = inputlist(iilow);
  iilow2nd=1;
  for ii=iilow+1:sizelist
    if(inputlist(ii)>inputlist(ii-1))
      lowest2nd = inputlist(ii);
      iilow2nd=ii;
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % find highest value beyond first non-monotonic region
  highest = inputlist(iilow);
  iihigh=1;
  for ii=iilow+1:sizelist
    if(inputlist(ii)>inputlist(ii-1))
      highest = inputlist(ii);
      iihigh=ii;
    end
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % check that can become monotonic (requires at least 2 points be different)
  if(highest==lowest)
    fprintf('ii=%d lowest=%g highest=%g\n',ii,lowest,highest);
    return;
  end
  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % process data so becomes monotonic

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % enforce up to lowest to be lower using 
  if lowest>0.0
    % then can do log extrapolation (avoids negatives in sspec and stot)
    for ii=1:iilow-1
      inputlist(ii) = 10.^(log10(lowest) + (ii-iilow)*(log10(lowest2nd)-log10(lowest))/(iilow2nd-iilow));
    end
  else
    for ii=1:iilow-1
      inputlist(ii) = lowest + (ii-iilow)*(lowest2nd-lowest)/(iilow2nd-iilow);
    end    
  end

  % DEBUG:
  %fprintf('%g %g %g %g %g %g : %g\n',lowest,iilow,lowest2nd,lowest,iilow2nd,iilow,(ii-iilow)*(lowest2nd-lowest)/(iilow2nd-iilow));
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % enforce that from lowest to highest that monotonic.  Fill-in with linear interpolation if not monotonically increasing
  for ii=iilow+1:sizelist
    
    if(inputlist(ii)>inputlist(ii-1))
      % use as new lowest so linear interpolation is most local
      iilow=ii;
      lowest=inputlist(ii);
    else
      % Then fix
      iihigh=-1;
      % First find next monotinic point
      for jj=ii+1:sizelist
        if(inputlist(jj)>lowest)
          highest=inputlist(jj);
          iihigh=jj;
          break; % stop then
        end
      end
      if iihigh==-1
        % then no more higher points, so just extrapolate from last 2 points that were monotonic
        for jj=ii:sizelist
          if lowest>0.0
            inputlist(jj) = 10.^(log10(inputlist(ii-2)) + ( jj - (ii-2) )*(log10(inputlist(ii-1))-log10(inputlist(ii-2)))/((ii-1) - (ii-2)));
          else
            inputlist(jj) = inputlist(ii-2) + ( jj - (ii-2) )*(inputlist(ii-1)-inputlist(ii-2))/((ii-1) - (ii-2));
          end
        end
      else
        % fix all between lowest and highest that are part of non-monotonic values
        for jj=ii:iihigh-1
          % Second use current lowest and current highest to do linear interpolation
          if lowest>0.0
            inputlist(ii) = 10.^(log10(lowest) + (ii-iilow)*(log10(highest)-log10(lowest))/(iihigh-iilow));
          else
            inputlist(ii) = lowest + (ii-iilow)*(highest-lowest)/(iihigh-iilow);
          end
        end
      end
    end
    
    
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




function out = monotonizeold(inputlist)
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
 lowest = inputlist(1);
 iilow=1;
 iihigh=iilow; % default
 sizelisttemp=size(inputlist(:));
 sizelist=sizelisttemp(1);
 
 
 % loop over data array
 for ii=2:sizelist-1

  
  if(inputlist(ii)<=inputlist(ii-1))  % not monotonically increasing if less or equal
    % then keep going until inputlist(ii) is finally > lowest
    for jj=ii:sizelist
		if(inputlist(jj)>lowest)  % only strictly greater is allowed (could put a machine threshold on this)
		   highest=inputlist(jj);
		   iihigh=jj;
		   %fprintf('got here0.3 ii=%d iilow=%d iihigh=%d\n',ii,iilow,iihigh);
		   break; % stop then
		end
		%fprintf('got here0.5 ii=%d iilow=%d iihigh=%d\n',ii,iilow,iihigh);
    end
    %   between the range of ii=iilow to iihigh the values are non-monotonic.  So replace them with a linear interpolation
    %fprintf('got here0.7 ii=%d iilow=%d iihigh=%d\n',ii,iilow,iihigh);
		if(highest==lowest)
		  fprintf('ii=%d lowest=%g highest=%g\n',ii,lowest,highest);
		end
    for jj=iilow:iihigh
		inputlist(jj) = lowest + (jj-iilow)*(highest-lowest)/(iihigh-iilow);
    end

    % jump ahead: set correct next ii as in the outer loop
    ii=iihigh;
    %fprintf('got here1 ii=%d iilow=%d iihigh=%d\n',ii,iilow,iihigh);
  else
    % replace lowest
    lowest = inputlist(ii);
    iilow=ii;
    iihigh=iilow;
    %fprintf('got here2 ii=%d iilow=%d iihigh=%d\n',ii,iilow,iihigh);
  end
		
 end

	%	sizelist

 % check last value and extrapolate instead of interpolate
 % presumes rest of list is monotonicically increasing now
 if(inputlist(sizelist)<=inputlist(sizelist-1))
   inputlist(sizelist) = inputlist(sizelist-2) + ((sizelist) - (sizelist-2))*(inputlist(sizelist-1)-inputlist(sizelist-2))/((sizelist-1) - (sizelist-2));
 end

if(inputlist(sizelist-1)==inputlist(sizelist-2))
  fprintf('finalpoint lowest=%g highest=%g\n',inputlist(sizelist-1),inputlist(sizelist-2));
end

		

 % output modified inputlist
 out = inputlist;

 % BEGIN DEBUG
 %inputlist
 %outputlist
 % END DEBUG
		

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
