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
