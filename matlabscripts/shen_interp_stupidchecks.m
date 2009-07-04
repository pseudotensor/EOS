  % Below needed for xcheck and yecheck
  if(1)
    % However, in general want to correct largest fractions that will
    % absorb smallest relative error
    ifxh=((ishenxh_row>ishenxneut_row) & (ishenxh_row>ishenxprot_row) & (ishenxh_row>ishenxalfa_row) );
    ifxneut=((ishenxneut_row>ishenxh_row) & (ishenxneut_row>ishenxprot_row) & (ishenxneut_row>ishenxalfa_row) );
    ifxprot=((ishenxprot_row>ishenxh_row) & (ishenxprot_row>ishenxneut_row) & (ishenxprot_row>ishenxalfa_row) );
    ifxalfa=((ishenxalfa_row>ishenxh_row) & (ishenxalfa_row>ishenxneut_row) & (ishenxalfa_row>ishenxprot_row) );
    ifnonelargest=(ifxh==0)&(ifxneut==0)&(ifxprot==0)&(ifxalfa==0);
  end

  
  
  
  if(1)
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % xcheck
    %
    % consistency of fractions relation #2 in Shen EOS guide.ps
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    xcheckbefore = 1.0 - (ishenxneut_row + ishenxprot_row + ishenxalfa_row + ishenxh_row);
    % likely that error is mostly in xh since reaches large values
    %
    % now correct
    ishenxh_row    = ishenxh_row    .*(1-ifxh)    + ifxh   .*(1.0 - (ishenxneut_row + ishenxprot_row + ishenxalfa_row));
    ishenxneut_row = ishenxneut_row .*(1-ifxneut) + ifxneut.*(1.0 - (ishenxh_row    + ishenxprot_row + ishenxalfa_row));
    ishenxprot_row = ishenxprot_row .*(1-ifxprot) + ifxprot.*(1.0 - (ishenxneut_row + ishenxh_row + ishenxalfa_row));
    ishenxalfa_row = ishenxalfa_row .*(1-ifxalfa) + ifxalfa.*(1.0 - (ishenxneut_row + ishenxprot_row + ishenxh_row));

    % if no largest (i.e. all equal), then just correct xh
    ishenxh_row    = ishenxh_row    .*(1-ifnonelargest)    + ifnonelargest   .*(1.0 - (ishenxneut_row + ishenxprot_row + ishenxalfa_row));

    
    %

    %%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Enforce minimums after correction
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%
    ishenxneut_row(ishenxneut_row<xmin)=xmin;
    ishenxprot_row(ishenxprot_row<xmin)=xmin;
    ishenxalfa_row(ishenxalfa_row<xmin)=xmin;
    ishenxh_row(ishenxh_row<xmin)=xmin;

    
    xcheckafter = 1.0 - (ishenxneut_row + ishenxprot_row + ishenxalfa_row + ishenxh_row);
    numarray = xcheckafter*0.0 + 1.0;
    totalxcheckafter=sum(sum(sum(abs(xcheckafter))))./sum(sum(sum(abs(numarray))));

    mynotnan = isfinite(xcheckafter);
    notnanxcheckafter(mynotnan)=xcheckafter(mynotnan);
    problems(notnanxcheckafter>0.1)=1;
    sumproblems=sum(sum(sum(problems)));
    if(sumproblems>=1)
      fprintf('Problems with xcheck: %d\n',sumproblems);
    end
    
    if(totalxcheckafter>0.1)
      fprintf('Bad totalxcheckafter: %g\n',totalxcheckafter);
    end
    
    xcheckaftercopy=xcheckafter;
    xcheckaftercopy(isnan(xcheckaftercopy))=0;
    problems=(xcheckaftercopy>0.1);
    sumproblems=sum(sum(sum(problems)));
    if(sumproblems>=1)
      fprintf('ver2: Problems with xcheck: %d\n',sumproblems);
    end
  end
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Y_e check
  %
  % Now enforce Y_e condition
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%
  %
  if(1)

    
    % Code below is from jon_lsbox.f
    ytot1temp   = ishenxneut_row  + ishenxprot_row + 0.25d0.*ishenxalfa_row;
    trustheavy=(ishenxh_row>xhtrust) & (ishenaheav_row>aheavtrust);
    % the 1E-30 below just avoids division by 0
    ytot1 = ytot1temp + (ishenxh_row./(ishenaheav_row+1E-30)).*(trustheavy);
    zbarxxtemp  = ishenxprot_row + 0.5d0*ishenxalfa_row;
    x = ishenzheav_row./ishenaheav_row;
    zbarxx = zbarxxtemp + x.*ishenxh_row.*(trustheavy);

    abarnum    = 1.0d0./ytot1;
    abar=abarnum;
    zbar    = zbarxx .* abar;
    % yelocal should be equal to ishenyp_row, but don't check, just fix below
    yelocalbefore = zbar./abarnum;
    
    % again, least trustable quantities are related to interpolations of
    % heavy nucleus quantities.  Already fixed xh, now choose to put error
    % equally in aheav and zheav since similar type of quantities
    % for now try only putting into aheav
    % see mathematica Sheneos.nb
    %
    % correct either a or z, whatever is largest and put minimum on answer
    % below returns back good (interpolated) yp
    
    % only makes sense to change aheav or zheav if xh dominates (and note can change aheav and zheav if xh was changed by xcheck)
    % otherwise fix largest species (xprot,xneut,xalfa) not already changed for xcheck
   
    xhdominates=(ifxh | ifnonelargest);
    ifaheavlarger=(ishenaheav_row>ishenzheav_row);
    changeaheav = (xhdominates & ifaheavlarger==1 & trustheavy==1);
    changezheav = (xhdominates & ifaheavlarger==0 & trustheavy==1);
    changeazheav = (changeaheav | changezheav);
    ishenaheav_row = ishenaheav_row.*(1-changeaheav) + changeaheav.*((ishenxh_row.*ishenzheav_row)./(ishenyp_row-0.5*ishenxalfa_row-ishenxprot_row));
    ishenzheav_row = ishenzheav_row.*(1-changezheav) + changezheav.*(ishenaheav_row.*(ishenyp_row-0.5*ishenxalfa_row-ishenxprot_row)./ishenxh_row);

    % xhange xprot if not chosen by xcheck, xh didn't dominate (so could change aheav or zheav), but larger than one of the other species
    % Figure out what x is second (i.e. between charged species)
    changexprot=(ifxprot==0 & changeazheav==0 & (ifxalfa==1 | ishenxprot_row>ishenxalfa_row ) );
    changexalfa=(ifxalfa==0 & changeazheav==0 & (ifxprot==1 | ishenxalfa_row>ishenxprot_row ) );
    
    ishenxprot_row = ishenxprot_row.*(1-changexprot) ...
        + changexprot.*(ishenyp_row - 0.5.*ishenxalfa_row - ishenxh_row.*ishenzheav_row./ishenaheav_row.*(trustheavy));

    ishenxalfa_row = ishenxalfa_row.*(1-changexalfa) ...
        + changexalfa.*(2.0.*(ishenyp_row - ishenxprot_row - ishenxh_row.*ishenzheav_row./ishenaheav_row.*(trustheavy)));
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Enforce minimums after corrections
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%
    ishenaheav_row(ishenaheav_row<azmin)=azmin;
    ishenzheav_row(ishenzheav_row<azmin)=azmin;
    % also may have changed these:
    ishenxneut_row(ishenxneut_row<xmin)=xmin;
    ishenxprot_row(ishenxprot_row<xmin)=xmin;
    ishenxalfa_row(ishenxalfa_row<xmin)=xmin;
    ishenxh_row(ishenxh_row<xmin)=xmin;

    
    
    
    
    % below would return same answer
    %ishenaheav_row = (ishenxh_row.*ishenzheav_row)./(yelocalbefore-0.5*ishenxalfa_row-ishenxprot_row);

    % Code below is from jon_lsbox.f
    ytot1temp   = ishenxneut_row  + ishenxprot_row + 0.25d0.*ishenxalfa_row;
    trustheavy=(ishenxh_row>xhtrust) & (ishenaheav_row>aheavtrust);
    % the 1E-30 below just avoids division by 0
    ytot1 = ytot1temp + ishenxh_row./(ishenaheav_row+1E-30).*(trustheavy);
    zbarxxtemp  = ishenxprot_row + 0.5d0*ishenxalfa_row;
    x = ishenzheav_row./ishenaheav_row;
    zbarxx = zbarxxtemp + x.*ishenxh_row.*(trustheavy);

    abarnum    = 1.0d0./ytot1;
    abar=abarnum;
    zbar    = zbarxx .* abar;
    % yelocal should be equal to ishenyp_row, but don't check, just fix below
    yelocalafter = zbar./abarnum;
    
    % GODMARK
    % could compare yelocalafter with ishenyp_row to make sure same to
    % machine accuracy
    errorye = abs(yelocalafter-ishenyp_row)./(abs(yelocalafter)+abs(ishenyp_row));
    mynotnan = isfinite(errorye);
    notnanerrorye(mynotnan)=errorye(mynotnan);
    problems(notnanerrorye>0.1)=1;
    sumproblems=sum(sum(sum(problems)));
    if(sumproblems>=1)
      fprintf('Problems with ye: %d\n',sumproblems);
    end

    % another method to preserve array structure so can debug where
    % problem is
    erroryecopy=errorye;
    erroryecopy(isnan(erroryecopy))=0;
    problems=(erroryecopy>0.1);
    sumproblems=sum(sum(sum(problems)));
    if(sumproblems>=1)
      fprintf('ver2: Problems with ye: %d\n',sumproblems);
    end
    
  
  end  

