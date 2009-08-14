function eos_extract()

%TODO:
% 1) In principle, should use fakelsoffset and fakeentropylsoffset in functions that don't involve independent variables, but corrections extremely small.
  
  
  
% TODO:
% 1) Perhaps use 2D instead of 1D interpolation : No, because not changing density, just remapping T or U,P,CHI,S
% 2) Figure out how to use higher-than-linear order (what's wrong with v5cubic or spline) : Done, just need to take more care with crazy extrapolations and how to avoid them using faketkof? and tkof? quantities.
% 3) For temp's, can still use linear or whatever required for good NaN behavior to define boundary of valid EOS : Done as per #2 above.
% 4) Sure Stot(T) should be monotonic?  Appears highly non-monotonic in a large portion of rho,T space.  Perhaps reconsider monotonization at end rather than at beginning.
%    Then would ensure monotonic functions (e.g.) for P(rho0,chi), so final result would be reasonable for inversion.  Monotonizing in T-space might lead (after all interpolations) to non-monotonic result in chi-space.
% 5) gradient should be fine as long as done on super-sampled higher-order interpolated quantities! : Yes.
% 6) Use 2D interpolation for downsampling as well as supersampling : No, as #1.
% 7) Perhaps fill-in cs2 drop-outs with ideal gas versions? : drop-outs probably related to entropy issue -- how does cs2helm compare to my cs2?
% 8) In HARM and Matlab, define ideal gas as ideal but gamma=4/3 or 5/3 depending upon if beyond where nuclear density ensures definitely 5/3.
%    Switch in HARM based upon rhonuccode = rhonuccgs*c*c/rhounit

  
% OLD TODO (that isn't important for now):
% 1) Need to specify base and log-offset for each variable
% 2) Then need to output a table for each (rho(X)Tdynorye(X)Hcm) the umin, utotmax, ptotmin, ptotmax, chimin, chimax.  This way I can get high accuracy for T~0Kelvin by starting the table near T~0 and resolving the space near T~0.
%  ACTUALLY: data is already really there.  I just need to use MATLAB to find min/max of each and make sure to have that as first and last data point for each rho,H,T.  Then in HARM I use 0 and N-1 values as min/max.  Only issue is amount of computation using logs and pow's to process that raw data.
% Should probably still store it separately, but can be done in HARM instead of MATLAB
% In matlab each N points should range from the min/max found for EACH rho,H,T instead of globally
%%%%%%%%%%%%%%%%%
% Aassume for now offset is 0 and base is 10
% For now use functional offset

      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SOME PARAMETERS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % linear ruins benefit of supersampling and truncates values near edges (i.e. doesn't extrapolate)
%  interptype='linear';
%  interptype2='linear';

%    interptype='cubic';
%    interptype2='cubic';

    interptype='pchip';
    interptype2='pchip';

  % Can't use spline since too oscillatory and introduces sometimes crazy deviant values (e.g. 28 orders of magnitude oscillation near temperature edges)
%  interptype='spline';
%  interptype2='spline';
  % force linear for fake temperature variable interpolation so can put in NaN in correct place inside spline interpolation temperature variables
  interptypefaketemp='linear';
  
  
  % reaches "out of bounds" and uses to get local values
  %  interptype='v5cubic';
  %  interptype2='v5cubic';
  % http://www.mathworks.com/access/helpdesk/help/techdoc/index.html?/access/helpdesk/help/techdoc/ref/interp2.html&http://www.google.com/search?q=interp2&ie=utf-8&oe=utf-8&aq=t&rls=org.mozilla:en-US:official&client=firefox-a

  
  % whether to smooth input quantities that will be used as independent variables
  % makes utot and so utotdiff less accurate (weights to larger values, so final HARM estimate of temperature is systematiclaly lower)
  smoothinputs=0; 
  
  % whether to smooth before doing derivatives
  smoothbeforederivativestep1=1;
  smoothbeforederivativestep2=1;


  % whether to use consolidate to ensure unique functions (no longer needed if monotonic enforced)
  consolid=1;
  
  % whether to use cleanvar method (1) or simple isfinite method (0)
  usecleanvar=0;
  
  
  % set whether on windows laptop (interactive) or in linux (matlab script)
  matlabscript=1;

  
  
  % whether to fix utot if using old Kaz code result where I didn't subtract
  % off rhob c^2
  utotfix = 0;

  % whether to fix stot if using old HELM/TIMMES code where forgot to add entropy of rest-mass of electrons
  % or general fix
  stotfix = 1;

  % whether to use analytical fit or numerical values to set degenerate (offset) values
  utotdegenanalytic=0;

  % 0 = log but "0" is u=0.0
  % 1 = log but "0" is utotoffset
  % 2 = log with "0" as utotoffset, but choosen lower and upper "u" as well
  % whether to set "0" for log output as degeneracy fitting formula (to be added by in when inside HARM)
  % trying to get better temperature resolution for the low-temperature domain
  utotdegencut = 2;


  % whether to force functions to be monotonic as functions of tk
  forcemonotk = 1;

  % whether to "clean" solution if out of bounds
  doclean=1;
  
  % whether to fix-up derivatives at end if detect bad
  finaldoclean=1;

  % whether to specify a limit log(u) range (for zoom tables or if user knows that automatically determined range is too large for some reason)
  %specifylurange=1

  % whether to compute sound speed before or after interpolation
  % 0=compute after using entropy (HELM's entropy is messed up, but LS and Shen are fine)
  % 1=compute before using only non-entropy data (BEST)
  % 2=compute both and use smaller of two versions if non-negative and use larger if one is negative (avoids difference errors) (shouldn't use for HELM since HELM's entropy is messed up)
  preinterpsoundspeed=0;




  %CONTOL=1E-13;
  % chosen so small that only exact equality counts
  CONTOL=1E-17;
  CONTOL2=1E-17;

  OUTBOUNDSVALUE=1E-20;

  
  


  % nx
  % ny
  % nc=number of columns
  %fid=fopen('/home/jondata/data.head');
  %dir='/home/jondata/grmhd-a.9375-rout400new/';
  %dir='f:\\matlabscripts\\matlabscripts\\';
  %dir='C:\\Documents and Settings\\jon\\My Documents\\';
  if matlabscript==0
    % BELOW USED ON laptop
    %dir='C:\\Documents and Settings\\jon\\My Documents\\grbjet\\eoslarge\\';
    dir='C:\\Documents and Settings\\jon\\My Documents\\grbjet\\eoslargesimple\\';
  end
  if matlabscript==1
    % BELOW USED in linux
    dir='./';
  end
  %dir='C:\\Documents and Settings\\jon\\My Documents\\eossmall\\';
  prefix='eos';
  prefixo='eosother';
  prefix2='eosnew';
  prefix3='eosdegennew';
  prefix4='eosparms';
  prefix5='eosmonodegen';
  suf1='.head';
  suf2='.dat';
  suf3='.debug';

  % BELOW are from const.dek -- should really read these in
  % speed of light (cm/s)
  c = 2.99792458E10;
  % Boltzmann's constant in erg/K
  kb = 1.380658E-16;
  % Planck's constant
  % hbar = 1.054592483915517E-27;
  % below 4 only used for stotfix==1
  me      = 9.10938215E-28;
  avoreal = 6.02214179d23;
  mb      = (1.0D0/avoreal);
  ergPmev = 1.782661758E-30*1.0E3*c*c;
  
  
  file1=strcat(dir,prefix,suf1);
  file2=strcat(dir,prefix,suf2);
  file3=strcat(dir,prefix2,suf2);
  file4=strcat(dir,prefixo,suf2);
  file5=strcat(dir,prefix2,suf1);
  file6=strcat(dir,prefix3,suf2);
  file7=strcat(dir,prefix4,suf1);
  file8=strcat(dir,prefix5,suf2);
  filedebug=strcat(dir,prefix,suf3);
  
  % open debug file
  fiddebug=fopen(filedebug,'w');
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %
  % READ HEADER
  %
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  

  
  fid=fopen(file1,'rt');
  % 3+5+15 = 23 total header entries
  [myhead,count]=fscanf(fid,'%d',[3]);
  % 3 header things
  ii=1;
  whichrnpmethod=myhead(ii); ii=ii+1;
  whichynumethod=myhead(ii); ii=ii+1;
  whichhcmmethod=myhead(ii); ii=ii+1;
  
  % 5 header things
  % for input
  %          First row indicates:
  %     1) which data output type
  %     2) number independent variables
  %     3) number of primary data columns
  %     4) number of auxillary data columns
  %     5) out of primary data, how many "extra" variables
  %     8 original values (5 indeps and 3 vars)
  [myhead,count]=fscanf(fid,'%d',[5]);
  ii=1;
  whichdatatype=myhead(ii); ii=ii+1;
  ndim=myhead(ii); ii=ii+1;
  nc=myhead(ii); ii=ii+1;
  nco=myhead(ii); ii=ii+1;
  numextras=myhead(ii); ii=ii+1;

  % 5 indeps and 3 things each indep = 15 total + lsoffset + fakelsoffset
  [myhead,count]=fscanf(fid,'%d %g %g',[3]);
  ii=1;
  nrhob=myhead(ii); ii=ii+1;
  rhobmin=myhead(ii); ii=ii+1;
  rhobmax=myhead(ii); ii=ii+1;

  [myhead,count]=fscanf(fid,'%d %g %g',[3]);
  ii=1;
  ntk=myhead(ii); ii=ii+1;
  tkmin=myhead(ii); ii=ii+1;
  tkmax=myhead(ii); ii=ii+1;
  
  [myhead,count]=fscanf(fid,'%d %g %g',[3]);
  ii=1;
  ntdynorye=myhead(ii); ii=ii+1;
  tdynoryemin=myhead(ii); ii=ii+1;
  tdynoryemax=myhead(ii); ii=ii+1;

  [myhead,count]=fscanf(fid,'%d %g %g',[3]);
  ii=1;
  ntdynorynu=myhead(ii); ii=ii+1;
  tdynorynumin=myhead(ii); ii=ii+1;
  tdynorynumax=myhead(ii); ii=ii+1;

  [myhead,count]=fscanf(fid,'%d %g %g',[3]);
  ii=1;
  nhcm=myhead(ii); ii=ii+1;
  hcmmin=myhead(ii); ii=ii+1;
  hcmmax=myhead(ii); ii=ii+1;

  % lsoffset stuff
  [myhead,count]=fscanf(fid,'%g %g',[2]);
  ii=1;
  lsoffset=myhead(ii); ii=ii+1;
  fakelsoffset=myhead(ii); ii=ii+1;


  fclose(fid);
  
  
  % set log min/max
  lrhobmin=log10(rhobmin);
  lrhobmax=log10(rhobmax);
  ltkmin=log10(tkmin);
  ltkmax=log10(tkmax);
  ltdynoryemin=log10(tdynoryemin);
  ltdynoryemax=log10(tdynoryemax);
  ltdynorynumin=log10(tdynorynumin);
  ltdynorynumax=log10(tdynorynumax);
  lhcmmin=log10(hcmmin);
  lhcmmax=log10(hcmmax);

  
  % Set steps (should be consistent with kazloopfunctions.f)
  stepltk=(ltkmax-ltkmin)/(ntk-1);
  steplrhob=(lrhobmax-lrhobmin)/(nrhob-1);

  if(whichrnpmethod==0 || ntdynorye<=1)
    stepltdynorye=0;
  else
    stepltdynorye=(ltdynoryemax-ltdynoryemin)/(ntdynorye-1);
  end

  if(whichynumethod==0 || ntdynorynu<=1)
    stepltdynorynu=0;
  else
    stepltdynorynu=(ltdynorynumax-ltdynorynumin)/(ntdynorynu-1);
  end

  if(whichhcmmethod==0 || nhcm<=1)
    steplhcm=0;
  else
    steplhcm=(lhcmmax-lhcmmin)/(nhcm-1);
  end







  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Read parameter information
  %
  % E.g., eosparms.head contains:
  % 1E-8 1E-8 1E-8 1E-3 1E-13 1E-13 1E-13 1E-6
  % 1E-16 1E-16 1E-16 1E-16 1E-1 1E-1 1E-1 1E-1
  %
  %
  numparms=16;
  
  fid=fopen(file7);
  [myhead,count]=fscanf(fid,'%g',[numparms]);
  fclose(fid);

  ii=1;
  % for degenfit
  UTOT0=myhead(ii); ii=ii+1;
  PTOT0=myhead(ii); ii=ii+1;
  CHI0=myhead(ii); ii=ii+1;
  STOT0=myhead(ii); ii=ii+1;
  
  UTOTF=myhead(ii); ii=ii+1;
  PTOTF=myhead(ii); ii=ii+1;
  CHIF=myhead(ii); ii=ii+1;
  STOTF=myhead(ii); ii=ii+1;

  % for Rin
  UTOTIN0=myhead(ii); ii=ii+1;
  PTOTIN0=myhead(ii); ii=ii+1;
  CHIIN0=myhead(ii); ii=ii+1;
  STOTIN0=myhead(ii); ii=ii+1;
  
  UTOTINF=myhead(ii); ii=ii+1;
  PTOTINF=myhead(ii); ii=ii+1;
  CHIINF=myhead(ii); ii=ii+1;
  STOTINF=myhead(ii); ii=ii+1;
  
  
  
  
  
  
  
  
  
  
  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 
  % BIG LOOP OVER EXTRA INDEPENDENT VARIABLES not used for complicated tasks like derivatives or interpolation
  %
  
  truentdynorye=ntdynorye;
  truentdynorynu=ntdynorynu;
  truenhcm=nhcm;
  
  % fake so don't have to change interior code from old way where did full memory at once
  % assuming never doing big memory method with ntdynorynu>1
  ntdynorye=1;
  ntdynorynu=1; % actually not used
  nhcm=1;
  

  
  
  if utotdegencut==2
    % this method doesn't require 2 passes since min/max of grid are setup for *each* \rho_0, Y_e, Y_\nu, H
    numpasses = 1;
    require2passes = 0;
  end
  
  if utotdegencut==1 || utotdegencut==0
    % number of total passes, used to get min/max of independent variables
    numpasses = 2;
    require2passes = 1;

  % initialize for first pass for every hcm and tdynorye
    % initially choose very non-min and non-max to always gets overwritten
    lutotdiffoutmin=1E50;
    lutotdiffoutmax=-1E50;
    lptotdiffoutmin=1E50;
    lptotdiffoutmax=-1E50;
    lchidiffoutmin=1E50;
    lchidiffoutmax=-1E50;
    lstotdiffoutmin=1E50;
    lstotdiffoutmax=-1E50;
    lsspecdiffoutmin=1E50;
    lsspecdiffoutmax=-1E50;
  end

  
  
  

  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % loop over 2 passes
  % first pass to obtain correct most general min/max
  % second pass to do normal calculation and outputs
  for passiter=1:numpasses

    
    fprintf(fiddebug,'pass: %d require2passes: %d\n',passiter,require2passes);

    
    
    % Open all files
    fid2=fopen(file2,'rt');
    fid8=fopen(file8,'w');
    if 0
      % don't read eosother.dat
      fid4=fopen(file4,'rt');
    end
    fid3=fopen(file3,'w');
    fid6=fopen(file6,'w');


    
    
    
    % consistent order for loops as output and read-in into HARM
    % must also be consistent with how jon_helm.f writes order so read-in in correct order
    for hiter=1:truenhcm
    for ynuiter=1:truentdynorynu
      for titer=1:truentdynorye

        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % open eos.dat
        %
        %
        [mydata,count]=fscanf(fid2,'%g',[nc*nrhob*ntk*ntdynorye*nhcm]);

        % make vector into 6-D array
        temp=reshape(mydata,[nc,nrhob,ntk,ntdynorye,nhcm]);
        clear mydata;
        % set primary matricies to be each column(field) read in from SM data (i.e.
        % read data is setup with column / nx / ny / nz order
        mynewdata=permute(temp,[2,3,4,5,1]);
        clear temp;
        %
        % so now all functions are F(rhob,tk)
        %

        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % clean-up the input data
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for p=1:nrhob
          for q=1:ntk
            for r=1:ntdynorye
              for s=1:nhcm
                nanonline=0;
                for t=1:nc
                  nanonline=nanonline+isnan(mynewdata(p,q,r,s,t));
                end
                if nanonline>0
                  fprintf(fiddebug,'Found nan at p=%d q=%d r=%d s=%d\n',p,q,r,s);
                  % assume NaN's occur at higher densities and just copy from lower densities
                  % first 4 quantities are independent variables that shouldn't change
                  for t=5:nc
                    mynewdata(p,q,r,s,t) = mynewdata(p-1,q,r,s,t);
                  end
                  %A(isnan(A))=0
                end
              end
            end
          end
        end

        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %
        % Assign labels to physical quantities that need special attention
        %
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
        % {tdynorye, hcm, rhob, tk, etae, npratio, p_tot, rho_tot, s_tot, p_photon,
        % p_eleposi,p_N,p_nu}, extra1,etc.
        ii=1;
        rhob(:,:,:,:) = mynewdata(:,:,:,:,ii); ii=ii+1;
        tk = mynewdata(:,:,:,:,ii); ii=ii+1;
        tdynorye = mynewdata(:,:,:,:,ii); ii=ii+1;
        tdynorynu = mynewdata(:,:,:,:,ii); ii=ii+1;
        hcm = mynewdata(:,:,:,:,ii); ii=ii+1;
        % below ptot,utot,stot actually don't include neutrino contribution if
        % using whichrnpmethod=0,1, but are true totals if nhcm!=1
        % which will be reconstructed within HARM
        ptot = mynewdata(:,:,:,:,ii); ii=ii+1;
        utot = mynewdata(:,:,:,:,ii); ii=ii+1;
        stot = mynewdata(:,:,:,:,ii); ii=ii+1;

        
        
        
        if stotfix==1

          % BEGIN DEBUG:
          % will use sspec to check stot
          %sspecdimless = stot./(rhob./mb);
%          for p=1:nhcm
%            for o=1:ntdynorye
%              for n=1:ntk
%                for m=1:nrhob
%                  blob=(lsoffset*ergPmev/mb)/(kb);
%                  blob2=(lsoffset*ergPmev*rhob(m,n,o,p)/mb)/(kb*tk(m,n,o,p));
%                  fprintf(fiddebug,'BEFORE: %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n',stot(m,n,o,p),sspecdimless(m,n,o,p),rhob(m,n,o,p),tk(m,n,o,p),blob,blob2);
%                end
%              end
%            end
%          end
          % END DEBUG:

          % Temporary fix for stot since in HELM/TIMMES forgot to add rest-mass of e+/e- to entropy
          % Recall stot in [1/cc]
          % no longer so this since was stupid:
          %stot = stot + ((me.*c.*c./mb)*tdynorye.*rhob)./(kb.*tk);
          
          % Remove lsoffset term since shouldn't have been there in way it was introduced
          stot = stot - (lsoffset.*ergPmev.*rhob./mb)./(kb.*tk);
          
          % add-in fake entropy offset since not yet in table.  Ensure use in HARM still!
          %fakeentropylsoffset=8.921;
          % minimum in data file is -12.94 for the dimensionless entropy
          %fakeentropylsoffset=14.0;
          % machine-error occurs with subtraction
          fakeentropylsoffset=1E2;
          stot = stot + fakeentropylsoffset.*(rhob./mb);
          
          % recompute sspecdimless
          % will use sspec to check stot
          sspecdimless = stot./(rhob./mb);

          % Note that entropy offset at nuclear-to-nonnuclear EOS boundary is unchanged, so above accounts for all new EOS changes

          % BEGIN DEBUG:
%          for p=1:nhcm
%            for o=1:ntdynorye
%              for n=1:ntk
%                for m=1:nrhob
%                  fprintf(fiddebug,'AFTER: %21.15g %21.15g\n',stot(m,n,o,p),sspecdimless(m,n,o,p));
%                end
%              end
%            end
%          end
          % END DEBUG:

          stotmin=min(min(min(min(stot))));
          stotmax=max(max(max(max(stot))));
          if(stotmin<=0.0)
            fprintf(fiddebug,'stot must be positive %21.15g %21.15g\n',stotmin,stotmax);
            quit
          end


          sspecdimlessmin=min(min(min(min(sspecdimless))));
          sspecdimlessmax=max(max(max(max(sspecdimless))));
          if(sspecdimlessmin<=0.0)
            fprintf(fiddebug,'sspecdimless must be positive %21.15g %21.15g\n',sspecdimlessmin,sspecdimlessmax);
            quit
          end

          %stot = stot/kb; % TEMPORARILY CONVERT ENTROPY TO 1/cc form for table that was in erg/K/cc form

        end
        
        
        

        % rest of below don't come into detailed calculations and are just
        % other quantities to interpolate from T to UTOT
        extraii=ii;
        % instead of using extra1,extra2, etc., just directly access mynewdata
        %  extra = mynewdata(:,:,:,:,extraii+extraiter); ii=ii+1;

        % COMPUTE CHI (from now on have to deal with chi separately)
        chi = utot + ptot;

        % COMPUTE sspec (from now on have to deal with sspec separately from stot) [as an independent variable only]
        rhobcsq = rhob.*c.*c;
        sspec = stot./(rhobcsq);


        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %
        % Perform some smoothing of native input functions
        %
        % Helps with accuracy on derivatives esp. d/dlrho0 at high temperatures and low densities
        %
        % But, problem since the average systematically lowers value average log-log averaging
        %
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if smoothinputs
          
          for q=1:ntdynorye
            for r=1:nhcm
              
              % assumes p,u,s >0
              roughlptot(:,:)=log10(ptot(:,:,q,r));
              roughlutot(:,:)=log10(utot(:,:,q,r));
              roughlstot(:,:)=log10(stot(:,:,q,r));

              % Smooth functions since sharp boundaries cause speed to be artificially large
              ptot(:,:,q,r) = 10.0.^(moving_average2(roughptot(:,:),1,1));
              utot(:,:,q,r) = 10.0.^(moving_average2(roughutot(:,:),1,1));
              stot(:,:,q,r) = 10.0.^(moving_average2(roughstot(:,:),1,1));

            end
          end
          
          % recomptue for independents still
          chi = utot+ptot;
          sspec = stot./(rhobcsq);
          
        end
        

        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %
        % Force monotonicity of the U(T), P(T), and S(T) otherwise multiple
        % solutions for T for each U,P,S
        %
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % whether to force monotonicity for functions as functions of temperature
        % otherwise EOS is non-convex and multivalued for a single internal energy
        % for the HARM code, all that really matters is that chi=u+p is single valued
        %
        %
        % http://www.sciencedirect.com/science?_ob=ArticleURL&_udi=B6WHY-4TP49TB-2&_user=145269&_rdoc=1&_fmt=&_orig=search&_sort=d&_docanchor=&view=c&_acct=C000012078&_version=1&_urlVersion=0&_userid=145269&md5=35910ec6d3376be6c25aa7d25e3259f4
        %
        % Need EOS to be convex.
        %
        if forcemonotk==1
          
          % will use sspec to process stot's montonicity

          for p=1:nrhob
            for q=1:ntdynorye
              for r=1:nhcm
                % only need to force these to be monotonic so chi is monotonic so can have single-valued inversion
                % ptot and utot (and so chi)  need to be monotonically increasing functions of temperature in order to obtain a single temperature for a single ptot,utot,chi
                utot(p,:,q,r)=monotonize(utot(p,:,q,r));

                %p
                %q
                %r
                %stot(p,:,q,r)

                % entropy should be monotonic so that can be used as independent variable, but only needs to be monotonic for P(S) and U(S), not T(S).
                % montonizing stot means sspec will be monotonic as required for entropy evolution
                sspec(p,:,q,r)=monotonize(sspec(p,:,q,r));


                if 0
                  % monotonize chi and recompute ptot from monotonized chi and utot
                  %chi(p,:,q,r)=monotonize(chi(p,:,q,r));
                  % recompute ptot from monotonized utot and chi since it's more important to be monotonic in chi and ptot
                  % leaves ptot<0 if chi chopped alot, so use old method of fixing utot and ptot
                  %ptot(p,:,q,r)=chi(p,:,q,r)-utot(p,:,q,r);
                end

                if 1
                  ptot(p,:,q,r)=monotonize(ptot(p,:,q,r));
                  % Re-compute chi using monotonized versions of u and p
                  %chi(p,:,q,r) = utot(p,:,q,r) + ptot(p,:,q,r);
                end
              end
            end
          end

          % recompute chi
          chi=utot+ptot;
          % recover stot as correctly monotonic (should really see why stot/rhob so non-monotonic due to nuclear terms and how I extrapolate it)
          stot=(sspec).*(rhobcsq);

          
        end %%%%% end if forcemonotk==1




        
        

        if utotfix==1
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %
          % Temporary fix for utot since I was subtracting off rhob instead of rhob
          % c^2
          %
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%

          utotnew= utot + rhob - rhob.*c.*c;

          % by itself the above generates negative utot's, so must put a lower limit
          % this wouldn't occur if did things correctly in Kaz's code

          isneg = (utotnew<=0);
          myfloor = utotnew.*0 + 1E-20;
          utot(isneg) = myfloor(isneg);
          
          % recompute chi
          chi=utot+ptot;

        end

        

        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % set "0" of utot, ptot, and chi so that outputted grid resolves
        % lowest temperatures when u,p,chi vary very little
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
        if utotdegenanalytic==1

          %%%%%%%%%%%%%%%%%%%%%%%%%
          % use analytical 0-point
          %%%%%%%%%%%%%%%%%%%%%%%%%

          %utot = utot - udegenfit_fun(rhob);
          %ptot = ptot - pdegenfit_fun(rhob);
          % chi is not read-in, and since we just modified utot and ptot, now anything (such as chi) computed from these will already be offset

          % instead of above subtraction, treat as general offset, so that output is uniform when including offset, but things are computed with standard utot and ptot
          utotoffset = udegenfit_fun(rhob);
          ptotoffset = pdegenfit_fun(rhob);
          % essentially, this means offset's only affect gridding of chosen u,p,chi as dependent variables
          % that is, we need to modify lu

          % set chidegenfit
          chidegenfit=udegenfit_fun(rhob) + pdegenfit_fun(rhob);
          chioffset  = utotoffset+ptotoffset;
          
          % no stot version for this analytic method
          % fake:
          % note that stot is in 1/cc form
          stotoffset = (utotoffset+ptotoffset)./(kb.*tk)
          
          sspecoffset = stotoffset./rhob;

        end %%%%%if utotdegenanalytic==1


        
        
        if  utotdegenanalytic==0

          %%%%%%%%%%%%%%%%%%%%%%%%%
          % use numerical 0-point
          %%%%%%%%%%%%%%%%%%%%%%%%%

          % pick lowest temperature as bottom
          % Duplicate this for other temperatures
          if(0)
            for q=1:ntk

              utotdegenfit(:,q,:,:) = utot(:,1,:,:);
              ptotdegenfit(:,q,:,:) = ptot(:,1,:,:);
              chidegenfit(:,q,:,:) = chi(:,1,:,:); % chi treated independently
            end
          end

          if(1)
            % don't assume monotonicity
            for p=1:nhcm
              for o=1:ntdynorye
                for m=1:nrhob

                  %%%%%%%%%% U
                  for n=1:ntk
                    mytemporary(n) = utot(m,n,o,p);
                  end
                  mytempmin = min(mytemporary);
                  mytempmax = max(mytemporary);
                  for n=1:ntk
                    utotdegenfit(m,n,o,p) = mytempmin;
                    utotmax(m,n,o,p) = mytempmax;
                  end

                  %%%%%%%%%% P
                  for n=1:ntk
                    mytemporary(n) = ptot(m,n,o,p);
                  end
                  mytempmin = min(mytemporary);
                  mytempmax = max(mytemporary);
                  for n=1:ntk
                    ptotdegenfit(m,n,o,p) = mytempmin;
                    ptotmax(m,n,o,p) = mytempmax;
                  end

                  %%%%%%%%%% CHI
                  for n=1:ntk
                    mytemporary(n) = chi(m,n,o,p);
                  end
                  mytempmin = min(mytemporary);
                  mytempmax = max(mytemporary);
                  for n=1:ntk
                    % chi treated independently
                    chidegenfit(m,n,o,p) = mytempmin;
                    chimax(m,n,o,p) = mytempmax;
                  end

                  %%%%%%%%%% S
                  for n=1:ntk
                    mytemporary(n) = stot(m,n,o,p);
                  end
                  mytempmin = min(mytemporary);
                  mytempmax = max(mytemporary);
                  for n=1:ntk
                    stotdegenfit(m,n,o,p) = mytempmin;
                    stotmax(m,n,o,p) = mytempmax;
                  end

                  %%%%%%%%%% SS
                  for n=1:ntk
                    mytemporary(n) = sspec(m,n,o,p);
                  end
                  mytempmin = min(mytemporary);
                  mytempmax = max(mytemporary);
                  for n=1:ntk
                    sspecdegenfit(m,n,o,p) = mytempmin;
                    sspecmax(m,n,o,p) = mytempmax;
                  end
                
                end
              end
            end
          end %%%%end if 1

          % now include multiplicative offset, where offset was chosen ahead of time so that utotdiff,ptotdiff are small and positive for all temperatures for each rhob,hcm,tdynorye,tdynorynu
          % This procedure doesn't make sense if the degen value is near 0, then need to offset by some fraction of (max-min)
          % chi treated independently

          %utotoffset = utotdegenfit   - abs(utotdegenfit)*(UTOT0);
          %ptotoffset = ptotdegenfit   - abs(utotdegenfit)*(PTOT0);
          %chioffset  = chidegenfit - abs(utotdegenfit)*(CHI0);
          %stotoffset = stotdegenfit - abs(stotdegenfit)*(STOT0);

          % general approach so log-intervals always equally resolved in log of temperature
          % concept is that minimum sets no scale.  Scale only set by maximum (assumes maximum is not too far from where quantity becomes as degenerate as not degenerate

          
          % also now vary offset since at higher density du/dT is much smaller
          numshift = ntk;
          
          %UTOTF=1.0-1E-14;
          %UTOTF=UTOT0;
          UTOTD=10.^(log10(UTOT0) + (log10(rhob)-lrhobmin)./(lrhobmax-lrhobmin).*(log10(UTOTF)-log10(UTOT0)));
          %PTOTF=1.0-1E-14;
          %PTOTF=PTOT0;
          PTOTD=10.^(log10(PTOT0) + (log10(rhob)-lrhobmin)./(lrhobmax-lrhobmin).*(log10(PTOTF)-log10(PTOT0)));
          %CHIF=1.0-1E-14;
          %CHIF=CHI0;
          CHID=10.^(log10(CHI0) + (log10(rhob)-lrhobmin)./(lrhobmax-lrhobmin).*(log10(CHIF)-log10(CHI0)));
          %STOTF=1.0-1E-14;
          %STOTF=STOT0;
          STOTD=10.^(log10(STOT0) + (log10(rhob)-lrhobmin)./(lrhobmax-lrhobmin).*(log10(STOTF)-log10(STOT0)));
          
          utotoffset = utotdegenfit   - max(10.^(log10(abs(utotmax))-numshift),abs(utotdegenfit).*(UTOTD));
          ptotoffset = ptotdegenfit   - max(10.^(log10(abs(ptotmax))-numshift),abs(ptotdegenfit).*(PTOTD));
          chioffset  = chidegenfit - max(10.^(log10(abs(chimax))-numshift),abs(chidegenfit).*(CHID));
          stotoffset = stotdegenfit   - max(10.^(log10(abs(stotmax))-numshift),abs(stotdegenfit).*(STOTD));

          % could have its own SSTOTD, but should be ok without one.  Below is correct since linear in stot.
          sspecoffset = stotoffset./(rhobcsq);

        end %%%%% end if utotdegenanalytic==0



        
        
        if utotdegencut==2
          % More advanced case where grid is chosen such that grid is r = R0 + exp(x1), with r=utot and R0=utotoffset with r-R0 = utotdiff,
          % but where Rin=chosen separately and Rout=chosen separately
          % Here grid is U = U0 + (Uin-U0)*pow( (Uout-U0)/(Uin-U0) , i/N ) for N points with "i" as the index
          % Only have to ensure that Uin>U0 and Uout>Uin.
          % That is, do NOT have to have all U's in original table such that U>U0.  Those U's missed will deresolve away as NaN's in terms of the new index fraction "i/N"
          % This means utotoffset does not have to be exactly offsetted from lowest values like done above when setting utotoffset, etc.

          % used for grid of data
          % meaning of utotdiff or ptotdiff or chidiff or stotdiff is "i/N" or the fraction of the grid from 0 to 1

          % reasonable upper limits -- upper limits haven't been much of a problem except for a 1 zone drop-out at the top.
          utotout = utotmax;
          ptotout = ptotmax;
          chiout = chimax;
          stotout = stotmax;
          sspecout = sspecmax;

          UTOTIND=10.^(log10(UTOTIN0) + (log10(rhob)-lrhobmin)./(lrhobmax-lrhobmin).*(log10(UTOTINF)-log10(UTOTIN0)));
          PTOTIND=10.^(log10(PTOTIN0) + (log10(rhob)-lrhobmin)./(lrhobmax-lrhobmin).*(log10(PTOTINF)-log10(PTOTIN0)));
          CHIIND=10.^(log10(CHIIN0) + (log10(rhob)-lrhobmin)./(lrhobmax-lrhobmin).*(log10(CHIINF)-log10(CHIIN0)));
          STOTIND=10.^(log10(STOTIN0) + (log10(rhob)-lrhobmin)./(lrhobmax-lrhobmin).*(log10(STOTINF)-log10(STOTIN0)));

          % utotdegencut==1 effectively starts with utotin=utotdegenfit
          % Start with utotdegenfit, but allow offset from that in proportion to utotdegenfit-utotoffset = Rinold - R0
          utotin = utotdegenfit   + abs(utotoffset).*(UTOTIND);
          ptotin = ptotdegenfit   + abs(ptotoffset).*(PTOTIND);
          chiin  = chidegenfit    + abs(chioffset).*(CHIIND);
          stotin = stotdegenfit   + abs(stotoffset).*(STOTIND);
          sspecin = stotin./(rhobcsq); % forces consistency to stotin since dont' have SSIND.  Above is linear in stot, so it may be ok to do directly with same STOTIND, but not an issue.

          % finally set "i/N" that always goes from 0 to 1
          lutotdiff = log10( (utot - utotoffset)./(utotin - utotoffset))./log10( (utotout - utotoffset)./(utotin - utotoffset));
          lptotdiff = log10( (ptot - ptotoffset)./(ptotin - ptotoffset))./log10( (ptotout - ptotoffset)./(ptotin - ptotoffset));
          lchidiff = log10( (chi - chioffset)./(chiin - chioffset))./log10( (chiout - chioffset)./(chiin - chioffset));
          lstotdiff = log10( (stot - stotoffset)./(stotin - stotoffset))./log10( (stotout - stotoffset)./(stotin - stotoffset));
          % Easier to compute sspecdiff independently rather than inverting stotdiff and then appying relationship and then reobtaining sspecdiff
          lsspecdiff = log10( (sspec - sspecoffset)./(sspecin - sspecoffset))./log10( (sspecout - sspecoffset)./(sspecin - sspecoffset));

        end

        
        if utotdegencut==1 || utotdegencut==0
          % Not used, just filler
          utotout = utotmax;
          ptotout = ptotmax;
          chiout = chimax;
          stotout = stotmax;
          sspecout = sspecmax;
          
          % Not used, just filler
          utotin = utotdegenfit;
          ptotin = ptotdegenfit;
          chiin  = chidegenfit;
          stotin = stotdegenfit;          
          sspecin = sspecdegenfit;
        end

        if utotdegencut==1
          % Simpler case where grid is chosen such that grid is r = R0 + exp(x1), with r=utot and R0=utotoffset with r-R0 = utotdiff, where Rin=min(r-R0) and Rout=max(r-R0) for all existing r's.

          % used for grid of data:
          lutotdiff  = log10(utot - utotoffset);
          lptotdiff  = log10(ptot - ptotoffset);
          lchidiff   = log10(chi  - chioffset);
          lstotdiff  = log10(stot - stotoffset);
          lsspecdiff = log10(sspec - sspecoffset);
        end

        if utotdegencut==0
          % Simplest case where grid is chosen as r=exp(x1) with Rin=min(r) and Rout=max(r-R0) for all existing r's.
          
          lutotdiff = log10(utot); % used for grid of data
          lptotdiff = log10(ptot); % used for grid of data
          lchidiff  = log10(chi);  % used for grid of data
          lstotdiff = log10(stot); % used for grid of data
          % sspecdiff just linear in stot, so below is correct
          lsspecdiff = log10(sspec);

        end



        
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % If passiter==2, then really do computations.  Otherwise, if passiter==1, then just reading in and setting up independent variables so can get global min/max over all independent varaibles
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
        if(passiter==2 || require2passes==0)
          


          
          
          %%%%%%%%%%%%%%%%%%%%%%%%
          %
          % Compute sound speed from temperature-based data.  This is taken from HELM TM EOS
          % TM took it from:: references: cox & giuli chapter 24 ; timmes & swesty apj 1999
          %
          %%%%%%%%%%%%%%%%%%%%%%%

          % 1 below means always set this so can compare to other sound speed
          if 1 || preinterpsoundspeed==1 ||  preinterpsoundspeed==2

            %..the temperature and density exponents (c&g 9.81 9.82) 
            %..the specific heat at constant volume (c&g 9.92)
            %..the third adiabatic exponent (c&g 9.93)
            %..the first adiabatic exponent (c&g 9.97) 
            %..the second adiabatic exponent (c&g 9.105)
            %..the specific heat at constant pressure (c&g 9.98) 
            %..and relativistic formula for the sound speed (c&g 14.29)

            %      JCM: my definitions from HELM TM EOS
            den=rhob.*c.*c;
            pres=ptot;
            temp=tk.*kb;
            % specific energy
            ener=utot./den;
            % speed of light
            %       clight=c;

            
            % get derivatives
            for q=1:ntdynorye
              for r=1:nhcm


                if smoothbeforederivativestep1==1
                  roughptot(:,:)=ptot(:,:,q,r);
                  roughener(:,:)=ener(:,:,q,r);
                  roughtemp(:,:)=temp(:,:,q,r);
                  roughden(:,:)=den(:,:,q,r);

                  % Smooth functions since sharp boundaries cause speed to be artificially large
                  myptot(:,:)=moving_average2(roughptot(:,:),1,1);
                  myener(:,:)=moving_average2(roughener(:,:),1,1);
                  mytemp(:,:)=moving_average2(roughtemp(:,:),1,1);
                  myden(:,:)=moving_average2(roughden(:,:),1,1);
                else
                  myptot(:,:)=ptot(:,:,q,r);
                  myener(:,:)=ener(:,:,q,r);
                  mytemp(:,:)=temp(:,:,q,r);
                  myden(:,:)=den(:,:,q,r);
                end
                
                % use rho c^2 so result is dimensionless
                rhocsqi=den(:,1,q,r)';
                %		rhocsqi2d=den(:,:,q,r);
                %		tsize=size(rhocsqi);
                %		sizerho=tsize(1);
                %		irho=1:sizerho;
                %		irhooffset=0.5:sizerho-0.5;
                lrhocsqi = log10(rhocsqi);
                %		dlrhocsqi = gradient(lrhocsqi);

                % use kb*T so result is dimensionless
                Ti=temp(1,:,q,r)';
                %		Ti2d=temp(:,:,1,1);
                %		tsize=size(Ti);
                %		sizeT=tsize(1);
                %		iT=1:sizeT;
                %		iToffset=0.5:sizeT-0.5;
                lTi = log10(Ti);
                %		dlTi = gradient(lTi);
                
                % given PofU(rho0,U,H) then derivative is dU, drho0, dH for independent variable order
                %

                % create things for change of variable
                [dtdlt(:,:), dtdld(:,:)] = gradient(mytemp(:,:),lTi,lrhocsqi);
                [dddlt(:,:), dddld(:,:)] = gradient(myden(:,:),lTi,lrhocsqi);
                
                  
                [dpresdlt(:,:), dpresdld(:,:)] = gradient(myptot(:,:),lTi,lrhocsqi);
                [denerdlt(:,:), denerdld(:,:)] = gradient(myener(:,:),lTi,lrhocsqi);

                % old comment about gradient:
                % correct for gradient position offset
                % That is, gradient takes simple difference divided by simple difference of positions
                % It does not interpolate back to central location
                % GODMARK: Apparently only at edges does it reduce to simple difference.
                % Otherwise it does put gradient at function location
                %dpresdlt(:)
                %dpresdld(:)
                %denerdlt(:)
                %denerdld(:)
                
                
                % Get actual derivative (assume rho and T are independent -- which is true since just took log, didn't mix them)
                dpresdt(:,:,q,r)=dpresdlt(:,:)./dtdlt(:,:);
                dpresdd(:,:,q,r)=dpresdld(:,:)./dddld(:,:);
                
                denerdt(:,:,q,r)=denerdlt(:,:)./dtdlt(:,:);
                denerdd(:,:,q,r)=denerdld(:,:)./dddld(:,:);
                

              end
            end

            
            %      dP/dT|rhob
            %       dpresdt

            %      dP/drhob|T
            %       dpresdd

            %      dener/dT|rhob
            %       denerdt

            %      HELM TM EOS extra definitions
            deni=1./den;

            %      copied from HELM TM EOS code
            zz    = pres .* deni;
            zzi   = den ./ pres;
            chit  = (temp ./ pres) .* dpresdt;
            chid  = (abs(dpresdd)+1E-20) .* zzi;
            % don't allow negative cv
            cv    = abs(denerdt)+1E-20;
            x     = (zz .* chit) ./ (temp .* cv);
            gam3  = x + 1.0d0;
            gam1  = (chit .* x) + chid;
            nabad = x ./ gam1;
            gam2  = 1.0d0 ./ (1.0d0 - nabad);
            cp    = (cv .* gam1) ./ chid;
            % GODMARK: Because "ener" appears here, means can't just offset ener
            % by any constant
            z     = 1.0d0 + (ener + 1.0d0) .* zzi;
            
            %      Finally quantity is sound speed squared in dimensionless units
            cs2rhoT=(gam1 ./ z);
            
            % Because of how computed, this may be negative or >1.0, so at least
            % enforce not negative
            cs2rhoT=abs(cs2rhoT)+1E-20;

            
          
            % adjust speed of sound to be no smaller than 0 and no larger than 1
            for p=1:nhcm
              for o=1:ntdynorye
                for n=1:ntk % really ntk!  Still in T-space
                  for m=1:nrhob
                    
                    if cs2rhoT(m,n,o,p)>1.0
                      cs2rhoT(m,n,o,p)=1.0-CONTOL;
                    end
                    if cs2rhoT(m,n,o,p)<0.0
                      cs2rhoT(m,n,o,p)=0.0;
                    end

                  end
                end
              end
            end

            
            
          end  %%%%%% end if 1


          

          %%%%%%%%%%%%%%%%%%%%%%%%
          %
          % Output monotonized EOS so can check in SM
          %
          % Assumes nutotdiffout is same size as nptotout and nchiout
          %
          % 22 + extras
          %
          % Note that stot (and so its children quantities) may be modified if stotfix==1
          %
          % Treat chi as independent
          %
          %
          %%%%%%%%%%%%%%%%%%%%%%%
          
          fprintf(fiddebug,'Outputting monotinized EOS %d %d %d\n',hiter,titer,ynuiter);
          
          % for large memory method
          %m-1, n-1, o-1, p-1, ...
          % for small memory method
          
          for o=1:nhcm
            for p=1:ntdynorye
              for n=1:ntk
                for m=1:nrhob %+0                 +5                                      +5                                      +5                                      +5                                       +5                                      +5                                     +5                                      +5      +1
                  fprintf(fid8,'%3d %3d %3d %3d %3d %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g ', ...
                          m-1, n-1, titer-1, ynuiter-1, hiter-1, ...
                          rhob(m,n,o,p), tk(m,n,o,p), tdynorye(m,n,o,p), tdynorynu(m,n,o,p), hcm(m,n,o,p), ...
                          utot(m,n,o,p), ptot(m,n,o,p), chi(m,n,o,p), stot(m,n,o,p), sspec(m,n,o,p), ...
                          utotdegenfit(m,n,o,p), ptotdegenfit(m,n,o,p), chidegenfit(m,n,o,p), stotdegenfit(m,n,o,p), sspecdegenfit(m,n,o,p), ...
                          utotoffset(m,n,o,p), ptotoffset(m,n,o,p), chioffset(m,n,o,p), stotoffset(m,n,o,p), sspecoffset(m,n,o,p), ...
                          utotin(m,n,o,p), ptotin(m,n,o,p), chiin(m,n,o,p), stotin(m,n,o,p), sspecin(m,n,o,p), ...
                          utotout(m,n,o,p), ptotout(m,n,o,p), chiout(m,n,o,p), stotout(m,n,o,p), sspecout(m,n,o,p), ...
                          lutotdiff(m,n,o,p), lptotdiff(m,n,o,p), lchidiff(m,n,o,p), lstotdiff(m,n,o,p), lsspecdiff(m,n,o,p), ...
                          cs2rhoT(m,n,o,p) ...
                          );
                  for ei=extraii:extraii+numextras-1
                    fprintf(fid8,'%21.15g ', mynewdata(m,n,o,p,ei));
                  end
                  %                fprintf(fid8,'%3d %3d %3d %3d', ...
                  %                m-1, n-1, o-1, p-1 ...
                  %		);
                  
                  fprintf(fid8,'\n');
                end
              end
            end
          end



        end % end passiter==2 || require2passes==0 section



        % DEBUG
        % for ii=1:48 fprintf(fiddebug,'%21.15g %21.15g %21.15g\n',log10(utot(28,ii,1,1)/(c*c)),log10(utotnew(28,ii,1,1)/(c*c)),log10(tk(28,ii,1,1))); end


        

        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % some computed things as related to also temporary independent variables needed for passiter==1
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %        rhobcsq		= rhob.*c.*c;

        wtot = (rhobcsq+chi);

        % hspec is dimensionless
        hspec = wtot./rhobcsq;
        

        % stot is entropy density (erg/K/cc)
        % so compute sound speed using specific entropy (entropy per baryon)
        % below makes no sense with pure photons
        %sspec = stot./(rhobcsq);
        % can do below, but ideal gas case gets complicated for U[rho0,Sden]
        %sspec = stot./(rhobcsq+utot);

        %x=log10(rhob(:,:,:,8));
        %y=log10(tk(:,:,:,8));
        %z=log10(npratio(:,:,:,8)+1);
        %figure; contour(x,y,z);


        
        
        

        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Continue passiter==2 || require2passes==0 for real computations
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if(passiter==2 || require2passes==0)


          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % open q_volume.dat (GODMARK: outdated)
          if 0

            [mydata,count]=fscanf(fid4,'%g',[nco*nrhob*ntk*ntdynorye*nhcm]);

            % make vector into 5-D array
            temp=reshape(mydata,[nco,ntk,nrhob,ntdynorye,nhcm]);
            clear mydata;
            % set primary matricies to be each column(field) read in from SM data (i.e.
            % read data is setup with column / nx / ny / nz order
            mynewdata=permute(temp,[2,3,4,5,1]);
            clear temp;
            %
            %
            % so now all functions are F(rhob,tk,tdynorye,hcm)
            %

            %


            ii=1;
            etae = mynewdata(:,:,:,:,ii); ii=ii+1;
            xnuc = mynewdata(:,:,:,:,ii); ii=ii+1;
            npratio = mynewdata(:,:,:,:,ii); ii=ii+1;

            pphoton = mynewdata(:,:,:,:,ii); ii=ii+1;
            peleposi = mynewdata(:,:,:,:,ii); ii=ii+1;
            pN = mynewdata(:,:,:,:,ii); ii=ii+1;
            pnu = mynewdata(:,:,:,:,ii); ii=ii+1;

            rhophoton = mynewdata(:,:,:,:,ii); ii=ii+1;
            rhoeleposi = mynewdata(:,:,:,:,ii); ii=ii+1;
            rhoN = mynewdata(:,:,:,:,ii); ii=ii+1;
            rhonu = mynewdata(:,:,:,:,ii); ii=ii+1;

            sphoton = mynewdata(:,:,:,:,ii); ii=ii+1;
            seleposi = mynewdata(:,:,:,:,ii); ii=ii+1;
            sN = mynewdata(:,:,:,:,ii); ii=ii+1;
            snu = mynewdata(:,:,:,:,ii); ii=ii+1;

            tauael = mynewdata(:,:,:,:,ii); ii=ii+1;
            taus = mynewdata(:,:,:,:,ii); ii=ii+1;
            tautel = mynewdata(:,:,:,:,ii); ii=ii+1;
            tautmu = mynewdata(:,:,:,:,ii); ii=ii+1;
            tauamu = mynewdata(:,:,:,:,ii); ii=ii+1;

            Qmel = mynewdata(:,:,:,:,ii); ii=ii+1;
            Qmmu = mynewdata(:,:,:,:,ii); ii=ii+1;
            Qmtau = mynewdata(:,:,:,:,ii); ii=ii+1;

            qminusel = mynewdata(:,:,:,:,ii); ii=ii+1;
            qminusmu = mynewdata(:,:,:,:,ii); ii=ii+1;

          end


          
        end %%%%%%%%% end if passiter==2 || require2passes==0
        

        
        
        
        %%%%%%%%%%%%%%%%%
        %
        % Continue passiter=1 or 2 || require2passes==0
        %
        %%%%%%%%%%%%%%%%%

        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % First construct 1-D grid of new dependent quantities
        %
        % These are used for internal calculations only, not final output
        %
        % Generally, gridding is (e.g.) utotdiffgrid = 10^(lutotdiffgrid) where lutotdiffgrid is on a uniform grid from lutotdiffmin to lutotdiffmax.
        % But note that lutotdiffmin and lutotdiffmax do not have to be min or max of utotdiff.  That would only be desirable if really want to resolve logarithmically right next to the lowest temperature, which can overresolve it sometimes.
        % This is just like r = R_0 + exp(x1) HARM grid, where R_0~utotoffset and (r-R_0)~10^lutotdiff, but Rin and Rout can be chosen arbitrarily as long as they are larger than R0.  We normally default to choosing Rin=min(10^lutotdiff), but one does not have to.
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % factor by which to enhance resolution of internal values before downgrading to output resolution
        factor=10;
        %factor=1;


        nutotdiff=ntk*factor;
        nptotdiff=nutotdiff;
        nchidiff=nutotdiff;
        nstotdiff=nutotdiff;
        nsspecdiff=nutotdiff;
        
        % anything can be put here, but assume basically same as size of ntk
        nutotdiffout=ntk;
        nptotdiffout=nutotdiffout;
        nchidiffout=nutotdiffout;
        nstotdiffout=nutotdiffout;
        nsspecdiffout=nutotdiffout; % not used really
        
        
        if utotdegencut==2
          % utotdegencut==2 sets up lutotdiff to already vary from 0 to 1 in a log way
          
          [lutotdiffmin lutotdiffmax steplutotdiff lutotdiffgrid] = setupgrid(nutotdiff);
          [lptotdiffmin lptotdiffmax steplptotdiff lptotdiffgrid] = setupgrid(nptotdiff);
          [lchidiffmin lchidiffmax steplchidiff lchidiffgrid] = setupgrid(nchidiff);
          [lstotdiffmin lstotdiffmax steplstotdiff lstotdiffgrid] = setupgrid(nstotdiff);
          [lsspecdiffmin lsspecdiffmax steplsspecdiff lsspecdiffgrid] = setupgrid(nsspecdiff); % used to be sspecgrid alone without diff

          [lutotdiffoutmin lutotdiffoutmax steplutotdiffout lutotdiffoutgrid] = setupgrid(nutotdiffout);
          [lptotdiffoutmin lptotdiffoutmax steplptotdiffout lptotdiffoutgrid] = setupgrid(nptotdiffout);
          [lchidiffoutmin lchidiffoutmax steplchidiffout lchidiffoutgrid] = setupgrid(nchidiffout);
          [lstotdiffoutmin lstotdiffoutmax steplstotdiffout lstotdiffoutgrid] = setupgrid(nstotdiffout);
          [lsspecdiffoutmin lsspecdiffoutmax steplsspecdiffout lsspecdiffoutgrid] = setupgrid(nsspecdiffout); % not used really
          
          % should never need (e.g.) utotdiffgrid from lutotdiffgrid since there is no real meaning to utotdiffgrid!  Only loginterp stuff used it, probably incorrectly so.
          % If really need things as functions of utot again, have to change variable using u(utotdiff).
          % Consider (e.g.) lutotdiff grid as fundamental grid upon which all interpolations are done.
          
        end

        
        
        if utotdegencut==1 || utotdegencut==0

          % Consider (e.g.) lutotdiff grid as fundamental grid upon which all interpolations are done.

          % grid of internal energy (utot) (used to be also grid for pressure (ptot), and \chi=u+p)
          %lutotmin=log10(rhobmin*c*c);
          %lutotmax=log10(rhobmax*c*c);

          % This shift is ok since this defines new grid as will be used and
          % original values of utotdiff/ptotdiff/chidiff can extent further than that grid
          SHIFTU=0.99;
          %  SHIFTU=1.0;

          %%%%%%%%%%%%%%%%%%%%%%%%%%%
          % grid of lutotdiff
          if utotfix==0
            lutotdiffmin=min(min(min(min(lutotdiff))));
            lutotdiffmax=max(max(max(max(lutotdiff))))*SHIFTU;
          end
          if utotfix==1
            % presently utot is wrong after fixup when small due to rhob subraction
            % choose reasonable lower limit
            lutotdiffmin=15.0;
            lutotdiffmax=max(max(max(max(lutotdiff))))*SHIFTU;
          end

          steplutotdiff = (lutotdiffmax-lutotdiffmin)/(nutotdiff-1);
          lutotdiffgrid=lutotdiffmin:steplutotdiff:lutotdiffmax;
          %utotdiffgrid = 10.^lutotdiffgrid;


          %%%%%%%%%%%%%%%%%%%%%%%%%%%
          % grid of lptotdiff
          lptotdiffmin=min(min(min(min(lptotdiff))));
          lptotdiffmax=max(max(max(max(lptotdiff))))*SHIFTU;
          steplptotdiff = (lptotdiffmax-lptotdiffmin)/(nptotdiff-1);
          lptotdiffgrid=lptotdiffmin:steplptotdiff:lptotdiffmax;
          %ptotdiffgrid = 10.^lptotdiffgrid;



          %%%%%%%%%%%%%%%%%%%%%%%%%%%
          % grid of lchidiff=log10(chidiff)=log10(utotdiff + ptotdiff)
          lchidiffmin=min(min(min(min(lchidiff))));
          lchidiffmax=max(max(max(max(lchidiff))))*SHIFTU;
          steplchidiff = (lchidiffmax-lchidiffmin)/(nchidiff-1);
          lchidiffgrid=lchidiffmin:steplchidiff:lchidiffmax;
          %chidiffgrid = 10.^lchidiffgrid;

          %%%%%%%%%%%%%%%%%%%%%%%%%%%
          % grid of lstotdiff
          lstotdiffmin=min(min(min(min(lstotdiff))));
          lstotdiffmax=max(max(max(max(lstotdiff))))*SHIFTU;
          steplstotdiff = (lstotdiffmax-lstotdiffmin)/(nstotdiff-1);
          lstotdiffgrid=lstotdiffmin:steplstotdiff:lstotdiffmax;
          %stotdiffgrid = 10.^lstotdiffgrid;


          %%%%%%%%%%%%%%%%%%%%%%%%%%%
          % grid of lsspecdiff
          % specific entropy as a independent variable
          % Note that could use stot gridding above since if here, stotdiff simply related to sspecdiff, but not necessary since no loss of accuracy here (unlike possible in interpolations farther below)
          lsspecdiffmin=min(min(min(min(lsspecdiff))));
          lsspecdiffmax=max(max(max(max(lsspecdiff))))*SHIFTU;
          steplsspecdiff = (lsspecdiffmax-lsspecdiffmin)/(nsspecdiff-1);
          lsspecdiffgrid=lsspecdiffmin:steplsspecdiff:lsspecdiffmax;
          %sspecdiffgrid = 10.^lsspecdiffgrid;


          
          
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %
          % OUTPUT grid of functions of utotdiff/ptotdiff/chidiff
          %
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %lutotdiffmin=log10(rhobmin*c*c);
          %lutotdiffmax=log10(rhobmax*c*c);
          lutotdiffoutmin=min(lutotdiffoutmin,lutotdiffmin);
          lutotdiffoutmax=max(lutotdiffoutmax,lutotdiffmax);
          steplutotdiffout = (lutotdiffoutmax-lutotdiffoutmin)/(nutotdiffout-1);
          lutotdiffoutgrid=lutotdiffoutmin:steplutotdiffout:lutotdiffoutmax;
          %utotdiffoutgrid = 10.^lutotdiffoutgrid;
          
          %fprintf(fiddebug,'debug: %21.15g %21.15g\n',utotdiffoutgrid(1),10.0.^lutotdiffoutmin);
          
          %if(abs(utotdiffoutgrid(1)-10.0.^lutotdiffoutmin)>1E-13)
          %  fprintf(fiddebug,'Problem with utotdiffoutgrid %d %d: %21.15g %21.15g\n',hiter,titer,utotdiffoutgrid(1),10.0.^lutotdiffoutmin);
          %end

          % OUTPUT grid of functions of ptotdiff
          %lptotdiffmin=log10(rhobmin*c*c);
          %lptotdiffmax=log10(rhobmax*c*c);
          lptotdiffoutmin=min(lptotdiffoutmin,lptotdiffmin);
          lptotdiffoutmax=max(lptotdiffoutmax,lptotdiffmax);
          steplptotdiffout = (lptotdiffoutmax-lptotdiffoutmin)/(nptotdiffout-1);
          lptotdiffoutgrid=lptotdiffoutmin:steplptotdiffout:lptotdiffoutmax;
          %ptotdiffoutgrid = 10.^lptotdiffoutgrid;

          % OUTPUT grid of functions of chidiff
          %lchidiffmin=log10(rhobmin*c*c);
          %lchidiffmax=log10(rhobmax*c*c);
          lchidiffoutmin=min(lchidiffoutmin,lchidiffmin);
          lchidiffoutmax=max(lchidiffoutmax,lchidiffmax);
          steplchidiffout = (lchidiffoutmax-lchidiffoutmin)/(nchidiffout-1);
          lchidiffoutgrid=lchidiffoutmin:steplchidiffout:lchidiffoutmax;
          %chidiffoutgrid = 10.^lchidiffoutgrid;

          % OUTPUT grid of functions of stotdiff
          %lstotdiffmin=log10(rhobmin*c*c);
          %lstotdiffmax=log10(rhobmax*c*c);
          lstotdiffoutmin=min(lstotdiffoutmin,lstotdiffmin);
          lstotdiffoutmax=max(lstotdiffoutmax,lstotdiffmax);
          steplstotdiffout = (lstotdiffoutmax-lstotdiffoutmin)/(nstotdiffout-1);
          lstotdiffoutgrid=lstotdiffoutmin:steplstotdiffout:lstotdiffoutmax;
          %stotdiffoutgrid = 10.^lstotdiffoutgrid;

          % below not really used
          lsspecdiffoutmin=min(lsspecdiffoutmin,lsspecdiffmin);
          lsspecdiffoutmax=max(lsspecdiffoutmax,lsspecdiffmax);
          steplsspecdiffout = (lsspecdiffoutmax-lsspecdiffoutmin)/(nsspecdiffout-1);
          lsspecdiffoutgrid=lsspecdiffoutmin:steplsspecdiffout:lsspecdiffoutmax;
          %sspecdiffoutgrid = 10.^lsspecdiffoutgrid;
          
        
        end %%%%%% end if utotdegencut==1 || utotdegencut==0

        
        
        
        

        
        

        
        if(passiter==2 || require2passes==0)


          fprintf(fiddebug,'Begin Interpolate %d %d %d\n',hiter,titer,ynuiter);

          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %
          % Now interpolate to obtain desired functions
          %
          % interpolate in log10 if makes sense
          %
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

          % clear any old-defined quantities (only really needed if running at
          % command line where old values exist and mismatch in dimensions can occur)
          clear extraofUdiff;
          extraofUdiff=zeros(nrhob,nutotdiff,ntdynorye,nhcm,numextras);

          
          % say X, Y, Xi and get Yi
          for p=1:nrhob
            for q=1:ntdynorye
              for r=1:nhcm

                %p
                %q
                %r

                
                if consolid==1
                  
                  %%%%%%%%%%%%%%%%%%%%
                  %
                  % Ensure can later do interpolation by here enforcing values of independent variable vs. Tk to be unique and monotonically increasing
                  %
                  %%%%%%%%%%%%%%%%%%%%

                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  %
                  % BEGIN OLD CODE
                  %
                  % check for duplicates in utot
                  %[xc,yc] = consolidator(utot(p,:,q,r),[],'count');
                  %myone=ones(length(yc),1);
                  %mydiffutot=sum(abs(abs(yc)-myone))

                  %[xc,yc] = consolidator(ptot(p,:,q,r),[],'count');
                  %myone=ones(length(yc),1);
                  %mydiffptot=sum(abs(abs(yc)-myone))

                  %[xc,yc] = consolidator(chi(p,:,q,r),[],'count');
                  %myone=ones(length(yc),1);
                  %mydiffchi=sum(abs(abs(yc)-myone))

                  %[xc,yc] = consolidator(stot(p,:,q,r),[],'count',.05);
                  %myone=ones(length(yc),1);
                  %mydiffstot=sum(abs(abs(yc)-myone))

                  %[xc,yc] = consolidator(log10(stot(p,:,q,r)),[],'count',.05);
                  %myone=ones(length(yc),1);
                  %mydiffstot=sum(abs(abs(yc)-myone))


                  %fprintf(fiddebug,'%21.15g\n',log10(stot(p,:,q,r)))

                  %
                  % END OLD CODE
                  %
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                  
                  %Consolidator: http://www.mathworks.com/matlabcentral/fileexchange/8354

                  % make sure things interpolating are unique
                  % Here x = lutotdiff,lptotdiff,lchidiff,lstotdiff,lsspecdiff
                  % Here y = ptot, hspec, utot, ptot, stot, sspec, Tk, extra?, etc.
                  
                  % below 4 are used for change of variable
                  [lutotdiffutotx,lutotdiffutoty] = consolidator(lutotdiff(p,:,q,r),log10(utot(p,:,q,r)),'mean',CONTOL2);
                  [lptotdiffptotx,lptotdiffptoty] = consolidator(lptotdiff(p,:,q,r),log10(ptot(p,:,q,r)),'mean',CONTOL2);
                  [lchidiffchix,lchidiffchiy] = consolidator(lchidiff(p,:,q,r),log10(chi(p,:,q,r)),'mean',CONTOL2);
                  [lstotdiffstotx,lstotdiffstoty] = consolidator(lstotdiff(p,:,q,r),log10(stot(p,:,q,r)),'mean',CONTOL2);
                  [lsspecdiffsspecx,lsspecdiffsspecy] = consolidator(lsspecdiff(p,:,q,r),log10(sspec(p,:,q,r)),'mean',CONTOL2);

                  % P(lutotdiff)
                  [lutotdiffptotx,lutotdiffptoty] = consolidator(lutotdiff(p,:,q,r),log10(ptot(p,:,q,r)),'mean',CONTOL);

                  % h(lutotdiff)
                  [lhspecx,lhspecy] = consolidator(lutotdiff(p,:,q,r),log10(hspec(p,:,q,r)),'mean',CONTOL);

                  % cs2(lutotdiff)
                  [lcs2rhoTx,lcs2rhoTy] = consolidator(lutotdiff(p,:,q,r),log10(cs2rhoT(p,:,q,r)),'mean',CONTOL);

                  % U(lptotdiff)
                  [lptotdiffutotx,lptotdiffutoty] = consolidator(lptotdiff(p,:,q,r),log10(utot(p,:,q,r)),'mean',CONTOL);
                  % U(lstotdiff)
                  [lstotdiffutotx,lstotdiffutoty] = consolidator(lstotdiff(p,:,q,r),log10(utot(p,:,q,r)),'mean',CONTOL);

                  % P(lchidiff)
                  [lchidiffptotx,lchidiffptoty] = consolidator(lchidiff(p,:,q,r),log10(ptot(p,:,q,r)),'mean',CONTOL);
                  % Ss(lchidiff)
                  [lchidiffsspecx,lchidiffsspecy] = consolidator(lchidiff(p,:,q,r),log10(sspec(p,:,q,r)),'mean',CONTOL);

                  % P(lsspecdiff)
                  [lsspecdiffptotx,lsspecdiffptoty] = consolidator(lsspecdiff(p,:,q,r),log10(ptot(p,:,q,r)),'mean',CONTOL);
                  % U(lsspecdiff)
                  [lsspecdiffutotx,lsspecdiffutoty] = consolidator(lsspecdiff(p,:,q,r),log10(utot(p,:,q,r)),'mean',CONTOL);
                  % lutotdiff(lsspecdiff)
                  [lsspecdifflutotdiffx,lsspecdifflutotdiffy] = consolidator(lsspecdiff(p,:,q,r),lutotdiff(p,:,q,r),'mean',CONTOL);

                  % stot(lutotdiff)
                  [lutotdiffstotx,lutotdiffstoty] = consolidator(lutotdiff(p,:,q,r),log10(stot(p,:,q,r)),'mean',CONTOL);


                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  %
                  % Temperature is used to establish if really within valid EOS
                  % (there was a valid inversion)
                  %
                  % Using temperature as marker to tell if within original EOS or outside and within
                  % interpolated/extrapolted/filled-in region
                  %
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  % Tk(lutotdiff)
                  [lutotdifftkx,lutotdifftky] = consolidator(lutotdiff(p,:,q,r),log10(tk(p,:,q,r)),'mean',CONTOL);

                  % Tk(lptotdiff)
                  [lptotdifftkx,lptotdifftky] = consolidator(lptotdiff(p,:,q,r),log10(tk(p,:,q,r)),'mean',CONTOL);

                  % Tk(lchitotdiff)
                  [lchidifftkx,lchidifftky] = consolidator(lchidiff(p,:,q,r),log10(tk(p,:,q,r)),'mean',CONTOL);
                  
                  % Tk(lstotdiff)
                  [lstotdifftkx,lstotdifftky] = consolidator(lstotdiff(p,:,q,r),log10(tk(p,:,q,r)),'mean',CONTOL);

                  % EXTRAS
                  % mynewdata(:,:,:,:,ei)
                  % assume extras are all >0 (replace 0 with 1E-20)
                  %p
                  %  q
                  %  r
                  %  size(lutotextrax)
                  %  size(utotdiff)
                  %  size(mynewdata)
                  
                  %setup size of extras array
                  %lutotextrax = zeros(numextras,ntk);
                  %lutotextray = zeros(numextras,ntk);

                  
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  % Setup size of the extra matrix since otherwise Matlab complains
                  % about size
                  [tempx tempy] = consolidator(lutotdiff(p,:,q,r),log10(abs(mynewdata(p,:,q,r,1))+1E-20),'mean',CONTOL);
                  sizex=size(tempx(:,1));
                  sizexx=sizex(1);
                  lutotdiffextrax=zeros(numextras,sizexx);
                  sizey=size(tempy(:,1));
                  sizeyy=sizey(1);
                  lutotdiffextray=zeros(numextras,sizeyy);

                  for ei=extraii:extraii+numextras-1
                    % for extra : comes second for stupid reason that
                    % size(utotdiff(p,:,q,r)) is 1 200
                    [lutotdiffextrax(ei-extraii+1,:) lutotdiffextray(ei-extraii+1,:)] = consolidator(lutotdiff(p,:,q,r),log10(abs(mynewdata(p,:,q,r,ei))+1E-20),'mean',CONTOL);
                  end
                  % now transpose result so that 200 1
                  temp=lutotdiffextrax';
                  lutotdiffextrax=temp;
                  temp=lutotdiffextray';
                  lutotdiffextray=temp;

                  
                  
                  
                else % don't use consolidate method -- doesn't generally work due to accuracy of utot at high density
                  
                  

                  lutotdiffutotx=lutotdiff(p,:,q,r);
                  lutotdiffutoty=log10(utot(p,:,q,r));

                  lptotdiffptotx=lptotdiff(p,:,q,r);
                  lptotdiffptoty=log10(ptot(p,:,q,r));

                  lchidiffchix=lchidiff(p,:,q,r);
                  lchidiffchiy=log10(chi(p,:,q,r));

                  lstotdiffstotx=lstotdiff(p,:,q,r);
                  lstotdiffstoty=log10(stot(p,:,q,r));

                  lsspecdiffsspecx=lsspecdiff(p,:,q,r);
                  lsspecdiffsspecy=log10(sspec(p,:,q,r));
                  
                  
                  % P(lutotdiff)
                  lutotdiffptotx=lutotdiff(p,:,q,r);
                  lutotdiffptoty=log10(ptot(p,:,q,r));

                  % h(lutotdiff)
                  lhspecx=lutotdiff(p,:,q,r);
                  lhspecy=log10(hspec(p,:,q,r));

                  % cs2(lutotdiff)
                  lcs2rhoTx=lutotdiff(p,:,q,r);
                  lcs2rhoTy=log10(cs2rhoT(p,:,q,r));


                  % U(lptotdiff)
                  lptotdiffutotx=lptotdiff(p,:,q,r);
                  lptotdiffutoty=log10(utot(p,:,q,r));
                  % U(lstotdiff)
                  lstotdiffutotx=lstotdiff(p,:,q,r);
                  lstotdiffutoty=log10(utot(p,:,q,r));

                  
                  % P(lchidiff)
                  lchidiffptotx=lchidiff(p,:,q,r);
                  lchidiffptoty=log10(ptot(p,:,q,r));

                  % Ss(lchidiff)
                  lchidiffsspecx=lchidiff(p,:,q,r);
                  lchidiffsspecy=log10(sspec(p,:,q,r));

                  
                  % P(lsspecdiff)
                  lsspecdiffptotx=lsspecdiff(p,:,q,r);
                  lsspecdiffptoty=log10(ptot(p,:,q,r));
                  
                  % U(lsspecdiff)
                  lsspecdiffutotx=lsspecdiff(p,:,q,r);
                  lsspecdiffutoty=log10(utot(p,:,q,r));
                  
                  % lutotdiff(lsspecdiff)
                  lsspecdifflutotdiffx=lsspecdiff(p,:,q,r);
                  lsspecdifflutotdiffy=lutotdiff(p,:,q,r);

                  
                  % stot(lutotdiff)
                  lutotdiffstotx=lutotdiff(p,:,q,r);
                  lutotdiffstoty=log10(stot(p,:,q,r));

                  
                  % Tk(lutotdiff)
                  lutotdifftkx=lutotdiff(p,:,q,r);
                  lutotdifftky=log10(tk(p,:,q,r));

                  % Tk(lptotdiff)
                  lptotdifftkx=lptotdiff(p,:,q,r);
                  lptotdifftky=log10(tk(p,:,q,r));

                  % Tk(lchidiff)
                  lchidifftkx=lchidiff(p,:,q,r);
                  lchidifftky=log10(tk(p,:,q,r));

                  % Tk(lstotdiff)
                  lstotdifftkx=lstotdiff(p,:,q,r);
                  lstotdifftky=log10(tk(p,:,q,r));

                  
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  % Setup size of the extra matrix since otherwise Matlab complains
                  % about size
                  tempx=lutotdiff(p,:,q,r);
                  sizex=size(tempx(1,:));
                  sizexx=sizex(2);
                  lutotdiffextrax=zeros(numextras,sizexx);
                  tempy=log10(abs(mynewdata(p,:,q,r,:))+1E-20);
                  sizey=size(tempy(1,:,1,1,1));
                  sizeyy=sizey(2);
                  lutotdiffextray=zeros(numextras,sizeyy);

                  for ei=extraii:extraii+numextras-1
                    % for extra : comes second for stupid reason that
                    % size(utotdiff(p,:,q,r)) is 1 200
                    lutotdiffextrax(ei-extraii+1,:)=lutotdiff(p,:,q,r);
                    lutotdiffextray(ei-extraii+1,:)=log10(abs(mynewdata(p,:,q,r,ei))+1E-20);
                  end
                  % now transpose result so that 200 1
                  temp=lutotdiffextrax';
                  lutotdiffextrax=temp;
                  temp=lutotdiffextray';
                  lutotdiffextray=temp;

                  
                  
                  
                  
                  
                  
                  
                  
                end %%%%% end "else if" not using consididation method
                
                
                

                
                

                %%%%%%% BEGIN DEBUG
                %if		p==197
                
                %    lsspecptotx'
                %    lsspecutoty'


                % for mmm=1:200 fprintf(fiddebug,'%21.15g\n',sspec(p,mmm,q,r)); end
                %lutotdiffptotx'
                %lutotdiffptoty'

                %end

                %  f(x) = interp(x(t),f(t),x_i)
                %%%%%%% END DEBUG
                


                
                
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %
                % Interpolate F(rho0,?).  Note that if capitalized Udiff, Pdiff, CHIdiff, Sdiff, SSdiff, then really related to lutotdiff, lptotdiff, etc.  That is, really are log or some other arbitrary mapping.
                %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % below 4 are used for change of variable:
                % These 4 are actually stored now to check degen offset method in HARM tables (not actually read-in there, just a test)
                UofUdiff(p,:,q,r) =     10.^(myinterp1(200,lutotdiffutotx, lutotdiffutoty, lutotdiffgrid',interptype));
                PofPdiff(p,:,q,r) =     10.^(myinterp1(201,lptotdiffptotx, lptotdiffptoty, lptotdiffgrid',interptype));
                CHIofCHIdiff(p,:,q,r) = 10.^(myinterp1(201,lchidiffchix, lchidiffchiy, lchidiffgrid',interptype));
                SofSdiff(p,:,q,r) =     10.^(myinterp1(201,lstotdiffstotx, lstotdiffstoty, lstotdiffgrid',interptype));
                SSofSSdiff(p,:,q,r) =   10.^(myinterp1(201,lsspecdiffsspecx, lsspecdiffsspecy, lsspecdiffgrid',interptype));

                
                % BELOW dependent variables are actually differenced (offsetted) versions
                % below is P(rho0,u,H)
                % below directly used as P(rho0,u)
                % dP/du |rho0
                % dP/drho0 | u
                PofUdiff(p,:,q,r) = 10.^(myinterp1(1,lutotdiffptotx, lutotdiffptoty, lutotdiffgrid',interptype));
                HofUdiff(p,:,q,r) = 10.^(myinterp1(2,lhspecx, lhspecy, lutotdiffgrid',interptype));

                % c_s^2 in dimensionless units as function of rhob,U
                cs2ofUdiff(p,:,q,r) = 10.^(myinterp1(2,lcs2rhoTx, lcs2rhoTy, lutotdiffgrid',interptype));

                % below is u(rho0,P,H)
                % below directly used as u(rho0,p)
                UofPdiff(p,:,q,r) = 10.^(myinterp1(3,lptotdiffutotx, lptotdiffutoty, lptotdiffgrid',interptype));
                % below is u(rho0,Sden,H)
                % below directly used as u(rho0,Sden)
                UofSdiff(p,:,q,r) = 10.^(myinterp1(3,lstotdiffutotx, lstotdiffutoty, lstotdiffgrid',interptype));

                % below is P(rho0,\chi,H)
                % below used for P[rho0,chi]
                % 1/ (dchi/dp)|rho0
                % 1/(drho0/dp)|chi
                PofCHIdiff(p,:,q,r) = 10.^(myinterp1(4,lchidiffptotx, lchidiffptoty, lchidiffgrid',interptype));
                % below is specificentropy(rho0,\chi,H)
                % below used for Ss[rho0,chi]
                % 1/ (dchi/dSs)|rho0
                % 1/(drho0/dSs)|chi
                SSofCHIdiff(p,:,q,r) = 10.^(myinterp1(4,lchidiffsspecx, lchidiffsspecy, lchidiffgrid',interptype));
                
                % Below is PofSSdiff(rho0,SS,H)
                % below used for c_s^2 = 1/h dp/drho0|SS
                PofSSdiff(p,:,q,r) = 10.^(myinterp1(5,lsspecdiffptotx, lsspecdiffptoty, lsspecdiffgrid',interptype));
                UofSSdiff(p,:,q,r) = 10.^(myinterp1(6,lsspecdiffutotx, lsspecdiffutoty, lsspecdiffgrid',interptype));
                UdiffofSSdiff(p,:,q,r) = 10.^(myinterp1(6,lsspecdifflutotdiffx, lsspecdifflutotdiffy, lsspecdiffgrid',interptype));

                % Below is SofU(rho0,u,H) : S = Sden
                SofUdiff(p,:,q,r) = 10.^(myinterp1(8,lutotdiffstotx, lutotdiffstoty, lutotdiffgrid',interptype));

                % Below is T(rho0,u)
                tkofUdiff(p,:,q,r) = 10.^(myinterp1(29,lutotdifftkx, lutotdifftky, lutotdiffgrid',interptype));
                % Below is T(rho0,p)
                tkofPdiff(p,:,q,r) = 10.^(myinterp1(30,lptotdifftkx, lptotdifftky, lptotdiffgrid',interptype));
                % Below is T(rho0,chi)
                tkofCHIdiff(p,:,q,r) = 10.^(myinterp1(31,lchidifftkx, lchidifftky, lchidiffgrid',interptype));
                % Below is T(rho0,Sden)
                tkofSdiff(p,:,q,r) = 10.^(myinterp1(32,lstotdifftkx, lstotdifftky, lstotdiffgrid',interptype));

                % Below is T(rho0,u)
                faketkofUdiff(p,:,q,r) = 10.^(interp1(lutotdifftkx, lutotdifftky, lutotdiffgrid',interptypefaketemp));
                % Below is T(rho0,p)
                faketkofPdiff(p,:,q,r) = 10.^(interp1(lptotdifftkx, lptotdifftky, lptotdiffgrid',interptypefaketemp));
                % Below is T(rho0,chi)
                faketkofCHIdiff(p,:,q,r) = 10.^(interp1(lchidifftkx, lchidifftky, lchidiffgrid',interptypefaketemp));
                % Below is T(rho0,Sden)
                faketkofSdiff(p,:,q,r) = 10.^(interp1(lstotdifftkx, lstotdifftky, lstotdiffgrid',interptypefaketemp));

                % Below is extra1ofU(rho0,u,H)
                % mynewdata(:,:,:,:,ei)
                % so far all these things are positive definite, so log interpolate ok
                for ei=1:numextras
                  % ei comes first for issues of size()
                  extraofUdiff(p,:,q,r,ei) = 10.^(myinterp1(7,lutotdiffextrax(:,ei), lutotdiffextray(:,ei), lutotdiffgrid',interptype));
                end

                
                
              end
            end
          end


          
          

          fprintf(fiddebug,'End Interpolate %d %d %d\n',hiter,titer,ynuiter);


          %DEBUG
          %HofUdiff(1,14,1,1)
          %HofUdiff(1,:,1,1)
          %DEBUG



          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %
          % rhobgrid, utotgrid, ptotgrid, chigrid, and sspecgrid need to be turned into full 4-D quantity for cleanvar
          % these are super-sampled versions so far
          %
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
          for p=1:nutotdiff
            for q=1:ntdynorye
              for r=1:nhcm
                % using 1 instead of p because of size of p is different
                rhobgrid4d(:,p,q,r)    = rhob(:,1,q,r);
                rhobcsqgrid4d(:,p,q,r) = rhobcsq(:,1,q,r);
              end
            end
          end
          for p=1:nrhob
            for q=1:ntdynorye
              for r=1:nhcm
                % diff values
                lutotdiffgrid4d(p,:,q,r)  = lutotdiffgrid(:);
                lptotdiffgrid4d(p,:,q,r)  = lptotdiffgrid(:);
                lchidiffgrid4d(p,:,q,r)   = lchidiffgrid(:);
                lstotdiffgrid4d(p,:,q,r)  = lstotdiffgrid(:);
                lsspecdiffgrid4d(p,:,q,r) = lsspecdiffgrid(:);

                % actual values (just different label)
                utotgrid4d(p,:,q,r)  = UofUdiff(p,:,q,r);
                ptotgrid4d(p,:,q,r)  = PofPdiff(p,:,q,r);
                chigrid4d(p,:,q,r)   = CHIofCHIdiff(p,:,q,r);
                stotgrid4d(p,:,q,r)  = SofSdiff(p,:,q,r);
                sspecgrid4d(p,:,q,r)  = SSofSSdiff(p,:,q,r);
              end
            end
          end




          
          
          

          
          if doclean
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % Clean/fixup F(rho0,?)
            %
            % Adjust quantities to be physical when interpolated range is beyond
            % existing range
            %
            % cleanvar presumes input independent variables are real U,P,CHI, not differences
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            % say X, Y, Xi and get Yi


            %%%%%%%%%%% Clean F(rho0,?)
            
            
            if usecleanvar==1
              % GODMARK: something is apparently wrong with how I'm accessing memory in cleanvar
              
              % clean change of variables
              UofUdiff = cleanvar(200, UofUdiff, rhobcsqgrid4d, utotgrid4d);
              PofPdiff = cleanvar(201, PofPdiff, rhobcsqgrid4d, ptotgrid4d);
              CHIofCHIdiff = cleanvar(202, CHIofCHIdiff, rhobcsqgrid4d, chigrid4d);
              SofSdiff = cleanvar(203, SofSdiff, rhobcsqgrid4d, stotgrid4d);
              SSofSSdiff = cleanvar(204, SSofSSdiff, rhobcsqgrid4d, sspecgrid4d);

              PofUdiff = cleanvar(1, PofUdiff, rhobgrid4d, utotgrid4d);
              HofUdiff = cleanvar(2, HofUdiff, rhobcsqgrid4d, utotgrid4d); % specific enthalpy is dimensionless

              cs2ofUdiff = cleanvar(40, cs2ofUdiff, rhobcsqgrid4d, utotgrid4d); % here cs2 is dimensionless

              UofPdiff = cleanvar(3, UofPdiff, rhobgrid4d, ptotgrid4d);
              UofSdiff = cleanvar(65, UofSdiff, rhobgrid4d, stotgrid4d);

              PofCHIdiff = cleanvar(4, PofCHIdiff, rhobgrid4d, chigrid4d);
              SSofCHIdiff = cleanvar(400, SSofCHIdiff, rhobgrid4d, chigrid4d);

              PofSSdiff = cleanvar(5, PofSSdiff, rhobcsqgrid4d, sspecgrid4d); % here SS is specific entropy and uses rhobcsq
              UofSSdiff = cleanvar(6, UofSSdiff, rhobcsqgrid4d, sspecgrid4d); % here SS is specific entropy and uses rhobcsq
              UdiffofSSdiff = cleanvar(6, UdiffofSSdiff, rhobcsqgrid4d, sspecgrid4d); % here SS is specific entropy and uses rhobcsq
              
              SofUdiff = cleanvar(8, SofUdiff, rhobcsqgrid4d, utotgrid4d); % here S is entropy density and uses rhobcsq

              tkofUdiff = cleanvar(29, tkofUdiff, rhobgrid4d, utotgrid4d);
              tkofPdiff = cleanvar(30, tkofPdiff, rhobgrid4d, ptotgrid4d);
              tkofCHIdiff = cleanvar(31, tkofCHIdiff, rhobgrid4d, chigrid4d);
              tkofSdiff = cleanvar(30, tkofSdiff, rhobgrid4d, stotgrid4d);
              
              % DO NOT clean faketk's

              %size(tkofU)
              %size(rhobgrid4d)
              %size(lutotdiffgrid4d)


              
              % GODMARK: Something wrong with how using memory?
              % extras:
              % mynewdata(:,:,:,:,ei)
              %extraofU(ei,:,:,:,:) = cleanvar(7, extraofU(ei,:,:,:,:),rhobgrid4d, utotgrid4d);
              % clean all at once
              for ei=1:numextras
                extraofUdiff(:,:,:,:,ei) = cleanvar(7, extraofUdiff(:,:,:,:,ei), rhobgrid4d(:,:,:,:), utotgrid4d(:,:,:,:));
              end
              
              
              

              
            else %%%%%% else if not using cleanvar() function
              
              

              
              if 1 % CURRENTLY USED METHOD:
                
                % go ahead and check finiteness of log10(var) for var's that should be >0.0 and never <=0.0
                % OK to get log10(0) warnings from the below -- just means extrapolation went a bit nuts -- won't use that data since faketkof? will save us.
                UofUdiff(~isfinite(log10(UofUdiff)))=OUTBOUNDSVALUE;
                PofPdiff(~isfinite(log10(PofPdiff)))=OUTBOUNDSVALUE;
                CHIofCHIdiff(~isfinite(log10(CHIofCHIdiff)))=OUTBOUNDSVALUE;
                SofSdiff(~isfinite(log10(SofSdiff)))=OUTBOUNDSVALUE;
                SSofSSdiff(~isfinite(log10(SSofSSdiff)))=OUTBOUNDSVALUE;

                PofUdiff(~isfinite(log10(PofUdiff)))=OUTBOUNDSVALUE;
                HofUdiff(~isfinite(log10(HofUdiff)))=OUTBOUNDSVALUE;

                cs2ofUdiff(~isfinite(log10(cs2ofUdiff)))=OUTBOUNDSVALUE;

                UofPdiff(~isfinite(log10(UofPdiff)))=OUTBOUNDSVALUE;
                UofSdiff(~isfinite(log10(UofSdiff)))=OUTBOUNDSVALUE;

                PofCHIdiff(~isfinite(log10(PofCHIdiff)))=OUTBOUNDSVALUE;
                SSofCHIdiff(~isfinite(log10(SSofCHIdiff)))=OUTBOUNDSVALUE;

                PofSSdiff(~isfinite(log10(PofSSdiff)))=OUTBOUNDSVALUE;
                UofSSdiff(~isfinite(log10(UofSSdiff)))=OUTBOUNDSVALUE;
                UdiffofSSdiff(~isfinite(log10(UdiffofSSdiff)))=OUTBOUNDSVALUE;

                SofUdiff(~isfinite(log10(SofUdiff)))=OUTBOUNDSVALUE;

                tkofUdiff(~isfinite(log10(tkofUdiff)))=OUTBOUNDSVALUE;
                tkofPdiff(~isfinite(log10(tkofPdiff)))=OUTBOUNDSVALUE;
                tkofCHIdiff(~isfinite(log10(tkofCHIdiff)))=OUTBOUNDSVALUE;
                tkofSdiff(~isfinite(log10(tkofSdiff)))=OUTBOUNDSVALUE;

                % DO NOT clean faketk's

                extraofUdiff(~isfinite(log10(extraofUdiff)))=OUTBOUNDSVALUE;
                
              
              
              else % not used currently
                
                for p=1:nrhob
                  for q=1:nutotdiff
                    for r=1:ntdynorye
                      for s=1:nhcm
                        %              p
                        %              q
                        %              r
                        %              s
                        % clean change of variables
                        if(~isfinite(UofUdiff(p,q,r,s)))
                          UofUdiff(p,q,r,s)=OUTBOUNDSVALUE;
                        end
                        if(~isfinite(PofPdiff(p,q,r,s)))
                          PofPdiff(p,q,r,s)=OUTBOUNDSVALUE;
                        end
                        if(~isfinite(CHIofCHIdiff(p,q,r,s)))
                          CHIofCHIdiff(p,q,r,s)=OUTBOUNDSVALUE;
                        end
                        if(~isfinite(SofSdiff(p,q,r,s)))
                          SofSdiff(p,q,r,s)=OUTBOUNDSVALUE;
                        end
                        if(~isfinite(SSofSSdiff(p,q,r,s)))
                          SSofSSdiff(p,q,r,s)=OUTBOUNDSVALUE;
                        end

                        if(~isfinite(PofUdiff(p,q,r,s)))
                          PofUdiff(p,q,r,s)=OUTBOUNDSVALUE;
                        end
                        if(~isfinite(HofUdiff(p,q,r,s)))
                          HofUdiff(p,q,r,s)=OUTBOUNDSVALUE;
                        end

                        if(~isfinite(cs2ofUdiff(p,q,r,s)))
                          cs2ofUdiff(p,q,r,s)=OUTBOUNDSVALUE;
                        end

                        if(~isfinite(UofPdiff(p,q,r,s)))
                          UofPdiff(p,q,r,s)=OUTBOUNDSVALUE;
                        end
                        if(~isfinite(UofSdiff(p,q,r,s)))
                          UofSdiff(p,q,r,s)=OUTBOUNDSVALUE;
                        end

                        if(~isfinite(PofCHIdiff(p,q,r,s)))
                          PofCHIdiff(p,q,r,s)=OUTBOUNDSVALUE;
                        end
                        if(~isfinite(SSofCHIdiff(p,q,r,s)))
                          SSofCHIdiff(p,q,r,s)=OUTBOUNDSVALUE;
                        end

                        if(~isfinite(PofSSdiff(p,q,r,s)))
                          PofSSdiff(p,q,r,s)=OUTBOUNDSVALUE;
                        end
                        if(~isfinite(UofSSdiff(p,q,r,s)))
                          UofSSdiff(p,q,r,s)=OUTBOUNDSVALUE;
                        end
                        if(~isfinite(UdiffofSSdiff(p,q,r,s)))
                          UdiffofSSdiff(p,q,r,s)=OUTBOUNDSVALUE;
                        end
                        
                        if(~isfinite(SofUdiff(p,q,r,s)))
                          SofUdiff(p,q,r,s)=OUTBOUNDSVALUE;
                        end

                        if(~isfinite(tkofUdiff(p,q,r,s)))
                          tkofUdiff(p,q,r,s)=OUTBOUNDSVALUE;
                        end
                        if(~isfinite(tkofPdiff(p,q,r,s)))
                          tkofPdiff(p,q,r,s)=OUTBOUNDSVALUE;
                        end
                        if(~isfinite(tkofCHIdiff(p,q,r,s)))
                          tkofCHIdiff(p,q,r,s)=OUTBOUNDSVALUE;
                        end
                        if(~isfinite(tkofSdiff(p,q,r,s)))
                          tkofSdiff(p,q,r,s)=OUTBOUNDSVALUE;
                        end

                        % DO NOT clean faketk's
                        
                        for ei=1:numextras
                          if(~isfinite(extraofUdiff(p,q,r,s,ei)))
                            extraofUdiff(p,q,r,s,ei)=OUTBOUNDSVALUE;
                          end
                        end
                        
                      end
                    end
                  end
                end
                
              end

            end


            %totalnan=sum(sum(sum(sum(isnan(UofS)))))
          end

          %HofUdiff(1,14,1,1)
          %HofUdiff(1,:,1,1)



          
          
          
          



          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %
          % Now compute derivatives using interpolated functions
          %
          % when using gradient, crucial that list of independent grid values is of size that's larger first if 1-D grid
          %
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

          % first get rho(i) and S(i)
          % assumes rhob is same for all other indicies (true)
          %rhobcsq=c.*c.*rhob;
          rhocsqi = rhobcsq(:,1,1,1);
          lrhocsqi = log10(rhocsqi);
          %drhocsqi = gradient(rhocsqi);

          %rhoi = rhob(:,1,1,1);
          %lrhoi = log10(rhoi);
          %drhoi = gradient(rhoi);

          % entropy density (stot) as independent quantity is setup to be same as ugrid
          % stoti = utotgrid;
          % dstoti=gradient(stoti);

          % specific entropy as independent quantity
          % sspeci = utotgrid./rhocsqi;
          %sspeci = sspecgrid';
          %lsspeci = lsspecgrid';
          %dsspeci=gradient(sspeci);

          
          %%%%%%%%%%%%%%%%%%%%%%
          %
          % Also setup differenced version of independent variables
          %
          %%%%%%%%%%%%%%%%%%%%%%%
          %utotdiffi = utotdiffgrid';
          lutotdiffi = lutotdiffgrid';
          %dutoti=gradient(utoti);

          %ptotdiffi = ptotdiffgrid';
          lptotdiffi = lptotdiffgrid';
          %dptoti=gradient(ptoti);

          %chidiffi = chidiffgrid';
          lchidiffi = lchidiffgrid';
          %dchii=gradient(chii);

          % note that entropy density not currently used as independent variable for derivatives, only for function U(rho0,S)
          %stotdiffi = stotdiffgrid';
          lstotdiffi = lstotdiffgrid';
          %dstoti=gradient(stoti);

          % note that entropy density not currently used as independent variable for derivatives, only for function U(rho0,S)
          %sspecdiffi = sspecdiffgrid';
          lsspecdiffi = lsspecdiffgrid';
          %dsspeci=gradient(sspeci);

          
          % say X, Y, Xi and get Yi
          for q=1:ntdynorye
            for r=1:nhcm
              % NOTE
              % The first output FX is always the gradient along the 2nd dimension of F, going across columns. The second output FY is always the gradient along the 1st dimension of F, going across rows. For the third output FZ and the outputs that follow, the Nth output is the gradient along the Nth dimension of F. 

              % if PofS(rho0,S,H) then derivative is dS, drho, dH
              % gradient gives df/di along each direction
              % so if want dP/dS = dP/di / dS/di

              % note order of dsspeci and drhocsqi is such that corresponds to output
              % derivatives that are flipped for 1 and 2 (row/col) as described above

              % when vector is given for second+ terms in gradient, assumed to be
              % positions instead of spacing as when scalar input

              % Need UofUdiff, etc. for change of variables since independent variables need
              % to be simple diff version and can't use wildly varying actual utot/etc.
              
              if smoothbeforederivativestep2==1
                
                roughUofUdiff(:,:)=UofUdiff(:,:,q,r);
                roughPofPdiff(:,:)=PofPdiff(:,:,q,r);
                roughCHIofCHIdiff(:,:)=CHIofCHIdiff(:,:,q,r);
                roughSofSdiff(:,:)=SofSdiff(:,:,q,r);
                roughSSofSSdiff(:,:)=SSofSSdiff(:,:,q,r);
                
                roughPofSSdiff(:,:)=PofSSdiff(:,:,q,r);
                roughPofUdiff(:,:)=PofUdiff(:,:,q,r);

                roughPofCHIdiff(:,:)=PofCHIdiff(:,:,q,r);
                roughSSofCHIdiff(:,:)=SSofCHIdiff(:,:,q,r);

                roughSofUdiff(:,:)=SofUdiff(:,:,q,r);
                
                % moving average uses 1,1 for size of averaging region
                myUofUdiff(:,:)=moving_average2(roughUofUdiff(:,:),1,1);
                myPofPdiff(:,:)=moving_average2(roughPofPdiff(:,:),1,1);
                myCHIofCHIdiff(:,:)=moving_average2(roughCHIofCHIdiff(:,:),1,1);
                mySofSdiff(:,:)=moving_average2(roughSofSdiff(:,:),1,1);
                mySSofSSdiff(:,:)=moving_average2(roughSSofSSdiff(:,:),1,1);

                myPofSSdiff(:,:)=moving_average2(roughPofSSdiff(:,:),1,1);
                myPofUdiff(:,:)=moving_average2(roughPofUdiff(:,:),1,1);

                myPofCHIdiff(:,:)=moving_average2(roughPofCHIdiff(:,:),1,1);
                mySSofCHIdiff(:,:)=moving_average2(roughSSofCHIdiff(:,:),1,1);

                mySofUdiff(:,:)=moving_average2(roughSofUdiff(:,:),1,1);
              else
                myUofUdiff(:,:)=UofUdiff(:,:,q,r);
                myPofPdiff(:,:)=PofPdiff(:,:,q,r);
                myCHIofCHIdiff(:,:)=CHIofCHIdiff(:,:,q,r);
                mySofSdiff(:,:)=SofSdiff(:,:,q,r);
                mySSofSSdiff(:,:)=SSofSSdiff(:,:,q,r);

                myPofSSdiff(:,:)=PofSSdiff(:,:,q,r);
                myPofUdiff(:,:)=PofUdiff(:,:,q,r);

                myPofCHIdiff(:,:)=PofCHIdiff(:,:,q,r);
                mySSofCHIdiff(:,:)=SSofCHIdiff(:,:,q,r);

                mySofUdiff(:,:)=SofUdiff(:,:,q,r);
              end
              

              % below are change of variable
              [dUofUdiffdludiff(:,:,q,r), dUofUdiffdlrho0(:,:,q,r)] = gradient(myUofUdiff(:,:),lutotdiffi,lrhocsqi);
              [dPofPdiffdlpdiff(:,:,q,r), dPofPdiffdlrho0(:,:,q,r)] = gradient(myPofPdiff(:,:),lptotdiffi,lrhocsqi);
              [dCHIofCHIdiffdlchidiff(:,:,q,r), dCHIofCHIdiffdlrho0(:,:,q,r)] = gradient(myCHIofCHIdiff(:,:),lchidiffi,lrhocsqi);
              [dSofSdiffdlstotdiff(:,:,q,r), dSofSdiffdlrho0(:,:,q,r)] = gradient(mySofSdiff(:,:),lstotdiffi,lrhocsqi);
              [dSSofSSdiffdlsspecdiff(:,:,q,r), dSSofSSdiffdlrho0(:,:,q,r)] = gradient(mySSofSSdiff(:,:),lsspecdiffi,lrhocsqi);

              % functions:
              [dPofSSdiffdlsspecdiff(:,:,q,r), dPofSSdiffdlrho0(:,:,q,r)] = gradient(myPofSSdiff(:,:),lsspecdiffi,lrhocsqi);
              
              % PofU(rho0,U,H) then derivative is dU, drho0, dH
              [dPofUdiffdludiff(:,:,q,r), dPofUdiffdlrho0(:,:,q,r)] = gradient(myPofUdiff(:,:),lutotdiffi,lrhocsqi);
              
              % PofCHI(rho0, chi ,H) then derivative is dchi, drho0, dH
              [dPofCHIdiffdlchidiff(:,:,q,r), dPofCHIdiffdlrho0(:,:,q,r)] = gradient(myPofCHIdiff(:,:),lchidiffi,lrhocsqi);

              % SSofCHI(rho0, chi ,H) then derivative is dchi, drho0, dH
              [dSSofCHIdiffdlchidiff(:,:,q,r), dSSofCHIdiffdlrho0(:,:,q,r)] = gradient(mySSofCHIdiff(:,:),lchidiffi,lrhocsqi);

              % SofU(rho0,U,H) then derivative is dU, drho0, dH
              [dSofUdiffdludiff(:,:,q,r), dSofUdiffdlrho0(:,:,q,r)] = gradient(mySofUdiff(:,:), lutotdiffi,lrhocsqi);
              

              
              
              % Now obtain true derivative for those derivatives that used diff instead of explicit variable
              % versions as independent variables
              % Note that we assume no mixing between rho0 and T, so change of variable just relies on one derivative and not the other
              % The below change of variables seem to cause division by 0, so trap below those things making the value out of bounds
              dPofSSdiffdsspec(:,:,q,r) = dPofSSdiffdlsspecdiff(:,:,q,r)./dSSofSSdiffdlsspecdiff(:,:,q,r);
              dPofUdiffdu(:,:,q,r) = dPofUdiffdludiff(:,:,q,r)./dUofUdiffdludiff(:,:,q,r);
              dPofCHIdiffdchi(:,:,q,r) = dPofCHIdiffdlchidiff(:,:,q,r)./dCHIofCHIdiffdlchidiff(:,:,q,r);
              dSSofCHIdiffdchi(:,:,q,r) = dSSofCHIdiffdlchidiff(:,:,q,r)./dCHIofCHIdiffdlchidiff(:,:,q,r);
              dSofUdiffdu(:,:,q,r) = dSofUdiffdludiff(:,:,q,r)./dUofUdiffdludiff(:,:,q,r);

              % since rho is constant vs. SSdiff, Udiff, CHIdiff, then no need to compute derivative using gradient since know it analytically
              % drho0/dlrho0 = rho0
              drho0ofSSdiffdlrho0(:,:,q,r) = rhobcsqgrid4d(:,:,q,r); % must use super-sampled version
              drho0ofUdiffdlrho0(:,:,q,r) = rhobcsqgrid4d(:,:,q,r);
              drho0ofCHIdiffdlrho0(:,:,q,r) = rhobcsqgrid4d(:,:,q,r);
              
              % now fix derivatives on density
              dPofSSdiffdrho0(:,:,q,r) = dPofSSdiffdlrho0(:,:,q,r)./drho0ofSSdiffdlrho0(:,:,q,r);
              dPofUdiffdrho0(:,:,q,r) = dPofUdiffdlrho0(:,:,q,r)./drho0ofUdiffdlrho0(:,:,q,r);
              dPofCHIdiffdrho0(:,:,q,r) = dPofCHIdiffdlrho0(:,:,q,r)./drho0ofCHIdiffdlrho0(:,:,q,r);
              dSSofCHIdiffdrho0(:,:,q,r) = dSSofCHIdiffdlrho0(:,:,q,r)./drho0ofCHIdiffdlrho0(:,:,q,r);
              dSofUdiffdrho0(:,:,q,r) = dSofUdiffdlrho0(:,:,q,r)./drho0ofUdiffdlrho0(:,:,q,r);

              % check if out of bounds
              dPofSSdiffdsspec(~isfinite(dPofSSdiffdsspec)) = OUTBOUNDSVALUE;
              dPofUdiffdu(~isfinite(dPofUdiffdu)) = OUTBOUNDSVALUE;
              dPofCHIdiffdchi(~isfinite(dPofCHIdiffdchi)) = OUTBOUNDSVALUE;
              dSSofCHIdiffdchi(~isfinite(dSSofCHIdiffdchi)) = OUTBOUNDSVALUE;
              dSofUdiffdu(~isfinite(dSofUdiffdu)) = OUTBOUNDSVALUE;

              dPofSSdiffdrho0(~isfinite(dPofSSdiffdrho0)) = OUTBOUNDSVALUE;
              dPofUdiffdrho0(~isfinite(dPofUdiffdrho0)) = OUTBOUNDSVALUE;
              dPofCHIdiffdrho0(~isfinite(dPofCHIdiffdrho0)) = OUTBOUNDSVALUE;
              dSSofCHIdiffdrho0(~isfinite(dSSofCHIdiffdrho0)) = OUTBOUNDSVALUE;
              dSofUdiffdrho0(~isfinite(dSofUdiffdrho0)) = OUTBOUNDSVALUE;

              
            end
          end

          % DEBUG CODE:
          %gradPofCHIdiff=gradient(ptot(1,:,1,1),chidiff(1,:,1,1));
          %figure; plot(chidiff(1,:,1,1),gradPofCHIdiff)





          %%%%%%%%%%%%%%%%%%
          %
          % Done with temporary variables used for gradient operator
          %
          %%%%%%%%%%%%%%%%%%

          %clear rhocsqi drhocsqi rhoi drhoi sspeci dsspeci utoti dutoti ptoti
          %dptoti chii dchii
          clear rhocsqi
          %rhoi
          %sspeci
          %utoti ptoti chii




          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %
          % Some derivatives need to be have a full change of variable on the grid itself
          %
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






          % say X, Y, Xi and get Yi
          for p=1:nrhob
            for q=1:ntdynorye
              for r=1:nhcm
                
                
                % HACK
                %UdiffofS(p,:,q,r) = UdiffofS(p,:,q,r) + OUTBOUNDSVALUE;
                
                
                % BEGIN DEBUG:
                %hiter
                %truenhcm

                %ynuiter
                %truentdynorynu

                %titer
                %truentdynorye

                %x(1,:)=log10(UdiffofS(p,:,q,r));
                %sumit = sum(log10(UdiffofS(p,:,q,r)));
                %if ~isfinite(sumit)

                %  nsspec
                %  lsspecmin
                %  lsspecmax
                %  steplsspec
                %  lsspecgrid
                %  sspecgrid

                %p
                %  q
                %  r
                  
                  
                %  sspecgrid4d(p,:,q,r)
                
                %  log10(UdiffofS(p,:,q,r))
                  
                %  dPofSdrho0(p,:,q,r)
                %end
                % END DEBUG:
                
                % below is presently dP/drho(rho0,S)|S and need dP/drho(rho0,u)|SS so change SS->u
                % dPofSdrho0(:,:,q,r)

                %[lutotdpx,lutotdpy] = consolidator(log10(UofS(p,:,q,r)),dPofSSdiffdrho0(p,:,q,r),'mean',CONTOL);
                [lutotdiffdpx,lutotdiffdpy] = consolidator(UdiffofSSdiff(p,:,q,r),dPofSSdiffdrho0(p,:,q,r),'mean',CONTOL);

                % Below is dP/drho0(U,rho0,H)|SS
                % below used for c_s^2 = 1/h dp/drho0|SS
                
                % Using UdiffofSSdiff since interpolation needs to use consistent quantity : diff version
                % Below changes from SSdiff to Udiff dependence
                % This may cause problems since may be more accurate (and avoid chopping off grid sections) by moving from T to U direction for sound speed (i.e. avoid entropy as in other sound speed method)
                dPofUdiffdrho0cSS(p,:,q,r) = myinterp1(9,lutotdiffdpx, lutotdiffdpy, lutotdiffgrid',interptype);

              end
            end
          end




          if doclean
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % Need to clean dPofUdrho0cS since reinterpolated from F(rho0,S) back to F(rho0,U)
            %
            % Notice that this uses rhobcsq instead of rhob because derivatives are
            % assumed to in dimensionless ratios
            %
            % Again, cleaning assumes input independents are actual utot/ptot/chi
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if(usecleanvar)
              dPofUdiffdrho0cSS = cleanvar(9, dPofUdiffdrho0cSS, rhobcsqgrid4d, utotgrid4d);
            else
              % this derivative should be checked for being positive definite since comes into sound speed that has otherwise positive definite contributions to a final positive definite quantity
              dPofUdiffdrho0cSS(~isfinite(log10(dPofUdiffdrho0cSS)))=OUTBOUNDSVALUE;
            end

          end




          
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %
          % Now compute things that don't involve taking more derivatives
          %
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

          if preinterpsoundspeed==0 || preinterpsoundspeed==2 

            % c_s^2 in (cm/s)^2
            % dPofUdrho0cS already has 1/c^2 in it yso below is per unit c^2
            % GODMARK: Note that explicit energy/baryon comes in here, unlike wanted
            cs2ofUdiffpost = (1.0./HofUdiff).*(dPofUdiffdrho0cSS);

          end

          if preinterpsoundspeed==0
            cs2ofUdiff=cs2ofUdiffpost;
          end

          if preinterpsoundspeed==2
            % Use cs2post to fix cs2
            for p=1:nhcm
              for o=1:ntdynorye
                for n=1:nutotdiff % supersampled so far
                  for m=1:nrhob
                    
                    if cs2ofUdiff(m,n,o,p)>2.0*cs2ofUdiffpost(m,n,o,p) && cs2ofUdiffpost(m,n,o,p)>=0.0
                      cs2ofUdiff(m,n,o,p) = cs2ofUdiffpost(m,n,o,p);
                    end
                    if cs2ofUdiff(m,n,o,p)<=0.0 && cs2ofUdiffpost(m,n,o,p)>0.0
                      if cs2ofUdiffpost(m,n,o,p)<1.0
                        cs2ofUdiff(m,n,o,p) = cs2ofUdiffpost(m,n,o,p);
                      else
                        cs2ofUdiff(m,n,o,p) = OUTBOUNDSVALUE;
                      end
                    end

                  end
                end
              end
            end



          end


          %%%%%%%%%%%%%%%%%%%%%%%%%%
          % do in any case
          %%%%%%%%%%%%%%%%%%%%%%%%%%
          cs2ofUdiffcgs = cs2ofUdiff.*c.*c;

          
          
          
          
          

          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %
          %		list of things to output
          %
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %
          % PofUdiff(rho0,U)
          % UofPdiff(rho0,P)

          % dPofUdiffdrho0(rho0,U,H)
          % dPofUdiffdu(rho0,U,H)

          % cs2ofUdiff(rho0,U,H)

          % sofUdiff(rho0,U,H)
          % dSofUdiffdrho0(rho0,U,H)
          % dSofUdiffdu(rho0,U,H)

          % PofCHIdiff(rho0,CHI,H)
          % dPofCHIdiffdrho0(rho0,CHI,H)
          % dPofCHIdiffdchi(rho0,CHI,H)


          % extraofUdiff(rho0,u,H,ei)

          % tkofUdiff(rho0,u,H)


          
          
          
          
          
          

          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %
          % Downsample for writing file
          %
          % downsample: ugrid since all F(:,Udiff,:,:) are of larger-than-desired size
          %
          % for quantities that have large dynamic range use log interpolation
          % unless quantities can naturally be 0 or negative such as derivatives or
          % things from derivatives (e.g. cs2ofUdiff)
          %
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          fprintf(fiddebug,'Begin Downsample %d %d %d\n',hiter,titer,ynuiter);

          for p=1:nrhob
            for q=1:ntdynorye
              for r=1:nhcm

                % for degen checks in HARM
                UofUdiffout(p,:,q,r)          = 10.^(myinterp1(1,lutotdiffgrid, log10(UofUdiff(p,:,q,r)), lutotdiffoutgrid',interptype));
                PofPdiffout(p,:,q,r)          = 10.^(myinterp1(2,lptotdiffgrid, log10(PofPdiff(p,:,q,r)), lptotdiffoutgrid',interptype));
                CHIofCHIdiffout(p,:,q,r)      = 10.^(myinterp1(3,lchidiffgrid, log10(CHIofCHIdiff(p,:,q,r)), lchidiffoutgrid',interptype));
                SofSdiffout(p,:,q,r)      = 10.^(myinterp1(3,lstotdiffgrid, log10(SofSdiff(p,:,q,r)), lstotdiffoutgrid',interptype));
                %SSofSSdiffout(p,:,q,r)      = 10.^(myinterp1(3,lsspecdiffgrid, log10(SSofSSdiff(p,:,q,r)), lsspecdiffoutgrid',interptype)); % not needed in output
                
                
                %            [lutotdpx,lutotdpy] = consolidator(log10(UofS(p,:,q,r)),log10(dPofSdrho0(p,:,q,r)),'mean',CONTOL);

                PofUdiffout(p,:,q,r)          = 10.^(myinterp1(1,lutotdiffgrid, log10(PofUdiff(p,:,q,r)), lutotdiffoutgrid',interptype));
                HofUdiffout(p,:,q,r)          = 10.^(myinterp1(2,lutotdiffgrid, log10(HofUdiff(p,:,q,r)), lutotdiffoutgrid',interptype));

                UofPdiffout(p,:,q,r)          = 10.^(myinterp1(3,lptotdiffgrid, log10(UofPdiff(p,:,q,r)), lptotdiffoutgrid',interptype));
                UofSdiffout(p,:,q,r)          = 10.^(myinterp1(3,lstotdiffgrid, log10(UofSdiff(p,:,q,r)), lstotdiffoutgrid',interptype));

                % dPofUdiffdrho0out(p,:,q,r)   = 10.^(myinterp1(20,lutotdiffgrid, log10(dPofUdiffdrho0(p,:,q,r)), lutotdiffoutgrid',interptype));
                % dPofUdiffduout(p,:,q,r)      = 10.^(myinterp1(21,lutotdiffgrid, log10(dPofUdiffdu(p,:,q,r)), lutotdiffoutgrid',interptype));
                dPofUdiffdrho0out(p,:,q,r)    = myinterp1(20,lutotdiffgrid, dPofUdiffdrho0(p,:,q,r), lutotdiffoutgrid',interptype);
                dPofUdiffduout(p,:,q,r)       = myinterp1(21,lutotdiffgrid, dPofUdiffdu(p,:,q,r), lutotdiffoutgrid',interptype);

                cs2ofUdiffout(p,:,q,r)           = myinterp1(22,lutotdiffgrid, cs2ofUdiff(p,:,q,r), lutotdiffoutgrid',interptype2);
                cs2ofUdiffcgsout(p,:,q,r)        = myinterp1(23,lutotdiffgrid, cs2ofUdiffcgs(p,:,q,r), lutotdiffoutgrid',interptype2);

                SofUdiffout(p,:,q,r)          = 10.^(myinterp1(8,lutotdiffgrid, log10(SofUdiff(p,:,q,r)), lutotdiffoutgrid',interptype));
                % dSofUdiffdrho0out(p,:,q,r)   = 10.^(myinterp1(24,lutotdiffgrid, log10(dSofUdiffdrho0(p,:,q,r)), lutotdiffoutgrid',interptype));
                % dSofUdiffduout(p,:,q,r)      = 10.^(myinterp1(25,lutotdiffgrid, log10(dSofUdiffdu(p,:,q,r)), lutotdiffoutgrid',interptype));
                dSofUdiffdrho0out(p,:,q,r)    = myinterp1(24,lutotdiffgrid, dSofUdiffdrho0(p,:,q,r), lutotdiffoutgrid',interptype);
                dSofUdiffduout(p,:,q,r)       = myinterp1(25,lutotdiffgrid, dSofUdiffdu(p,:,q,r), lutotdiffoutgrid',interptype);

                PofCHIdiffout(p,:,q,r)        = 10.^(myinterp1(4,lchidiffgrid, log10(PofCHIdiff(p,:,q,r)), lchidiffoutgrid',interptype));
                % dPofCHIdiffdrho0out(p,:,q,r) = 10.^(myinterp1(26,lchidiffgrid, log10(dPofCHIdiffdrho0(p,:,q,r)), lchidiffoutgrid',interptype));
                % dPofCHIdiffdchiout(p,:,q,r)  = 10.^(myinterp1(27,lchidiffgrid, log10(dPofCHIdiffdchi(p,:,q,r)), lchidiffoutgrid',interptype));
                dPofCHIdiffdrho0out(p,:,q,r)  = myinterp1(26,lchidiffgrid, dPofCHIdiffdrho0(p,:,q,r), lchidiffoutgrid',interptype);
                dPofCHIdiffdchiout(p,:,q,r)   = myinterp1(27,lchidiffgrid, dPofCHIdiffdchi(p,:,q,r), lchidiffoutgrid',interptype);

                SSofCHIdiffout(p,:,q,r)        = 10.^(myinterp1(4,lchidiffgrid, log10(SSofCHIdiff(p,:,q,r)), lchidiffoutgrid',interptype));
                % dSSofCHIdiffdrho0out(p,:,q,r) = 10.^(myinterp1(26,lchidiffgrid, log10(dSSofCHIdiffdrho0(p,:,q,r)), lchidiffoutgrid',interptype));
                % dSSofCHIdiffdchiout(p,:,q,r)  = 10.^(myinterp1(27,lchidiffgrid, log10(dSSofCHIdiffdchi(p,:,q,r)), lchidiffoutgrid',interptype));
                dSSofCHIdiffdrho0out(p,:,q,r)  = myinterp1(26,lchidiffgrid, dSSofCHIdiffdrho0(p,:,q,r), lchidiffoutgrid',interptype);
                dSSofCHIdiffdchiout(p,:,q,r)   = myinterp1(27,lchidiffgrid, dSSofCHIdiffdchi(p,:,q,r), lchidiffoutgrid',interptype);

                tkofUdiffout(p,:,q,r)         = 10.^(myinterp1(29,lutotdiffgrid, log10(tkofUdiff(p,:,q,r)), lutotdiffoutgrid',interptype));
                tkofPdiffout(p,:,q,r)         = 10.^(myinterp1(30,lptotdiffgrid, log10(tkofPdiff(p,:,q,r)), lptotdiffoutgrid',interptype));
                tkofCHIdiffout(p,:,q,r)       = 10.^(myinterp1(31,lchidiffgrid, log10(tkofCHIdiff(p,:,q,r)), lchidiffoutgrid',interptype));
                tkofSdiffout(p,:,q,r)         = 10.^(myinterp1(30,lstotdiffgrid, log10(tkofSdiff(p,:,q,r)), lstotdiffoutgrid',interptype));

                % NOTE: Below will produce warning's about NaN found in interp1, but that's correct behavior since we use those NaN's to decide if really within reasonable part of table rather than using the extrapolation
                faketkofUdiffout(p,:,q,r)         = 10.^(interp1(lutotdiffgrid, log10(faketkofUdiff(p,:,q,r)), lutotdiffoutgrid',interptypefaketemp));
                faketkofPdiffout(p,:,q,r)         = 10.^(interp1(lptotdiffgrid, log10(faketkofPdiff(p,:,q,r)), lptotdiffoutgrid',interptypefaketemp));
                faketkofCHIdiffout(p,:,q,r)       = 10.^(interp1(lchidiffgrid, log10(faketkofCHIdiff(p,:,q,r)), lchidiffoutgrid',interptypefaketemp));
                faketkofSdiffout(p,:,q,r)         = 10.^(interp1(lstotdiffgrid, log10(faketkofSdiff(p,:,q,r)), lstotdiffoutgrid',interptypefaketemp));

                % extras:
                % mynewdata(:,:,:,:,ei)
                for ei=1:numextras
                  extraofUdiffout(p,:,q,r,ei) = 10.^(myinterp1(28,lutotdiffgrid, log10(extraofUdiff(p,:,q,r,ei)), lutotdiffoutgrid',interptype));
                end
                
                
              end
            end
          end

          %HofUout(1,:,1,1)


          fprintf(fiddebug,'End Downsample %d %d %d\n',hiter,titer,ynuiter);

          
          
          
          
          
          
          
          
          
          
          %%%%%%%%%%%%%%%%%%
          %
          % Force higher-order temperature variable to be NaN when linear (non-extrapolated) temperature variable is a NaN
          %
          %%%%%%%%%%%%%%%%%%
          myisnan=isnan(faketkofUdiffout);
          tkofUdiffout(myisnan)=NaN;

          myisnan=isnan(faketkofPdiffout);
          tkofPdiffout(myisnan)=NaN;

          myisnan=isnan(faketkofCHIdiffout);
          tkofCHIdiffout(myisnan)=NaN;

          myisnan=isnan(faketkofSdiffout);
          tkofSdiffout(myisnan)=NaN;

          % generate down-sampled rho
          for p=1:nutotdiffout
            for q=1:ntdynorye
              for r=1:nhcm
                % using 1 instead of p because of size of p is different
                rhobgrid4dout(:,p,q,r)    = rhob(:,1,q,r);
                rhobcsqgrid4dout(:,p,q,r) = rhobcsq(:,1,q,r);
              end
            end
          end
          for p=1:nrhob
            for q=1:ntdynorye
              for r=1:nhcm
                % diff values
                lutotdiffgrid4dout(p,:,q,r)  = lutotdiffoutgrid(:);
                lptotdiffgrid4dout(p,:,q,r)  = lptotdiffoutgrid(:);
                lchidiffgrid4dout(p,:,q,r)   = lchidiffoutgrid(:);
                lstotdiffgrid4dout(p,:,q,r)  = lstotdiffoutgrid(:);
                lsspecdiffgrid4dout(p,:,q,r)     = lsspecdiffoutgrid(:); % not really used

                % actual values (just different label)
                utotgrid4dout(p,:,q,r)  = UofUdiffout(p,:,q,r);
                ptotgrid4dout(p,:,q,r)  = PofPdiffout(p,:,q,r);
                chigrid4dout(p,:,q,r)   = CHIofCHIdiffout(p,:,q,r);
                stotgrid4dout(p,:,q,r)  = SofSdiffout(p,:,q,r);
              end
            end
          end

          
          if(usecleanvar)
            tkofUdiffout = cleanvar(29, tkofUdiffout, rhobgrid4dout, utotgrid4dout);
            tkofPdiffout = cleanvar(30, tkofPdiffout, rhobgrid4dout, ptotgrid4dout);
            tkofCHIdiffout = cleanvar(31, tkofCHIdiffout, rhobgrid4dout, chigrid4dout);
            tkofSdiffout = cleanvar(30, tkofSdiffout, rhobgrid4dout, stotgrid4dout);
          else
            tkofUdiffout(~isfinite(log10(tkofUdiffout)))=OUTBOUNDSVALUE;
            tkofPdiffout(~isfinite(log10(tkofPdiffout)))=OUTBOUNDSVALUE;
            tkofCHIdiffout(~isfinite(log10(tkofCHIdiffout)))=OUTBOUNDSVALUE;
            tkofSdiffout(~isfinite(log10(tkofSdiffout)))=OUTBOUNDSVALUE;
          end
          
          
          

          
          
          
          
          
          
          
          %          stot(p,:,q,r)=monotonize(stot(p,:,q,r));
          % Below is done because S(T) need not be monotonic, but U(S) and U(T) should be.  So can't force monotonicity of S(T) like did at top of this file with U(T) and P(T) and CHI(T)
          % NOT sure if necessary since no derivatives needed of it, so don't require monotonicity.
          % Below used in HARM for entropy-tracking of same fluid element.  Entropy lookup can be arbitrary and non-monotonic and that's fine.
%          fprintf(fiddebug,'Begin monotonize of final output for U(S) %d %d %d\n',hiter,titer,ynuiter);
%
%          for p=1:nrhob
%            for q=1:ntdynorye
%              for r=1:nhcm
%
%                UofSdiffout(p,:,q,r) = monotonize(UofSdiffout(p,:,q,r));
%
%              end
%            end
%          end
%          
%          fprintf(fiddebug,'End monotonize of final output for U(S) %d %d %d\n',hiter,titer,ynuiter);
          
          
          
          
          
          
     

          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %
          % Clear interpolated but pre-output version of variables
          %
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

          clear UofUdiff PofPdiff CHIofCHIdiff SofSdiff SSofSSdiff
          clear PofUdiff HofUdiff
          clear UofPdiff UofSdiff
          clear dPofUdiffdrho0 dPofUdiffdu cs2ofUdiff cs2ofUdiffcgs SofUdiff
          clear dSofUdiffdrho0 dSofUdiffdu PofCHIdiff dPofCHIdiffdrho0 dPofCHIdiffdchi SSofCHIdiff dSSofCHIdiffdrho0 dSSofCHIdiffdchi tkofUdiff tkofPdiff tkofCHIdiff tkofSdiff
          clear HofU cs2 PofU UofP dPofUdrho0 dPofUdu cs2cgs SofU dSofUdrho0 dSofUdu PofCHI dPofCHIdrho0 dPofCHIdchi SSofCHI dSSofCHIdrho0 dSSofCHIdchi
          %tkofU tkofP tkofCHI tkofS;
          clear extraofU;



          
          
          
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %
          % Adjust quantities to be physical in case of numerical error
          %
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


          [mincs2out,mincs2outI]=min(min(min(min(cs2ofUdiffout))));
          fprintf(fiddebug,'Old Min c_s^2/c^2: %21.15g\n',mincs2out);
          [maxcs2out,maxcs2outI]=max(max(max(max(cs2ofUdiffout))));
          fprintf(fiddebug,'Old Max c_s^2/c^2: %21.15g\n',maxcs2out);

          [mincs2cgsout,mincs2cgsoutI]=min(min(min(min(cs2ofUdiffcgsout))));
          fprintf(fiddebug,'Old Min c_s^2[cgs]: %21.15g\n',mincs2cgsout);
          [maxcs2cgsout,maxcs2cgsoutI]=max(max(max(max(cs2ofUdiffcgsout))));
          fprintf(fiddebug,'Old Max c_s^2[cgs]: %21.15g\n',maxcs2cgsout);

          
          % adjust speed of sound to be no smaller than 0 and no larger than 1

          for p=1:nhcm
            for o=1:ntdynorye
              for n=1:nutotdiffout
                for m=1:nrhob
                  
                  if cs2ofUdiffout(m,n,o,p)>1.0
                    cs2ofUdiffout(m,n,o,p)=1.0-CONTOL;
                  end
                  if cs2ofUdiffout(m,n,o,p)<0.0
                    cs2ofUdiffout(m,n,o,p)=0.0;
                  end

                  if cs2ofUdiffcgsout(m,n,o,p)>c.*c
                    cs2ofUdiffcgsout(m,n,o,p)=(1.0-CONTOL).*c.*c;
                  end
                  if cs2ofUdiffcgsout(m,n,o,p)<0.0
                    cs2ofUdiffcgsout(m,n,o,p)=0.0;
                  end

                end
              end
            end
          end

          [mincs2out,mincs2outI]=min(min(min(min(cs2ofUdiffout))));
          fprintf(fiddebug,'New Min c_s^2/c^2: %21.15g\n',mincs2out);
          [maxcs2out,maxcs2outI]=max(max(max(max(cs2ofUdiffout))));
          fprintf(fiddebug,'New Max c_s^2/c^2: %21.15g\n',maxcs2out);

          [mincs2cgsout,mincs2cgsoutI]=min(min(min(min(cs2ofUdiffcgsout))));
          fprintf(fiddebug,'New Min c_s^2[cgs]: %21.15g\n',mincs2cgsout);
          [maxcs2cgsout,maxcs2cgsoutI]=max(max(max(max(cs2ofUdiffcgsout))));
          fprintf(fiddebug,'New Max c_s^2[cgs]: %21.15g\n',maxcs2cgsout);





          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %
          % write data to file
          %
          % Note that all ratios (or derivatives) are dimensionless
          %
          % all other quantities, including sound speed, have dimensions
          %
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




          % GODMARK: Check units of entropy stuff

          %                HofUout(m,n,o,p) , cs2out(m,n,o,p), ...
          % %21.15g %21.15g 
          

          
          
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %
          % DO FINAL CLEANING since interpolation above still might lead to NaN if thinks out of bounds
          %
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%

          if doclean

            UofUdiffout(~isfinite(UofUdiffout))=OUTBOUNDSVALUE;
            PofPdiffout(~isfinite(PofPdiffout))=OUTBOUNDSVALUE;
            CHIofCHIdiffout(~isfinite(CHIofCHIdiffout))=OUTBOUNDSVALUE;
            SofSdiffout(~isfinite(SofSdiffout))=OUTBOUNDSVALUE;

            PofUdiffout(~isfinite(PofUdiffout))=OUTBOUNDSVALUE;

            UofPdiffout(~isfinite(UofPdiffout))=OUTBOUNDSVALUE;
            UofSdiffout(~isfinite(UofSdiffout))=OUTBOUNDSVALUE;

            dPofUdiffdrho0out(~isfinite(dPofUdiffdrho0out))=OUTBOUNDSVALUE;
            dPofUdiffduout(~isfinite(dPofUdiffduout))=OUTBOUNDSVALUE;

            cs2ofUdiffcgsout(~isfinite(cs2ofUdiffcgsout))=OUTBOUNDSVALUE;

            SofUdiffout(~isfinite(SofUdiffout))=OUTBOUNDSVALUE;
            dSofUdiffdrho0out(~isfinite(dSofUdiffdrho0out))=OUTBOUNDSVALUE;
            dSofUdiffduout(~isfinite(dSofUdiffduout))=OUTBOUNDSVALUE;


            PofCHIdiffout(~isfinite(PofCHIdiffout))=OUTBOUNDSVALUE;
            dPofCHIdiffdrho0out(~isfinite(dPofCHIdiffdrho0out))=OUTBOUNDSVALUE;
            dPofCHIdiffdchiout(~isfinite(dPofCHIdiffdchiout))=OUTBOUNDSVALUE;

            SSofCHIdiffout(~isfinite(SSofCHIdiffout))=OUTBOUNDSVALUE;
            dSSofCHIdiffdrho0out(~isfinite(dSSofCHIdiffdrho0out))=OUTBOUNDSVALUE;
            dSSofCHIdiffdchiout(~isfinite(dSSofCHIdiffdchiout))=OUTBOUNDSVALUE;

            tkofUdiffout(~isfinite(tkofUdiffout))=OUTBOUNDSVALUE;
            tkofPdiffout(~isfinite(tkofPdiffout))=OUTBOUNDSVALUE;
            tkofCHIdiffout(~isfinite(tkofCHIdiffout))=OUTBOUNDSVALUE;
            tkofSdiffout(~isfinite(tkofSdiffout))=OUTBOUNDSVALUE;


            extraofUdiffout(~isfinite(extraofUdiffout))=OUTBOUNDSVALUE;
            
            
          end
          
          
          
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %
          % Clean derivatives because of numerical roundoff error in d/drho0 calculations at high T and low rho0
          %
          % d/drho0 can be negative and correctly so, so only correct region where we know calculation fails to give correct answer
          %
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
          if finaldoclean
            
            for p=1:nhcm
              for o=1:ntdynorye
                for n=1:nutotdiffout

                  for m=1:nrhob
                    % based upon P=(\gamma-1)u with \gamma~2 being normal highest and assuming can have up to \gamma=3.  Trying to keep inversion method stable.
                    if(dPofUdiffduout(m,n,o,p)>2.0)
                      dPofUdiffduout(m,n,o,p)=2.0;
                    end
                  end

                  for m=1:nrhob
                    % based upon P=(\gamma-1)u assuming P=(\gamma-1)/\gamma chi so maximum of dP/dchi is 1.0 no matter what gamma.   Trying to keep inversion method stable.
                    if(dPofCHIdiffdchiout(m,n,o,p)>1.0)
                      dPofCHIdiffdchiout(m,n,o,p)=1.0;
                    end
                  end

                  
                  
                  m0=-1;
                  for m=nrhob:-1:1
                    %            if(dPofUdiffdrho0out(m,n,o,p)<0.0 || dPofUdiffdrho0out(m,n,o,p)>5.0)
                    if((dPofUdiffdrho0out(m,n,o,p)<0.0 || dPofUdiffdrho0out(m,n,o,p)>5.0)&&(tkofUdiffout(m,n,o,p)>1E11)&&(rhob(m,n,o,p)<1E8)  )
                      m0=m;
                    end
                    if(m0>=1)
                      %fprintf(fiddebug,'dP1fix %d %d %d %d\n',p,o,n,m0);
                      for mm=m0:-1:1
                        dPofUdiffdrho0out(mm,n,o,p)=OUTBOUNDSVALUE;
                      end
                      break
                    end
                  end
                  
                  m0=-1;
                  for m=nrhob:-1:1
                    %            if(dSofUdiffdrho0out(m,n,o,p)<0.0 || dSofUdiffdrho0out(m,n,o,p)>5.0)
                    % not sure if entropy version can be negative
                    if((dSofUdiffdrho0out(m,n,o,p)<0.0 || dSofUdiffdrho0out(m,n,o,p)>5.0)&&(tkofUdiffout(m,n,o,p)>1E11)&&(rhob(m,n,o,p)<1E8)  )
                      %            if(dSofUdiffdrho0out(m,n,o,p)<0.0 && 0)
                      m0=m;
                    end
                    if(m0>=1)
                      %fprintf(fiddebug,'dP2fix %d %d %d %d\n',p,o,n,m0);
                      for mm=m0:-1:1
                        dSofUdiffdrho0out(mm,n,o,p)=OUTBOUNDSVALUE;
                      end
                      break
                    end
                  end
                  
                  m0=-1;
                  for m=nrhob:-1:1
                    if((dPofCHIdiffdrho0out(m,n,o,p)<0.0 || dPofCHIdiffdrho0out(m,n,o,p)>5.0)&&(tkofCHIdiffout(m,n,o,p)>1E11)&&(rhob(m,n,o,p)<1E8)  )
                      %            if(dPofCHIdiffdrho0out(m,n,o,p)<0.0 || dPofCHIdiffdrho0out(m,n,o,p)>5.0)
                      m0=m;
                    end
                    if(m0>=1)
                      %fprintf(fiddebug,'dP3fix %d %d %d %d\n',p,o,n,m0);
                      for mm=m0:-1:1
                        dPofCHIdiffdrho0out(mm,n,o,p)=OUTBOUNDSVALUE;
                      end
                      break
                    end
                  end

                  m0=-1;
                  for m=nrhob:-1:1
                    if((dSSofCHIdiffdrho0out(m,n,o,p)<0.0 || dSSofCHIdiffdrho0out(m,n,o,p)>5.0)&&(tkofCHIdiffout(m,n,o,p)>1E11)&&(rhob(m,n,o,p)<1E8)  )
                      %            if(dSSofCHIdiffdrho0out(m,n,o,p)<0.0 || dSSofCHIdiffdrho0out(m,n,o,p)>5.0)
                      m0=m;
                    end
                    if(m0>=1)
                      %fprintf(fiddebug,'dP3fix %d %d %d %d\n',p,o,n,m0);
                      for mm=m0:-1:1
                        dSSofCHIdiffdrho0out(mm,n,o,p)=OUTBOUNDSVALUE;
                      end
                      break
                    end
                  end

                  
                end
              end
            end
            
            
          end


          
          
          
          
          
          
          
          
          
          %%%%%%%%%%%%%%%%%%%%%%%%
          %
          % Notice that below indicies are in C-order and C-style (start with 0 instead of 1)
          %
          % Assumes nutotdiffout is same size as nptotout and nchiout
          %
          %%%%%%%%%%%%%%%%%%%%%%%
          fprintf(fiddebug,'Begin final output: %d %d %d',hiter,titer,ynuiter);


          
          % large-memory
          %m-1, n-1, o-1, p-1, ...
          % small-memory
          
          for p=1:nhcm
            for o=1:ntdynorye
              for n=1:nutotdiffout
                for m=1:nrhob
                  %            0                  +5                                                      +8                              +4      +1               +2              +2      +1                     +3                      +3                      +3                              +4 
                  fprintf(fid3,'%3d %3d %3d %3d %3d %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g ', ...
                          m-1, n-1, titer-1, ynuiter-1, hiter-1,...
                          rhob(m,n,o,p), lutotdiffoutgrid(n), lptotdiffoutgrid(n), lchidiffoutgrid(n), lstotdiffoutgrid(n), tdynorye(m,n,o,p), tdynorynu(m,n,o,p), hcm(m,n,o,p), ...
                          UofUdiffout(m,n,o,p), PofPdiffout(m,n,o,p), CHIofCHIdiffout(m,n,o,p), SofSdiffout(m,n,o,p), ...
                          PofUdiffout(m,n,o,p), ...
                          UofPdiffout(m,n,o,p), UofSdiffout(m,n,o,p), ...
                          dPofUdiffdrho0out(m,n,o,p), dPofUdiffduout(m,n,o,p), ...
                          cs2ofUdiffcgsout(m,n,o,p), ...
                          SofUdiffout(m,n,o,p), dSofUdiffdrho0out(m,n,o,p), dSofUdiffduout(m,n,o,p), ...
                          SSofCHIdiffout(m,n,o,p), dSSofCHIdiffdrho0out(m,n,o,p), dSSofCHIdiffdchiout(m,n,o,p), ...
                          PofCHIdiffout(m,n,o,p), dPofCHIdiffdrho0out(m,n,o,p), dPofCHIdiffdchiout(m,n,o,p), ...
                          tkofUdiffout(m,n,o,p),tkofPdiffout(m,n,o,p),tkofCHIdiffout(m,n,o,p),tkofSdiffout(m,n,o,p) ...
                          );
                  % mynewdata(:,:,:,:,ei)
                  for ei=1:numextras
                    fprintf(fid3,'%21.15g ',extraofUdiffout(m,n,o,p,ei));
                  end
                  fprintf(fid3,'\n');
                end
              end
            end
          end








          
          
          %%%%%%%%%%%%%%%%%%%%%%%%
          %
          %
          %  Output utot,ptot,chi as functions of rhob for T=0
          %
          %
          %%%%%%%%%%%%%%%%%%%%%%%

          % large-memory
          %m-1, o-1, p-1, ...
          % small-memory

          % corresponds to n=1 solution since u~0 implies T~0.  Reduced to degenerate solution independent of temperature
          
          
          
          n=1;
          for p=1:nhcm
            for o=1:ntdynorye
              for m=1:nrhob
                
                if 1
                  %utotoffset(m,n,o,p) = UofUdiffout(m,n,o,p) - utotdiffoutgrid(n) NO!
                  
                end
                
                %           +0              +4                              +4                              +4                              +4                              +4
                fprintf(fid6,'%3d %3d %3d %3d %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g ', ...
                        m-1, titer-1, ynuiter-1, hiter-1, ...
                        rhob(m,n,o,p), tdynorye(m,n,o,p), tdynorynu(m,n,o,p), hcm(m,n,o,p), ...
                        utotoffset(m,n,o,p), ptotoffset(m,n,o,p), chioffset(m,n,o,p), stotoffset(m,n,o,p), ...
                        utotin(m,n,o,p), ptotin(m,n,o,p), chiin(m,n,o,p), stotin(m,n,o,p), ...
                        utotout(m,n,o,p), ptotout(m,n,o,p), chiout(m,n,o,p), stotout(m,n,o,p) ...
                        );
                fprintf(fid6,'\n');
              end
            end
          end
            

          %  Should have utotdiffoutgrid + utotoffset = UofUdiffout





          % some debug stuff:
          %
          % for ii=1:48 fprintf(fiddebug,'%21.15g %21.15g\n',log10(utot(28,ii,1,1)),log10(tk(28,ii,1,1)));end




          fprintf(fiddebug,'End final output %d %d %d\n',hiter,titer,ynuiter);




        end

        
        
        
        % end loops over truehcm and truetdynorye and truetdynorynu
      end
    end
    end
    
    
    
    
    % close files
    fclose(fid2);
    fclose(fid8);
    if 0
      % didn't open eosother.dat
      fclose(fid4);
    end
    fclose(fid3);
    fclose(fid6);

    

  
    % go ahead and output header after first pass is done
    if(passiter==1 || require2passes==0)
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %
      % write header to file
      %
      % Note that steps and all such things are computed same as in Kaz's code
      %
      % We are outside loop here, so all these quantities are defined over ALL rhob,(u,p,chi),tdynorye,tdynornynu,hcm
      % notice that the "steps" will be redefined correctly during second pass
      %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      fid5=fopen(file5,'w');

      
      %number of columns outputted (including numbers indicating position of element in table)
      %NUMOUTCOLUMNS=25 % was 23 before new degen offset method
      % was (maybe) 27 + numextras before ynu
      NUMINDEPDIMENS=5; % same as in kazfulleos.global.h (for checking HARM expectation with Matlab output)
      NUMEOSINDEPS=8; % same as in kazfulleos.global.h (for checking HARM expectation with Matlab output)
      NUMVAR1=4; % utotdiff,ptotdiff,chidiff,stotdiff (for checking HARM expectation with Matlab output)
      % below begins PofRHOU, etc. as in kazfulleos.c
      NUMFUN1=1+2+2;
      NUMCS=1;
      NUMFUN2=3+3+3;
      NUMTEMP=4;
      % 29 + numextras
      NUMOUTCOLUMNS=NUMINDEPDIMENS+NUMEOSINDEPS+NUMVAR1+NUMFUN1+NUMCS+NUMFUN2+NUMTEMP+numextras;

      
      % first is table size, then as read-in by HARM code:
      % second [4] : 0 = lower log_base limit, 1 = upper log_base limit, 2=step, 3 = divisor of grid position 4=base of log, 5 = linear value of offset for log_base stepping so can control how resolved


      % set base and linear offset here since not using this in eos_extract.m yet
      baselrhob=10.0;
      linearoffsetlrhob=0.0;

      baselutotdiffout=10.0;
      linearoffsetlutotdiffout=0.0;
      baselptotdiffout=10.0;
      linearoffsetlptotdiffout=0.0;
      baselchidiffout=10.0;
      linearoffsetlchidiffout=0.0;
      baselstotdiffout=10.0;
      linearoffsetlstotdiffout=0.0;

      baseltdynorye=10.0;
      linearoffsetltdynorye=0.0;
      baseltdynorynu=10.0;
      linearoffsetltdynorynu=0.0;
      baselhcm=10.0;
      linearoffsetlhcm=0.0;

      baseltk=10.0;
      linearoffsetltk=0.0;

      % method, number of output colums, num extras
      fprintf(fid5,'%d %d %d\n',whichrnpmethod,whichynumethod,whichhcmmethod);
      fprintf(fid5,'%d %d %d\n',whichdatatype,NUMOUTCOLUMNS,numextras);
      fprintf(fid5,'%d %21.15g %21.15g %21.15g %21.15g %21.15g\n',nrhob,lrhobmin,lrhobmax,steplrhob,baselrhob,linearoffsetlrhob);

      fprintf(fid5,'%d %21.15g %21.15g %21.15g %21.15g %21.15g\n',nutotdiffout,lutotdiffoutmin,lutotdiffoutmax,steplutotdiffout,baselutotdiffout,linearoffsetlutotdiffout);
      fprintf(fid5,'%d %21.15g %21.15g %21.15g %21.15g %21.15g\n',nptotdiffout,lptotdiffoutmin,lptotdiffoutmax,steplptotdiffout,baselptotdiffout,linearoffsetlptotdiffout);
      fprintf(fid5,'%d %21.15g %21.15g %21.15g %21.15g %21.15g\n',nchidiffout,lchidiffoutmin,lchidiffoutmax,steplchidiffout,baselchidiffout,linearoffsetlchidiffout);
      fprintf(fid5,'%d %21.15g %21.15g %21.15g %21.15g %21.15g\n',nstotdiffout,lstotdiffoutmin,lstotdiffoutmax,steplstotdiffout,baselstotdiffout,linearoffsetlstotdiffout);

      fprintf(fid5,'%d %21.15g %21.15g %21.15g %21.15g %21.15g\n',truentdynorye,ltdynoryemin,ltdynoryemax,stepltdynorye,baseltdynorye,linearoffsetltdynorye);
      fprintf(fid5,'%d %21.15g %21.15g %21.15g %21.15g %21.15g\n',truentdynorynu,ltdynorynumin,ltdynorynumax,stepltdynorynu,baseltdynorynu,linearoffsetltdynorynu);
      fprintf(fid5,'%d %21.15g %21.15g %21.15g %21.15g %21.15g\n',truenhcm,lhcmmin,lhcmmax,steplhcm,baselhcm,linearoffsetlhcm);

      fprintf(fid5,'%21.15g %21.15g\n',lsoffset,fakelsoffset);

      fprintf(fid5,'%d %21.15g %21.15g %21.15g %21.15g %21.15g\n',ntk,ltkmin,ltkmax,stepltk,baseltk,linearoffsetltk);

      fclose(fid5);
    end

  
    
  end   %%%%%%%% end over passes
  


  
  %%%%%%%%%%%%%%%%%%
  %
  % Close files
  %
  %%%%%%%%%%%%%%%%%%
  
  fclose(fiddebug);

  fclose('all');

  % BELOW FOR in linux
  if matlabscript
    quit;
  end


end






















% See phivsr.m
% should not be function of temperature (tk)
function out = udegenfitOLD(rho0)

K = 1.22929981645031E15;
A = 1.51663;
B = 4.6E6;
C = 0.9999;
D = 1.0;
Ye = yefit(rho0,tk);

out = A.*K.*((rho0+B).*Ye).^(C/3).*rho0.^D;

% from SM:
% set udegenfit1=1.51663*K*((rhob+4.6E6)*yefit)**(0.9999/3)*rhob**(1.0)

end


% See phivsr.m
% should not be function of temperature (tk)
function out = udegenfit_fun(rhob)

Ye = yefit(rhob);
K = 1.22929981645031E15;

A = 1.516;
B = 0.0E6;
C = 4.0/3.0;
D = 1.0/3.0;
udegenfit1 = A.*K.*(rhob+B).^(C).*Ye.^(D);

A = 1.66*1E2*1.516;
B = 0.0E6;
C = 1.0;
D = 1.0/3.0;
udegenfit2 = A.*K.*(rhob+B).^(C).*Ye.^(D);

npow=1.95;

% larger function wins
out = (udegenfit1.^npow + udegenfit2.^npow).^(1/npow);


%		FROM SM:
% try log with break
%		set npow=1.95
%		set udegenfit1=1.516*K*((rhob+0.0E6))**(4.0/3)* yefit**(1.0/3)
%		set udegenfit2=1.66*1E2*1.516*K*((rhob+0.0E6))**(1.0)* yefit**(1.0/3)
%		#set udegenfit = (1.0/(1.0/udegenfit1**npow+1.0/udegenfit2**npow))**(1/npow)
%		set udegenfit = (udegenfit1**npow + udegenfit2**npow)**(1/npow)

end








% See phivsr.m
function out = pdegenfit_fun(rhob)

K = 1.22929981645031E15;
A1 = 3.0310;
A2 = 0.0235;
B = 0.0;
Ye = yefit(rhob);

pdegenfit1=(A1.*K.*((rhob+B).*yefit).^(1/3).*rhob).*(4.0/3.0-1.0).*yefit;
pdegenfit2=(A2.*K.*((rhob+B).*yefit).^(2/3).*rhob).*(4.0/3.0-1.0).*yefit;
npow=1.94;

out = (1.0./(1.0./pdegenfit1.^npow+1.0./pdegenfit2.^npow)).^(1./npow);

% from SM:
% below best fit so far for ptot
%	set pdegenfit1=(3.0310*K*((rhob+0E6)*yefit)**(1/3)*rhob)*(4.0/3.0-1.0)*yefit
%		set pdegenfit2=(0.0235*K*((rhob+0E6)*yefit)**(2.00/3)*rhob)*(4.0/3.0-1.0)*yefit
%		set npow=1.94
%		set pdegenfit = (1.0/(1.0/pdegenfit1**npow+1.0/pdegenfit2**npow))**(1/npow)


end


% \chi\equiv u+p
function out = chidegenfit_fun(rhob)

out = udegenfit_fun(rhob) + pdegenfit_fun(rhob);

end






% From jon_helm.f
%     From SM, compute Ye for presupernova-like abundances
function out = yefit(rhob)

tk=0.0;

if(tk<10.^(9.532))
  yefit = 0.5;
else
  if(tk>10.^(9.92))
    yefit = 0.43;
  else
    yefit = (0.43-0.5)/(10^(9.92)-10^(9.532))*(tk-10.^(9.532))+0.5;
  end
end

out = yefit;

end


