function shen_interp()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SOME PARAMETERS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set whether on laptop (interactive) or relativity (matlab script)
  onlinux=1

  
  %  For mex stuff:
  %  1) modified /opt/matlab/bin/f90opts.sh to use ifort with HELM FFLAGS
  %  2) removed -libstdc++ and added -cpp to FFLAGS so can handle #include's
  %  3) ran mex -setup and chose #1 F90 option and say "y" to save
  %  4) run below: mex namefunction.mex.f [extra files needed to satify dependencies]

  fprintf('Begin compile mex functions\n');
  % Use -largeArrayDims on 64-bit machines.  Not on relativity.
  % I have some problem with 64-bit machines and mex processing: segfaults on mxgetpr on LHS data
%  mex -largeArrayDims -v -g -output enforce_x_consistency /home/matlab/matlab/enforce_x_consistency.mex.f /home/matlab/matlab/general_compute_functions.f
%  mex -largeArrayDims -v -g  -output enforce_ye_consistency /home/matlab/matlab/enforce_ye_consistency.mex.f /home/matlab/matlab/general_compute_functions.f
%  mex -largeArrayDims -v -g  -output compute_nuclear_azbar /home/matlab/matlab/compute_nuclear_azbar.mex.f /home/matlab/matlab/general_compute_functions.f
  % RELATIVITY (32-bit machines)
  mex -v -g -output enforce_x_consistency /home/matlab/matlab/enforce_x_consistency.mex.f /home/matlab/matlab/general_compute_functions.f
  mex -v -g  -output enforce_ye_consistency /home/matlab/matlab/enforce_ye_consistency.mex.f /home/matlab/matlab/general_compute_functions.f
  mex -v -g  -output compute_nuclear_azbar /home/matlab/matlab/compute_nuclear_azbar.mex.f /home/matlab/matlab/general_compute_functions.f
  fprintf('End compile mex functions\n');

  % test:
  fprintf('Begin mex test1\n');
  [abarnum abar zbar yelocal] = compute_nuclear_azbar(0.9,0.1,0.0,0.1,1,1)
  fprintf('Begin mex test2\n');
  [xconsistent xnut xprot xalfa xh a x] = enforce_x_consistency(0.5,0.9,0.1,0.0,0.1,1,1)
  fprintf('Begin mex test3\n');
  [yeconsistent  xnut xprot xalfa xh a x] = enforce_ye_consistency(0.5,0.9,0.1,0.0,0.1,1,1)
  fprintf('End mex tests\n');
  
  % STILL:  yetrue,xnut,xprot,xalfa,xh,a,x);
  % func input: (ye_inp, xnut, xprot, xalfa, xh, a, x)
  %   [xconsistent xnut xprot xalfa xh a x] = enforce_x_consistency(0.13744,6.86564842562921,2.28653862719783,2.000000000000000E-010,1.000000013351432E-010,0.388093,5.44071e-12)
  
  
  
  
  
  
  if onlinux==0
    % BELOW USED ON laptop
    %dir='C:\\Documents and Settings\\jon\\My Documents\\grbjet\\eoslarge\\';
    dir='C:\\Documents and Settings\\jon\\My Documents\\grbjet\\eoslargesimple\\';
  end
  if onlinux==1
    % BELOW USED ON relativity
    dir='./';
  end
  %dir='C:\\Documents and Settings\\jon\\My Documents\\eossmall\\';
  prefix='sheneos.tab';
  prefix2='sheneos.t00';
  prefix3='sheneos.yp0';
  prefix4='sheneos.dat';
  prefix5='sheneos.head';
  
  file1=strcat(dir,prefix);
  file2=strcat(dir,prefix2);
  file3=strcat(dir,prefix3);
  file4=strcat(dir,prefix4);
  file5=strcat(dir,prefix5);

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % from const.dek (mb used since output is consistent with definition of
  % mb instead of Shen definition of \rho_b=m_b/amu so that there is
  % alignment between rho_b outputted and used in Fortran code)
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % 1MeV = mev2K  K
  mev2K=1.1605d10;
  %avoreal = 6.02214179d23;
  %amu    = 1.0/avoreal;
  C=2.99792458E10;
  ergPmev = 1.782661758E-30*1.0E3*C*C;
  % below is $m_{u}$ in MeV from Shen's guide.tex
  % can't change this!  If mumev different today, then below stays same and m_b will be changed elsewhere to compensate
  mumev = 931.49432;
  % below amu is just from mumev
  amu=mumev*ergPmev/C/C;
  avoreal = 1.0/amu;
  me=9.10938215E-28;
  mn=1.674927211E-24;
  mp=1.672621637E-24;
  %mn      = 1.6749286d-24;
  %mp      = 1.6726231d-24;
  % below was for Kaz EOS
  %   mb      = (0.5*(mp+mn));
  % below now for consistency across all EOSs
  % below should be JON's code value of mb, which is amu but amu might be newer today
  modernamu=1.660538782E-24;
  mb      = modernamu;

  
  % these numbers should be consistent with vector_eos.dek numbers
  azmin = 1E-10;
  xmin=1E-10;
  aheavtrust=1.0;
  zheavtrust=azmin;
  xhtrust=xmin;
  yetolerance=1E-10;
  xchecktolerance=1E-10;

  % these should be consistent with jon_sheneos.f
  nbmin=1E-30;
  mstarfloor=1e-10;


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % interpolate each function with some method
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % whether to use interpn (doesn't allow pchip and doesn't extrapolate) or to use multiple 1D interpolations
  useinterpn=0;

  % funout = interpn(xin,yin,zin,funin,xout,yout,zout,'spline');
  % nc items
  fakeinterpntype='linear';
  %interpntype='spline';
  interpntype='cubic';
  %interpntypemu='v5cubic';
  % spline too crazy near boundaries
  %interpntypemu='spline';
  interpntypemu='cubic';

  fakeinterp1type='linear';
  %interpntype='spline';
  interp1type='pchip';
  % spline too crazy near boundaries
  %interp1typemu='spline';
  interp1typemu='pchip';


    % pchip is cubic -- keeps monotonicity and otherwise higher-order

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Input table sizes
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  nc=17; % number of columns
  nrhob=104; % fastest index
  nyp=72; % 1 is Yp=0 % next fastest
  ntk=32; % 1 is T=0 % slowest index

  % From my fortran code
  % Objects:
  %
  %  1) log10(rhob) [g/cm^3] : log10 Baryon rest-mass density
  %  2) n_b [fm^{-3}]        : Baryon number density (rhob = n_b*amu)
  %  3) log10(Yp)            : log10 Proton number
  %  4) Y_p                  : Proton fraction
  %                            Y_p = N_p/N_B in general, or Y_p = (n_p+2n_\alpha)/(n_n+n_p+4n_\alpha)
  %  5) F [MeV]              : Free energy per baryon
  %                            F = f/n_b - M  where f=total free energy and M=938MeV
  %  6) E [MeV]              : Internal energy per baryon
  %                            E = \ep/n_b - m_u where \ep is total energy
  %                            m_u = 931.49432MeV
  %  7) S [k_b]              : Entropy per baryon
  %                            S = s/n_b where s=entropy density
  %  8) aheav                : Mass number of heavy nucleus
  %  9) zheav                : Charge number of heavy nucleus
  % 10) M* [MeV]             : Effective mass
  % 11) xneut                : Free neutron mass fraction
  % 12) xprot                : Free proton mass fraction
  % 13) xalfa                : Alpha mass fraction
  % 14) xheav                : Heavy nucleus mass fraction
  % 15) P [MeV/fm^3]         : Pressure
  % 16) \tilde{\mu_n} [MeV]  : Chemical potential of neutrons relative to free nucleon mass M
  %                            n_n = (1-Y_p)n_b
  % 17) \tilde{\mu_p} [MeV]  : Chemical potential of protons relative to free nucleon mass M
  %                            n_p = Y_p n_b
  %
  


  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Input data ranges
  % x = rho_b (converted from Shen definition to Kaz definition)
  % y = Yp
  % z = Temp
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



  % (1) eos.tab (main EOS table, size: 54.5MB)
  % \item temperature $T (MeV)$: 
  % $ -1.0 \leq \log_{10}(T) \leq 2.0$
  % mesh of $\log_{10}(T) \sim 0.1      $
  % \item proton fraction $Y_p$:
  % $ -2.00 \leq \log_{10}(Y_p) \leq -0.25 $
  % mesh of $\log_{10}(Y_p) \sim 0.025 $
  % \item baryon mass density $\rho_B (g/cm^3)$:
  % $ 5.1 \leq \log_{10}(\rho_B) \leq 15.4$
  % mesh of $\log_{10}(\rho_B) \sim 0.1 $

  
  
  % in g/cc [cgs]
  lrhobminin=log10(10.^5.1/amu*mb); % converted to mb from amu!
  lrhobmaxin=log10(10.^15.4/amu*mb); % converted to mb from amu!
  lrhobstepin=(lrhobmaxin-lrhobminin)./(nrhob-1);
  
  lypminin=-2.0;
  lypmaxin=-0.25;
  lypstepin=(lypmaxin-lypminin)./(nyp-1-1);

  % in MeV [nuclear units]
  ltkminin=-1.0;
  ltkmaxin=2.0;
  ltkstepin=(ltkmaxin-ltkminin)./(ntk-1-1);

  % lowest effective values of Yp and T[K]

  % for easy interplation purposes
  % Yp=10^(-9) is used instead of Yp=0 for Shen's sheneos.yp0 table (i.e. far
  % from minimum yp on sheneos.tab table
  lYpfloor=lypminin-7;
  % T[MeV]=10^(-9) MeV is used instead of T= 0 MeV (i.e. far from minimum
  % T in sheneos.t00 table
  ltempfloor=ltkminin-7;

  Ypfloor=10.^(lYpfloor);
  tempfloor=10.^(ltempfloor);


  
  
  xv = [            lrhobminin:lrhobstepin:lrhobmaxin];
  yv = [lYpfloor    lypminin:lypstepin:lypmaxin];
  zv = [ltempfloor  ltkminin:ltkstepin:ltkmaxin];
  [xin,yin,zin] = ndgrid(xv,yv,zv);



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %
  % open Shen EOS files
  %
  %
  fprintf('Open File\n');

  
  
  fid=fopen(file1);
  [mydata,count]=fscanf(fid,'%g',[nc*(ntk-1)*nrhob*(nyp-1)]);
  fclose(fid);
  % make vector into 4-D array
  temp=reshape(mydata,[nc,nrhob,nyp-1,ntk-1]);
  clear mydata;
  % set primary matricies to be each column(field) read in from SM data (i.e.
  % read data is setup with column / nx / ny / nz order
  mynewdata1=permute(temp,[2,3,4,1]);
  clear temp;

  fid=fopen(file2);
  [mydata,count]=fscanf(fid,'%g',[nc*(1)*nrhob*(nyp-1)]);
  fclose(fid);
  % make vector into 4-D array
  temp=reshape(mydata,[nc,nrhob,nyp-1,1]);
  clear mydata;
  % set primary matricies to be each column(field) read in from SM data (i.e.
  % read data is setup with column / nx / ny / nz order
  mynewdata2=permute(temp,[2,3,4,1]);
  clear temp;

  fid=fopen(file3);
  [mydata,count]=fscanf(fid,'%g',[nc*(ntk-1)*nrhob*(1)]);
  fclose(fid);
  % make vector into 4-D array
  temp=reshape(mydata,[nc,nrhob,1,ntk-1]);
  clear mydata;
  % set primary matricies to be each column(field) read in from SM data (i.e.
  % read data is setup with column / nx / ny / nz order
  mynewdata3=permute(temp,[2,3,4,1]);
  clear temp;

  %
  %		so now all functions are F(rhob,yp,tk,column)
  %

  % error reporting file
  fid=fopen('error.txt','w');

  fprintf('Distribute Data\n');

  for r=1:ntk
    for q=1:nyp
      for p=1:nrhob

        if((r==1)&(q>1))
          
          % T=0 case
          
          % nc items
          ii=1;
          shenlrhob_row(p,q,r) =mynewdata2(p,q-1,r,ii); ii=ii+1;
          shennb_row(p,q,r) =mynewdata2(p,q-1,r,ii); ii=ii+1;
          shenlyp_row(p,q,r) =mynewdata2(p,q-1,r,ii); ii=ii+1;
          shenyp_row(p,q,r) =mynewdata2(p,q-1,r,ii); ii=ii+1;
          shenf_row(p,q,r) =mynewdata2(p,q-1,r,ii); ii=ii+1;
          shenebulk_row(p,q,r) =mynewdata2(p,q-1,r,ii); ii=ii+1;
          shensbulk_row(p,q,r) =mynewdata2(p,q-1,r,ii); ii=ii+1;
          shenaheav_row(p,q,r) =mynewdata2(p,q-1,r,ii); ii=ii+1;
          shenzheav_row(p,q,r) =mynewdata2(p,q-1,r,ii); ii=ii+1;
          shenmstar_row(p,q,r) =mynewdata2(p,q-1,r,ii); ii=ii+1;
          shenxneut_row(p,q,r) =mynewdata2(p,q-1,r,ii); ii=ii+1;
          shenxprot_row(p,q,r) =mynewdata2(p,q-1,r,ii); ii=ii+1;
          shenxalfa_row(p,q,r) =mynewdata2(p,q-1,r,ii); ii=ii+1;
          shenxh_row(p,q,r) =mynewdata2(p,q-1,r,ii); ii=ii+1;
          shenpbulk_row(p,q,r) =mynewdata2(p,q-1,r,ii); ii=ii+1;
          shenmun_row(p,q,r) =mynewdata2(p,q-1,r,ii); ii=ii+1;
          shenmup_row(p,q,r) =mynewdata2(p,q-1,r,ii); ii=ii+1;

          shenltemp_row(p,q,r) = ltempfloor;
          shentemp_row(p,q,r) = 10.^(shenltemp_row(p,q,r));

        elseif((q==1)&(r>1))
          
          % nc items
          ii=1;
          shenlrhob_row(p,q,r) =mynewdata3(p,q,r-1,ii); ii=ii+1;
          shennb_row(p,q,r) =mynewdata3(p,q,r-1,ii); ii=ii+1;
          shenlyp_row(p,q,r) =mynewdata3(p,q,r-1,ii); ii=ii+1;
          shenyp_row(p,q,r) =mynewdata3(p,q,r-1,ii); ii=ii+1;
          shenf_row(p,q,r) =mynewdata3(p,q,r-1,ii); ii=ii+1;
          shenebulk_row(p,q,r) =mynewdata3(p,q,r-1,ii); ii=ii+1;
          shensbulk_row(p,q,r) =mynewdata3(p,q,r-1,ii); ii=ii+1;
          shenaheav_row(p,q,r) =mynewdata3(p,q,r-1,ii); ii=ii+1;
          shenzheav_row(p,q,r) =mynewdata3(p,q,r-1,ii); ii=ii+1;
          shenmstar_row(p,q,r) =mynewdata3(p,q,r-1,ii); ii=ii+1;
          shenxneut_row(p,q,r) =mynewdata3(p,q,r-1,ii); ii=ii+1;
          shenxprot_row(p,q,r) =mynewdata3(p,q,r-1,ii); ii=ii+1;
          shenxalfa_row(p,q,r) =mynewdata3(p,q,r-1,ii); ii=ii+1;
          shenxh_row(p,q,r) =mynewdata3(p,q,r-1,ii); ii=ii+1;
          shenpbulk_row(p,q,r) =mynewdata3(p,q,r-1,ii); ii=ii+1;
          shenmun_row(p,q,r) =mynewdata3(p,q,r-1,ii); ii=ii+1;
          shenmup_row(p,q,r) =mynewdata3(p,q,r-1,ii); ii=ii+1;

          %       reset yp
          shenlyp_row(p,q,r) = lYpfloor;
          shenyp_row(p,q,r) = 10.^(shenlyp_row(p,q,r));

          shenltemp_row(p,q,r) = -1 + (r-2).*(2+1)./(ntk-2);
          shentemp_row(p,q,r) = 10.^(shenltemp_row(p,q,r));
          
        elseif((q>1)&(r>1))

          % nc items
          ii=1;
          shenlrhob_row(p,q,r) =mynewdata1(p,q-1,r-1,ii); ii=ii+1;
          shennb_row(p,q,r) =mynewdata1(p,q-1,r-1,ii); ii=ii+1;
          shenlyp_row(p,q,r) =mynewdata1(p,q-1,r-1,ii); ii=ii+1;
          shenyp_row(p,q,r) =mynewdata1(p,q-1,r-1,ii); ii=ii+1;
          shenf_row(p,q,r) =mynewdata1(p,q-1,r-1,ii); ii=ii+1;
          shenebulk_row(p,q,r) =mynewdata1(p,q-1,r-1,ii); ii=ii+1;
          shensbulk_row(p,q,r) =mynewdata1(p,q-1,r-1,ii); ii=ii+1;
          shenaheav_row(p,q,r) =mynewdata1(p,q-1,r-1,ii); ii=ii+1;
          shenzheav_row(p,q,r) =mynewdata1(p,q-1,r-1,ii); ii=ii+1;
          shenmstar_row(p,q,r) =mynewdata1(p,q-1,r-1,ii); ii=ii+1;
          shenxneut_row(p,q,r) =mynewdata1(p,q-1,r-1,ii); ii=ii+1;
          shenxprot_row(p,q,r) =mynewdata1(p,q-1,r-1,ii); ii=ii+1;
          shenxalfa_row(p,q,r) =mynewdata1(p,q-1,r-1,ii); ii=ii+1;
          shenxh_row(p,q,r) =mynewdata1(p,q-1,r-1,ii); ii=ii+1;
          shenpbulk_row(p,q,r) =mynewdata1(p,q-1,r-1,ii); ii=ii+1;
          shenmun_row(p,q,r) =mynewdata1(p,q-1,r-1,ii); ii=ii+1;
          shenmup_row(p,q,r) =mynewdata1(p,q-1,r-1,ii); ii=ii+1;

          shenltemp_row(p,q,r) = -1 + (r-2).*(2+1)./(ntk-2);
          shentemp_row(p,q,r) = 10.^(shenltemp_row(p,q,r));

          
        else

          % average to get corner

          % nc items
          ii=1;
          shenlrhob_row(p,q,r) =0.5.*(mynewdata3(p,q,r,ii) + mynewdata2(p,q,r,ii)); ii=ii+1;
          shennb_row(p,q,r) =0.5.*(mynewdata3(p,q,r,ii) + mynewdata2(p,q,r,ii)); ii=ii+1;
          shenlyp_row(p,q,r) =0.5.*(mynewdata3(p,q,r,ii) + mynewdata2(p,q,r,ii)); ii=ii+1;
          shenyp_row(p,q,r) =0.5.*(mynewdata3(p,q,r,ii) + mynewdata2(p,q,r,ii)); ii=ii+1;
          shenf_row(p,q,r) =0.5.*(mynewdata3(p,q,r,ii) + mynewdata2(p,q,r,ii)); ii=ii+1;
          shenebulk_row(p,q,r) =0.5.*(mynewdata3(p,q,r,ii) + mynewdata2(p,q,r,ii)); ii=ii+1;
          shensbulk_row(p,q,r) =0.5.*(mynewdata3(p,q,r,ii) + mynewdata2(p,q,r,ii)); ii=ii+1;
          shenaheav_row(p,q,r) =0.5.*(mynewdata3(p,q,r,ii) + mynewdata2(p,q,r,ii)); ii=ii+1;
          shenzheav_row(p,q,r) =0.5.*(mynewdata3(p,q,r,ii) + mynewdata2(p,q,r,ii)); ii=ii+1;
          shenmstar_row(p,q,r) =0.5.*(mynewdata3(p,q,r,ii) + mynewdata2(p,q,r,ii)); ii=ii+1;
          shenxneut_row(p,q,r) =0.5.*(mynewdata3(p,q,r,ii) + mynewdata2(p,q,r,ii)); ii=ii+1;
          shenxprot_row(p,q,r) =0.5.*(mynewdata3(p,q,r,ii) + mynewdata2(p,q,r,ii)); ii=ii+1;
          shenxalfa_row(p,q,r) =0.5.*(mynewdata3(p,q,r,ii) + mynewdata2(p,q,r,ii)); ii=ii+1;
          shenxh_row(p,q,r) =0.5.*(mynewdata3(p,q,r,ii) + mynewdata2(p,q,r,ii)); ii=ii+1;
          shenpbulk_row(p,q,r) =0.5.*(mynewdata3(p,q,r,ii) + mynewdata2(p,q,r,ii)); ii=ii+1;
          shenmun_row(p,q,r) =0.5.*(mynewdata3(p,q,r,ii) + mynewdata2(p,q,r,ii)); ii=ii+1;
          shenmup_row(p,q,r) =0.5.*(mynewdata3(p,q,r,ii) + mynewdata2(p,q,r,ii)); ii=ii+1;

          shenlyp_row(p,q,r) = lYpfloor;
          shenyp_row(p,q,r) = 10.^(shenlyp_row(p,q,r));

          shenltemp_row(p,q,r) = ltempfloor;
          shentemp_row(p,q,r) = 10.^(shenltemp_row(p,q,r));

        end

        
        if(shenaheav_row(p,q,r)<0.0)
          fprintf(fid,'Problem with aheav upon input: %d %d %d\n',p,q,r);
        end
        if(shenzheav_row(p,q,r)<0.0)
          fprintf(fid,'Problem with zheav upon input: %d %d %d\n',p,q,r);
        end
        if(shenxneut_row(p,q,r)<0.0)
          fprintf(fid,'Problem with xneut upon input: %d %d %d\n',p,q,r);
        end
        if(shenxprot_row(p,q,r)<0.0)
          fprintf(fid,'Problem with xprot upon input: %d %d %d\n',p,q,r);
        end
        if(shenxalfa_row(p,q,r)<0.0)
          fprintf(fid,'Problem with xalfa upon input: %d %d %d\n',p,q,r);
        end
        if(shenxh_row(p,q,r)<0.0)
          fprintf(fid,'Problem with xh upon input: %d %d %d\n',p,q,r);
        end

        
      end
    end
  end


  

  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Output data ranges and sizes (nrhobout, nypout, and ntkout should be
  % consistent with vector_sheneos.dek)
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  ncout=nc+2;
  
 

  if(useinterpn==1)
    % then no point in expanding physical grid size
    
    % since using linear interpolation in jon_shen.f in helmeoscode, user more points to more accurate
    %  nrhobout=200;
    nrhobout=nrhob*2;
    %  rhobminout=1.e2; % this is Kaz-definition of \rho_b using mb
    rhobminout=10.^(lrhobminin); % this is Kaz-definition of \rho_b using mb
    rhobmaxout=10.^(lrhobmaxin); % this is Kaz-definition of \rho_b using mb
    
    lrhobminout=log10(rhobminout);
    lrhobmaxout=log10(rhobmaxout);
    lrhobstepout=(lrhobmaxout-lrhobminout)./(nrhobout-1);

    %  nypout=20;
    nypout=nyp*2;
    %  ypminout=1.e-2;
    %  ypmaxout=10.^(-.25);
    ypminout=10.^(lypminin);
    ypmaxout=10.^(lypmaxin);
    
    lypminout=log10(ypminout);
    lypmaxout=log10(ypmaxout);
    lypstepout=(lypmaxout-lypminout)./(nypout-1);
    
    %  ntkout=200;
    %  ntkout=ntk;
    ntkout=2*2*ntk; % -1 to 2 in logMeV originally, extended 2 more orders, to
                    % 5 orders total, so 5/3 bigger grid needed to resolve same
                    % way, just choose 2X bigger then
    tkminout=1.e7/mev2K; % go down a bit since we do have T=0 solution to interpolate to
                         %  tkmaxout=1.e13/mev2K;
    tkminout=10.^(ltkminin);
    tkmaxout=10.^(ltkmaxin);
    
    ltkminout=log10(tkminout);
    ltkmaxout=log10(tkmaxout);
    ltkstepout=(ltkmaxout-ltkminout)./(ntkout-1);
  else
    % then expand grid
    
    % since using linear interpolation in jon_shen.f in helmeoscode, user more points to more accurate
    %  nrhobout=200;
    nrhobout=nrhob*2;
    rhobminout=10.^(lrhobminin-0.5); % this is Kaz-definition of \rho_b using mb
    rhobmaxout=10.^(lrhobmaxin+0.5); % this is Kaz-definition of \rho_b using mb
    
    lrhobminout=log10(rhobminout);
    lrhobmaxout=log10(rhobmaxout);
    lrhobstepout=(lrhobmaxout-lrhobminout)./(nrhobout-1);

    %  nypout=20;
    nypout=nyp*2;
    %  ypminout=1.e-2;
    %  ypmaxout=10.^(-.25);
    ypminout=10.^(lypminin-0.1);
    ypmaxout=10.^(lypmaxin+0.1);
    
    lypminout=log10(ypminout);
    lypmaxout=log10(ypmaxout);
    lypstepout=(lypmaxout-lypminout)./(nypout-1);
    
    %  ntkout=200;
    %  ntkout=ntk;
    ntkout=2*2*ntk; % -1 to 2 in logMeV originally, extended 2 more orders, to
                    % 5 orders total, so 5/3 bigger grid needed to resolve same
                    % way, just choose 2X bigger then
    tkminout=1.e7/mev2K; % go down a bit since we do have T=0 solution to interpolate to
                         %  tkmaxout=1.e13/mev2K;
    tkminout=10.^(ltkminin-0.5);
    tkmaxout=10.^(ltkmaxin+0.5);
    
    ltkminout=log10(tkminout);
    ltkmaxout=log10(tkmaxout);
    ltkstepout=(ltkmaxout-ltkminout)./(ntkout-1);    
  end
  
  

  [xout,yout,zout] = ndgrid(lrhobminout:lrhobstepout:lrhobmaxout,lypminout:lypstepout:lypmaxout,ltkminout:ltkstepout:ltkmaxout);

  % test of fortran code
  % fortran test
  %lrhob = log10(rhob)
  %   rhobi = 1 + (lrhob - lrhobminout).*(nrhobout-1)./(lrhobmaxout-lrhobminout)


  %%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Output header for Shen EOS to be read-in by HELM code
  %
  %%%%%%%%%%%%%%%%%%%%%%%
  fid=fopen(file5,'w');
  fprintf(fid,'%d\n%d\n%d\n%d %d %d\n%d %d %d\n%21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n%21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n', ...
          useinterpn, ...
          nc, ...
          ncout, ...
          nrhob, nyp, ntk, ...
          nrhobout, nypout, ntkout, ...
          lrhobminin,lrhobmaxin,lypminin,lypmaxin,ltkminin,ltkmaxin, ...
          lrhobminout,lrhobmaxout,lypminout,lypmaxout,ltkminout,ltkmaxout ...
          );
  fprintf(fid,'\n');
  fclose(fid);

  



  
  %%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Enforce minimums on a,z,and x's
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%
  shenaheav_row(shenaheav_row<aheavtrust)=aheavtrust;
  shenzheav_row(shenzheav_row<zheavtrust)=zheavtrust;

  shenxneut_row(shenxneut_row<xmin)=xmin;
  shenxprot_row(shenxprot_row<xmin)=xmin;
  shenxalfa_row(shenxalfa_row<xmin)=xmin;
  shenxh_row(shenxh_row<xhtrust)=xhtrust;

  shennb_row(shennb_row<nbmin)=nbmin;
  shenmstar_row(shenmstar_row<mstarfloor)=mstarfloor;

  
  fprintf('Begin Interpolate Data\n');
  
  %%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Interpolate
  %
  % use ifake_row to indicate if really within table by using 'linear' that generates NaN outside range of data
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%
  
  fake_row=shenlrhob_row;


  if(useinterpn==1)
    
    fprintf('Begin Interpolate Data using interpn -- all dimensions.\n');

    iiifake_row  =      interpn(xin,yin,zin,fake_row              ,xout,yout,zout,fakeinterpntype);
    
    iiishenlrhob_row  =      interpn(xin,yin,zin,shenlrhob_row              ,xout,yout,zout,interpntype);
    iiishennb_row     = 10.^(interpn(xin,yin,zin,log10(shennb_row+nbmin)    ,xout,yout,zout,interpntype));
    iiishenlyp_row    =      interpn(xin,yin,zin,shenlyp_row                ,xout,yout,zout,interpntype);
    iiishenyp_row     = 10.^(interpn(xin,yin,zin,log10(shenyp_row+Ypfloor)  ,xout,yout,zout,interpntype));
    iiishenf_row      =      interpn(xin,yin,zin,shenf_row                  ,xout,yout,zout,interpntype);
    iiishenebulk_row  =      interpn(xin,yin,zin,shenebulk_row              ,xout,yout,zout,interpntype);
    iiishensbulk_row  =      interpn(xin,yin,zin,shensbulk_row              ,xout,yout,zout,interpntype);
    iiishenaheav_row  = 10.^(interpn(xin,yin,zin,log10(shenaheav_row) ,xout,yout,zout,interpntype));
    %  iiishenaheav_row  =      interpn(xin,yin,zin,shenaheav_row              ,xout,yout,zout,interpntype);
    iiishenzheav_row  = 10.^(interpn(xin,yin,zin,log10(shenzheav_row) ,xout,yout,zout,interpntype));
    %  iiishenzheav_row  =      interpn(xin,yin,zin,shenzheav_row              ,xout,yout,zout,interpntype);
    iiishenmstar_row  = 10.^(interpn(xin,yin,zin,log10(shenmstar_row+mstarfloor) ,xout,yout,zout,interpntype));
    iiishenxneut_row  = 10.^(interpn(xin,yin,zin,log10(shenxneut_row+xmin)  ,xout,yout,zout,interpntype));
    %  iiishenxneut_row  =      interpn(xin,yin,zin,shenxneut_row              ,xout,yout,zout,interpntype);
    iiishenxprot_row  = 10.^(interpn(xin,yin,zin,log10(shenxprot_row+xmin)  ,xout,yout,zout,interpntype));
    %  iiishenxprot_row  =      interpn(xin,yin,zin,shenxprot_row              ,xout,yout,zout,interpntype);
    iiishenxalfa_row  = 10.^(interpn(xin,yin,zin,log10(shenxalfa_row+xmin)  ,xout,yout,zout,interpntype));
    %  iiishenxalfa_row  =      interpn(xin,yin,zin,shenxalfa_row)             ,xout,yout,zout,interpntype);
    iiishenxh_row     = 10.^(interpn(xin,yin,zin,log10(shenxh_row+xmin)     ,xout,yout,zout,interpntype));
    %  iiishenxh_row     =      interpn(xin,yin,zin,shenxh_row                 ,xout,yout,zout,interpntype);
    iiishenpbulk_row  =      interpn(xin,yin,zin,shenpbulk_row              ,xout,yout,zout,interpntype);
    iiishenmunminusmup_row =  interpn(xin,yin,zin,shenmun_row-shenmup_row    ,xout,yout,zout,interpntypemu);
    iiishenmup_row    =      interpn(xin,yin,zin,shenmup_row                ,xout,yout,zout,interpntypemu);

    iiishenltemp_row  =      interpn(xin,yin,zin,shenltemp_row              ,xout,yout,zout,interpntype);
    iiishentemp_row   = 10.^(interpn(xin,yin,zin,log10(shentemp_row+tempfloor)  ,xout,yout,zout,interpntype));

    fprintf('End Interpolate Data using interpn -- all dimensions.\n');

  
  else

    
    % do multiple 1D interpolations so can use pchip to avoid crazy oscillations with cubic or spline but get higher-order interpolations if monotonic
    
    
    %%% idea is that we are expanding/interpolating in each direction separately step-by-step in what could be an order-dependent process
    
    fprintf('Begin Interpolate Data using interp1 - x-dimension.\n');
    % 1D's along rho (i.e. "x" variable and "p" iterator)
    % expanding in rho-direction.
    % xin(:,1,1) and xout(:,1,1) used since x doesn't vary in q,r.  Not required here, but required for next expansion for sure.
    % Stupid matlab needs ' on xout(:,1,1) due to need to be in list as single row, while other yout(1,:,1) and zout(1,1,:) cannot have the ' .
    
    myxin=squeeze(xin(:,1,1));
    myxout=squeeze(xout(:,1,1));
 
    for r=1:ntk
      for q=1:nyp

        ifake_row(:,q,r)  =  interp1(myxin(:),squeeze(fake_row(:,q,r)),myxout(:)',fakeinterp1type);

        ishenlrhob_row(:,q,r)  =      extrap1(myxin(:),squeeze(shenlrhob_row(:,q,r)              ),myxout(:),interp1type);
        ishennb_row(:,q,r)     = 10.^(extrap1(myxin(:),squeeze(log10(shennb_row(:,q,r)+nbmin)    ),myxout(:)',interp1type));
        ishenlyp_row(:,q,r)    =      extrap1(myxin(:),squeeze(shenlyp_row(:,q,r)                ),myxout(:)',interp1type);
        ishenyp_row(:,q,r)     = 10.^(extrap1(myxin(:),squeeze(log10(shenyp_row(:,q,r)+Ypfloor)  ),myxout(:)',interp1type));
        ishenf_row(:,q,r)      =      extrap1(myxin(:),squeeze(shenf_row(:,q,r)                  ),myxout(:)',interp1type);
        ishenebulk_row(:,q,r)  =      extrap1(myxin(:),squeeze(shenebulk_row(:,q,r)              ),myxout(:)',interp1type);
        ishensbulk_row(:,q,r)  =      extrap1(myxin(:),squeeze(shensbulk_row(:,q,r)              ),myxout(:)',interp1type);
        ishenaheav_row(:,q,r)  = 10.^(extrap1(myxin(:),squeeze(log10(shenaheav_row(:,q,r)) ),myxout(:)',interp1type));
        ishenzheav_row(:,q,r)  = 10.^(extrap1(myxin(:),squeeze(log10(shenzheav_row(:,q,r)) ),myxout(:)',interp1type));
        ishenmstar_row(:,q,r)  = 10.^(extrap1(myxin(:),squeeze(log10(shenmstar_row(:,q,r)+mstarfloor) ),myxout(:)',interp1type));
        ishenxneut_row(:,q,r)  = 10.^(extrap1(myxin(:),squeeze(log10(shenxneut_row(:,q,r)+xmin)  ),myxout(:)',interp1type));
        ishenxprot_row(:,q,r)  = 10.^(extrap1(myxin(:),squeeze(log10(shenxprot_row(:,q,r)+xmin)  ),myxout(:)',interp1type));
        ishenxalfa_row(:,q,r)  = 10.^(extrap1(myxin(:),squeeze(log10(shenxalfa_row(:,q,r)+xmin)  ),myxout(:)',interp1type));
        ishenxh_row(:,q,r)     = 10.^(extrap1(myxin(:),squeeze(log10(shenxh_row(:,q,r)+xmin)     ),myxout(:)',interp1type));
        ishenpbulk_row(:,q,r)  =      extrap1(myxin(:),squeeze(shenpbulk_row(:,q,r)              ),myxout(:)',interp1type);
        ishenmunminusmup_row(:,q,r) =  extrap1(myxin(:),squeeze(shenmun_row(:,q,r)-shenmup_row(:,q,r)    ),myxout(:)',interp1typemu);
        ishenmup_row(:,q,r)    =      extrap1(myxin(:),squeeze(shenmup_row(:,q,r)                ),myxout(:)',interp1typemu);

        ishenltemp_row(:,q,r)  =      extrap1(myxin(:),squeeze(shenltemp_row(:,q,r)              ),myxout(:)',interp1type);
        ishentemp_row(:,q,r)   = 10.^(extrap1(myxin(:),squeeze(log10(shentemp_row(:,q,r)+tempfloor)  ),myxout(:)',interp1type));

        
      end
    end

    fprintf('Begin Interpolate Data using interp1 - y-dimension.\n');
    % 1D's along yp (i.e. "y" variable and "q" iterator)
    % now use nrhobout since we already expanded that direction.  Still using old ntk.
    % Use yin(1,:,1) and yout(1,:,1) since other directions don't change value
    
    myyin=squeeze(yin(1,:,1));
    myyout=squeeze(yout(1,:,1));
 
    for r=1:ntk
      for p=1:nrhobout

        iifake_row(p,:,r)  =  interp1(myyin(:),squeeze(ifake_row(p,:,r)),myyout(:),fakeinterp1type);

        iishenlrhob_row(p,:,r)  =      extrap1(myyin(:),squeeze(ishenlrhob_row(p,:,r)              ),myyout(:),interp1type);
        iishennb_row(p,:,r)     = 10.^(extrap1(myyin(:),squeeze(log10(ishennb_row(p,:,r)+nbmin)    ),myyout(:),interp1type));
        iishenlyp_row(p,:,r)    =      extrap1(myyin(:),squeeze(ishenlyp_row(p,:,r)                ),myyout(:),interp1type);
        iishenyp_row(p,:,r)     = 10.^(extrap1(myyin(:),squeeze(log10(ishenyp_row(p,:,r)+Ypfloor)  ),myyout(:),interp1type));
        iishenf_row(p,:,r)      =      extrap1(myyin(:),squeeze(ishenf_row(p,:,r)                  ),myyout(:),interp1type);
        iishenebulk_row(p,:,r)  =      extrap1(myyin(:),squeeze(ishenebulk_row(p,:,r)              ),myyout(:),interp1type);
        iishensbulk_row(p,:,r)  =      extrap1(myyin(:),squeeze(ishensbulk_row(p,:,r)              ),myyout(:),interp1type);
        iishenaheav_row(p,:,r)  = 10.^(extrap1(myyin(:),squeeze(log10(ishenaheav_row(p,:,r)) ),myyout(:),interp1type));
        iishenzheav_row(p,:,r)  = 10.^(extrap1(myyin(:),squeeze(log10(ishenzheav_row(p,:,r)) ),myyout(:),interp1type));
        iishenmstar_row(p,:,r)  = 10.^(extrap1(myyin(:),squeeze(log10(ishenmstar_row(p,:,r)+mstarfloor) ),myyout(:),interp1type));
        iishenxneut_row(p,:,r)  = 10.^(extrap1(myyin(:),squeeze(log10(ishenxneut_row(p,:,r)+xmin)  ),myyout(:),interp1type));
        iishenxprot_row(p,:,r)  = 10.^(extrap1(myyin(:),squeeze(log10(ishenxprot_row(p,:,r)+xmin)  ),myyout(:),interp1type));
        iishenxalfa_row(p,:,r)  = 10.^(extrap1(myyin(:),squeeze(log10(ishenxalfa_row(p,:,r)+xmin)  ),myyout(:),interp1type));
        iishenxh_row(p,:,r)     = 10.^(extrap1(myyin(:),squeeze(log10(ishenxh_row(p,:,r)+xmin)     ),myyout(:),interp1type));
        iishenpbulk_row(p,:,r)  =      extrap1(myyin(:),squeeze(ishenpbulk_row(p,:,r)              ),myyout(:),interp1type);
        iishenmunminusmup_row(p,:,r) =  extrap1(myyin(:),squeeze(ishenmunminusmup_row(p,:,r)       ),myyout(:),interp1typemu);
        iishenmup_row(p,:,r)    =      extrap1(myyin(:),squeeze(ishenmup_row(p,:,r)                ),myyout(:),interp1typemu);

        iishenltemp_row(p,:,r)  =      extrap1(myyin(:),squeeze(ishenltemp_row(p,:,r)              ),myyout(:),interp1type);
        iishentemp_row(p,:,r)   = 10.^(extrap1(myyin(:),squeeze(log10(ishentemp_row(p,:,r)+tempfloor)  ),myyout(:),interp1type));

        
      end
    end

    fprintf('Begin Interpolate Data using interp1 - z-dimension.\n');
    % 1D's along tk (i.e. "z" variable and "r" iterator)
    % now use nrhobout and nypout, but internally ":" for the inputs is still old ntk in size
    % Use zin(1,:,1) and zout(1,:,1) since other directions don't change value
    
    myzin=squeeze(zin(1,1,:));
    myzout=squeeze(zout(1,1,:));
    
    % Stupid Matlab doesn't like certain "implied formats" in interp1.

    for q=1:nypout
      for p=1:nrhobout

        iiifake_row(p,q,:)  =  interp1(myzin(:),squeeze(iifake_row(p,q,:)),myzout(:),fakeinterp1type);

        iiishenlrhob_row(p,q,:)  =      extrap1(myzin(:),squeeze(iishenlrhob_row(p,q,:)              ),myzout(:),interp1type);
        iiishennb_row(p,q,:)     = 10.^(extrap1(myzin(:),squeeze(log10(iishennb_row(p,q,:)+nbmin)    ),myzout(:),interp1type));
        iiishenlyp_row(p,q,:)    =      extrap1(myzin(:),squeeze(iishenlyp_row(p,q,:)                ),myzout(:),interp1type);
        iiishenyp_row(p,q,:)     = 10.^(extrap1(myzin(:),squeeze(log10(iishenyp_row(p,q,:)+Ypfloor)  ),myzout(:),interp1type));
        iiishenf_row(p,q,:)      =      extrap1(myzin(:),squeeze(iishenf_row(p,q,:)                  ),myzout(:),interp1type);
        iiishenebulk_row(p,q,:)  =      extrap1(myzin(:),squeeze(iishenebulk_row(p,q,:)              ),myzout(:),interp1type);
        iiishensbulk_row(p,q,:)  =      extrap1(myzin(:),squeeze(iishensbulk_row(p,q,:)              ),myzout(:),interp1type);
        iiishenaheav_row(p,q,:)  = 10.^(extrap1(myzin(:),squeeze(log10(iishenaheav_row(p,q,:)) ),myzout(:),interp1type));
        iiishenzheav_row(p,q,:)  = 10.^(extrap1(myzin(:),squeeze(log10(iishenzheav_row(p,q,:)) ),myzout(:),interp1type));
        iiishenmstar_row(p,q,:)  = 10.^(extrap1(myzin(:),squeeze(log10(iishenmstar_row(p,q,:)+mstarfloor) ),myzout(:),interp1type));
        iiishenxneut_row(p,q,:)  = 10.^(extrap1(myzin(:),squeeze(log10(iishenxneut_row(p,q,:)+xmin)  ),myzout(:),interp1type));
        iiishenxprot_row(p,q,:)  = 10.^(extrap1(myzin(:),squeeze(log10(iishenxprot_row(p,q,:)+xmin)  ),myzout(:),interp1type));
        iiishenxalfa_row(p,q,:)  = 10.^(extrap1(myzin(:),squeeze(log10(iishenxalfa_row(p,q,:)+xmin)  ),myzout(:),interp1type));
        iiishenxh_row(p,q,:)     = 10.^(extrap1(myzin(:),squeeze(log10(iishenxh_row(p,q,:)+xmin)     ),myzout(:),interp1type));
        iiishenpbulk_row(p,q,:)  =      extrap1(myzin(:),squeeze(iishenpbulk_row(p,q,:)              ),myzout(:),interp1type);
        iiishenmunminusmup_row(p,q,:) =  extrap1(myzin(:),squeeze(iishenmunminusmup_row(p,q,:)       ),myzout(:),interp1typemu);
        iiishenmup_row(p,q,:)    =      extrap1(myzin(:),squeeze(iishenmup_row(p,q,:)                ),myzout(:),interp1typemu);

        iiishenltemp_row(p,q,:)  =      extrap1(myzin(:),squeeze(iishenltemp_row(p,q,:)              ),myzout(:),interp1type);
        iiishentemp_row(p,q,:)   = 10.^(extrap1(myzin(:),squeeze(log10(iishentemp_row(p,q,:)+tempfloor)  ),myzout(:),interp1type));

        
      end
    end

    fprintf('End Interpolate Data using interp1 - x,y,z-dimensions.\n');

  end
  
  fprintf('End Interpolate Data\n');
    
  %%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Enforce minimums on a,z,and x's
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%
  iiishenaheav_row(iiishenaheav_row<aheavtrust)=aheavtrust;
  iiishenzheav_row(iiishenzheav_row<zheavtrust)=zheavtrust;

  iiishenxneut_row(iiishenxneut_row<xmin)=xmin;
  iiishenxprot_row(iiishenxprot_row<xmin)=xmin;
  iiishenxalfa_row(iiishenxalfa_row<xmin)=xmin;
  iiishenxh_row(iiishenxh_row<xhtrust)=xhtrust;

  iiishennb_row(iiishennb_row<nbmin)=nbmin;
  iiishenmstar_row(iiishenmstar_row<mstarfloor)=mstarfloor;
  

  
  %%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Enforce conditions on interpolation since interpolation of a and z and x's separately is unconstrained
  %
  % Note that "xout" contains non-interpolated actually-desired value of lrhob, so use this as reference value
  % Note that "yout" contains non-interpolated actually-desired value of lyp, so use this as reference value
  % Note that "yout" contains non-interpolated actually-desired value of ltk, so use this as reference value
  %
  % Use these instead of ishenlrhob_row, ishenyp_row, and ishenltemp_row that have interpolation errors
  % Use true n_b instead of ishennb_row
  %
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%
  fprintf('Table Indicator\n');

  truelrhob = xout;
  truerhob = 10.^(truelrhob);
  % below consistent with how Kaz defines \rho_b
  truenb = truerhob./mb;
  shenrhob = truenb*amu;
  truelyp = yout;
  trueyp = 10.^(truelyp);
  trueltk = zout;
  truetk = 10.^(trueltk);

  
  %  Create table indicator
  withintable_row = ...
      ~(~isfinite(iiifake_row) ...
        +~isfinite(iiishenlrhob_row)+~isfinite(iiishennb_row) ...
        +~isfinite(iiishennb_row)+~isfinite(iiishenlyp_row) ...
        +~isfinite(iiishenf_row)+~isfinite(iiishenebulk_row) ...
        +~isfinite(iiishensbulk_row)+~isfinite(iiishenaheav_row) ...
        +~isfinite(iiishenzheav_row)+~isfinite(iiishenmstar_row) ...
        +~isfinite(iiishenxneut_row)+~isfinite(iiishenxprot_row) ...
        +~isfinite(iiishenxalfa_row)+~isfinite(iiishenxh_row) ...
        +~isfinite(iiishenpbulk_row)+~isfinite(iiishenmunminusmup_row) ...
        +~isfinite(iiishenmup_row)+~isfinite(iiishenltemp_row) ...
        +~isfinite(iiishentemp_row));


  %%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % 
  %
  % For each state enforce reasonable consistency on unconstrained interpolated values
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%
  fprintf('Enforce\n');

  xcheckproblems=0;
  yeconsistentproblems=0;
  totalpoints=0;
  totalwithintable=0;
  
  for r=1:ntkout
    for q=1:nypout
      for p=1:nrhobout

        totalpoints=totalpoints+1;

        % only check inside normal table
        if(withintable_row(p,q,r)==1)
          
          totalwithintable=totalwithintable+1;
          
          %%%%%%%%%%%%%%%%%%%%
          % Setup variables for functions
          yetrue = trueyp(p,q,r);
          xnut = iiishenxneut_row(p,q,r);
          xprot = iiishenxprot_row(p,q,r);
          xalfa = iiishenxalfa_row(p,q,r);
          xh = iiishenxh_row(p,q,r);
          a     = iiishenaheav_row(p,q,r);
          x =  iiishenzheav_row(p,q,r)./iiishenaheav_row(p,q,r);

          
          
          % output overwrites input
          [xconsistent xnut xprot xalfa xh a x] = enforce_x_consistency(yetrue,xnut,xprot,xalfa,xh,a,x);

          if(xconsistent==0)
            fprintf('STILL problem with x(1): %d %d %d\n',p,q,r);
            fprintf('STILL problem with x(2): %g %g %g %g %g %g %g\n',yetrue,xnut,xprot,xalfa,xh,a,x);
            xcheckproblems=xcheckproblems+1;
          end

          
          % output overwrites input
          [yeconsistent  xnut xprot xalfa xh a x] = enforce_ye_consistency(yetrue,xnut,xprot,xalfa,xh,a,x);

          [abarnum abar zbar yelocal] = compute_nuclear_azbar(xnut,xprot,xalfa,xh,a,x);
          errorye = abs(yelocal-yetrue)/(abs(yelocal)+abs(yetrue));

          % same 100.0*yetolerance as in jon_shen.f, although for
          % different reasons.  Above we *always* do enforce of ye, while
          % in Fortran code we only do if errorye<yetolerance, hence
          % don't want to have machine error get triggered below.  Here
          % we just duplicate code section rather than idea
          if(errorye>100.0*yetolerance)
            fprintf('STILL problem with Y_e(1): (yeconsistent=%d, errorye=%g): %d %d %d\n',yeconsistent,errorye,p,q,r);
            fprintf('STILL problem with Y_e(2): %g %g %g %g %g %g %g\n',yetrue,xnut,xprot,xalfa,xh,a,x);
            fprintf('STILL problem with Y_e(3): %g %g %g %g\n',abarnum,abar,zbar,yelocal);
            yeconsistentproblems=yeconsistentproblems+1;
          end
          
          
          %%%%%%%%%%%%%%%%
          % save results
          iiishenxneut_row(p,q,r)=xnut;
          iiishenxprot_row(p,q,r)=xprot;
          iiishenxalfa_row(p,q,r)=xalfa;
          iiishenxh_row(p,q,r)=xh;
          iiishenaheav_row(p,q,r)=a;
          iiishenzheav_row(p,q,r)=a*x;
          
        end
        
        
      end
    end
  end
  
  fprintf('Problem count: BAD xcheck=%d (frac=%g) NOT yeconsistent=%d (frac=%g) totalwithintable=%d totalpoints=%d\n',xcheckproblems,xcheckproblems/totalwithintable,yeconsistentproblems,yeconsistentproblems/totalwithintable,totalwithintable,totalpoints);
  

  
  %%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Enforce minimums on a,z,and x's
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%
  shenaheav_row(shenaheav_row<aheavtrust)=aheavtrust;
  shenzheav_row(shenzheav_row<zheavtrust)=zheavtrust;

  shenxneut_row(shenxneut_row<xmin)=xmin;
  shenxprot_row(shenxprot_row<xmin)=xmin;
  shenxalfa_row(shenxalfa_row<xmin)=xmin;
  shenxh_row(shenxh_row<xhtrust)=xhtrust;

  shennb_row(shennb_row<nbmin)=nbmin;
  shenmstar_row(shenmstar_row<mstarfloor)=mstarfloor;

  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % 
  %
  % Thermodynamic consistency relation #3 in Shen EOS guide.ps
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%

  % assume error is in F
  iiishenf_row = iiishenebulk_row - truetk.*iiishensbulk_row + mumev - iiishenmstar_row;

  
  %%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % 
  %
  % Thermodynamic consistency relation #1 in Shen EOS guide.ps
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%

  if 0
    % causes problems at high T and low \rho for some unknown reason (assume f above absorbs all error)
    % assume error is in \mu_n and \mu_p equally
    % for now assume only in \mu_n for simplicity
    iiishenmunminusmup_row = (iiishenf_row + iiishenpbulk_row./truenb - iiishenmup_row.*trueyp)./(1.0-trueyp) - iiishenmup_row;
  end
  
  
  %
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%


  
  fprintf('2nd table indicator\n');
  
  %  re-create table indicator in case some quantities became NaN or
  %  infinity after checks
  %  Create table indicator
  withintable_row = ...
      ~(~isfinite(iiifake_row) ...
        +~isfinite(iiishenlrhob_row)+~isfinite(iiishennb_row) ...
        +~isfinite(iiishennb_row)+~isfinite(iiishenlyp_row) ...
        +~isfinite(iiishenf_row)+~isfinite(iiishenebulk_row) ...
        +~isfinite(iiishensbulk_row)+~isfinite(iiishenaheav_row) ...
        +~isfinite(iiishenzheav_row)+~isfinite(iiishenmstar_row) ...
        +~isfinite(iiishenxneut_row)+~isfinite(iiishenxprot_row) ...
        +~isfinite(iiishenxalfa_row)+~isfinite(iiishenxh_row) ...
        +~isfinite(iiishenpbulk_row)+~isfinite(iiishenmunminusmup_row) ...
        +~isfinite(iiishenmup_row)+~isfinite(iiishenltemp_row) ...
        +~isfinite(iiishentemp_row));
  

  finalwithintable=sum(sum(sum(withintable_row)));
  fprintf('finalwithintable=%d\n',finalwithintable);



  %%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Output interpolated Shen EOS (ncout=19 things)
  % use true values for grid positions
  %
  %%%%%%%%%%%%%%%%%%%%%%%
  fprintf('Write to File\n');
  fid=fopen(file4,'w');
  for r=1:ntkout
    for q=1:nypout
      for p=1:nrhobout
        fprintf(fid,'%3d %3d %3d %3d %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g', ...
                p, q, r, withintable_row(p,q,r), ...
        truelrhob(p,q,r)  ,truenb(p,q,r)  ,truelyp(p,q,r)  ,trueyp(p,q,r)  ,iiishenf_row(p,q,r)  ,iiishenebulk_row(p,q,r)  ,iiishensbulk_row(p,q,r)  ,iiishenaheav_row(p,q,r)  ,iiishenzheav_row(p,q,r)  ,iiishenmstar_row(p,q,r)  ,iiishenxneut_row(p,q,r)  ,iiishenxprot_row(p,q,r)  ,iiishenxalfa_row(p,q,r)  ,iiishenxh_row(p,q,r)  ,iiishenpbulk_row(p,q,r)  ,iiishenmunminusmup_row(p,q,r)  ,iiishenmup_row(p,q,r)  ,  trueltk(p,q,r) ,truetk(p,q,r) ...
                );
        fprintf(fid,'\n');
      end
    end
  end
  fclose(fid);

  %                iiishenlrhob_row(p,q,r)  ,iiishennb_row(p,q,r)  ,iiishenlyp_row(p,q,r)  ,iiishenyp_row(p,q,r)  ,iiishenf_row(p,q,r)  ,iiishenebulk_row(p,q,r)  ,iiishensbulk_row(p,q,r)  ,iiishenaheav_row(p,q,r)  ,iiishenzheav_row(p,q,r)  ,iiishenmstar_row(p,q,r)  ,iiishenxneut_row(p,q,r)  ,iiishenxprot_row(p,q,r)  ,iiishenxalfa_row(p,q,r)  ,iiishenxh_row(p,q,r)  ,iiishenpbulk_row(p,q,r)  ,iiishenmun_row(p,q,r)  ,iiishenmup_row(p,q,r)  ,  iiishenltemp_row(p,q,r) ,iiishentemp_row(p,q,r) ...

  
  

  

  fprintf('Done\n');


  fclose('all');

  % BELOW FOR relativity
  if onlinux
    quit;
  end


end





