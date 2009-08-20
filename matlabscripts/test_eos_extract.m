%%%%%%%%%%% test of consolidator() + monotonize()

x = round(rand(1000,1)*100);
y = x+randn(size(x));
[xg,yg] = consolidator(x,y,@mean);

xmono = monotonize(x);
ymono = monotonize(y);

xgmono = monotonize(xg);
ygmono = monotonize(yg);

figure;
plot(x,y);
figure;
plot(xg,yg);
figure;
plot(xmono,ymono);
figure;
plot(xgmono,ygmono);

issorted(x)
issorted(xg)
issorted(xmono)
issorted(xgmono)

issorted(y)
issorted(yg)
issorted(ymono)
issorted(ygmono)



%%%%%%% Now test interpolation

% assume already monotonic, deal with later
x=0:0.01:3.14;
y1=x.^(1/2); % like U(T)
y2=cos(x);   % like P(T)

figure;
plot(y1,y2);

% chosen so small that only exact equality counts
CONTOL=1E-17;
CONTOL2=1E-17;

[y1x y2y] = consolidator(y1,y2,'mean',CONTOL2);

figure;
plot(y1x,y2y);

%newy1=-2:0.01:2; % like U
newy1=-2:0.1:2; % like U

% linear produces NaN at edges where plotting shows good fit!
y2ofy1 = interp1(y1x,y2y,newy1','linear');
figure;
plot(newy1,y2ofy1);

% spline is fine as long as using resolved-enough newy1
%y2ofy1 = interp1(y1x,y2y,newy1','spline');
%figure;
%plot(newy1,y2ofy1);

% cubic goes arbitrarily nuts at edges sometimes
%y2ofy1 = interp1(y1x,y2y,newy1','cubic');
%figure;
%plot(newy1,y2ofy1);

% cubic goes arbitrarily nuts at edges sometimes
%y2ofy1 = interp1(y1x,y2y,newy1','v5cubic');
%figure;
%plot(newy1,y2ofy1);

% spline is fine as long as using resolved-enough newy1
y2ofy1 = extrap1(y1x,y2y,newy1','spline');
figure;
plot(newy1,y2ofy1);


% spline is fine as long as using resolved-enough newy1
y2ofy1 = interpn(y1x,y2y,newy1','linear');
figure;
plot(newy1,y2ofy1);


%%END


f = @(x,y,z,t) t.*exp(-x.^2 - y.^2 - z.^2);
[x,y,z,t] = ndgrid(-1:0.2:1,-1:0.2:1,-1:0.2:1,0:2:10);
v = f(x,y,z,t);
[xi,yi,zi,ti] = ndgrid(-1:0.05:1,-1:0.08:1,-1:0.05:1,0:0.5:10);
vi = interpn(x,y,z,t,v,xi,yi,zi,ti,'spline');
nframes = size(ti, 4);
for j = 1:nframes
  slice(yi(:,:,:,j), xi(:,:,:,j), zi(:,:,:,j), vi(:,:,:,j),0,0,0);
  caxis([0 10]);
  M(j) = getframe;
end
movie(M); 



