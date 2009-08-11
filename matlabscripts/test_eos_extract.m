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
y1=x.^(1/2);
y2=cos(x);

figure;
plot(y1,y2);

% chosen so small that only exact equality counts
CONTOL=1E-17;
CONTOL2=1E-17;

[y1x y2y] = consolidator(y1,y2,'mean',CONTOL2);

figure;
plot(y1x,y2y);

newy1=-2:0.01:2;

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


%%END


