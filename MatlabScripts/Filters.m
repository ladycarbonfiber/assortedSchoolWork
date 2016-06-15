%10.5
%High pass filter with specs
%passband .8pi - pi
%stopband(0,.75pi
%1% ripple in passband = 
%at least 46dB stopband attenuation

%butterworth
WS = .75;
WP = .8;
R = .086;
S = 46;
[n,Wn] = buttord(WP,WS,3,S) %n = 22, Wn = .7996
[z,p,k] = butter(n, Wn,'high');
sos = zp2sos(z,p,k);
fvtool(sos,'Analysis','freq')

%Chebyshev 1
[n,Wp] = cheb1ord(WP,WS,R,S)
[z,p,k]=cheby1(n,R,Wp,'high');%n = 11
sos = zp2sos(z,p,k);
fvtool(sos,'Analysis','freq')
%Elliptic
[n, Wp] = ellipord(WP, WS,R, S);
[z,p,k] = ellip(n,R,S,Wp,'high');%n=22
sos = zp2sos(z,p,k);
fvtool(sos,'Analysis','freq')