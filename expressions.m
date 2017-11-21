% Author: Nikhil Chandra Admal
%
% This matlab script is used to generate three .txt files: 
% 1) parameters.txt, 2) variables.txt and 3) pde.txt
% The above files subsequently used to run the comsol implementation of 
% the two-dimensional simulation of the polycrystal with grain boundary
% energy as described in the paper:
%
%   N. C. Admal, G. Po, J. Marian. ?A unified framework for polycrystal 
%   plasticity with grain boundary evolution.? 
%   arXiv preprint arXiv:1709.10176, 2017. 
%   Submitted to International Journal of Plasticity

% For comparison purpose, an equivalent kwc model is also included.
% Additional details on the usage of the above three files for the
% implementation of COMSOL are given in the README file. 

clear all
close all
clc

dim=2;
numSlips=2;

% kwc parameters
bphi     = 1;
mvmin    = 1e-9;
mvmax    = 1.0;
thetac   = 30;
alphasq  = 0.0053;
epsisq   = 2.1333e-4;
gamma    = 300;
sv       = 0.0017;
e0       = 0.0021;
delrhosq = 1e-9;

% parameter related to dissipative microstress
g       = 1;

% Elastic constants
mu      = 0.044776;
lam     = 0.09515;

fID=fopen('parameters.txt','w');
lineFormat='%s\n';

% slip 1
fprintf(fID,lineFormat,'s1_1 1' );
fprintf(fID,lineFormat,'s1_2 0' );
fprintf(fID,lineFormat,'n1_1 0' );
fprintf(fID,lineFormat,'n1_2 1' );

% slip 2
fprintf(fID,lineFormat,'s1_1 0' );
fprintf(fID,lineFormat,'s1_2 1' );
fprintf(fID,lineFormat,'n1_1 -1');
fprintf(fID,lineFormat,'n1_2 0' );

fprintf(fID,lineFormat,['bphi'     ' ' num2str(bphi)     '[fJ*s/(nm)^3]' ]);
fprintf(fID,lineFormat,['mvmax'    ' ' num2str(mvmax)    '[(nm)^3/(fJ*s)]'   ]);
fprintf(fID,lineFormat,['mvmin'    ' ' num2str(mvmin)    '[(nm)^3/(fJ*s)]'   ]);
fprintf(fID,lineFormat,['thetac'   ' ' num2str(thetac)   '[deg]'         ]);
fprintf(fID,lineFormat,['alphasq'  ' ' num2str(alphasq)  '[fJ/nm]'       ]);
fprintf(fID,lineFormat,['epsisq'   ' ' num2str(epsisq)   '[fJ/nm]'       ]);
fprintf(fID,lineFormat,['gamma'    ' ' num2str(gamma)    '[nm]'          ]);
fprintf(fID,lineFormat,['sv'       ' ' num2str(sv)       '[fJ/(nm)^2]'   ]);
fprintf(fID,lineFormat,['e0'       ' ' num2str(e0)       '[fJ/(nm)^3]'   ]);
fprintf(fID,lineFormat,['delrhosq' ' ' num2str(delrhosq) '[1/(nm)^2]'    ]);

fprintf(fID,lineFormat,['g'        ' ' num2str(g)        '[fJ*s/nm]'     ]);
fprintf(fID,lineFormat,['lam'      ' ' num2str(lam)      '[fJ/(nm)^3]'   ]);
fprintf(fID,lineFormat,['mu'       ' ' num2str(mu)       '[fJ/(nm)^3]'   ]);


fID=fopen('variables.txt','w');
lineFormat='%s\n';

%% u
u={'u1','u2'};
X={'x','y'};
fp={{'theta','Up11'},{'Up12','Up22'}};

%% w(rho)
syms w;
syms lam;
syms mu;
syms gphi;
syms phi;
syms phix;
syms phiy;
syms sv;
syms e0;
syms epsisq;
syms alphasq;
syms gamma;
syms delrhosq;


syms absrho;
syms rhosquare;

gphi = sv*phi^2;
fphi = e0*(1-phi)^2;

%% Define symbolic components of Ee
for bi=1:dim
    for bj=1:dim
        eval(['syms Ee' num2str(bi) num2str(bj) ' real']);
    end
end

%% Define symbolic components of G
for bi=1:3
    for bj=1:3
        eval(['syms G' num2str(bi) num2str(bj)]);
    end
end

%% Assemble the symbolic matrices Ee,dEe,G


Ee = [Ee11, Ee12, 0; 
      Ee21, Ee22, 0;
         0,    0, 0];
 



G = [0, 0, 0;
     0, 0, 0;
     G31, G32, 0];



%% H denotes lattice rotation which is a function of G, Ee and its derivative.
% But for this implementation, we have H=G

syms H31;
syms H32;

H=sym(zeros(3,3));
for bi=1:3
    for bj=1:3
        H(bi,bj) = G(bi,bj);
    end
end

varName='H31';
description='31 component of H';
expression = char(eval('H(3,1)'));
expression=expression(setdiff([1:length(expression)],strfind(expression,' ')));
fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);

varName='H32';
description='32 component of H';
expression = char(eval('H(3,2)'));
expression=expression(setdiff([1:length(expression)],strfind(expression,' ')));
fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);

% Derivative of H31 and H32 w.r.t G. Due to the assumption H=G these
% deriavatives turn out to be trivial


for bi=1:3
    for bj=1:3
        varName=['dH31dG' num2str(bi) num2str(bj)];
        description=['Derivative dH31/dG' num2str(bi) num2str(bj)];
        expression=[];
        expression = char(eval(['diff(H(3,1),G' num2str(bi) num2str(bj) ')']));
        expression=expression(setdiff([1:length(expression)],strfind(expression,' ')));
        fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);
    end
end

for bi=1:3
    for bj=1:3
        varName=['dH32dG' num2str(bi) num2str(bj)];
        description=['Derivative dH32/dG' num2str(bi) num2str(bj)];
        expression=[];
        expression = char(eval(['diff(H(3,2),G' num2str(bi) num2str(bj) ')']));
        expression=expression(setdiff([1:length(expression)],strfind(expression,' ')));
        fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);
    end
end

%% Energy w=w_el+w_theta+w_phi has three components

% w_theta is a function of H31 and H32

rhosquare = H31^2+H32^2+delrhosq;
absrho = sqrt(rhosquare);
w_theta=0.5*epsisq*rhosquare + gphi*(log(cosh(gamma*absrho))/gamma-log(cosh(gamma*sqrt(delrhosq)))/gamma);

% w_el is a function of Ee
a_el = mu + 0.5*(-mu-lam/4);
b_el = -0.5*mu-0.5*(-mu-lam/4);
c_el = 0.25*(-mu-lam/4 + mu + 3*lam/4);
d_el = 0.5*(mu + 3*lam/4 + mu+lam/4);

w_el=a_el*(3+2*trace(Ee))+...
     b_el*(3+4*trace(Ee)+2*trace(Ee)^2-2*trace(Ee*Ee))+...
     c_el*det(2*Ee+eye(3))-0.5*d_el*log(det(2*Ee+eye(3)))-...
     (3*a_el+3*b_el+c_el);
 
% w_phi is a function of phi and its gradient
w_phi=0.5*alphasq*(phix^2+phiy^2)+fphi;

w=w_theta+w_el+w_phi;

varName='absrho';
description=['Absoloute value of the gradient of rotation'];
expression = char(eval(absrho));
expression=expression(setdiff([1:length(expression)],strfind(expression,' ')));
fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);

varName='rhosquare';
description=['square of the absoloute value of the gradient of rotation'];
expression = char(eval(rhosquare));
expression=expression(setdiff([1:length(expression)],strfind(expression,' ')));
fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);


varName='w_el';
description=['Elastic energy'];
expression = char(eval(w_el));
expression=expression(setdiff([1:length(expression)],strfind(expression,' ')));
fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);

varName='w_phi';
description=['phi energy'];
expression = char(eval(w_phi));
expression=expression(setdiff([1:length(expression)],strfind(expression,' ')));
fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);

varName='w_theta';
description=['theta energy'];
expression = char(eval(w_theta));
expression=expression(setdiff([1:length(expression)],strfind(expression,' ')));
fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);

varName='w';
description=['total energy'];
expression = char(eval(w));
expression=expression(setdiff([1:length(expression)],strfind(expression,' ')));
fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);


%% mv: slip mobility

varName=['mv'];
description=['slip mobility'];
expression=['(mvmin+(1-phi^3*(10-15*phi+6*phi^2))*(mvmax-mvmin))'];
fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);


%% bv :  inverse phi mobility

varName=['bv'];
description=['inverse slip mobility'];
expression=['1/mv'];
fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);



%% Derivatives of energy
% dw/dH31, dw/dH32, dw/dEe, dw/dG, dw/dphi

varName='dwdH31';
description='dw/dH31';
expression=[];
expression = char(eval('diff(w,H31)'));
expression=expression(setdiff([1:length(expression)],strfind(expression,' ')));
fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);

varName='dwdH32';
description='dw/dH32';
expression=[];
expression = char(eval('diff(w,H32)'));
expression=expression(setdiff([1:length(expression)],strfind(expression,' ')));
fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);


% W (dw/dEe)
% W_{bi,bj} = d(w,Ee_{bi,bj})
for bi=1:dim
    for bj=1:dim
        varName=['W' num2str(bi) num2str(bj)];
        description=['intermediate stress, component ' num2str(bi) num2str(bj)];
        expression=[];
        expression = char(eval(['diff(w,Ee' num2str(bi) num2str(bj) ')']));
        expression=expression(setdiff([1:length(expression)],strfind(expression,' ')));

        fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);
    end
end


% tt=dw/dG
% tt_{bi,bj} = d(w,G_{bi,bj})
for bi=1:3
    for bj=1:3
        varName=['tt' num2str(bi) num2str(bj)];
        description=[' derivative  dw/dG' num2str(bi) num2str(bj)];
        expression=[];
        expression = ['dwdH31*dH31dG' num2str(bi) num2str(bj)...
                     '+dwdH32*dH32dG' num2str(bi) num2str(bj)];  
        fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);
    end
end


% dwdphi = dw/dphi
varName='dwdphi';
description='dw/dphi';
expression = char(eval(['diff(w,phi)']));
expression=expression(setdiff([1:length(expression)],strfind(expression,' ')));
fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);



%% Plastic rotation Rp

syms theta;


Rp = [cos(theta), -sin(theta);
      sin(theta), cos(theta)];

for bi=1:dim
    for I=1:dim
        varName=['Rp' num2str(bi) num2str(I)];
        description=['Rp, component ' num2str(bi) num2str(I)];
        expression = [char(Rp(bi,I))];
        fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);
    end
end

%% Plastic stretch Up
syms Up11;
syms Up12;
syms Up22;
Up = [Up11, Up12;
      Up12, Up22];


%% \dot Rp
dotRp = diff(Rp);

%% Fp = Rp Up
Fp = Rp*Up;
for bi=1:dim
    for I=1:dim
        varName=['Fp' num2str(bi) num2str(I)];
        description=['Fp, component ' num2str(bi) num2str(I)];
        expression = char(Fp(bi,I));
        expression=expression(setdiff([1:length(expression)],strfind(expression,' ')));
        fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);
    end
end
fprintf(fID,lineFormat,['Fp13' ' 0 ' 'Fp component 13']);
fprintf(fID,lineFormat,['Fp31' ' 0 ' 'Fp component 31']);
fprintf(fID,lineFormat,['Fp23' ' 0 ' 'Fp component 23']);
fprintf(fID,lineFormat,['Fp32' ' 0 ' 'Fp component 32']);
fprintf(fID,lineFormat,['Fp33' ' 1 ' 'Fp component 33']);

%%%%%%%%%

%% da
% The l.h.s of the flow rule \dot \Fp = Lp Fp is implemented as the matrix da 
% times the time-derivative of the unknown array [theta, Up11, Up12, Up22]^T
syms FpU11;
syms FpU12;
syms FpU22;
syms FpT;
da = [-sin(theta)*Up11-cos(theta)*Up12, cos(theta), -sin(theta), 0;
    -sin(theta)*Up12-cos(theta)*Up22, 0, cos(theta), -sin(theta);
    cos(theta)*Up11-sin(theta)*Up12, sin(theta), cos(theta),0;
    cos(theta)*Up12-sin(theta)*Up22, 0, sin(theta), cos(theta)];


for i=1:dim*dim
    for j=1:dim*dim
        varName=['da' num2str(i) num2str(j)];
        description=['da, component ' num2str(i) num2str(j)];
        expression=char(da(i,j));
        expression=expression(setdiff([1:length(expression)],findstr(expression,' ')));
        fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);
    end
end


%% Gp
for bi=1:dim
    for I=1:dim
        eval(['syms Fp' num2str(bi) num2str(I) ' real;']);
        eval(['Fps(' num2str(bi) ',' num2str(I) ')= Fp' num2str(bi) num2str(I) ';']);
    end
end
Jp=det(Fps);
expression=char(Jp);
expression=expression(setdiff([1:length(expression)],strfind(expression,' ')));
fprintf(fID,lineFormat,['Jp' ' ' expression ' ' 'determinant of Fp']);

Gp=inv(Fps);
for bi=1:dim
    for I=1:dim
        varName=['Gp' num2str(I) num2str(bi)];
        description=['inverse Fp, component ' num2str(I) num2str(bi)];
        expression=char(Gp(I,bi)*Jp);
        expression=expression(setdiff([1:length(expression)],strfind(expression,' ')));
        fprintf(fID,lineFormat,[varName ' (' expression ')/Jp ' description]);
    end
end

%% F
% F_{i,J} = d(u_i,X{J}) + delta_{i,J}
for i=1:dim
    for J=1:dim
        varName=['F' num2str(i) num2str(J)];
        description=['deformation gradient, component ' num2str(i) num2str(J)];
        expression=[u{i} X{J} '+' num2str(i==J)];
        fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);
    end
end


%% Fe = F*Gp
% Fe_{i,bi} = F_{i,J}*Gp_{J,bi}
for i=1:dim
    for bi=1:dim
        varName=['Fe' num2str(i) num2str(bi)];
        description=['elastic distortion, component ' num2str(i) num2str(bi)];
        expression=[];
        for J=1:dim
            expression=[expression '+F' num2str(i) num2str(J) '*Gp' num2str(J) num2str(bi)];
        end
        fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);
    end
end


%% Ce_{bi,bj}=Fe_{k,bi} Fe_{k,bj}
for bi=1:dim
    for bj=1:dim
        varName=['Ce' num2str(bi) num2str(bj)];
        description=['left CG stretch tensor, component ' num2str(bi) num2str(bj)];
        expression=[];
        for k=1:dim
            expression=[expression '+Fe' num2str(k) num2str(bi) '*Fe' num2str(k) num2str(bj)];
        end
        fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);
    end
end

%% Ee_{bi,bj} = 0.5*Ce_{bi,bj} - 0.5 delta_{bi,bj}

for bi=1:dim
    for bj=1:dim
        varName=['Ee' num2str(bi) num2str(bj)];
        description=['elastic strain, component ' num2str(bi) num2str(bj)];
        expression= ['0.5*Ce' num2str(bi) num2str(bj) '-0.5*' num2str(bi==bj)];
        fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);
    end
end


%% Piola P
% P_{k,I}=Fe_{k,bi}*W_{bi,bj}*Gp_{I,bj}
%        
for k=1:dim
    for I=1:dim
        varName=['P' num2str(k) num2str(I)];
        description=['first Piola stress, component ' num2str(k) num2str(I)];
        expression=[];
        for bi=1:dim
            for bj=1:dim
                expression=[expression '+Fe' num2str(k) num2str(bi) '*W' num2str(bi) num2str(bj) '*Gp' num2str(I) num2str(bj)];
            end
        end
        
        fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);
    end
end


%% Resolved shear stress
% rss = W_{bi,bj}*n_{bj}*Ce_{bi,bk}*s_{bk}
for slipID=1:numSlips
    varName=['rss' num2str(slipID)];
    description=['Resolved shear stress ' num2str(slipID)];
    expression=[];
    % W_{bi,bj}*n_{bj}*Ce_{bi,bk}*s_{bk}
    for bi=1:dim
        for bj=1:dim
            for bk=1:dim
                expression=[expression '+W' num2str(bi) num2str(bj) '*n' num2str(slipID) '_' num2str(bj) '*Ce' num2str(bi) num2str(bk) '*s' num2str(slipID) '_' num2str(bk)];
            end
        end
    end
    fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);
end




%% Lp
for bi=1:dim
    for bj=1:dim
        varName=['Lp' num2str(bi) num2str(bj)];
        description=['Lp, component ' num2str(bi) num2str(bj)];
        expression=[];
        for slipID=1:numSlips
            expression=[expression '+v' num2str(slipID) '*s' num2str(slipID) '_' num2str(bi) '*n' num2str(slipID) '_' num2str(bj)];
        end
        fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);
    end
end


%% LpFp
for bi=1:dim
    for I=1:dim
        varName=['LpFp' num2str(bi) num2str(I)];
        description=['LpFp, component ' num2str(bi) num2str(I)];
        expression=[];
        for bj=1:dim
            expression=[expression '+Lp' num2str(bi) num2str(bj) '*Fp' num2str(bj) num2str(I)];
        end
        fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);
    end
end

%% Xi_phi = dw/ddphi

for I=1:dim
    varName=['Xi_phi_' num2str(I)];
    description=[' component ' num2str(I) ' of the back stress Xi_phi_' num2str(slipID)];
    expression = char(eval(['diff(w,phi' X{I} ')']));
    expression=expression(setdiff([1:length(expression)],strfind(expression,' ')));
    fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);
end




%% Xi_en
% Xi_en_I= Jp*(Gp_{I,bj} * eps(bj,bp,bq)*n_{bp}*tt_{bq,br}*s_{br} 
for slipID=1:numSlips
    for I=1:dim
        varName=['Xi_en' num2str(slipID) '_' num2str(I)];
        description=[' component ' num2str(I) ' of the energetic back stress Xi_en' num2str(slipID)];
        expression=[];
        for bj=1:dim
            for bp=1:dim
                for bq=1:3
                    for br=1:dim
                        expression=[expression '+Gp' num2str(I) num2str(bj) ...
                                    '*' num2str(levCiv(bj,bp,bq))...
                                    '*n' num2str(slipID) '_' num2str(bp) ...
                                    '*tt' num2str(bq) num2str(br) ...
                                    '*s' num2str(slipID) '_' num2str(br)];
                    end
                end
            end
        end
        % Reverting to Gurtin to maintain consistency
        expression=['Jp*(' expression ')'];
        
        fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);
    end
end
  
%% Xi_diss
for slipID=1:numSlips
    for i=1:dim
        varName=['Xi_diss' num2str(slipID) '_' num2str(i)];
        description=[' component ' num2str(slipID) ' of the dissipative back stress Xi_diss' num2str(slipID)];
        expression=[];
  
        expression=['g*v' num2str(slipID) X{i}];

        fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);
    end
end

%% Xi
for slipID=1:numSlips
    for I=1:dim
        varName=['Xi' num2str(slipID) '_' num2str(I)];
        description=[' component ' num2str(I) ' of the total back stress Xi' num2str(slipID)];
        expression=['Xi_en' num2str(slipID) '_' num2str(I) '+Xi_diss'  num2str(slipID) '_' num2str(I)];
        fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);
    end
end

%% dotFp = \dot Rp * Up + Rp * \dot Up

for bi=1:dim
    for J=1:dim
        varName=['dotFp' num2str(bi) num2str(J)];
        description=['dotFp, component ' num2str(bi) num2str(J)];
        expression=[];
        for K=1:dim
            if (K<J)
                expression=[expression '+thetat*(' char(dotRp(bi,K)) ')*Up' num2str(K) num2str(J) ...
                                       '+Rp' num2str(bi) num2str(K) '*Up' num2str(K) num2str(J) 't'];
            else
                expression=[expression '+thetat*(' char(dotRp(bi,K)) ')*Up' num2str(J) num2str(K) ...
                                       '+Rp' num2str(bi) num2str(K) '*Up' num2str(J) num2str(K) 't'];
            end
        end
        fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);
    end
end


%% G 
% G = Fp_{bi,K} * eps(K,R,S)*d(Fp_{bj,S},R)
for bi=1:3
    for bj=1:3
        varName=['G' num2str(bi) num2str(bj)];
        description=['G, component ' num2str(bi) num2str(bj)];
        expression=[];
        for K=1:3
            for R=1:2
                for S=1:3
                    expression=[expression '+' num2str(levCiv(K,R,S)) '*d(Fp' num2str(bj) num2str(S) ',' X{R} ')*Fp' num2str(bi) num2str(K)];
                end
            end
        end
            fprintf(fID,lineFormat,[varName ' (' expression ')/Jp ' description]);
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%% kwc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

syms theta0x;
syms theta0y;
syms phi0x;
syms phi0y;
syms phi0;

gphi0 = sv*phi0^2;
fphi0 = e0*(1-phi0)^2;
ndtheta = sqrt(theta0x^2+theta0y^2+delrhosq);

% Approximation: abs(x) ~ ln(cosh(gamma*x)/gamma) 
w_kwc=gphi0*(log(cosh(gamma*ndtheta))/gamma-log(cosh(gamma*sqrt(delrhosq)))/gamma) + ...
      0.5*epsisq*(theta0x^2+theta0y^2)+...
      0.5*alphasq*(phi0x^2+phi0y^2)+...
      fphi0;

varName='w_kwc';
description=['total kwc energy'];
expression = char(eval(w_kwc));
expression=expression(setdiff([1:length(expression)],strfind(expression,' ')));
fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);

%% Xi_phi0 = dw_kwc/ddphi0

for I=1:dim
    varName=['Xi_phi0_' num2str(I)];
    description=[' component ' num2str(I) ' of the back stress Xi_phi0_' num2str(slipID)];
    expression = char(eval(['diff(w_kwc,phi0' X{I} ')']));
    expression=expression(setdiff([1:length(expression)],strfind(expression,' ')));
    fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);
end


%% Xi_theta0 = dw_kwc/ddtheta
for I=1:dim
    varName=['Xi_theta0_' num2str(I)];
    description=[' component ' num2str(I) ' of the back stress Xi_theta0_' num2str(slipID)];
    expression = char(eval(['diff(w_kwc,theta0' X{I} ')']));
    expression=expression(setdiff([1:length(expression)],strfind(expression,' ')));
    fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);
end



%% dw_kwcdphi0 = dw/dphi0
varName='dw_kwcdphi0';
description='dw_kwc/dphi0';
expression = char(eval(['diff(w_kwc,phi0)']));
expression=expression(setdiff([1:length(expression)],strfind(expression,' ')));
fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);

%% mtheta0: slip mobility

varName=['mtheta0'];
description=['theta0 mobility'];
expression=['(mvmin+(1-phi0^3*(10-15*phi0+6*phi0^2))*(mvmax-mvmin))'];
fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);

%% btheta0 :  inverse phi mobility

varName=['btheta0'];
description=['inverse theta0 mobility'];
expression=['1/mtheta0'];
fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);

%% bphi0 :  inverse phi mobility

varName=['bphi0'];
description=['inverse phi0 mobility'];
expression=['bphi'];
fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);

%% Check elastic constants

varName='d2wdH312';
description='d2w/d(H31)2';
expression=[];
expression = char(eval('diff(diff(w,H31),H31)'));
expression=expression(setdiff([1:length(expression)],strfind(expression,' ')));
fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);


%% Close the file

fclose(fID);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PDE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PDE for crystal plasticity

FID=fopen('pde.txt','w');


% P_{iJ} * test(ui,X{J})
for i=1:dim
    expression=[];
    for J=1:dim
        expression=[expression '+P' num2str(i) num2str(J) '*test(u' num2str(i) X{J} ')'];
    end
    fprintf(FID,[expression '\n']);
end

% 
for slipID=1:numSlips
    expression=[];
    for J=1:dim
        expression=[expression '+Xi' num2str(slipID) '_' num2str(J) '*test(v' num2str(slipID) X{J} ')'];
    end
        expression=[expression '-rss' num2str(slipID) '*test(v' num2str(slipID) ')'];
        expression=[expression '+(tt31*(G31*s' num2str(slipID) '_1*n' num2str(slipID) '_1'...
                                      '+G32*s' num2str(slipID) '_1*n' num2str(slipID) '_2)'...
                               '+tt32*(G31*s' num2str(slipID) '_2*n' num2str(slipID) '_1'...
                                     '+G32*s' num2str(slipID) '_2*n' num2str(slipID) '_2))*test(v' num2str(slipID) ')'];

        expression=[expression '+bv*v' num2str(slipID) '*test(v' num2str(slipID) ')'];
        

        fprintf(FID,[expression '\n']);
end

expression=[];
for j=1:dim
    expression=[expression '+Xi_phi_' num2str(j) '*test(phi' X{j} ')'];
end
expression=[expression '+(dwdphi+bphi*phit)*test(phi)'];
fprintf(FID,[expression '\n']);


% pde for kwc 

expression=[];
for j=1:dim
    expression=[expression '+Xi_theta0_' num2str(j) '*test(theta0' X{j} ')'];
end
expression=[expression '+btheta0*theta0t*test(theta0)'];
fprintf(FID,[expression '\n']);

expression=[];
for j=1:dim
    expression=[expression '+Xi_phi0_' num2str(j) '*test(phi0' X{j} ')'];
end
expression=[expression '+(dw_kwcdphi0+bphi0*phi0t)*test(phi0)'];
fprintf(FID,[expression '\n']);

fclose(fID);


