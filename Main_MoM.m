format long;
clear all;
% Constant
eps0 = 8.85E-12;      % Vacuum permittivity (F/m)
mu0 = 4*pi*1E-7;      % Vacuum permeability (H/m)
c = 1/sqrt(eps0*mu0); % vacuum light speed (m/s)

% Define incident wave properties
freq = 1E9;           % 1GHz frequency 
lambda = c/freq;      % wavelength (m)
k0 = 2*pi/lambda;     % wavevector (1/m)
amp = 1;              % amplitude
lenBox = 10 * lambda; % size of the simulation domain 

% Define the cross-sectional shape of the cylinderical scatterer
Z0 = 1;               % impedence
sides = 4;            % number of sides of the scatterer (4 stands for a square)   
sideN = 10;           % # of elements per side 
N = sideN * sides;    % total number of elements
len = 3 * lambda * sind(180/sides);  % side length  
dlen = len / sideN;   % length of each element
NBox = floor(lenBox/dlen)+1;  % total number of elements in the simulation domain
rot = 0;              % rotation (degree)

deg=[0:1:361];		  % generate a circle 
circlex=cosd(deg);
circley=sind(deg);
step=360/sides;       % sample points on the circle for the polygons
indiceX = circlex(1:step:361);
indiceY = circley(1:step:361);
indiceXY = [indiceX; indiceY;];
indiceXY = [ cosd(rot) -sind(rot) ; sind(rot) cosd(rot) ] * indiceXY;  % rotate the polygons
indiceX = indiceXY(1,:);
indiceY = indiceXY(2,:);


% Matrix Initialization 
coord = struct('xi', 0,'yi', 0,'xe', 0, 'ye', 0); % initialize a coordinate parameter for the scatterer
coordBox = struct('x', 0,'y', 0);                 % initialize a coordinate parameter for the simulation box
flagReg = zeros(NBox,NBox); 		% initialize a 2D matrix to distinguished the regions between the scatterer and vacuum
V = zeros(N,1); 				    % initialize a 1D array for the incident fields (E or H) 
Z = zeros(N,N);                     % initialize a 2D array for the interacting impedences in EFIE
J = zeros(N,1);						% initialize a 1D array for the current densities in EFIE			  
Einc = zeros(NBox,NBox);            % initialize a 2D array for the incident E-fields in EFIE
Esc  = zeros(NBox,NBox); 			% initialize a 2D array for the scattered E-fields in EFIE
Etot = zeros(NBox,NBox); 			% initialize a 2D array for the total E-fields in EFIE 
ZH = zeros(N,N);                    % initialize a 2D array for the interacting impedences in MFIE
JH   = zeros(N,1);                  % initialize a 1D array for the current densities in MFIE		
Hinc = zeros(NBox,NBox);            % initialize a 2D array for the incident H-fields in MFIE 
Hsc  = zeros(NBox,NBox); 			% initialize a 2D array for the scattered H-fields in MFIE 
Htot = zeros(NBox,NBox); 			% initialize a 2D array for the total H-fields in MFIE  

% Connect the sampling points on the perimeter of the scatterer into elements (2 points for 1 element)
global_index = 1; % the index of each element
for i = 1:sides
    xini = indiceX(i);
    yini = indiceY(i);
    xend = indiceX(i+1);
    yend = indiceY(i+1);
    len_ori = sqrt(abs(xini-xend)^2 + abs(yini-yend)^2 );
    
    tmpx = linspace(xini,xend,sideN+1) * len / len_ori;
    tmpy = linspace(yini,yend,sideN+1) * len / len_ori;
    
    for j = 1:sideN
        coord(global_index).xi = tmpx(j);
        coord(global_index).yi = tmpy(j);
        coord(global_index).xe = tmpx(j+1);
        coord(global_index).ye = tmpy(j+1);
        global_index = global_index + 1;
    end
end


% Initialize the coordinates in the simulation box
xcenter = ceil(NBox/2) * dlen  - dlen /2; 
ycenter = ceil(NBox/2) * dlen  - dlen /2;
for i = 1:NBox
    for j = 1:NBox
        x = j * dlen  - dlen /2; 
        y = i * dlen  - dlen /2;
        coordBox(i,j).x = x - xcenter; 
        coordBox(i,j).y = y - ycenter; 
    end
end

% A general algorithm to find out the region enclosed by the scatterer
for s = 1:N
    a = coord(s).yi - coord(s).ye;
    b = coord(s).xe - coord(s).xi;
    c = -a*coord(s).xi-b*coord(s).yi; 
    for i = 1:NBox
        for j = 1:NBox
            val = a*coordBox(i,j).x + b*coordBox(i,j).y + c; 
            if val >=0 
                flagReg(i,j) = flagReg(i,j) + 1 ;    % The enclosed regions should have the largest value
            end
        end
    end
end
maxval = max(max(flagReg));                         
for i = 1:NBox
        for j = 1:NBox
            if flagReg(i,j) == maxval
                flagReg(i,j) = 1 ;                   % The region is inside the scatterer
            else
                flagReg(i,j) = 0 ;                   % The region is outside the scatterer
            end
        end
end


% Matrix constructions for eqs. 2 and 3 (EFIE) and 7 and 8 (MFIE) 
coef = 1j * k0 * Z0;  
for m = 1:N       % construct the matrix element-by-element 
    xm = (coord(m).xi + coord(m).xe)/2;
    ym = (coord(m).yi + coord(m).ye)/2;
    for n = 1:N   % the interaction between the m-th and the n-th elements
        xn = (coord(n).xi + coord(n).xe)/2;
        yn = (coord(n).yi + coord(n).ye)/2;
        nvec = -[coord(n).yi-coord(n).ye,coord(n).xe-coord(n).xi];
        nvec = nvec/norm(nvec);
        
        % EFIE
        if m==n  % Small argument approximation is used 
            Z(m,n) = k0*Z0*dlen/4 * (1-1j*(2/pi)*log(k0*1.7180724*dlen/(4*exp(1))));
        else
            Z(m,n) = k0*Z0*dlen/4 * besselh(0,2, k0 * sqrt( (xn-xm).^2 + (yn-ym).^2 ));
        end  
        % MFIE
        if m==n  % Small argument approximation is used 
            ZH(m,n) = -1/2; 
        else
            distance = sqrt( (xn-xm).^2 + (yn-ym).^2 ); 
            disvec = [xm-xn, ym-yn];
            ZH(m,n) = k0*dlen/4j * besselh(1,2, k0 * distance) * dot(nvec,disvec)/distance;
        end  

    end
    V(m) = amp * exp(-1j*k0*xm);
end
% Solve the current densities on the PEC for EFIE and MFIE
J  = inv(Z) * V; 
JH = inv(ZH) * V;


% Once the current densities are determined
% E and H at arbitary points can be dobtained with minor efforts 
% Far-Field Calculations 
for i = 1:NBox
    for j = 1:NBox
        
        xm = coordBox(i,j).x;  % the (x,y) coordinate of a observation point  
        ym = coordBox(i,j).y;
        cum_intG0  = 0;
        cum_intG0H = 0;
        for n = 1:N       % calculate the contribution from each element
            xn = (coord(n).xi + coord(n).xe)/2;
            yn = (coord(n).yi + coord(n).ye)/2;
            nvec = -[coord(n).yi-coord(n).ye,coord(n).xe-coord(n).xi];
            nvec = nvec/norm(nvec);
            %EFIE (associated with eq. 4)
            intG0 = -k0*Z0*dlen/4 * besselh(0,2, k0 * sqrt( (xn-xm).^2 + (yn-ym).^2 ))*J(n);
            %MFIE (associated with eq. 9)
            distance = sqrt( (xn-xm).^2 + (yn-ym).^2 ); 
            disvec = [xm-xn, ym-yn];
            intG0H = -k0*dlen/4j * besselh(1,2, k0 * distance) * dot(nvec,disvec)/distance *JH(n);
            cum_intG0  = cum_intG0 + intG0; 
            cum_intG0H = cum_intG0H + intG0H; 
        end
		
        if flagReg(i,j) == 1  % If the observation point is inside PEC, then both E and H are 0
            Esc(i,j)  = 0; 
            Einc(i,j) = 0;
            Etot(i,j) = 0; 
            Hsc(i,j)  = 0; 
            Hinc(i,j) = 0;
            Htot(i,j) = 0;            
        else 
            Esc(i,j)  = cum_intG0; 
            Einc(i,j) = amp * exp(-1j*k0*xm);
            Etot(i,j) = Einc(i,j) + Esc(i,j); 
            Hsc(i,j)  = cum_intG0H; 
            Hinc(i,j) = amp * exp(-1j*k0*xm);
            Htot(i,j) = Hinc(i,j) + Hsc(i,j); 
        end
    end
end


% normalize the x and y scales with respect to the incident wavelength 
normxyrange = [ coordBox(1,1).x/lambda : dlen/lambda : coordBox(1,NBox).x/lambda]; 

% plot the shape of the scatterer enclosed by a circle
figure(1);
subplot(1,2,1);
plot(indiceX,indiceY);
hold on;
plot(circlex,circley);
axis([-1 1 -1 1]);
axis equal;
hold off;

% plot the shape of the scatterer
subplot(1,2,2);
title('(b) Region');
contourf(normxyrange,normxyrange,flagReg,1,'edgecolor','none');
cobj = colorbar('Ticks',[0,0.5],'TickLabels',{'Vacuum','Scatterer'});
hold on;
for i = 1:N
    plot([coord(i).xi/lambda coord(i).xe/lambda], [coord(i).yi/lambda coord(i).ye/lambda]);
end
axis([-5 5 -5 5]);
xticks([-5 : 1 : 5]);
yticks([-5 : 1 : 5]);
xlabel('x (\lambda)');
ylabel('y (\lambda)');
axis equal;
hold off;


% plot the current densities in the TM mode along the edge of the shape
fig = figure('Position',[0 0 1000 2000]);
subplot(3,2,1);
plot([dlen/lambda:dlen/lambda:len/lambda * sides],abs(J));
title('(a) J_{TM}');
axis([0 len/lambda * sides 0 max(abs(J))]);
xticks([0 : 2 : len/lambda * sides]);
yticks([0 : 0.5 : max(abs(J)) ]);
xlabel('s (\lambda)');
ylabel('J_{s} (A/m)');

% plot the current densities in the TE mode along the edge of the shape
subplot(3,2,2);
plot([dlen/lambda:dlen/lambda:len/lambda * sides],abs(JH));
title('(b) J_{TE}');
axis([0 len/lambda * sides 0 max(abs(JH))]);
xticks([0 : 2 : len/lambda * sides]);
yticks([0 : 0.5 : max(abs(JH)) ]);
xlabel('s (\lambda)');
ylabel('J_{s} (A/m)');

% plot the scattered E-fields in the TM mode
subplot(3,2,3);
title('(c) E_{sc}');
hold on;
contourf(normxyrange,normxyrange,abs(Esc),40,'edgecolor','none');
colormap jet;
cobj = colorbar;
for i = 1:N
    plot([coord(i).xi/lambda coord(i).xe/lambda], [coord(i).yi/lambda coord(i).ye/lambda]);
end
axis([-5 5 -5 5]);
caxis([0 max(max(abs(Esc)))]);
xticks([-5 : 1 : 5]);
yticks([-5 : 1 : 5]);
set(cobj,'YTick',[0:0.2:max(max(abs(Esc)))]);
xlabel('x (\lambda)');
ylabel('y (\lambda)');
cobj.Label.String = 'V/m';
axis equal;
hold off;

% plot the scattered H-fields in the TE mode
subplot(3,2,4);
title('(d) H_{sc}');
hold on;
contourf(normxyrange,normxyrange,abs(Hsc),40,'edgecolor','none');
colormap jet;
cobj = colorbar;
for i = 1:N
    plot([coord(i).xi/lambda coord(i).xe/lambda], [coord(i).yi/lambda coord(i).ye/lambda]);
end
axis([-5 5 -5 5]);
caxis([0 max(max(abs(Hsc)))]);
xticks([-5 : 1 : 5]);
yticks([-5 : 1 : 5]);
set(cobj,'YTick',[0: 0.5 :max(max(abs(Hsc)))]);
xlabel('x (\lambda)');
ylabel('y (\lambda)');
cobj.Label.String = 'A/m';
axis equal;
hold off;


% plot the total E-fields in the TM mode
subplot(3,2,5);
title('(e) E_{tot}');
hold on;
contourf(normxyrange,normxyrange,abs(Etot),40,'edgecolor','none');
colormap jet;
cobj = colorbar;
for i = 1:N
    plot([coord(i).xi/lambda coord(i).xe/lambda], [coord(i).yi/lambda coord(i).ye/lambda]);
end
axis([-5 5 -5 5]);
caxis([0 max(max(abs(Etot)))]);
xticks([-5 : 1 : 5]);
yticks([-5 : 1 : 5]);
set(cobj,'YTick',[0:0.5:max(max(abs(Etot)))]);
xlabel('x (\lambda)');
ylabel('y (\lambda)');
cobj.Label.String = 'V/m';
axis equal;
hold off;

% plot the total H-fields in the TE mode
subplot(3,2,6);
title('(f) H_{tot}');
hold on;
contourf(normxyrange,normxyrange,abs(Htot),40,'edgecolor','none');
colormap jet;
cobj = colorbar;
for i = 1:N
    plot([coord(i).xi/lambda coord(i).xe/lambda], [coord(i).yi/lambda coord(i).ye/lambda]);
end
axis([-5 5 -5 5]);
caxis([0 max(max(abs(Htot)))]);
xticks([-5 : 1 : 5]);
yticks([-5 : 1 : 5]);
set(cobj,'YTick',[0:0.5:max(max(abs(Htot)))]);
xlabel('x (\lambda)');
ylabel('y (\lambda)');
cobj.Label.String = 'A/m';
axis equal;
hold off;
print(fig,'MySavedPlot','-r300' ,'-dpng');








