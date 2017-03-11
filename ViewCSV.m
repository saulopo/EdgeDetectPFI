% ViewCSV - a matlab/octave script to view I/O data
%
% ViewCSV display contour plots of the input or output data
% involved in the EdgeDetectPFI algorithm (EdgeDetectPFI.f90)
%
% The user chooses whether to plot the input or output. 
% input:  contours of the signum transforms of Mz and Mz+-Mh
% output: contours of depths (h) and widths (2*a)
%
% Ref: SP Oliveira, FJF Ferreira & J de Souza, "EdgeDetectPFI: an 
% algorithm for automatic edge detection from magnetic anomalies", 
% published in Computers & Geosciences
%
% http://dx.doi.org/10.1016/j.cageo.2017.02.006

io = input('Enter 0 to view input data and 1 to view output data: ');

if(io==0)
  M      = dlmread('input.csv');
  nplot  = 2;
  titles = [{'ST Mz'},{'ST Mz-Mh'}];
else
  ho     = input(['Enter flying height (0 if none): ']);
  M      = dlmread('output.csv'); 
  nplot  = 2;
  titles = [{'Depth'},{'Width'}];
  M(:,3) = M(:,3) - ho;
end

tol = 1.0e-10;

% create 2D grid

x = M(:,1); y = M(:,2); 

xo    = min(x); xf = max(x);               
dxv   = x(2:end)-x(1:end-1);
dx    = min( abs(dxv) + 10^10*(dxv<tol) ); 
xgrid = xo:dx:xf;  Nx = length(xgrid);

yo    = min(y); yf = max(y);               
dyv   = y(2:end)-y(1:end-1);
dy    = min( abs(dyv) + 10^10*(dyv<tol) );
ygrid = yo:dy:yf;  Ny = length(ygrid);

[xx,yy] = meshgrid(xgrid,ygrid);        

% gather and plot dependent variables

z = zeros(Ny,Nx,nplot) + NaN;
for i = 1:Ny
  for j = 1:Nx
    [dist,k] = min( (x-xx(i,j)).^2 + (y-yy(i,j)).^2 ); 
    if(dist<tol)
      for l = 1:nplot
        z(i,j,l) = M(k,l+2);
      end
    end
  end
end

for l = 1:nplot
  figure(l)
  contourf(xx,yy,z(:,:,l),'LineColor','none')
  set(gca,'FontSize',16); colorbar; axis equal; axis([xo,xf,yo,yf]);
  title(titles(l));
end
