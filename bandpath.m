
function ret = bandpath(nx,ny,Estates,Evals,npoints,nE,Emin, Emax, lambda, theta, SLBZ)
    
    a = 1; % nearest neighbor distance in angstroms
    
    %List of xy pairs to target on path
    if(SLBZ == false)
        Ktargets = ((2/3)*pi)*[[0 0];[0 1];[tand(30) 1];[0 0]];
    else
        n1 = lambda*10^9;
        n2 = lambda*10^9;
        a1 = [1 0];
        a2 = [1/2 sqrt(3)/2];
        b1 = n1*a1 + n2*a2;
        b2 = -n2*a1 + (n1+n2)*a2;
        Ktargets = [];
    end
    %Generates path from destination points
    Kvals = interpolate(Ktargets,npoints);
    nKvals = length(Kvals);
    Rvals = a*generateRealSpaceBasis(nx,ny); %gives (x,y) pos per index
    bz1PWs = kToPW(Kvals,Rvals,nx,ny); %PWs for 1st brillouin zone
    bz2PWs = kToPW(Kvals + [0 (4/3)*pi],Rvals,nx,ny); %PWs for 2nd brillouin zone
    DOS = zeros(nKvals, nE);
    %size(H)
    %[Estates, Evals] = eig(H);

    labels = {'\Gamma','M','K','\Gamma'};
    % move along the xpath
    for i = 1:nKvals
        %env = envelope(Rvals,Xvals(i,:)); %envelope to project wfc onto site
        % start dissecting the hamiltonian
        PW1 = bz1PWs(:,i);
        PW2 = bz2PWs(:,i);
        for j = 1:length(Estates)
            psi = Estates(:,j);
            %size(psi)
            %size(PW)
            overlap = abs((psi')*PW2)^2 + abs((psi')*PW1)^2; % performs |<psi|planewave(k)>|^2 == <PW|psi><psi|PW> for changing basis, adds 1st and 2nd BZ
            
            % get DOS index and bin the delta(E)*overlap
            iDOS = round(nE*(Evals(j,j) - Emin)/(Emax - Emin));
            if (iDOS > 0 && iDOS <= nE)
                DOS(i,iDOS) = DOS(i,iDOS) + overlap;
            end
        end
    end
    
    %Begin the plotting
    fdos = figure;
    figure(fdos);
    dE = (Emax-Emin)/nE + 0.00000001; %small number helps align plot axes
    ax = axes;
    [X,Y] = meshgrid(0:nKvals-1,Emin:dE:Emax);
    s = surf(X,Y,zeros(nE, nKvals),transpose(DOS));
    s.EdgeColor = 'none';
    xticks(ax,[0:npoints:nKvals]);
    ylim(ax,[Emin+0.01 Emax-0.01]);
    xlim(ax,[0 nKvals-1]);
    ax.FontSize = 10;
    xticklabels(ax,labels);
    hold on;
    for i = 0:(length(Ktargets)-1)
        plot(ax,i*[npoints, npoints],[-10, 10],'k','LineWidth',0.8)
    end
    view(2);
    hold on;
    xlabel('K', 'FontSize',10)
    ylabel('Energy (eV)','FontSize',10)
    %set(ax, 'FontSize', 20);
    colormap(ax,fireice)
    hcb=colorbar;
    caxis(ax, [0 (30*nx*ny)/(nE)]);
    title(hcb,'DOS (arb. units)');
    ret = fdos;
end

function ret = interpolate(xlist,npts) % returns a list with npts interpolated between target points
    d = 1/npts;
    nTargets = length(xlist);
    Xvals = [];
    x = 0;
    for x1 = 1:(nTargets-1)
		x2 = x1 + 1;
		for f = 0:d:(1-d) %fraction of the 2nd xpt
			x = (1-f)*xlist(x1,:) + f*xlist(x2,:);
			Xvals = [Xvals;x];
		end
	end
	Xvals = [Xvals;xlist(end,:)];
    ret = Xvals;
end

function ret = envelope(grid,x) %Will return a Gaussian or delta function centered at given xpt
    n = length(grid); %number of Carbons
    env = zeros(n,1);
    sum = 0;
    spread = 10.5;
    gaussEnv = true;
    if(gaussEnv)
    for i = 1:n
        val = exp(-norm(grid(i,:)- x)^2/(2*spread^2));
        env(i,1) = val;
        sum = sum + val;
    end
    else
        index = floor(x(1) + (nx-1)*x(2)); %trash, disregard
        env(index) = 1;
        sum = 1;
    end
    ret = (1/sum)*env; %normalize so Integral(env) = one
end

function ret = generateRealSpaceBasis(nx,ny)
    Rvals = zeros(nx*ny,2);
    supercell = true;
    if(supercell)
        for ix = 1:nx
            y = 0; % reset yval back to 0
            x = (ix-1)*sqrt(3); %xval back to starting x position before y loop
            for iy = 1:ny
                Rvals(iy+ny*(ix-1),:) = [x y];
                if mod(iy,4) == 1
                    y = y+0.5; 
                    x = x + 0.5*sqrt(3);
                elseif mod(iy,4) == 2
                    y = y+1; 
                    x = x;
                elseif mod(iy,4) == 3
                    y = y+0.5; 
                    x = x - 0.5*sqrt(3);
                else
                    y = y+1.0; 
                    x = x;   
                end
            end
        end
    else
        for ix = 1:nx
            y = 0; % reset yval back to 0
            x = 0; %xval back to starting x position before y loop
            for iy = 1:ny
                Rvals(iy+ny*(ix-1),:) = [x y];
                if mod(iy,4) == 1
                    y = y+0.5; 
                    x = x + 0.5*sqrt(3);
                elseif mod(iy,4) == 2
                    y = y+1; 
                    x = x;
                elseif mod(iy,4) == 3
                    y = y+0.5; 
                    x = x - 0.5*sqrt(3);
                else
                    y = 0; 
                    x = 0;   
                end
            end
        end
    end
    ret = Rvals;
end

function ret = kToPW(Kvals, Rvals, nx, ny)
    nkpts = length(Kvals);
    PWs = zeros(nx*ny, nkpts);
    for ik = 1:nkpts
        for iR = 1:(nx*ny)
            k = Kvals(ik,:);
            R = Rvals(iR,:);
            %size(k)
            %size(R)
            PWs(iR,ik) = exp(-1i*dot(k,R));
        end
    end
    ret = PWs;
end