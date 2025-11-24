steps = [1 2 3 4 5];
loadonly = 0;

addpath(genpath('/projects/domino/KNS/'));
%Organizer
org = organizer('repository','Models/','prefix','LIA_crashfix-','steps',steps);

if perform(org,'Mesh') % {{{

	md=model();
	md.miscellaneous.name='KNS';
	md.miscellaneous.notes='Setting up initialization';
	md.mesh.epsg=3413;

	% Creating initial mesh
	md=triangle(model,'/projects/domino/KNS/Exp/domain.exp',500);

	disp(' -- Refining Mesh');
	ncdatavx = '/projects/domino/KNS/Trunk/datasets/Velocity/Clipped/vx_vel-CL.nc';
	ncdatavy = '/projects/domino/KNS/Trunk/datasets/Velocity/Clipped/vy_vel-CL.nc' ;

	velx= ncread(ncdatavx,'Band1')';
	xx= ncread(ncdatavx,'x');
	yx= ncread(ncdatavx,'y');

	vely= ncread(ncdatavy,'Band1')';
	xy= ncread(ncdatavx,'x');
	yy= ncread(ncdatavx,'y');

	vx	= InterpFromGrid(xx,yx,velx,md.mesh.x,md.mesh.y);
	vy	= InterpFromGrid(xy,yy,vely,md.mesh.x,md.mesh.y);
	vel	= sqrt(vx.^2+vy.^2);
	in=ContourToNodes(md.mesh.x,md.mesh.y,'/projects/domino/KNS/Exp/no-ice-mask_0603.exp',1);
	vx(find(in))=0;
	vy(find(in))=0;
	vel(find(in))=0;

	%Refine TW
	h=NaN*ones(md.mesh.numberofvertices,1);
	in=ContourToNodes(md.mesh.x,md.mesh.y,'/projects/domino/KNS/Exp/refinement.exp',1);
	h(find(in))=500; % may be too fine

	% Refine mesh using surface velocities as metric
	md=bamg(md,'hmin',500,'hmax',1000000,'field',vel,'err',5,'hVertices',h); 
	[md.mesh.lat,md.mesh.long] = xy2ll(md.mesh.x,md.mesh.y,+1,39,71);
	clear h in ncdatavx ncdatavy vel velx vely vx vy xx xy yx yy %tidy

	savemodel(org, md);
end % }}}
if perform(org,'Parameterization') % {{{
	md = loadmodel(org, 'Mesh');

	disp('Interpolate geometry');
	ncsurf='/projects/domino/KNS/Trunk/datasets/BedMachine/Clipped/surface-bedmachine-CL.nc';
	nctopg='/projects/domino/KNS/Trunk/datasets/BedMachine/Clipped/bed-bedmachine-CL.nc'; 
	ncthick='/projects/domino/KNS/Trunk/datasets/BedMachine/Clipped/thickness-bedmachine-CL.nc'; 
	x1= ncread(ncsurf,'x'); y1=ncread(ncsurf,'y'); surf=ncread(ncsurf,'surface')';
	x2= ncread(nctopg,'x'); y2=ncread(nctopg,'y'); topg=ncread(nctopg,'bed')';
	x3= ncread(ncthick,'x'); y3=ncread(ncthick,'y'); thick=ncread(ncthick,'thickness')';
	md.geometry.surface=InterpFromGrid(x1,y1,surf,md.mesh.x,md.mesh.y,'linear');
	md.geometry.bed = InterpFromGrid(x2,y2,topg,md.mesh.x,md.mesh.y, 'linear');
	md.geometry.thickness = InterpFromGrid(x3,y3,thick,md.mesh.x,md.mesh.y, 'linear');

	md.geometry.base=md.geometry.bed; %because i don't have any floating ice
	pos=find(md.geometry.bed<0 & ContourToNodes(md.mesh.x,md.mesh.y,'/projects/domino/KNS/Exp/no-ice-mask_0603.exp',1));
	md.geometry.base(pos) = 0;
	md.geometry.surface(pos)=md.geometry.base(pos);
	loc=find(md.geometry.thickness<=10 & ~(ContourToNodes(md.mesh.x,md.mesh.y,'/projects/domino/KNS/Exp/no-ice-mask_0603.exp',1)));
	spt=find(md.geometry.bed>0 & ContourToNodes(md.mesh.x,md.mesh.y,'/projects/domino/KNS/Exp/no-ice-mask_0603.exp',1));
	md.geometry.surface(loc)=md.geometry.bed(loc)+10;
	md.geometry.surface(spt)=md.geometry.bed(spt)+10;

	md.geometry.thickness=md.geometry.surface-md.geometry.base;

	disp('   Interpolating velocities');
	ncdatavx = '/projects/domino/KNS/Trunk/datasets/Velocity/Clipped/vx_vel-CL.nc';
	ncdatavy = '/projects/domino/KNS/Trunk/datasets/Velocity/Clipped/vy_vel-CL.nc' ;
	ncdatavv = '/projects/domino/KNS/Trunk/datasets/Velocity/Clipped/cal_vel-CL.nc';
	xx= ncread(ncdatavx,'x'); yx= ncread(ncdatavx,'y'); velx= ncread(ncdatavx,'Band1')';
	xy= ncread(ncdatavx,'x'); yy= ncread(ncdatavx,'y'); velv= ncread(ncdatavv,'Band1')';
	xv= ncread(ncdatavv,'x'); yv= ncread(ncdatavv,'y'); vely= ncread(ncdatavy,'Band1')';
	md.inversion.vel_obs = InterpFromGrid(xv,yv,velv,md.mesh.x,md.mesh.y,'linear');
	md.inversion.vx_obs	= InterpFromGrid(xx,yx,velx,md.mesh.x,md.mesh.y,'linear');
	md.inversion.vy_obs= InterpFromGrid(xy,yy,vely,md.mesh.x,md.mesh.y,'linear');
	in = ContourToNodes(md.mesh.x, md.mesh.y, '/projects/domino/KNS/Exp/no-ice-mask_0603.exp', 1);
	pos = find(in | isnan(md.inversion.vel_obs));
	md.inversion.vel_obs(pos)=0; md.initialization.vel= md.inversion.vel_obs;
	md.inversion.vx_obs(pos)=0; md.initialization.vx = md.inversion.vx_obs;
	md.inversion.vy_obs(pos)=0; md.initialization.vy = md.inversion.vy_obs;
	md.initialization.vz=zeros(md.mesh.numberofvertices,1);
	
	in = ContourToNodes(md.mesh.x,md.mesh.y,'/projects/domino/KNS/Exp/no-ice-mask_0603.exp',1) & md.geometry.bed<0;
	md.levelset.spclevelset = -1*ones(md.mesh.numberofvertices,1);
	md.levelset.spclevelset(find(in)) = 1;
	md.levelset.spclevelset = reinitializelevelset(md, md.levelset.spclevelset);
	md.mask.ice_levelset    = md.levelset.spclevelset;
	md.mask.ocean_levelset  = sethydrostaticmask(md);

	disp('   Interpolating temperatures');
	ncdatatemp='/projects/domino/KNS/Trunk/datasets/RACMO23p2/Clipped/t2m_clipped.nc';
	xtemp= ncread(ncdatatemp,'x'); ytemp= ncread(ncdatatemp,'y');
	temp_month = [];
	for i = 25:385 % 30 year SMB average from 1960 to 1990, note data is given in bands each month since 1958-01-15
		varName = sprintf('Band%d', i); 
		tempBand = ncread(ncdatatemp, varName); 
		md.initialization.temperature = InterpFromGrid(xtemp, ytemp, tempBand', md.mesh.x, md.mesh.y,'linear');
		temp_month = cat(3, temp_month, md.initialization.temperature); 
	end
	md.initialization.temperature = mean(temp_month, 3);

	disp('   Construct ice rheological properties (assume -10C)');
	md.materials.rheology_n=3*ones(md.mesh.numberofelements,1);
	md.materials.rheology_B=cuffey(273.15-10)*ones(md.mesh.numberofvertices, 1); 
	md.damage.D=zeros(md.mesh.numberofvertices,1);
	md.materials.rheology_law='None';

	disp('   Set other boundary conditions');
	md=SetIceSheetBC(md);
	pos = find(md.mesh.vertexonboundary);
	md.masstransport.spcthickness(pos) = md.geometry.thickness(pos);
	md.masstransport.stabilization = 1; %5: SUPG, 2: SU, 1: art diff
%	in=ContourToNodes(md.mesh.x,md.mesh.y,'/projects/domino/KNS/Exp/adjustBC.exp',1);
%	md.stressbalance.spcvx(find(in))=NaN; md.stressbalance.spcvy(find(in))=NaN; md.stressbalance.spcvz(find(in))=NaN;
	md.stressbalance.spcvx(pos) = md.initialization.vx(pos); md.stressbalance.spcvy(pos) = md.initialization.vy(pos); md.stressbalance.spcvz(pos) = md.initialization.vz(pos);
	md = setflowequation(md,'SSA','all');
	md.toolkits.DefaultAnalysis=bcgslbjacobioptions();
	md.masstransport.min_thickness=10;
	savemodel(org, md);
end % }}}
if perform(org, 'Inversion') % {{{
	md = loadmodel(org, 'Parameterization');

	disp('   Initialize basal friction using driving stress');
	disp('      -- Compute surface slopes and use 10 L2 projections');
	[sx,sy,s]=slope(md); sslope=averaging(md,s,10);
	disp('      -- Process surface velocity data');
	vel = md.inversion.vel_obs;
	flags=(vel==0); pos1=find(flags); pos2=find(~flags);
	vel(pos1) = griddata(md.mesh.x(pos2),md.mesh.y(pos2),vel(pos2),md.mesh.x(pos1),md.mesh.y(pos1));
	vel=max(vel,0.1);
	disp('      -- Calculate effective pressure');
	Neff = md.materials.rho_ice*md.geometry.thickness+md.materials.rho_water*md.geometry.base;
	Neff(find(Neff<=0))=1;
	disp('      -- Deduce friction coefficient');
	md.friction.coefficient=sqrt(md.materials.rho_ice*md.geometry.thickness.*(sslope)./(Neff.*vel/md.constants.yts));
	md.friction.coefficient=min(md.friction.coefficient,150);
	md.friction.p = 1.0 * ones(md.mesh.numberofelements,1);
	md.friction.q = 1.0 * ones(md.mesh.numberofelements,1);
	disp('      -- Extrapolate on ice free and floating ice regions');
	flags=(md.mask.ice_levelset>0) | (md.mask.ocean_levelset<0); pos1=find(flags); pos2=find(~flags);
	md.friction.coefficient(pos1) = 1;
	pos=find(isnan(md.friction.coefficient));
	md.friction.coefficient(pos)  = 1;
	pos=find(md.mask.ice_levelset<0 & md.initialization.vel==0);
	md.friction.coefficient(pos)=150;

	md.inversion=m1qn3inversion();
	md.inversion.vx_obs=md.initialization.vx;
	md.inversion.vy_obs=md.initialization.vy;
	md.inversion.vel_obs=md.initialization.vel;
	md.inversion.iscontrol=1;
	md.inversion.cost_functions=[101 103 501];
	md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,3);
	md.inversion.cost_functions_coefficients(:,1)=400;
	md.inversion.cost_functions_coefficients(:,2)=0.01;
	md.inversion.cost_functions_coefficients(:,3)=1e-7;

	pos = find(md.geometry.thickness<150 & md.geometry.thickness>10);
	md.inversion.cost_functions_coefficients(pos,1)=300;

	md.inversion.maxsteps = 80;
	md.inversion.maxiter  = 80;
	md.inversion.dxmin=0.001;
	%no cost where no ice
	pos = find(md.mask.ice_levelset>0 | md.inversion.vel_obs==0);
	md.inversion.cost_functions_coefficients(pos,[1]) = 0;

	md.inversion.control_parameters={'FrictionCoefficient'};
	md.inversion.min_parameters=0.1*ones(md.mesh.numberofvertices,1);
	md.inversion.max_parameters=150*ones(md.mesh.numberofvertices,1);

	%	pos = find(md.geometry.thickness<150 & md.geometry.thickness>10);
	%	md.inversion.max_parameters(pos)=50;

	md.stressbalance.restol=0.01;
	md.stressbalance.reltol=0.1;
	md.stressbalance.abstol=NaN;

	disp('Running inversion');
	% Solve

		md.verbose=verbose(0);
		md.verbose.control=1;
		md.settings.waitonlock=1;
	%md.verbose=verbose('all');
	md.verbose.solution=0;
	md.cluster=generic('name',oshostname(),'np',2);
	md.miscellaneous.name='invert';
	md=solve(md,'Stressbalance');

	md.friction.coefficient=md.results.StressbalanceSolution.FrictionCoefficient;

	disp(' --Running regression');
	in = md.geometry.bed < 0; %below sea-level
	indepbed = md.geometry.bed(find(in));
	depfriction = md.friction.coefficient(find(in));
	p = polyfit(indepbed, depfriction, 1);
	fittedValues = polyval(p, indepbed);
	md.friction.coefficient(find(in))= p(1)*md.geometry.bed(find(in)) + p(2);
	md.friction.coefficient(find(in & md.friction.coefficient<0.1)) = 0.1;
	% Plot the data and the regression line
	figure;
	scatter(indepbed, depfriction, 'b'); % Scatter plot of the filtered data
	hold on;
	plot(indepbed, fittedValues, 'r', 'LineWidth', 2); % Regression line
	hold off;
	%pos=find(~coefficientontourToNodes(md.mesh.x,md.mesh.y,'/projects/domino/KNS/Exp/frontforcing/1920.exp',1) & md.mask.ice_levelset>0);
	%md.friction.coefficient(pos) = 0.1; %This is a test for the forced front. Depending on outcome, should consider a timeseries friction coefficient given 1946 and 1948 shallow sills are causing unrealistic stability in the advance.
	plotmodel(md,'data',md.friction.coefficient);
	savemodel(org, md);

end % }}}
if perform(org,'RACMO-forcing') % {{{
	md=loadmodel(org,'Inversion');
   md.smb=SMBforcing();
	ncdata='/projects/domino/KNS/Trunk/datasets/RACMO23p2/Clipped/mean_smb_clipped_3413_1960-1990.nc';
	% finfo = ncinfo(ncdatasmb);
	xsmb= ncread(ncdata,'x');
	ysmb= ncread(ncdata,'y');
%	smb_month = [];
%	for i = 1:30 % 30 year SMB average from 1960 to 1990, note data is given in bands each month since 1958-01-15
%		varName = sprintf('Band%d', i); 
	smbBand = ncread(ncdata, 'Band1'); 
	md.smb.mass_balance = InterpFromGridToMesh(xsmb, ysmb, smbBand', md.mesh.x, md.mesh.y, 0);
%		smb_month = [smb_month md.smb.mass_balance]; 
%	end
	% Calculate the average value for each cell across all time steps
%	md.smb.mass_balance = mean(smb_month, 3);
	md.smb.mass_balance=md.smb.mass_balance/1000*(1000/917);	
	pos = find(md.smb.mass_balance<0);

	md.smb.mass_balance(pos) = md.smb.mass_balance(pos)*2; %artifical forcing, so that dH/dt = -(grad)uH + b (where dH is change in thickness, dt is change in time, (grad)uH is short-hand vector notation for mass transport, and b is SMB). When reaching an initialized steady state, dH/dt needs to be approx. 0
	savemodel(org,md);

end  % }}}
if perform(org, 'Solve-Relaxation') % {{{
	md = loadmodel(org, 'RACMO-forcing');
	begin=300; % Year CE, start of forcing period
	finish=900; %YearCE, end of forcing period

%	md.friction.Cmax = 0.17 * ones(md.mesh.numberofvertices,1);


	name='artracmo-relaxation_spcthick-schoof';
%	md.initialization.vel = md.results.StressbalanceSolution.Vel;
%	md.initialization.vx = md.results.StressbalanceSolution.Vx;
%	md.initialization.vy = md.results.StressbalanceSolution.Vy;

	pos = find(md.mesh.vertexonboundary);
	md.masstransport.spcthickness(pos) = md.geometry.thickness(pos);

	md.transient.issmb=0;
	md.transient.ismasstransport=1;
	md.transient.isstressbalance=1;
	md.transient.ismovingfront=0;
	md.transient.isgroundingline=1;
	% md.transient.isgia=0;
	md.transient.isthermal=0;
	md.transient.isdamageevolution=0;
	md.transient.ishydrology=0;
	% md.transient.isslr=0;
	md.transient.isesa=0;
	md.transient.isoceancoupling=0;
	% md.transient.iscoupler=0;

	md.calving=calvinghab();     %This is a placeholder, not used while no moving front.
	md.calving.flotation_fraction=0;
	md.frontalforcings.meltingrate=zeros(md.mesh.numberofvertices,1);
	md.groundingline.migration='SubelementMigration';
	md.groundingline.friction_interpolation='SubelementFriction1';
	md.groundingline.melt_interpolation='NoMeltOnPartiallyFloating';

	md.timestepping = timesteppingadaptive(md.timestepping);
	md.timestepping.start_time=begin;
	md.timestepping.final_time=finish;
	md.timestepping.time_step_min=0.000001; 
	md.timestepping.time_step_max=0.1;
	md.settings.output_frequency=5;

	md.inversion.iscontrol=0;
	md.transient.requested_outputs={'default','IceVolume'}; % other options?

	md.verbose=verbose('all');

	loadonly = 0;

	md.settings.waitonlock = 0;
	md.cluster.interactive = 0; %only needed if you are using the generic cluster
	md.miscellaneous.name = name; 

	%		pos=find(md.initialization.vel>6000);
	%		md.initialization.vel(pos)=6000;
	md=solve(md,'Transient','runtimename',false);


	savemodel(org,md);
end % }}}

if perform(org, 'Load-Relaxation') % {{{
	md=loadmodel(org,'Solve-Relaxation');
	md=loadresultsfromcluster(md);
	% Check the CFL criteria and adjust adaptive timestep accordingly
	CFL = cfl_step(md,md.results.TransientSolution(end).Vx,md.results.TransientSolution(end).Vy)
	md.timestepping.time_step_min = 0.01; % 3.6 days -- This might be too big given CFL is 0.015
	md.results.TransientSolution(end).time

			md = transientrestart(md); % reinitialise model from end of the relaxation
			md.results = rmfield(md.results,'StressbalanceSolution');
			md.results = rmfield(md.results,'TransientSolution3');

	%		disp('Reset geometry to match mask after restart');
			md.geometry.base=md.geometry.bed; %because i don't have any floating ice
			pos=find(md.geometry.bed<0 & md.mask.ice_levelset>0);
			md.geometry.thickness(pos)=0;
			md.geometry.base(pos) = 0;
			loc=find(md.geometry.thickness<=10 & md.mask.ice_levelset<0);
			pos=find(md.geometry.bed>0 & md.mask.ice_levelset>0);
			md.geometry.thickness(loc)=10;
			md.geometry.thickness(pos)=10;
			md.geometry.surface=md.geometry.thickness+md.geometry.base;

	savemodel(org,md);

end %}}}
	if perform(org, 'PDD-Present') % {{{
		md = loadmodel(org, 'Load-Relaxation');
		[md.mesh.lat,md.mesh.long] = xy2ll(md.mesh.x,md.mesh.y,+1,39,71);

		% Set up SMB D180 class and set PDD factors

		md.smb = SMBd18opdd();  % Turn on Class
		md.smb.isd18opd=1;	% 1 means we want to load our monthly climatologies of temperature and precipitation 

		md.smb.delta18o = [-40.0110;0];  % This is meaningless, but needs to be set.  Basically, if you set md.smb.isd18opd=0, that would mean you want to take a modern climatology and scale it back through time based upon an ice core record.  For now, we just need to set this, but it will not be used since we set md.smb.isd18opd=1.

		md.smb.rlaps=6.0; % This is the Lapse rate you want to use
		md.smb.desfac = 1; % This is the elevation desertification effect.  0 is no change, 1 is maximum change.  Which is a reduction by a factor of 2 for every 1km change in surface elevation
		md.smb.rlapslgm = 5.5; % This just needs to be set, but is not used for your purposes

		% Set PDD Factors
		md.smb.issetpddfac = 1;
		md.smb.pddfac_snow = 4.3*ones(md.mesh.numberofvertices,1); % Notice these are spatially constant, but if you wanted could be spatially varying since they are on the vertices.
		md.smb.pddfac_ice = 8.3*ones(md.mesh.numberofvertices,1);

		disp('   ACC');
		ncdataacc='/projects/domino/KNS/Trunk/datasets/Box1840thr2012/accbox.nc';
		finfo = ncinfo(ncdataacc);
		xacc= ncread(ncdataacc,'x');
		yacc= ncread(ncdataacc,'y');

		md.smb.precipitations_presentday = NaN(md.mesh.numberofvertices,12);
		q=1;
		for i = 133:144 % 12 months
			index = i:12:((i+11)+12*(149)); %index across 1850 to 2000 for each month
			acc_month=[];
			for w= 1:150
				varName = sprintf('Band%d', index(w)); 
				accBand = ncread(ncdataacc, varName); 
				acc = InterpFromGridToMesh(xacc, yacc, accBand', md.mesh.x, md.mesh.y, 0);
				acc_month = cat(3, acc_month, acc); 
			end
			md.smb.precipitations_presentday(:,q)=mean(acc_month,3);
			q=q+1;
		end

		disp('   DEM');
		ncdatadem='/projects/domino/KNS/Trunk/datasets/Box1840thr2012/dembox.nc';
		% finfo = ncinfo(ncdatadem);
		xdem= ncread(ncdatadem,'x');
		ydem= ncread(ncdatadem,'y');
		dem=ncread(ncdatadem,'Band1');
		box_dem=InterpFromGridToMesh(xdem, ydem, dem', md.mesh.x, md.mesh.y, 0);
		md.smb.s0p = (max(box_dem,0));
		md.smb.s0t = (max(box_dem,0));

		disp('   Temp');
		ncdatatemp='/projects/domino/KNS/Trunk/datasets/Box1840thr2012/tempbox.nc';
		%finfo = ncinfo(ncdatatemp);
		xtemp= ncread(ncdatatemp,'x');
		ytemp= ncread(ncdatatemp,'y');

		md.smb.temperatures_presentday = NaN(md.mesh.numberofvertices,12);
		q=1;
		for i = 133:144 % 12 months
			index = i:12:((i+11)+12*(149)); %index across 1850 to 2000 for each month
			temp_month=[];
			for w= 1:150
				varName = sprintf('Band%d', index(w)); 
				tempBand = ncread(ncdatatemp, varName); 
				temp = InterpFromGridToMesh(xtemp, ytemp, tempBand', md.mesh.x, md.mesh.y, 0);
				temp_month = cat(3, temp_month, temp); 
			end
			md.smb.temperatures_presentday(:,q)=mean(temp_month,3);
			q=q+1;
		end
		temp_present = md.smb.temperatures_presentday;
		precip_present = md.smb.precipitations_presentday;
		save('PDD-presentday.mat', 'temp_present', 'precip_present')	

		md.smb.temperatures_presentday = md.smb.temperatures_presentday+274.15;
		md.smb.precipitations_presentday = md.smb.precipitations_presentday * (1/1000) * 12;
		savemodel(org,md);
	end % }}}	
if perform(org,'PDD-Paleo1100') % {{{
	md = loadmodel(org, 'PDD-Present');
	begin=1100; % Year CE, start of forcing period
	finish=1150; %Year CE, end of forcing period
	n=401-((1950-begin)/50); %Convert from year CE to index
	m=401-((1950-finish)/50);
	span=m-n;		

	load('PDD-presentday.mat');

	ncbadg='/projects/domino/KNS/Trunk/datasets/paleoran/briner2020recons.nc';
	pre1 = ncread(ncbadg,'P_moderate'); %Precipiation, options: P_low, P_moderate, P_high
	pre1 = permute(pre1,[4 3 2 1]);
	pre1 = pre1(n:m,:,:,:);
	latB  = ncread(ncbadg,'lat');
	lon  = ncread(ncbadg,'lon');
	x = find(lon>180);
	lonB(x) = lon(x)-360;
	pre_jessica = NaN(md.mesh.numberofvertices+1,span,12);
	disp('Interpolating precipitation output onto grid');
	for i = 1:span
		for k = 1:12
			pre_jessica(1:md.mesh.numberofvertices,i,k)=(InterpFromGrid(latB,lonB',squeeze(pre1(i,k,:,:))',md.mesh.lat,md.mesh.long,'linear'));
		end
	end

	precip = [];	% For precipitation, we multiply the Box climatology by the fraction from the Badgeley Reconstruction
	for i = 1:span
		for j = 1:12
			precip = [precip;precip_present(:,j)'.*pre_jessica(1:end-1,i,j)'];
		end
	end
	precip = permute(precip,[2 1]);
	precip=precip*(1/1000)*12; % Convert from mm/month WE to m/yr WE

	temp1 = ncread(ncbadg,'T_moderate'); %Temperature, options: T_low, T_moderate, T_high
	temp1 = permute(temp1,[4 3 2 1]);
	temp1 = temp1(n:m,:,:,:);
	latB  = ncread(ncbadg,'lat');
	lon  = ncread(ncbadg,'lon');
	x = find(lon>180);
	lonB(x) = lon(x)-360;
	temp_jessica = NaN(md.mesh.numberofvertices+1,span,12);
	for i = 1:span
		for k = 1:12
			temp_jessica(1:md.mesh.numberofvertices,i,k)=(InterpFromGrid(latB,lonB',squeeze(temp1(i,k,:,:))',md.mesh.lat,md.mesh.long,'linear'));
		end
	end

	tmp = [];  % For temperature, we add the Box climatology to the anomalies from the Badgeley Reconstruction
	for i = 1:span
		for j = 1:12
			tmp = [tmp;temp_present(:,j)'+temp_jessica(1:end-1,i,j)'];
		end
	end
	tmp = permute(tmp,[2 1]);
	tmp=tmp+274.15; %Convert from Celsius to Kelvin		

	x = [0+1/12:1/12:500+1/12];

	% Stride every 50 years - Time Resolution of UW Product
	a = x(1:600:end);b = x(2:600:end);c = x(3:600:end);d = x(4:600:end);e = x(5:600:end); f = x(6:600:end); g = x(7:600:end); h = x(8:600:end); ii = x(9:600:end);j = x(10:600:end); k = x(11:600:end); l = x(12:600:end);

	xx = [];
	for i =1:10;
		xx = [xx;a(i);b(i);c(i);d(i);e(i);f(i);g(i);h(i);ii(i);j(i);k(i);l(i)];
	end

	repprecip = [];
	reptemp = [];
	for i = 1:10
		repprecip = [repprecip precip];
		reptemp = [reptemp tmp];
	end

	md.smb.isprecipscaled=0; %  This allows us to set the precip ourselves to what we have done above.
	md.smb.precipitations_reconstructed = NaN(md.mesh.numberofvertices+1,10*12); % Needs to be size numberofvertices+1, so the last row can hold the timestep
	md.smb.precipitations_reconstructed(1:end-1,:) = repprecip;
	md.smb.precipitations_reconstructed(end,:) = xx;
	md.smb.precipitations_reconstructed(md.smb.precipitations_reconstructed<0) = 0.01; % There is a possibility for negative precipitation. values.  Set anything less than 0 to low value.

	md.smb.istemperaturescaled=0;
	md.smb.temperatures_reconstructed= NaN(md.mesh.numberofvertices+1,10*12);
	md.smb.temperatures_reconstructed(1:end-1,:) = reptemp;
	md.smb.temperatures_reconstructed(end,:) = xx;
	md.smb.rlaps=0;

	savemodel(org,md);
end % }}}
	if perform(org,'Relaxation-Calving-VonMISES') % {{{
		md = loadmodel(org, 'Relax-Inversion');
		md1 = loadmodel('/projects/domino/KNS/Models/KNS_Rheology10convergence-calving-vonmises.mat');
		md.smb=md1.smb;
		name='Comparing-withmaskfix';
		md.levelset.spclevelset=md1.levelset.spclevelset;
		md.mask.ocean_levelset=md1.mask.ocean_levelset;

%		in = ContourToNodes(md.mesh.x,md.mesh.y,'/projects/domino/KNS/Exp/frontforcing/900.exp',1) & md.geometry.bed<0;
%		md.levelset.spclevelset = -1*ones(md.mesh.numberofvertices,1);
%		md.levelset.spclevelset(find(in)) = 1;
%		pos = find(md.mesh.vertexonboundary);
%		md.levelset.spclevelset(:) = NaN;
%		md.levelset.spclevelset(pos) = md.mask.ice_levelset(pos);	
%		%md.levelset.spclevelset = reinitializelevelset(md, md.levelset.spclevelset);
%		md.mask.ocean_levelset  = sethydrostaticmask(md);
		md.transient.issmb=1;

%		md.calving=calvingvonmises();
		md.calving=md1.calving;
		
%		md.calving.min_thickness = md1.calving.min_thickness;
		md.calving.stress_threshold_groundedice=646000;
		md.calving.stress_threshold_floatingice=384000;
		md.basalforcings.floatingice_melting_rate=md1.basalforcings.floatingice_melting_rate;
		md.timestepping.start_time=0;
		md.timestepping.final_time=500;
		md.settings.output_frequency=5;
		md.transient.ismovingfront = 1;
		md.transient.requested_outputs={'default','IceVolume','TotalSmb','SmbMassBalance','CalvingCalvingrate'}; % other options?

		md.miscellaneous.name = name; 
		loadonly = 0;
      md.inversion.iscontrol=0;
      md.verbose=verbose('all');
      md.settings.waitonlock = 0;
      md.cluster.interactive = 0; %only needed if you are using the generic cluster


		md=solve(md,'Transient','runtimename',false,'loadonly',loadonly);

		savemodel(org,md);
	end	% }}}
if perform(org, 'Relaxation-Load-VonMISES') % {{{
	md=loadmodel(org,'Relaxation-Calving-VonMISES');
	md=loadresultsfromcluster(md);

%	md = transientrestart(md); % reinitialise model from end of the relaxation
%	md.results = rmfield(md.results,'StressbalanceSolution');
%	md.results = rmfield(md.results,'TransientSolution3');

	disp('Reset geometry to match mask after restart');
	md.geometry.base=md.geometry.bed; %because i don't have any floating ice
	pos=find(md.geometry.bed<0 & md.mask.ice_levelset>0);
	md.geometry.thickness(pos)=0;
	md.geometry.base(pos) = 0;
	loc=find(md.geometry.thickness<=10 & md.mask.ice_levelset<0);
	pos=find(md.geometry.bed>0 & md.mask.ice_levelset>0);
	md.geometry.thickness(loc)=10;
	md.geometry.thickness(pos)=10;
	md.geometry.surface=md.geometry.thickness+md.geometry.base;

	savemodel(org,md);
end %}}}


	return; %%%%%%%%%%%%%%%%%%%%%%
	md1 = loadmodel('/projects/domino/KNS/Models/KNS_1100-Load-VonMISES.mat'); 
	if perform(org,'ReadvanceToLIA') % {{{
		md = loadmodel(org, 'PDD-Paleo900-1900');

		%Setup moving boundary
		md.calving = calving();   %This is a placeholder, will not be used while forcing the front with spclevelset  
		md.calving.calvingrate = zeros(md.mesh.numberofvertices,1);
		%	md.calving=calvinghab();
		%	md.calving.flotation_fraction=0.8;	
		md.transient.ismovingfront = 1;

		pos = find(md.levelset.spclevelset(end,:)==1781);
		LIAlevelset = md.levelset.spclevelset(1:end-1, pos);
		md.levelset.spclevelset = LIAlevelset;
		pos = find(LIAlevelset<0 & ~md.mesh.vertexonboundary);
		md.levelset.spclevelset(pos) = NaN;

		in = md.geometry.bed < 0; %below sea-level
		indepbed = md.geometry.bed(find(in));
		depfriction = md.friction.coefficient(find(in));
		p = polyfit(indepbed, depfriction, 1);
		fittedValues = polyval(p, indepbed);
		md.friction.coefficient(find(in))= p(1)*md.geometry.bed(find(in)) + p(2);
		md.friction.coefficient(find(in & md.friction.coefficient<0.1)) = 0.1;

		md.timestepping = timesteppingadaptive(md.timestepping);
		md.timestepping.start_time=900;
		md.timestepping.final_time=1700;
		md.timestepping.time_step_min=0.01; 
		md.timestepping.time_step_max=1;
		md.settings.output_frequency=20;

		md.inversion.iscontrol=0;
		md.transient.requested_outputs={'default','IceVolume','TotalSmb','SmbMassBalance'}; % other options?

		name='LIA_initialized';
		md.miscellaneous.name = name; 

		md=solve(md,'Transient','runtimename',false,'loadonly',loadonly);

		savemodel(org,md);
	end% }}}
if perform(org, 'Load-ReadvanceToLIA') % {{{
	md=loadmodel(org,'ReadvanceToLIA');
	md=loadresultsfromcluster(md);
	alias
	savemodel(org,md);
end %}}}
if perform(org,'Calving-VonMISES') % {{{
	md = loadmodel(org, 'Load-ReadvanceToLIA');
	md = transientrestart(md); % reinitialise model from end of the relaxation
	md.results = rmfield(md.results,'TransientSolution2');

	md.calving=calvingvonmises();
	md.calving.min_thickness = 0;
	md.calving.stress_threshold_groundedice=7e5;
	md.calving.stress_threshold_floatingice=400e3;	
	md.timestepping.start_time=1760;
	md.timestepping.final_time=1900;
	md.settings.output_frequency=5;

	md.transient.requested_outputs={'default','IceVolume','TotalSmb','SmbMassBalance','CalvingCalvingrate'}; % other options?

	name='VonMises_7_0';
	md.miscellaneous.name = name; 

	md=solve(md,'Transient','runtimename',false,'loadonly',loadonly);

	savemodel(org,md);
end	% }}}
if perform(org,'PDD-Paleo') % {{{
	md = loadmodel(org, 'PDD-Present');
	begin=300; % Year CE, start of forcing period
	finish=900; %Year CE, end of forcing period
	n=401-((1950-begin)/50); %Convert from year CE to index
	m=401-((1950-finish)/50);
	span=m-n;		
	x = [begin+1/12:1/12:finish+1/12];

	load('PDD-presentday.mat');
	% Stride every 50 years - Time Resolution of UW Product
	a = x(1:600:end);b = x(2:600:end);c = x(3:600:end);d = x(4:600:end);e = x(5:600:end); f = x(6:600:end); g = x(7:600:end); h = x(8:600:end); ii = x(9:600:end);j = x(10:600:end); k = x(11:600:end); l = x(12:600:end);

	xx = [];
	for i =1:span;
		xx = [xx;a(i);b(i);c(i);d(i);e(i);f(i);g(i);h(i);ii(i);j(i);k(i);l(i)];
	end

	load ('/projects/domino/KNS/PDD-presentday.mat');
	ncbadg='/projects/domino/KNS/Trunk/datasets/paleoran/briner2020recons.nc';
	pre1 = ncread(ncbadg,'P_moderate'); %Precipiation, options: P_low, P_moderate, P_high
	pre1 = permute(pre1,[4 3 2 1]);
	pre1 = pre1(n:m,:,:,:);
	latB  = ncread(ncbadg,'lat');
	lon  = ncread(ncbadg,'lon');
	x = find(lon>180);
	lonB(x) = lon(x)-360;
	pre_jessica = NaN(md.mesh.numberofvertices+1,span,12);
	disp('Interpolating precipitation output onto grid');
	for i = 1:span
		for k = 1:12
			pre_jessica(1:md.mesh.numberofvertices,i,k)=(InterpFromGrid(latB,lonB',squeeze(pre1(i,k,:,:))',md.mesh.lat,md.mesh.long,'linear'));
		end
	end

	precip = [];	% For precipitation, we multiply the Box climatology by the fraction from the Badgeley Reconstruction
	for i = 1:span
		for j = 1:12
			precip = [precip;precip_present(:,j)'.*pre_jessica(1:end-1,i,j)'];
		end
	end
	precip = permute(precip,[2 1]);
	precip=precip*(1/1000)*12; % Convert from mm/month WE to m/yr WE

	monthly_precip = [];
	meanmonthly_precip = zeros(md.mesh.numberofvertices, 12); % Preallocate for efficiency

	for i = 1:12
		monthly_precip = precip(:, i:12:(12 * span)); % Extract all i-th months
		meanmonthly_precip(:, i) = mean(monthly_precip, 2); % Compute mean along rows
	end
	
	contmean_precip=[];
	for i=1:span
		contmean_precip= [contmean_precip, meanmonthly_precip];
	end

	md.smb.isprecipscaled=0; %  This allows us to set the precip ourselves to what we have done above.
	md.smb.precipitations_reconstructed = NaN(md.mesh.numberofvertices+1,span*12); % Needs to be size numberofvertices+1, so the last row can hold the timestep
	md.smb.precipitations_reconstructed(1:end-1,:) = contmean_precip;
	md.smb.precipitations_reconstructed(end,:) = xx;
	md.smb.precipitations_reconstructed(md.smb.precipitations_reconstructed<0) = 0.01; % There is a possibility for negative precipitation. values.  Set anything less than 0 to low value.

	temp1 = ncread(ncbadg,'T_moderate'); %Temperature, options: T_low, T_moderate, T_high
	temp1 = permute(temp1,[4 3 2 1]);
	temp1 = temp1(n:m,:,:,:);
	latB  = ncread(ncbadg,'lat');
	lon  = ncread(ncbadg,'lon');
	x = find(lon>180);
	lonB(x) = lon(x)-360;
	temp_jessica = NaN(md.mesh.numberofvertices+1,span,12);
	for i = 1:span
		for k = 1:12
			temp_jessica(1:md.mesh.numberofvertices,i,k)=(InterpFromGrid(latB,lonB',squeeze(temp1(i,k,:,:))',md.mesh.lat,md.mesh.long,'linear'));
		end
	end

	tmp = [];  % For temperature, we add the Box climatology to the anomalies from the Badgeley Reconstruction
	for i = 1:span
		for j = 1:12
			tmp = [tmp;temp_present(:,j)'+temp_jessica(1:end-1,i,j)'];
		end
	end
	tmp = permute(tmp,[2 1]);
	tmp=tmp+274.15; %Convert from Celsius to Kelvin

	monthly_temp = [];
	meanmonthly_temp = zeros(md.mesh.numberofvertices, 12); % Preallocate for efficiency

	for i = 1:12
		monthly_temp = tmp(:, i:12:(12 * span)); % Extract all i-th months
		meanmonthly_temp(:, i) = mean(monthly_temp, 2); % Compute mean along rows
	end

	contmean_temp=[];
	for i=1:span
		contmean_temp= [contmean_temp, meanmonthly_temp];
	end

	md.smb.istemperaturescaled=0;
	md.smb.temperatures_reconstructed= NaN(md.mesh.numberofvertices+1,span*12);
	md.smb.temperatures_reconstructed(1:end-1,:) = contmean_temp;
	md.smb.temperatures_reconstructed(end,:) = xx;
	md.smb.rlaps=0;

	savemodel(org,md);
end % }}}
if perform(org, 'Inversion-Schoof') % {{{
	md2 = loadmodel(org, 'Parameterization');
	md2.friction = frictionschoof();
	md2=EstimateFric_Schoof(md2); % initial guess from Driving Stress

	md2.inversion=m1qn3inversion();
	md2.inversion.vx_obs=md2.initialization.vx;
	md2.inversion.vy_obs=md2.initialization.vy;
	md2.inversion.vel_obs=md2.initialization.vel;
	md2.inversion.iscontrol=1;
	md2.inversion.cost_functions=[101 103 501];
	md2.inversion.cost_functions_coefficients=ones(md2.mesh.numberofvertices,3);

	md2.inversion.cost_functions_coefficients(:,1)=600;
	md2.inversion.cost_functions_coefficients(:,2)=.01;
	md2.inversion.cost_functions_coefficients(:,end)=1e-11;
	md2.inversion.control_parameters={'FrictionC'};
	md2.inversion.min_parameters=0*ones(md2.mesh.numberofvertices,1);
	md2.inversion.max_parameters=12000*ones(md2.mesh.numberofvertices,1);

	md2.inversion.maxsteps = 80;
	md2.inversion.maxiter  = 80;
	md2.inversion.dxmin=0.001;
		md2.stressbalance.restol=0.01;
	md2.stressbalance.reltol=0.1;
	md2.stressbalance.abstol=NaN;

	disp('Running inversion');
	% Solve
	%md2.verbose=verbose('all');
	md2.verbose.solution=0;
	md2.cluster=generic('name',oshostname(),'np',2);
	md2.miscellaneous.name='invert';
	md2=solve(md2,'Stressbalance');

	md2.friction.C=md2.results.StressbalanceSolution.FrictionC;

	disp(' --Running regression');
	in = md2.geometry.bed < 0; %below sea-level
	indepbed = md2.geometry.bed(find(in));
	depfriction = md2.friction.C(find(in));
	p = polyfit(indepbed, depfriction, 1);
	fittedValues = polyval(p, indepbed);
	md2.friction.C(find(in))= p(1)*md2.geometry.bed(find(in)) + p(2);
	md2.friction.C(find(in & md2.friction.C<0.1)) = 0.1;
	% Plot the data and the regression line
	figure;
	scatter(indepbed, depfriction, 'b'); % Scatter plot of the filtered data
	hold on;
	plot(indepbed, fittedValues, 'r', 'LineWidth', 2); % Regression line
	hold off;
%	pos=find(~ContourToNodes(md2.mesh.x,md2.mesh.y,'/projects/domino/KNS/Exp/frontforcing/1920.exp',1) & md2.mask.ice_levelset>0);
%	md2.friction.C(pos) = 0.1; %This is a test for the forced front. Depending on outcome, should consider a timeseries friction coefficient given 1946 and 1948 shallow sills are causing unrealistic stability in the advance.
	plotmodel(md2,'data',md2.friction.C);
	savemodel(org, md2);

end % }}}

if perform(org, 'Load-VonMISES') % {{{
	md=loadmodel(org,'Calving-VonMISES');
	md=loadresultsfromcluster(md);
	alias
	savemodel(org,md);
end %}}}
	if perform(org, 'Solve-Contemp-Relaxation') % {{{
		md = loadmodel(org, 'RACMO-forcing');
		begin=0; % Year CE, start of forcing period
		finish=50; %Year CE, end of forcing period

		md.transient.issmb=1;
		md.transient.ismasstransport=1;
		md.transient.isstressbalance=1;
		md.transient.ismovingfront=0;
		md.transient.isgroundingline=1;
		% md.transient.isgia=0;
		md.transient.isthermal=0;
		md.transient.isdamageevolution=0;
		md.transient.ishydrology=0;
		% md.transient.isslr=0;
		md.transient.isesa=0;
		md.transient.isoceancoupling=0;
		% md.transient.iscoupler=0;

		md.calving=calvinghab();     %This is a placeholder, not used while no moving front.
		md.calving.flotation_fraction=0;
		md.frontalforcings.meltingrate=zeros(md.mesh.numberofvertices,1);
		md.groundingline.migration='SubelementMigration';
		md.groundingline.friction_interpolation='SubelementFriction1';
		md.groundingline.melt_interpolation='NoMeltOnPartiallyFloating';

		md.timestepping = timesteppingadaptive(md.timestepping);
		md.timestepping.start_time=begin;
		md.timestepping.final_time=finish;
		md.timestepping.time_step_min=0.0001; 
		md.timestepping.time_step_max=1;
		md.settings.output_frequency=5;

		md.inversion.iscontrol=0;
		md.transient.requested_outputs={'default','IceVolume','TotalSmb','SmbMassBalance'}; % other options?

		md.verbose=verbose('all');

		name='relaxation_contemp';
		loadonly = 0;

		md.settings.waitonlock = 0;
		md.cluster.interactive = 0; %only needed if you are using the generic cluster
		md.miscellaneous.name = name; 

%		pos=find(md.initialization.vel>6000);
%		md.initialization.vel(pos)=6000;
		md=solve(md,'Transient','runtimename',false);


		savemodel(org,md);
	end % }}}
if perform(org, 'Load-Contemp-Relaxation') % {{{
	md=loadmodel(org,'Solve-Contemp-Relaxation');
	md=loadresultsfromcluster(md);
	savemodel(org,md);

end %}}}
if perform(org,'Contemp-RelaxationRestart') % {{{
	md = loadmodel(org,'Load-Contemp-Relaxation');
	%	md.geometry.thickness=md.results.TransientSolution(end).Thickness;
	%	md.geometry.base=md.results.TransientSolution(end).Base;
	%	pos=find(md.geometry.base<md.geometry.bed);
	%	md.geometry.base(pos)=md.geometry.bed(pos);
	%	md.geometry.surface=md.results.TransientSolution(end).Surface;
	%	md.initialization.vel=md.results.TransientSolution(end).Vel;
	%	md.initialization.vx=md.results.TransientSolution(end).Vx;
	%	md.initialization.vy=md.results.TransientSolution(end).Vy;

	% Check the CFL criteria and adjust adaptive timestep accordingly
	CFL = cfl_step(md,md.results.TransientSolution(end).Vx,md.results.TransientSolution(end).Vy);
	md.timestepping.time_step_min = 0.01; % 3.6 days -- This might be too big given CFL is 0.015

	md = transientrestart(md); % reinitialise model from end of the relaxation
	md.results = rmfield(md.results,'StressbalanceSolution');
	md.results = rmfield(md.results,'TransientSolution3');

	savemodel(org,md);
end %}}}

if perform(org, 'Relax-Inversion-Schoof') % {{{
	md = loadmodel(org, 'Load-Relaxation');
	md=EstimateFric_Schoof(md); % initial guess from Driving Stress

	md.inversion=m1qn3inversion();
	md.inversion.vx_obs=md.initialization.vx;
	md.inversion.vy_obs=md.initialization.vy;
	md.inversion.vel_obs=md.initialization.vel;
	md.inversion.iscontrol=1;
	md.inversion.cost_functions=[101 103 501];
	md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,3);

	md.inversion.cost_functions_coefficients(:,1)=500;
	md.inversion.cost_functions_coefficients(:,2)=.01;
	md.inversion.cost_functions_coefficients(:,end)=1e-11;
	md.inversion.control_parameters={'FrictionC'};
	md.inversion.min_parameters=0*ones(md.mesh.numberofvertices,1);
	md.inversion.max_parameters=12000*ones(md.mesh.numberofvertices,1);

	md.inversion.maxsteps = 80;
	md.inversion.maxiter  = 80;
	md.inversion.dxmin=0.001;
	md.stressbalance.restol=0.01;
	md.stressbalance.reltol=0.1;
	md.stressbalance.abstol=NaN;

	disp('Running inversion');
	% Solve
	md.verbose=verbose('000100000');
	md.settings.waitonlock=1;
	md.cluster=generic('name',oshostname(),'np',2);
	md.miscellaneous.name='relax-invert';
	md=solve(md,'Stressbalance');

	md.friction.C=md.results.StressbalanceSolution.FrictionC;

	disp(' --Running regression');
	in = md.geometry.bed < 0; %below sea-level
	indepbed = md.geometry.bed(find(in));
	depfriction = md.friction.C(find(in));
	p = polyfit(indepbed, depfriction, 1);
	fittedValues = polyval(p, indepbed);
	md.friction.C(find(in))= p(1)*md.geometry.bed(find(in)) + p(2);
	md.friction.C(find(in & md.friction.C<0.1)) = 0.1;
	% Plot the data and the regression line
	figure;
	scatter(indepbed, depfriction, 'b'); % Scatter plot of the filtered data
	hold on;
	plot(indepbed, fittedValues, 'r', 'LineWidth', 2); % Regression line
	hold off;
	%	pos=find(~ContourToNodes(md.mesh.x,md.mesh.y,'/projects/domino/KNS/Exp/frontforcing/1920.exp',1) & md.mask.ice_levelset>0);
	%	md.friction.C(pos) = 0.1; %This is a test for the forced front. Depending on outcome, should consider a timeseries friction coefficient given 1946 and 1948 shallow sills are causing unrealistic stability in the advance.
	plotmodel(md,'data',md.friction.C);
	savemodel(org, md);

end % }}}
if perform(org,'Contemp-Calving-VonMISES') % {{{
	md = loadmodel(org, 'Contemp-RelaxationRestart');
	md.calving=calvingvonmises();
	md.calving.min_thickness = 0;
	md.calving.stress_threshold_groundedice=600000;
	md.calving.stress_threshold_floatingice=400e3;	
	md.timestepping.start_time=1700;
	md.timestepping.final_time=1900;
	md.settings.output_frequency=5;

	md.transient.requested_outputs={'default','IceVolume','TotalSmb','SmbMassBalance'}; % other options?

	name='Contemp-VonMises_7_0';
	md.miscellaneous.name = name; 

	md=solve(md,'Transient','runtimename',false,'loadonly',loadonly);

	savemodel(org,md);
end	% }}}
if perform(org, 'Contemp-Load-VonMISES') % {{{
	md=loadmodel(org,'Contemp-Calving-VonMISES');
	md=loadresultsfromcluster(md);
	alias
	savemodel(org,md);
end %}}}
	if perform(org,'Transient') % {{{
		md = loadmodel(org, 'PDD-Paleo');
		md.transient.issmb=1;
		md.transient.ismasstransport=1;
		md.transient.isstressbalance=1;
		md.transient.isthermal=0;
		md.transient.isgroundingline=1;
		% md.transient.isgia=0;
		md.transient.isesa=0;
		md.transient.isdamageevolution=0;
		md.transient.ismovingfront=1;
		md.transient.ishydrology=0;
		% md.transient.isslr=0;
		md.transient.isoceancoupling=0;
		% md.transient.iscoupler=0;

		md.calving=calvinghab();     
		md.calving.flotation_fraction=2;
		md.frontalforcings.meltingrate=zeros(md.mesh.numberofvertices,1);
		% md.levelset.spclevelset=NaN(md.mesh.numberofvertices,1);
		md.groundingline.migration='None';

		md.timestepping = timesteppingadaptive(md.timestepping);
		md.timestepping.start_time=900;
		md.timestepping.final_time=1900;
		md.timestepping.time_step_min=0.05; 
		md.timestepping.time_step_max=1;
		md.settings.output_frequency=5;

		md.inversion.iscontrol=0;
		md.transient.requested_outputs={'default','IceVolume','TotalSmb','SmbMassBalance'}; % other options?

		md.toolkits.DefaultAnalysis=bcgslbjacobioptions();
		md.verbose=verbose('all');
		md.cluster=barkla();

		name='frontforce';
		loadonly = 0;

		md.settings.waitonlock = 0;
		md.cluster.interactive = 0; %only needed if you are using the generic cluster
		md.miscellaneous.name = name; 

		md=solve(md,'Transient','runtimename',false,'loadonly',loadonly);

		clear depfriction fittedValues in indepbase loadonly name p
		savemodel(org,md);
	end % }}}
if perform(org,'Calving-VonMISES') % {{{
	md = loadmodel(org, 'Load-ReadvanceToLIA');
	md = transientrestart(md); % reinitialise model from end of the relaxation
	md.results = rmfield(md.results,'TransientSolution2');

	md.calving=calvinghab();
	md.calving.min_thickness = 0;
	md.calving.stress_threshold_groundedice = 500;
	
	md.timestepping.start_time=1700;
	md.timestepping.final_time=1900;
	md.settings.output_frequency=5;

	md.transient.requested_outputs={'default','IceVolume','TotalSmb','SmbMassBalance'}; % other options?

	name='VonMises_test';
	md.miscellaneous.name = name; 

	md=solve(md,'Transient','runtimename',false,'loadonly',loadonly);

	savemodel(org,md);
end	% }}}
	if perform(org,'PDD-Paleo900-1900') % {{{
		md = loadmodel(org, 'LevelsetSeries');
		begin=900; % Year CE, start of forcing period
		finish=1900; %Year CE, end of forcing period
		n=401-((1950-begin)/50); %Convert from year CE to index
		m=401-((1950-finish)/50);
		span=m-n;		
		x = [begin+1/12:1/12:finish+1/12];

		load('PDD_presentday.mat');
		% Stride every 50 years - Time Resolution of UW Product
		a = x(1:600:end);b = x(2:600:end);c = x(3:600:end);d = x(4:600:end);e = x(5:600:end); f = x(6:600:end); g = x(7:600:end); h = x(8:600:end); ii = x(9:600:end);j = x(10:600:end); k = x(11:600:end); l = x(12:600:end);

		xx = [];
		for i =1:span;
			xx = [xx;a(i);b(i);c(i);d(i);e(i);f(i);g(i);h(i);ii(i);j(i);k(i);l(i)];
		end

		load ('/projects/domino/KNS/PDD_presentday.mat');
		ncbadg='/projects/domino/KNS/Trunk/datasets/paleoran/briner2020recons.nc';
		pre1 = ncread(ncbadg,'P_moderate'); %Precipiation, options: P_low, P_moderate, P_high
		pre1 = permute(pre1,[4 3 2 1]);
		pre1 = pre1(n:m,:,:,:);
		latB  = ncread(ncbadg,'lat');
		lon  = ncread(ncbadg,'lon');
		x = find(lon>180);
		lonB(x) = lon(x)-360;
		pre_jessica = NaN(md.mesh.numberofvertices+1,span,12);
		disp('Interpolating precipitation output onto grid');
		for i = 1:span
			for k = 1:12
				pre_jessica(1:md.mesh.numberofvertices,i,k)=(InterpFromGrid(latB,lonB',squeeze(pre1(i,k,:,:))',md.mesh.lat,md.mesh.long,'linear'));
			end
		end

		precip = [];	% For precipitation, we multiply the Box climatology by the fraction from the Badgeley Reconstruction
		for i = 1:span
			for j = 1:12
				precip = [precip;precip_present(:,j)'.*pre_jessica(1:end-1,i,j)'];
			end
		end
		precip = permute(precip,[2 1]);
		precip=precip*(1/1000)*12; % Convert from mm/month WE to m/yr WE

		md.smb.isprecipscaled=0; %  This allows us to set the precip ourselves to what we have done above.
		md.smb.precipitations_reconstructed = NaN(md.mesh.numberofvertices+1,span*12); % Needs to be size numberofvertices+1, so the last row can hold the timestep
		md.smb.precipitations_reconstructed(1:end-1,:) = precip;
		md.smb.precipitations_reconstructed(end,:) = xx;
		md.smb.precipitations_reconstructed(md.smb.precipitations_reconstructed<0) = 0.01; % There is a possibility for negative precipitation. values.  Set anything less than 0 to low value.

		temp1 = ncread(ncbadg,'T_moderate'); %Temperature, options: T_low, T_moderate, T_high
		temp1 = permute(temp1,[4 3 2 1]);
		temp1 = temp1(n:m,:,:,:);
		latB  = ncread(ncbadg,'lat');
		lon  = ncread(ncbadg,'lon');
		x = find(lon>180);
		lonB(x) = lon(x)-360;
		temp_jessica = NaN(md.mesh.numberofvertices+1,span,12);
		for i = 1:span
			for k = 1:12
				temp_jessica(1:md.mesh.numberofvertices,i,k)=(InterpFromGrid(latB,lonB',squeeze(temp1(i,k,:,:))',md.mesh.lat,md.mesh.long,'linear'));
			end
		end

		tmp = [];  % For temperature, we add the Box climatology to the anomalies from the Badgeley Reconstruction
		for i = 1:span
			for j = 1:12
				tmp = [tmp;temp_present(:,j)'+temp_jessica(1:end-1,i,j)'];
			end
		end
		tmp = permute(tmp,[2 1]);
		tmp=tmp+274.15; %Convert from Celsius to Kelvin
		md.smb.istemperaturescaled=0;
		md.smb.temperatures_reconstructed= NaN(md.mesh.numberofvertices+1,span*12);
		md.smb.temperatures_reconstructed(1:end-1,:) = tmp;
		md.smb.temperatures_reconstructed(end,:) = xx;
		md.smb.rlaps=0;

		savemodel(org,md);
	end % }}}

	if perform(org,'Solve-Forcefront') % {{{
		md = loadmodel(org, 'PDD-Paleo900-1900');
		md.transient.issmb=1;
		md.transient.ismasstransport=1;
		md.transient.isstressbalance=1;
		md.transient.ismovingfront=1;
		md.transient.isgroundingline=1;
		% md.transient.isgia=0;
		md.transient.isthermal=0;
		md.transient.isdamageevolution=0;
		md.transient.ishydrology=0;
		% md.transient.isslr=0;
		md.transient.isesa=0;
		md.transient.isoceancoupling=0;
		% md.transient.iscoupler=0;

		md.calving=calvinghab();   %This is a placeholder, will not be used while forcing the front with spclevelset  
		md.calving.flotation_fraction=2;
		md.frontalforcings.meltingrate=zeros(md.mesh.numberofvertices,1);
		% md.levelset.spclevelset=NaN(md.mesh.numberofvertices,1);
		md.groundingline.migration='None';

		begin=900;
		finish=1900;
		md.timestepping = timesteppingadaptive(md.timestepping);
		md.timestepping.start_time=begin;
		md.timestepping.final_time=finish-50;
		md.timestepping.time_step_min=0.05; 
		md.timestepping.time_step_max=1;
		md.settings.output_frequency=5;

		md.inversion.iscontrol=0;
		md.transient.requested_outputs={'default','IceVolume','TotalSmb','SmbMassBalance'}; % other options?

		md.verbose=verbose('all');

		name='frontforce_testlowfriction';
		loadonly = 0;

		md.settings.waitonlock = 0;
		md.cluster.interactive = 0; %only needed if you are using the generic cluster
		md.miscellaneous.name = name; 

		md=solve(md,'Transient','runtimename',false);

		clear depfriction fittedValues in indepbase loadonly name p
		savemodel(org,md);
	end % }}}
if perform(org, 'ResetTo900') % {{{
	md = loadmodel(org, 'Inversion');
	md2=loadmodel('/projects/domino/KNS/reference_model.mat');
	load('/projects/domino/KNS/GrIS_for_Domino.mat');
	pos = 3;
	data = [project2d(md2, GrIS_for_Domino(pos).Thickness, 4) project2d(md2, GrIS_for_Domino(pos).Surface, 4)];
	dataout = InterpFromMeshToMesh2d(md2.mesh.elements2d, md2.mesh.x2d, md2.mesh.y2d, data, md.mesh.x, md.mesh.y);
	% if 0
	T = dataout(:,1);
	md.geometry.thickness=T; %Thickness

	md.geometry.base=md.geometry.bed; %because i don't have any floating ice
	pos=find(md.geometry.bed<0 & ContourToNodes(md.mesh.x,md.mesh.y,'/projects/domino/KNS/Exp/no-ice-mask_0603.exp',1));
	md.geometry.thickness(pos)=0;
	md.geometry.base(pos) = 0;
	loc=find(md.geometry.thickness<=10 & ~(ContourToNodes(md.mesh.x,md.mesh.y,'/projects/domino/KNS/Exp/no-ice-mask_0603.exp',1)));
	pos=find(md.geometry.bed>0 & ContourToNodes(md.mesh.x,md.mesh.y,'/projects/domino/KNS/Exp/no-ice-mask_0603.exp',1));
	md.geometry.thickness(loc)=10;
	md.geometry.thickness(pos)=10;
	md.geometry.surface=md.geometry.thickness+md.geometry.base;

	%	else
	%		S  = dataout(:,2);
	%		md.geometry.surface=S;
	%		pos=find(md.geometry.bed<0 & ContourToNodes(md.mesh.x,md.mesh.y,'/projects/domino/KNS/Exp/no-ice-mask_0603.exp',1));
	%		md.geometry.surface(pos)=0;
	%		md.geometry.base(pos)=0;
	%		loc=find(md.geometry.surface<=10 & ~(ContourToNodes(md.mesh.x,md.mesh.y,'/projects/domino/KNS/Exp/no-ice-mask_0603.exp',1)));
	%		pos=find(md.geometry.bed>0 & ContourToNodes(md.mesh.x,md.mesh.y,'/projects/domino/KNS/Exp/no-ice-mask_0603.exp',1));
	%		md.geometry.surface(loc)=10;
	%		md.geometry.surface(pos)=10;
	%		md.geometry.thickness = md.geometry.surface - md.geometry.base;
	%	end
	pos = find(md.mesh.vertexonboundary);
	md.masstransport.spcthickness(pos) = md.geometry.thickness(pos);

	savemodel(org,md);
end % }}}	
	if perform(org,'PDD-PaleoRelaxation') % {{{
		md = loadmodel(org, 'PDD-Present');
		begin=850; % Year CE, start of forcing period
		finish=950; %Year CE, end of forcing period
		n=401-((1950-begin)/50); %Convert from year CE to index
		m=401-((1950-finish)/50);
		span=m-n;		
		x = [begin+1/12:1/12:finish+1/12];
		load('PDD_presentday.mat');
		% Stride every 50 years - Time Resolution of UW Product
		a = x(1:600:end);b = x(2:600:end);c = x(3:600:end);d = x(4:600:end);e = x(5:600:end); f = x(6:600:end); g = x(7:600:end); h = x(8:600:end); ii = x(9:600:end);j = x(10:600:end); k = x(11:600:end); l = x(12:600:end);

		xx = [];
		for i =1:span;
			xx = [xx;a(i);b(i);c(i);d(i);e(i);f(i);g(i);h(i);ii(i);j(i);k(i);l(i)];
		end

		ncbadg='/projects/domino/KNS/Trunk/datasets/paleoran/briner2020recons.nc';
		pre1 = ncread(ncbadg,'P_moderate'); %Precipiation
		pre1 = permute(pre1,[4 3 2 1]);
		pre1 = pre1(n:m,:,:,:);
		latB  = ncread(ncbadg,'lat');
		lon  = ncread(ncbadg,'lon');
		x = find(lon>180);
		lonB(x) = lon(x)-360;
		pre_jessica = NaN(md.mesh.numberofvertices+1,span,12);
		disp('Interpolating precipitation output onto grid');
		for i = 1:span
			for k = 1:12
				pre_jessica(1:md.mesh.numberofvertices,i,k)=(InterpFromGrid(latB,lonB',squeeze(pre1(i,k,:,:))',md.mesh.lat,md.mesh.long,'linear'));
			end
		end

		precip = [];	% For precipitation, we multiply the Box climatology by the fraction from the Badgeley Reconstruction
		for i = 1:span
			for j = 1:12
				precip = [precip;precip_present(:,j)'.*pre_jessica(1:end-1,i,j)'];
			end
		end
		precip = permute(precip,[2 1]);
		precip=precip*(1/1000)*12; % Convert from mm/month WE to m/yr WE

		md.smb.isprecipscaled=0; %  This allows us to set the precip ourselves to what we have done above.
		md.smb.precipitations_reconstructed = NaN(md.mesh.numberofvertices+1,span*12); % Needs to be size numberofvertices+1, so the last row can hold the timestep
		md.smb.precipitations_reconstructed(1:end-1,:) = precip;
		md.smb.precipitations_reconstructed(end,:) = xx;
		md.smb.precipitations_reconstructed(md.smb.precipitations_reconstructed<0) = 0.01; % There is a possibility for negative precipitation. values.  Set anything less than 0 to low value.

		temp1 = ncread(ncbadg,'T_moderate'); %Temperature
		temp1 = permute(temp1,[4 3 2 1]);
		temp1 = temp1(n:m,:,:,:);
		latB  = ncread(ncbadg,'lat');
		lon  = ncread(ncbadg,'lon');
		x = find(lon>180);
		lonB(x) = lon(x)-360;
		temp_jessica = NaN(md.mesh.numberofvertices+1,span,12);
		for i = 1:span
			for k = 1:12
				temp_jessica(1:md.mesh.numberofvertices,i,k)=(InterpFromGrid(latB,lonB',squeeze(temp1(i,k,:,:))',md.mesh.lat,md.mesh.long,'linear'));
			end
		end

		tmp = [];  % For temperature, we add the Box climatology to the anomalies from the Badgeley Reconstruction
		for i = 1:span
			for j = 1:12
				tmp = [tmp;temp_present(:,j)'+temp_jessica(1:end-1,i,j)'];
			end
		end
		tmp = permute(tmp,[2 1]);
		tmp=tmp+274.15; %Convert from Celsius to Kelvin
		md.smb.istemperaturescaled=0;
		md.smb.temperatures_reconstructed= NaN(md.mesh.numberofvertices+1,span*12);
		md.smb.temperatures_reconstructed(1:end-1,:) = tmp;
		md.smb.temperatures_reconstructed(end,:) = xx;
		md.smb.rlaps=0;

		savemodel(org,md);
	end % }}}
if perform(org, 'LevelsetSeries') % {{{
	md = loadmodel(org, 'Relax-Inversion');
	fronts = {'900','1200','LIA','1800','1920','1946','1948','2004'};
	spclevelset=zeros(md.mesh.numberofvertices+1,numel(fronts));
	for i=1:length(fronts);
		path=sprintf('/projects/domino/KNS/Exp/frontforcing/%s.exp',fronts{i});
		in = ContourToNodes(md.mesh.x,md.mesh.y,path,1) & md.geometry.bed<0;
		ice_levelset=-1*ones(md.mesh.numberofvertices,1);
		ice_levelset(find(in))=1;
		ice_levelset = reinitializelevelset(md, ice_levelset);
		if strcmp(fronts{i}, 'LIA')
			time = 1781;
		else
			time = str2num(fronts{i});
		end
		spclevelset(1:end-1,i) = ice_levelset;
		spclevelset(end,i)     = time;
	end
	md.levelset.spclevelset = spclevelset;
	savemodel(org, md);

end % }}}

if perform(org, 'EGU-model-rheo10inv') % {{{
	md = loadmodel('/projects/domino/KNS/Models/KNS_1100-Load-VonMISES.mat');
	md.materials.rheology_B = cuffey(273.15-10) * ones(md.mesh.numberofvertices, 1);
	
		disp('   Initialize basal friction using driving stress');
		disp('      -- Compute surface slopes and use 10 L2 projections');
		[sx,sy,s]=slope(md); sslope=averaging(md,s,10);
		disp('      -- Process surface velocity data');
		vel = md.inversion.vel_obs;
		flags=(vel==0); pos1=find(flags); pos2=find(~flags);
		vel(pos1) = griddata(md.mesh.x(pos2),md.mesh.y(pos2),vel(pos2),md.mesh.x(pos1),md.mesh.y(pos1));
		vel=max(vel,0.1);
		disp('      -- Calculate effective pressure');
		Neff = md.materials.rho_ice*md.geometry.thickness+md.materials.rho_water*md.geometry.base;
		Neff(find(Neff<=0))=1;
		disp('      -- Deduce friction coefficient');
		md.friction.coefficient=sqrt(md.materials.rho_ice*md.geometry.thickness.*(sslope)./(Neff.*vel/md.constants.yts));
		md.friction.coefficient=min(md.friction.coefficient,150);
		md.friction.p = 1.0 * ones(md.mesh.numberofelements,1);
		md.friction.q = 1.0 * ones(md.mesh.numberofelements,1);
		disp('      -- Extrapolate on ice free and floating ice regions');
		flags=(md.mask.ice_levelset>0) | (md.mask.ocean_levelset<0); pos1=find(flags); pos2=find(~flags);
		md.friction.coefficient(pos1) = 1;
		pos=find(isnan(md.friction.coefficient));
		md.friction.coefficient(pos)  = 1;
		pos=find(md.mask.ice_levelset<0 & md.initialization.vel==0);
		md.friction.coefficient(pos)=150;

		md.inversion=m1qn3inversion();
		md.inversion.vx_obs=md.initialization.vx;
		md.inversion.vy_obs=md.initialization.vy;
		md.inversion.vel_obs=md.initialization.vel;
		md.inversion.iscontrol=1;
		md.inversion.cost_functions=[101 103 501];
		md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,3);
		md.inversion.cost_functions_coefficients(:,1)=400;
		md.inversion.cost_functions_coefficients(:,2)=0.01;
		md.inversion.cost_functions_coefficients(:,3)=1e-7;

		pos = find(md.geometry.thickness<150 & md.geometry.thickness>10);
		md.inversion.cost_functions_coefficients(pos,1)=300;

		md.inversion.maxsteps = 80;
		md.inversion.maxiter  = 80;
		md.inversion.dxmin=0.001;
		%no cost where no ice
		pos = find(md.mask.ice_levelset>0 | md.inversion.vel_obs==0);
		md.inversion.cost_functions_coefficients(pos,[1]) = 0;

		md.inversion.control_parameters={'FrictionCoefficient'};
		md.inversion.min_parameters=0.1*ones(md.mesh.numberofvertices,1);
		md.inversion.max_parameters=150*ones(md.mesh.numberofvertices,1);

		%	pos = find(md.geometry.thickness<150 & md.geometry.thickness>10);
		%	md.inversion.max_parameters(pos)=50;

		md.stressbalance.restol=0.01;
		md.stressbalance.reltol=0.1;
		md.stressbalance.abstol=NaN;

		disp('Running inversion');
		md.verbose=verbose(0);
		md.verbose.control=1;
		md.settings.waitonlock=1;
		md.cluster=generic('name',oshostname(),'np',2);
		md.miscellaneous.name='invert';
		md=solve(md,'Stressbalance'); 

		md.friction.coefficient=md.results.StressbalanceSolution.FrictionCoefficient;

		disp(' --Running regression');
		in = md.geometry.bed < 0; %below sea-level
		indepbed = md.geometry.bed(find(in));
		depfriction = md.friction.coefficient(find(in));
		p = polyfit(indepbed, depfriction, 1);
		fittedValues = polyval(p, indepbed);
		md.friction.coefficient(find(in))= p(1)*md.geometry.bed(find(in)) + p(2);
		md.friction.coefficient(find(in & md.friction.coefficient<0.1)) = 0.1;
		% Plot the data and the regression line
		figure;
		scatter(indepbed, depfriction, 'b'); % Scatter plot of the filtered data
		hold on;
		plot(indepbed, fittedValues, 'r', 'LineWidth', 2); % Regression line
		hold off;
		%pos=find(~coefficientontourToNodes(md.mesh.x,md.mesh.y,'/projects/domino/KNS/Exp/frontforcing/1920.exp',1) & md.mask.ice_levelset>0);
		%md.friction.coefficient(pos) = 0.1; %This is a test for the forced front. Depending on outcome, should consider a timeseries friction coefficient given 1946 and 1948 shallow sills are causing unrealistic stability in the advance.
		plotmodel(md,'data',md.friction.coefficient);
		savemodel(org,md);

end % }}}

if perform(org,'EGU-model-RACMO-forcing') % {{{
	md=loadmodel(org,'EGU-model-rheo10inv');
   md.smb=SMBforcing();
	ncdata='/projects/domino/KNS/Trunk/datasets/RACMO23p2/Clipped/mean_smb_clipped_3413_1960-1990.nc';
	% finfo = ncinfo(ncdatasmb);
	xsmb= ncread(ncdata,'x');
	ysmb= ncread(ncdata,'y');
%	smb_month = [];
%	for i = 1:30 % 30 year SMB average from 1960 to 1990, note data is given in bands each month since 1958-01-15
%		varName = sprintf('Band%d', i); 
	smbBand = ncread(ncdata, 'Band1'); 
	md.smb.mass_balance = InterpFromGridToMesh(xsmb, ysmb, smbBand', md.mesh.x, md.mesh.y, 0);
%		smb_month = [smb_month md.smb.mass_balance]; 
%	end
	% Calculate the average value for each cell across all time steps
%	md.smb.mass_balance = mean(smb_month, 3);
	md.smb.mass_balance=md.smb.mass_balance/1000*(1000/917);	
	pos = find(md.smb.mass_balance<0);

	md.smb.mass_balance(pos) = md.smb.mass_balance(pos)*2; %artifical forcing, so that dH/dt = -(grad)uH + b (where dH is change in thickness, dt is change in time, (grad)uH is short-hand vector notation for mass transport, and b is SMB). When reaching an initialized steady state, dH/dt needs to be approx. 0
	savemodel(org,md);

end  % }}}
if perform(org, 'EGU-model_relax') % {{{
	md = loadmodel(org, 'EGU-model-RACMO-forcing');
	begin=300; % Year CE, start of forcing period
	finish=900; %YearCE, end of forcing period

%	md.friction.Cmax = 0.17 * ones(md.mesh.numberofvertices,1);


	name='egumod_rheo10-relax';
%	md.initialization.vel = md.results.StressbalanceSolution.Vel;
%	md.initialization.vx = md.results.StressbalanceSolution.Vx;
%	md.initialization.vy = md.results.StressbalanceSolution.Vy;

	pos = find(md.mesh.vertexonboundary);
	md.masstransport.spcthickness(pos) = md.geometry.thickness(pos);

	md.transient.issmb=0;
	md.transient.ismasstransport=1;
	md.transient.isstressbalance=1;
	md.transient.ismovingfront=0;
	md.transient.isgroundingline=1;
	% md.transient.isgia=0;
	md.transient.isthermal=0;
	md.transient.isdamageevolution=0;
	md.transient.ishydrology=0;
	% md.transient.isslr=0;
	md.transient.isesa=0;
	md.transient.isoceancoupling=0;
	% md.transient.iscoupler=0;

	md.calving=calvinghab();     %This is a placeholder, not used while no moving front.
	md.calving.flotation_fraction=0;
	md.frontalforcings.meltingrate=zeros(md.mesh.numberofvertices,1);
	md.groundingline.migration='SubelementMigration';
	md.groundingline.friction_interpolation='SubelementFriction1';
	md.groundingline.melt_interpolation='NoMeltOnPartiallyFloating';

	md.timestepping = timesteppingadaptive(md.timestepping);
	md.timestepping.start_time=begin;
	md.timestepping.final_time=finish;
	md.timestepping.time_step_min=0.000001; 
	md.timestepping.time_step_max=0.1;
	md.settings.output_frequency=5;

	md.inversion.iscontrol=0;
	md.transient.requested_outputs={'default','IceVolume'}; % other options?

	md.verbose=verbose('all');

	loadonly = 0;

	md.settings.waitonlock = 0;
	md.cluster.interactive = 0; %only needed if you are using the generic cluster
	md.miscellaneous.name = name; 

	%		pos=find(md.initialization.vel>6000);
	%		md.initialization.vel(pos)=6000;
	md=solve(md,'Transient','runtimename',false);


	savemodel(org,md);
end % }}}
if perform(org, 'EGU-model_relax-load') % {{{
	md=loadmodel(org,'EGU-model_relax');
	md=loadresultsfromcluster(md);
   % Check the CFL criteria and adjust adaptive timestep accordingly
   CFL = cfl_step(md,md.results.TransientSolution(end).Vx,md.results.TransientSolution(end).Vy)
   md.timestepping.time_step_min = 0.01; % 3.6 days -- This might be too big given CFL is 0.015

	md = transientrestart(md); % reinitialise model from end of the relaxation
	md.results = rmfield(md.results,'StressbalanceSolution');
	md.results = rmfield(md.results,'TransientSolution3');

	disp('Reset geometry to match mask after restart');
	md.geometry.base=md.geometry.bed; %because i don't have any floating ice
	pos=find(md.geometry.bed<0 & md.mask.ice_levelset>0);
	md.geometry.thickness(pos)=0;
	md.geometry.base(pos) = 0;
	loc=find(md.geometry.thickness<=10 & md.mask.ice_levelset<0);
	pos=find(md.geometry.bed>0 & md.mask.ice_levelset>0);
	md.geometry.thickness(loc)=10;
	md.geometry.thickness(pos)=10;
	md.geometry.surface=md.geometry.thickness+md.geometry.base;

	savemodel(org,md);
end %}}}
	if perform(org, 'EGU-model_PDD-Present') % {{{
		md = loadmodel(org, 'EGU-model_relax-load');
		[md.mesh.lat,md.mesh.long] = xy2ll(md.mesh.x,md.mesh.y,+1,39,71);

		% Set up SMB D180 class and set PDD factors

		md.smb = SMBd18opdd();  % Turn on Class
		md.smb.isd18opd=1;	% 1 means we want to load our monthly climatologies of temperature and precipitation 

		md.smb.delta18o = [-40.0110;0];  % This is meaningless, but needs to be set.  Basically, if you set md.smb.isd18opd=0, that would mean you want to take a modern climatology and scale it back through time based upon an ice core record.  For now, we just need to set this, but it will not be used since we set md.smb.isd18opd=1.

		md.smb.rlaps=6.0; % This is the Lapse rate you want to use
		md.smb.desfac = 1; % This is the elevation desertification effect.  0 is no change, 1 is maximum change.  Which is a reduction by a factor of 2 for every 1km change in surface elevation
		md.smb.rlapslgm = 5.5; % This just needs to be set, but is not used for your purposes

		% Set PDD Factors
		md.smb.issetpddfac = 1;
		md.smb.pddfac_snow = 4.3*ones(md.mesh.numberofvertices,1); % Notice these are spatially constant, but if you wanted could be spatially varying since they are on the vertices.
		md.smb.pddfac_ice = 8.3*ones(md.mesh.numberofvertices,1);

		disp('   ACC');
		ncdataacc='/projects/domino/KNS/Trunk/datasets/Box1840thr2012/accbox.nc';
		finfo = ncinfo(ncdataacc);
		xacc= ncread(ncdataacc,'x');
		yacc= ncread(ncdataacc,'y');

		md.smb.precipitations_presentday = NaN(md.mesh.numberofvertices,12);
		q=1;
		for i = 133:144 % 12 months
			index = i:12:((i+11)+12*(149)); %index across 1850 to 2000 for each month
			acc_month=[];
			for w= 1:150
				varName = sprintf('Band%d', index(w)); 
				accBand = ncread(ncdataacc, varName); 
				acc = InterpFromGridToMesh(xacc, yacc, accBand', md.mesh.x, md.mesh.y, 0);
				acc_month = cat(3, acc_month, acc); 
			end
			md.smb.precipitations_presentday(:,q)=mean(acc_month,3);
			q=q+1;
		end

		disp('   DEM');
		ncdatadem='/projects/domino/KNS/Trunk/datasets/Box1840thr2012/dembox.nc';
		% finfo = ncinfo(ncdatadem);
		xdem= ncread(ncdatadem,'x');
		ydem= ncread(ncdatadem,'y');
		dem=ncread(ncdatadem,'Band1');
		box_dem=InterpFromGridToMesh(xdem, ydem, dem', md.mesh.x, md.mesh.y, 0);
		md.smb.s0p = (max(box_dem,0));
		md.smb.s0t = (max(box_dem,0));

		disp('   Temp');
		ncdatatemp='/projects/domino/KNS/Trunk/datasets/Box1840thr2012/tempbox.nc';
		%finfo = ncinfo(ncdatatemp);
		xtemp= ncread(ncdatatemp,'x');
		ytemp= ncread(ncdatatemp,'y');

		md.smb.temperatures_presentday = NaN(md.mesh.numberofvertices,12);
		q=1;
		for i = 133:144 % 12 months
			index = i:12:((i+11)+12*(149)); %index across 1850 to 2000 for each month
			temp_month=[];
			for w= 1:150
				varName = sprintf('Band%d', index(w)); 
				tempBand = ncread(ncdatatemp, varName); 
				temp = InterpFromGridToMesh(xtemp, ytemp, tempBand', md.mesh.x, md.mesh.y, 0);
				temp_month = cat(3, temp_month, temp); 
			end
			md.smb.temperatures_presentday(:,q)=mean(temp_month,3);
			q=q+1;
		end
		temp_present = md.smb.temperatures_presentday;
		precip_present = md.smb.precipitations_presentday;
		save('PDD-presentday.mat', 'temp_present', 'precip_present')	

		md.smb.temperatures_presentday = md.smb.temperatures_presentday+274.15;
		md.smb.precipitations_presentday = md.smb.precipitations_presentday * (1/1000) * 12;
		savemodel(org,md);
	end % }}}	
if perform(org,'EGU-model_PDD-Paleo1100') % {{{
	md = loadmodel(org, 'EGU-model_PDD-Present');
	begin=1100; % Year CE, start of forcing period
	finish=1150; %Year CE, end of forcing period
	n=401-((1950-begin)/50); %Convert from year CE to index
	m=401-((1950-finish)/50);
	span=m-n;		

	load('PDD-presentday.mat');

	ncbadg='/projects/domino/KNS/Trunk/datasets/paleoran/briner2020recons.nc';
	pre1 = ncread(ncbadg,'P_moderate'); %Precipiation, options: P_low, P_moderate, P_high
	pre1 = permute(pre1,[4 3 2 1]);
	pre1 = pre1(n:m,:,:,:);
	latB  = ncread(ncbadg,'lat');
	lon  = ncread(ncbadg,'lon');
	x = find(lon>180);
	lonB(x) = lon(x)-360;
	pre_jessica = NaN(md.mesh.numberofvertices+1,span,12);
	disp('Interpolating precipitation output onto grid');
	for i = 1:span
		for k = 1:12
			pre_jessica(1:md.mesh.numberofvertices,i,k)=(InterpFromGrid(latB,lonB',squeeze(pre1(i,k,:,:))',md.mesh.lat,md.mesh.long,'linear'));
		end
	end

	precip = [];	% For precipitation, we multiply the Box climatology by the fraction from the Badgeley Reconstruction
	for i = 1:span
		for j = 1:12
			precip = [precip;precip_present(:,j)'.*pre_jessica(1:end-1,i,j)'];
		end
	end
	precip = permute(precip,[2 1]);
	precip=precip*(1/1000)*12; % Convert from mm/month WE to m/yr WE

	temp1 = ncread(ncbadg,'T_moderate'); %Temperature, options: T_low, T_moderate, T_high
	temp1 = permute(temp1,[4 3 2 1]);
	temp1 = temp1(n:m,:,:,:);
	latB  = ncread(ncbadg,'lat');
	lon  = ncread(ncbadg,'lon');
	x = find(lon>180);
	lonB(x) = lon(x)-360;
	temp_jessica = NaN(md.mesh.numberofvertices+1,span,12);
	for i = 1:span
		for k = 1:12
			temp_jessica(1:md.mesh.numberofvertices,i,k)=(InterpFromGrid(latB,lonB',squeeze(temp1(i,k,:,:))',md.mesh.lat,md.mesh.long,'linear'));
		end
	end

	tmp = [];  % For temperature, we add the Box climatology to the anomalies from the Badgeley Reconstruction
	for i = 1:span
		for j = 1:12
			tmp = [tmp;temp_present(:,j)'+temp_jessica(1:end-1,i,j)'];
		end
	end
	tmp = permute(tmp,[2 1]);
	tmp=tmp+274.15; %Convert from Celsius to Kelvin		

	x = [0+1/12:1/12:500+1/12];

	% Stride every 50 years - Time Resolution of UW Product
	a = x(1:600:end);b = x(2:600:end);c = x(3:600:end);d = x(4:600:end);e = x(5:600:end); f = x(6:600:end); g = x(7:600:end); h = x(8:600:end); ii = x(9:600:end);j = x(10:600:end); k = x(11:600:end); l = x(12:600:end);

	xx = [];
	for i =1:10;
		xx = [xx;a(i);b(i);c(i);d(i);e(i);f(i);g(i);h(i);ii(i);j(i);k(i);l(i)];
	end

	repprecip = [];
	reptemp = [];
	for i = 1:10
		repprecip = [repprecip precip];
		reptemp = [reptemp tmp];
	end

	md.smb.isprecipscaled=0; %  This allows us to set the precip ourselves to what we have done above.
	md.smb.precipitations_reconstructed = NaN(md.mesh.numberofvertices+1,10*12); % Needs to be size numberofvertices+1, so the last row can hold the timestep
	md.smb.precipitations_reconstructed(1:end-1,:) = repprecip;
	md.smb.precipitations_reconstructed(end,:) = xx;
	md.smb.precipitations_reconstructed(md.smb.precipitations_reconstructed<0) = 0.01; % There is a possibility for negative precipitation. values.  Set anything less than 0 to low value.

	md.smb.istemperaturescaled=0;
	md.smb.temperatures_reconstructed= NaN(md.mesh.numberofvertices+1,10*12);
	md.smb.temperatures_reconstructed(1:end-1,:) = reptemp;
	md.smb.temperatures_reconstructed(end,:) = xx;
	md.smb.rlaps=0;

	savemodel(org,md);
end % }}}

if perform(org, 'EGU-model_transient') % {{{
	md=loadmodel(org,'EGU-model_PDD-Paleo1100');
	in = ContourToNodes(md.mesh.x,md.mesh.y,'/projects/domino/KNS/Exp/frontforcing/900.exp',1) & md.geometry.bed<0;
	%	disp('Manually lowering PDD forcing');
	%	md.smb.temperatures_reconstructed(1:end-1,:) = md.smb.temperatures_reconstructed(1:end-1,:)+10; %10 degrees hotter
	%	md.smb.precipitations_reconstructed(1:end-1,:) = md.smb.precipitations_reconstructed(1:end-1,:)-0.3; %0.3 m/yr less precip
	md.transient.issmb=1;
	md.levelset.spclevelset = -1*ones(md.mesh.numberofvertices,1);
	md.levelset.spclevelset(find(in)) = 1;
	pos = find(md.mesh.vertexonboundary);
	md.levelset.spclevelset(:) = NaN;
	md.levelset.spclevelset(pos) = md.mask.ice_levelset(pos);	
	%md.levelset.spclevelset = reinitializelevelset(md, md.levelset.spclevelset);
	md.mask.ocean_levelset  = sethydrostaticmask(md);

   md.miscellaneous.name='1100-VonMises-g800kPa-f500kPa-icetemp10-contempiniti-meltoutput-bm15-fm0';
	md.calving=calvingvonmises();
	md.calving.stress_threshold_groundedice = 8e5 * ones(md.mesh.numberofvertices,1);
   md.calving.stress_threshold_floatingice = 5e5 * ones(md.mesh.numberofvertices,1);
   md.basalforcings.groundedice_melting_rate = 15 * ones(md.mesh.numberofvertices,1);
	md.calving.min_thickness = 0;
	%Notes: stable front for KNS g775kPa, f500kPa. It's insensitive to f at this point. Advance when g800kPa.
%	in = find(ContourToNodes(md.mesh.x, md.mesh.y, '/projects/domino/KNS/Exp/AS-calving.exp', 1));
%	md.calving.stress_threshold_groundedice(in) = 7e5;
%  md.calving.stress_threshold_floatingice(in) = 5e5;
	md.timestepping.start_time=0;
	md.timestepping.final_time=450;
	md.settings.output_frequency=5;
	md.transient.ismovingfront = 1;
	md.transient.requested_outputs={'default','IceVolume','TotalSmb','SmbSnowMelt','SmbMelt','SmbMassBalance','CalvingCalvingrate'}; % other options?

	loadonly = 0;
	md.inversion.iscontrol=0;
	md.verbose=verbose('all');
	md.settings.waitonlock = 0;
	md.cluster.interactive = 0; %only needed if you are using the generic cluster


	md=solve(md,'Transient','runtimename',false,'loadonly',loadonly);

   savemodel(org,md);
end % }}}
if perform(org, 'EGU-model_transient-load') % {{{
	md=loadmodel(org,'EGU-model_transient');
	md=loadresultsfromcluster(md);
	savemodel(org,md);
end %}}}

if perform(org,'1100-Calving-VonMISES') % {{{

%	md = loadmodel(org, 'Relaxation-Load-VonMISES');
	md = loadmodel(org, 'PDD-Paleo1100');
	
	in = ContourToNodes(md.mesh.x,md.mesh.y,'/projects/domino/KNS/Exp/frontforcing/900.exp',1) & md.geometry.bed<0;
	%	disp('Manually lowering PDD forcing');
	%	md.smb.temperatures_reconstructed(1:end-1,:) = md.smb.temperatures_reconstructed(1:end-1,:)+10; %10 degrees hotter
	%	md.smb.precipitations_reconstructed(1:end-1,:) = md.smb.precipitations_reconstructed(1:end-1,:)-0.3; %0.3 m/yr less precip
	md.transient.issmb=1;
	md.levelset.spclevelset = -1*ones(md.mesh.numberofvertices,1);
	md.levelset.spclevelset(find(in)) = 1;
	pos = find(md.mesh.vertexonboundary);
	md.levelset.spclevelset(:) = NaN;
	md.levelset.spclevelset(pos) = md.mask.ice_levelset(pos);	
	%md.levelset.spclevelset = reinitializelevelset(md, md.levelset.spclevelset);
	md.mask.ocean_levelset  = sethydrostaticmask(md);

	name='1100-VonMises-g800kPa-f600kPa-icetemp15-contempiniti-meltoutput';
	md.calving=calvingvonmises();
	md.calving.min_thickness = 0;
	md.calving.stress_threshold_groundedice = 8e5*ones(md.mesh.numberofvertices,1);
	md.calving.stress_threshold_floatingice = 6e5*ones(md.mesh.numberofvertices,1);	
	%Notes: stable front for KNS g775kPa, f500kPa. It's insensitive to f at this point. Advance when g800kPa.
%	in = find(ContourToNodes(md.mesh.x, md.mesh.y, '/projects/domino/KNS/Exp/AS-calving.exp', 1));
%	md.calving.stress_threshold_groundedice(in) = 7e5;
%  md.calving.stress_threshold_floatingice(in) = 5e5;
	md.timestepping.start_time=0;
	md.timestepping.final_time=450;
	md.settings.output_frequency=5;
	md.transient.ismovingfront = 1;
	md.transient.requested_outputs={'default','IceVolume','TotalSmb','SmbSnowMelt','SmbMelt','SmbMassBalance','CalvingCalvingrate'}; % other options?

	md.miscellaneous.name = name; 
	loadonly = 0;
	md.inversion.iscontrol=0;
	md.verbose=verbose('all');
	md.settings.waitonlock = 0;
	md.cluster.interactive = 0; %only needed if you are using the generic cluster


	md=solve(md,'Transient','runtimename',false,'loadonly',loadonly);

	savemodel(org,md);
end	% }}}
if perform(org, '1100-Load-VonMISES') % {{{
	md=loadmodel(org,'1100-Calving-VonMISES');
	md=loadresultsfromcluster(md);

%	cal = []; vel = [];
%	for i = 1:(md.results.TransientSolution(end).step/md.settings.output_frequency)
%		c = max(md.results.TransientSolution(i).CalvingCalvingrate);
%		v = max(md.results.TransientSolution(i).Vel);
%		cal = [cal c];
%		vel = [vel v];
%	end
%
%	time = [];
%	for i = 1:length(cal)
%		t =  md.results.TransientSolution(i).time;
%		time = [time t];
%	end
%	figure;
%	scatter(time, cal, 5,'r','filled');
%	hold on;
%	scatter(time, vel, 5, 'b','filled');
%	legend('Calving','Vel')
%	hold off

	savemodel(org,md);
end %}}}
if perform(org,'PDD-Paleo1800') % {{{
	md = loadmodel(org, '1100-Load-VonMISES');
	begin=1800; % Year CE, start of forcing period
	finish=1850; %Year CE, end of forcing period
	n=401-((1950-begin)/50); %Convert from year CE to index
	m=401-((1950-finish)/50);
	span=m-n;		

	load('PDD-presentday.mat');

	ncbadg='/projects/domino/KNS/Trunk/datasets/paleoran/briner2020recons.nc';
	pre1 = ncread(ncbadg,'P_moderate'); %Precipiation, options: P_low, P_moderate, P_high
	pre1 = permute(pre1,[4 3 2 1]);
	pre1 = pre1(n:m,:,:,:);
	latB  = ncread(ncbadg,'lat');
	lon  = ncread(ncbadg,'lon');
	x = find(lon>180);
	lonB(x) = lon(x)-360;
	pre_jessica = NaN(md.mesh.numberofvertices+1,span,12);
	disp('Interpolating precipitation output onto grid');
	for i = 1:span
		for k = 1:12
			pre_jessica(1:md.mesh.numberofvertices,i,k)=(InterpFromGrid(latB,lonB',squeeze(pre1(i,k,:,:))',md.mesh.lat,md.mesh.long,'linear'));
		end
	end

	precip = [];	% For precipitation, we multiply the Box climatology by the fraction from the Badgeley Reconstruction
	for i = 1:span
		for j = 1:12
			precip = [precip;precip_present(:,j)'.*pre_jessica(1:end-1,i,j)'];
		end
	end
	precip = permute(precip,[2 1]);
	precip=precip*(1/1000)*12; % Convert from mm/month WE to m/yr WE

	temp1 = ncread(ncbadg,'T_moderate'); %Temperature, options: T_low, T_moderate, T_high
	temp1 = permute(temp1,[4 3 2 1]);
	temp1 = temp1(n:m,:,:,:);
	latB  = ncread(ncbadg,'lat');
	lon  = ncread(ncbadg,'lon');
	x = find(lon>180);
	lonB(x) = lon(x)-360;
	temp_jessica = NaN(md.mesh.numberofvertices+1,span,12);
	for i = 1:span
		for k = 1:12
			temp_jessica(1:md.mesh.numberofvertices,i,k)=(InterpFromGrid(latB,lonB',squeeze(temp1(i,k,:,:))',md.mesh.lat,md.mesh.long,'linear'));
		end
	end

	tmp = [];  % For temperature, we add the Box climatology to the anomalies from the Badgeley Reconstruction
	for i = 1:span
		for j = 1:12
			tmp = [tmp;temp_present(:,j)'+temp_jessica(1:end-1,i,j)'];
		end
	end
	tmp = permute(tmp,[2 1]);
	tmp=tmp+274.15; %Convert from Celsius to Kelvin		

	x = [0+1/12:1/12:500+1/12];

	% Stride every 50 years - Time Resolution of UW Product
	a = x(1:600:end);b = x(2:600:end);c = x(3:600:end);d = x(4:600:end);e = x(5:600:end); f = x(6:600:end); g = x(7:600:end); h = x(8:600:end); ii = x(9:600:end);j = x(10:600:end); k = x(11:600:end); l = x(12:600:end);

	xx = [];
	for i =1:10;
		xx = [xx;a(i);b(i);c(i);d(i);e(i);f(i);g(i);h(i);ii(i);j(i);k(i);l(i)];
	end

	repprecip = [];
	reptemp = [];
	for i = 1:10
		repprecip = [repprecip precip];
		reptemp = [reptemp tmp];
	end

	md.smb.isprecipscaled=0; %  This allows us to set the precip ourselves to what we have done above.
	md.smb.precipitations_reconstructed = NaN(md.mesh.numberofvertices+1,10*12); % Needs to be size numberofvertices+1, so the last row can hold the timestep
	md.smb.precipitations_reconstructed(1:end-1,:) = repprecip;
	md.smb.precipitations_reconstructed(end,:) = xx;
	md.smb.precipitations_reconstructed(md.smb.precipitations_reconstructed<0) = 0.01; % There is a possibility for negative precipitation. values.  Set anything less than 0 to low value.

	md.smb.istemperaturescaled=0;
	md.smb.temperatures_reconstructed= NaN(md.mesh.numberofvertices+1,10*12);
	md.smb.temperatures_reconstructed(1:end-1,:) = reptemp;
	md.smb.temperatures_reconstructed(end,:) = xx;
	md.smb.rlaps=0;

	savemodel(org,md);
end % }}}
if perform(org,'1800-Calving-VonMISES') % {{{

%	md = loadmodel(org, 'Relaxation-Load-VonMISES');
	md = loadmodel(org, '1100-Load-VonMISES');

	md.smb=SMBforcing();
	md.smb.mass_balance = md.results.TransientSolution(end).SmbMassBalance;
	in = find(md.smb.mass_balance < 0);
	md.smb.mass_balance(in) = md.smb.mass_balance(in)*4;

	md = transientrestart(md); % reinitialise model from end of the relaxation
	md.results = rmfield(md.results,'StressbalanceSolution');
	md.results = rmfield(md.results,'TransientSolution3');
	
%	in = ContourToNodes(md.mesh.x,md.mesh.y,'/projects/domino/KNS/Exp/frontforcing/900.exp',1) & md.geometry.bed<0;
	%	disp('Manually lowering PDD forcing');
	%	md.smb.temperatures_reconstructed(1:end-1,:) = md.smb.temperatures_reconstructed(1:end-1,:)+10; %10 degrees hotter
	%	md.smb.precipitations_reconstructed(1:end-1,:) = md.smb.precipitations_reconstructed(1:end-1,:)-0.3; %0.3 m/yr less precip
	md.transient.issmb=0;
%	md.levelset.spclevelset = -1*ones(md.mesh.numberofvertices,1);
%	md.levelset.spclevelset(find(in)) = 1;
%	pos = find(md.mesh.vertexonboundary);
%	md.levelset.spclevelset(:) = NaN;
%	md.levelset.spclevelset(pos) = md.mask.ice_levelset(pos);	
	%md.levelset.spclevelset = reinitializelevelset(md, md.levelset.spclevelset);
%	md.mask.ocean_levelset  = sethydrostaticmask(md);

	name='1800-VonMises-g800kPa-f400kPa-Schoof-4x1100SMB-contempinit';
	md.calving=calvingvonmises();
	md.calving.min_thickness = 0;
	md.calving.stress_threshold_groundedice = 8e5*ones(md.mesh.numberofvertices,1);
	md.calving.stress_threshold_floatingice = 4e5*ones(md.mesh.numberofvertices,1);	
	%Notes: stable front for KNS g775kPa, f500kPa. It's insensitive to f at this point. Advance when g800kPa.
%	in = find(ContourToNodes(md.mesh.x, md.mesh.y, '/projects/domino/KNS/Exp/AS-calving.exp', 1));
%	md.calving.stress_threshold_groundedice(in) = 7e5;
%  md.calving.stress_threshold_floatingice(in) = 5e5;
	md.timestepping.start_time=0;
	md.timestepping.final_time=450;
	md.settings.output_frequency=5;
	md.transient.ismovingfront = 1;
%	md.transient.requested_outputs={'default','IceVolume','TotalSmb','SmbMassBalance','CalvingCalvingrate'}; % other options?
	md.transient.requested_outputs={'default','IceVolume','CalvingCalvingrate'}; % other options?

	md.miscellaneous.name = name; 
	loadonly = 0;
	md.inversion.iscontrol=0;
	md.verbose=verbose('all');
	md.settings.waitonlock = 0;
	md.cluster.interactive = 0; %only needed if you are using the generic cluster


	md=solve(md,'Transient','runtimename',false,'loadonly',loadonly);

	savemodel(org,md);
end	% }}}
if perform(org, '1800-Load-VonMISES') % {{{
	md=loadmodel(org,'1800-Calving-VonMISES');
	md=loadresultsfromcluster(md);

	cal = []; vel = [];
	for i = 1:(md.results.TransientSolution(end).step/md.settings.output_frequency)
		c = max(md.results.TransientSolution(i).CalvingCalvingrate);
		v = max(md.results.TransientSolution(i).Vel);
		cal = [cal c];
		vel = [vel v];
	end

	time = [];
	for i = 1:length(cal)
		t =  md.results.TransientSolution(i).time;
		time = [time t];
	end
%	figure;
%	scatter(time, cal, 5,'r','filled');
%	hold on;
%	scatter(time, vel, 5, 'b','filled');
%	legend('Calving','Vel')
%	hold off

	savemodel(org,md);
end %}}}

if perform(org,'test') % {{{
	md=loadmodel(org,'Inversion');

	md.smb = SMBd18opdd();  % Turn on Class
	ncdata='/projects/domino/KNS/air_MCruns_ensemble_mean_LMRv2.0_EPSG3413.nc';
	% finfo = ncinfo(ncdatasmb);
	xsmb= ncread(ncdata,'x');
	ysmb= ncread(ncdata,'y'); ysmb = flipud(ysmb);
	%	smb_month = [];
	%	for i = 1:30 % 30 year SMB average from 1960 to 1990, note data is given in bands each month since 1958-01-15
	%		varName = sprintf('Band%d', i); 
	smbBand = ncread(ncdata, 'air'); smbBand = smbBand(:,:,1); smbBand = flipud(smbBand); 
	md.smb.temperatures_reconstructed = InterpFromGridToMesh(xsmb, ysmb, smbBand', md.mesh.x, md.mesh.y, 0);
	%		smb_month = [smb_month md.smb.mass_balance]; 
	%	end
	% Calculate the average value for each cell across all time steps
	%	md.smb.mass_balance = mean(smb_month, 3);
	savemodel(org,md);

end  % }}}
if perform(org,'Meltrate') % {{{
	md = loadmodel(org, 'PDD-Paleo1800');
	md = transientrestart(md); % reinitialise model from end of the relaxation
	md.results = rmfield(md.results,'StressbalanceSolution');
	md.results = rmfield(md.results,'TransientSolution3');

	md.transient.issmb=1;

	md.frontalforcings.meltingrate = (0.5*365) * ones(md.mesh.numberofvertices,1);
	name='1800-VonMises-g800kPa-f400kPa-0.5mdmeltrate-contempinit-meltoutput';
	md.calving=calvingvonmises();
	md.calving.min_thickness = 0;
	md.calving.stress_threshold_groundedice = 8e5*ones(md.mesh.numberofvertices,1);
	md.calving.stress_threshold_floatingice = 4e5*ones(md.mesh.numberofvertices,1);	
	%Notes: stable front for KNS g775kPa, f500kPa. It's insensitive to f at this point. Advance when g800kPa.
%	in = find(ContourToNodes(md.mesh.x, md.mesh.y, '/projects/domino/KNS/Exp/AS-calving.exp', 1));
%	md.calving.stress_threshold_groundedice(in) = 7e5;
%  md.calving.stress_threshold_floatingice(in) = 5e5;
	md.timestepping.start_time=0;
	md.timestepping.final_time=450;
	md.settings.output_frequency=5;
	md.transient.ismovingfront = 1;
	md.transient.requested_outputs={'default','IceVolume','SmbSnowMelt','SmbMelt','TotalSmb','SmbMassBalance','CalvingCalvingrate'}; % other options?
%	md.transient.requested_outputs={'default','IceVolume','CalvingCalvingrate'}; % other options?

	md.miscellaneous.name = name; 
	loadonly = 0;
	md.inversion.iscontrol=0;
	md.verbose=verbose('all');
	md.settings.waitonlock = 0;
	md.cluster.interactive = 0; %only needed if you are using the generic cluster


	md=solve(md,'Transient','runtimename',false,'loadonly',loadonly);

	savemodel(org,md);
end % }}}
if perform(org, 'Load-Meltrate') % {{{
	md=loadmodel(org,'Meltrate');
	md=loadresultsfromcluster(md);

	cal = []; vel = [];
	for i = 1:(md.results.TransientSolution(end).step/md.settings.output_frequency)
		c = max(md.results.TransientSolution(i).CalvingCalvingrate);
		v = max(md.results.TransientSolution(i).Vel);
		cal = [cal c];
		vel = [vel v];
	end

	time = [];
	for i = 1:length(cal)
		t =  md.results.TransientSolution(i).time;
		time = [time t];
	end
%	figure;
%	scatter(time, cal, 5,'r','filled');
%	hold on;
%	scatter(time, vel, 5, 'b','filled');
%	legend('Calving','Vel')
%	hold off

	savemodel(org,md);
end %}}}
if perform(org,'Loadfromdisk') % {{{
	md1 = loadmodel('/projects/domino/KNS/Models/KNS_1100-Load-VonMISES.mat');
	md1=loadresultsfromdisk(md1,'/projects/domino/ISSM/execution/massbalance_senstest_seasonalbuttressing_stablels1__gr646kPa_fl384kPa_bm33_fm0/massbalance_senstest_seasonalbuttressing_stablels1__gr646kPa_fl384kPa_bm33_fm0.outbin');
%	md=loadmodel('/projects/domino/KNS/Models/KNS_RheologyLoad-Relaxation.mat');
%	md=loadresultsfromdisk(md1,'/projects/domino/ISSM/execution/10rheoB_stablels1_RACMOini_gr646kPa_fl384kPa_bm33_fm0/10rheoB_stablels1_RACMOini_gr646kPa_fl384kPa_bm33_fm0.outbin');
%	md1=loadmodel('/projects/domino/KNS/Models/KNS_Rheologyconvergence-calving-vonmises.mat');
%	md1=loadresultsfromdisk(md,'/projects/domino/ISSM/execution/15rheoB_stablels1_gr800kPa_fl400kPa_bm0_fm0/15rheoB_stablels1_gr800kPa_fl400kPa_bm0_fm0.outbin');

%	export2nc(md,md.geometry.base,'EGU-model.nc','base');
%	export2nc(md,md.geometry.thickness,'EGU-model.nc','thickness');
%	export2nc(md,md.geometry.surface,'EGU-model.nc','surface');
%	export2nc(md,md.initialization.temperature,'EGU-model.nc','temperature');
%	export2nc(md,md.initialization.vel,'EGU-model.nc','vel');
%	export2nc(md,md.initialization.vx,'EGU-model.nc','vx');
%	export2nc(md,md.initialization.vy,'EGU-model.nc','vy');
%	export2nc(md,md.initialization.vz,'EGU-model.nc','vz');
%	export2nc(md,md.levelset.spclevelset,'EGU-model.nc','spclevelset');
%	export2nc(md,md.mask.ice_levelset,'EGU-model.nc','ice_levelset');
%	export2nc(md,md.mask.ocean_levelset,'EGU-model.nc','ocean_levelset');
%	export2nc(md,md.masstransport.spcthickness,'EGU-model.nc','spcthickness');
%	export2nc(md,md.stressbalance.spcvx,'EGU-model.nc','spcvx');
%	export2nc(md,md.stressbalance.spcvz,'EGU-model.nc','spcvz');
%	export2nc(md,md.stressbalance.spcvy,'EGU-model.nc','spcvy');
%	export2nc(md,md.friction.coefficient,'EGU-model.nc','friccoef');

	savemodel(org,md);
end % }}}
if perform(org,'Joel-Parameterization') % {{{
	md = loadmodel('/projects/domino/KNS/Models/KNS_1100-Load-VonMISES.mat');
	md=loadresultsfromdisk(md,'/projects/domino/ISSM/execution/massbalance_senstest_seasonalbuttressing_stablels1__gr646kPa_fl384kPa_bm33_fm0/massbalance_senstest_seasonalbuttressing_stablels1__gr646kPa_fl384kPa_bm33_fm0.outbin');
%%	md=loadmodel('/projects/domino/KNS/Models/KNS_RheologyLoad-Relaxation.mat');
%	md=loadmodel(org,'Load-Relaxation');

%   ncEGU='/projects/domino/KNS/EGU-model.nc';
%   x= ncread(ncEGU,'x'); y=ncread(ncEGU,'y'); 
% 	thick=ncread(ncEGU,'thickness')'; base=ncread(ncEGU,'base')'; surf=ncread(ncEGU,'surface')';
%   md.geometry.surface=InterpFromGrid(x,y,surf,md.mesh.x,md.mesh.y,'linear');
%   md.geometry.base = InterpFromGrid(x,y,base,md.mesh.x,md.mesh.y, 'linear');
%   md.geometry.thickness = InterpFromGrid(x,y,thick,md.mesh.x,md.mesh.y, 'linear');

%	md.initialization.vel=md1.initialization.vel;	
%	md.initialization.vx=md1.initialization.vx; md.stressbalance.spcvx=md1.stressbalance.spcvx;
%	md.initialization.vy=md1.initialization.vy; md.stressbalance.spcvy=md1.stressbalance.spcvy;
%	md.initialization.vz=md1.initialization.vz; md.stressbalance.spcvz=md1.stressbalance.spcvz;

%	md.masstransport.spcthickness=md1.masstransport.spcthickness;
%	md.geometry.base=md1.geometry.base; md.geometry.surface=md1.geometry.surface; md.geometry.thickness=md1.geometry.thickness;

%	md.levelset.spclevelset=md1.levelset.spclevelset;
%	md.mask.ocean_levelset=md1.mask.ocean_levelset;
%	md.mask.ice_levelset=md1.mask.ice_levelset;
%	md.initialization.temperature=md1.initialization.temperature;
%	md.materials.rheology_B=md1.materials.rheology_B;

	md.materials.rheology_B=cuffey(273.15-10)*ones(md.mesh.numberofvertices, 1); 
%	md.calving=calvingvonmises();
%      md.calving.min_thickness = 0;
      md.calving.stress_threshold_groundedice=6.5e5;
      md.calving.stress_threshold_floatingice=4e5;
%      md.timestepping.start_time=0;
%      md.timestepping.final_time=500;
%      md.settings.output_frequency=5;
%      md.transient.ismovingfront = 1;
%      md.transient.requested_outputs={'default','IceVolume','TotalSmb','SmbMassBalance','CalvingCalvingrate'};;
 % other options?
		md.results = rmfield(md.results,'TransientSolution');

      md.miscellaneous.name = 'Joel-Model_test-10B';
      loadonly = 0;
      md.inversion.iscontrol=0;
      md.verbose=verbose('all');
      md.settings.waitonlock = 0;
      md.cluster.interactive = 0; %only needed if you are using the generic cluster


      md=solve(md,'Transient','runtimename',false,'loadonly',loadonly);

	savemodel(org, md);
end % }}}

if perform(org, 'Joel-Load') % {{{
	md=loadmodel(org,'Joel-Parameterization');
	md=loadresultsfromcluster(md);
	savemodel(org,md);
end %}}}

	if perform(org, 'Relax-Inversion') % {{{
		md = loadmodel(org, 'Load-Relaxation');
		disp('   Initialize basal friction using driving stress');
		disp('      -- Compute surface slopes and use 10 L2 projections');
		[sx,sy,s]=slope(md); sslope=averaging(md,s,10);
		disp('      -- Process surface velocity data');
		vel = md.inversion.vel_obs;
		flags=(vel==0); pos1=find(flags); pos2=find(~flags);
		vel(pos1) = griddata(md.mesh.x(pos2),md.mesh.y(pos2),vel(pos2),md.mesh.x(pos1),md.mesh.y(pos1));
		vel=max(vel,0.1);
		disp('      -- Calculate effective pressure');
		Neff = md.materials.rho_ice*md.geometry.thickness+md.materials.rho_water*md.geometry.base;
		Neff(find(Neff<=0))=1;
		disp('      -- Deduce friction coefficient');
		md.friction.coefficient=sqrt(md.materials.rho_ice*md.geometry.thickness.*(sslope)./(Neff.*vel/md.constants.yts));
		md.friction.coefficient=min(md.friction.coefficient,150);
		md.friction.p = 1.0 * ones(md.mesh.numberofelements,1);
		md.friction.q = 1.0 * ones(md.mesh.numberofelements,1);
		disp('      -- Extrapolate on ice free and floating ice regions');
		flags=(md.mask.ice_levelset>0) | (md.mask.ocean_levelset<0); pos1=find(flags); pos2=find(~flags);
		md.friction.coefficient(pos1) = 1;
		pos=find(isnan(md.friction.coefficient));
		md.friction.coefficient(pos)  = 1;
		pos=find(md.mask.ice_levelset<0 & md.initialization.vel==0);
		md.friction.coefficient(pos)=150;

		md.inversion=m1qn3inversion();
		md.inversion.vx_obs=md.initialization.vx;
		md.inversion.vy_obs=md.initialization.vy;
		md.inversion.vel_obs=md.initialization.vel;
		md.inversion.iscontrol=1;
		md.inversion.cost_functions=[101 103 501];
		md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,3);
		md.inversion.cost_functions_coefficients(:,1)=400;
		md.inversion.cost_functions_coefficients(:,2)=0.01;
		md.inversion.cost_functions_coefficients(:,3)=1e-7;

		md.inversion.maxsteps = 80;
		md.inversion.maxiter  = 80;
		md.inversion.dxmin=0.01;
		%no cost where no ice
		pos = find(md.mask.ice_levelset>0 | md.inversion.vel_obs==0);
		md.inversion.cost_functions_coefficients(pos,[1:2]) = 0;

		md.inversion.control_parameters={'FrictionCoefficient'};
		md.inversion.min_parameters=0.1*ones(md.mesh.numberofvertices,1);
		md.inversion.max_parameters=150*ones(md.mesh.numberofvertices,1);

		md.stressbalance.restol=0.01;
		md.stressbalance.reltol=0.1;
		md.stressbalance.abstol=NaN;

		plotmodel(md,'data',md.inversion.vel_obs);

		disp('Running inversion');
		% Solve
		md.verbose=verbose(0);
		md.verbose.control=1;
		md.settings.waitonlock=1;
		md.cluster=generic('name',oshostname(),'np',2);
		md.miscellaneous.name='invert';
		md=solve(md,'Stressbalance');

		md.friction.coefficient=md.results.StressbalanceSolution.FrictionCoefficient;

		disp(' --Running regression');
		%	in = md.geometry.bed < 0 & 5<md.friction.coefficient<145; disp('Based on basal elevation');
		in = md.initialization.vel > 100 & 5<md.friction.coefficient<145; disp('Based on velocity');
		indepbed = md.geometry.bed(find(in));
		depfriction = md.friction.coefficient(find(in));
		p = polyfit(indepbed, depfriction, 1);
		fittedValues = polyval(p, indepbed);
		pos = md.geometry.bed < 0;
		md.friction.coefficient(find(pos))= p(1)*md.geometry.bed(find(pos)) + p(2);
		md.friction.coefficient(find(pos & md.friction.coefficient<0.1)) = 0.1;
		% Plot the data and the regression line
		figure;
		scatter(indepbed, depfriction, 'b'); % Scatter plot of the filtered data
		hold on;
		plot(indepbed, fittedValues, 'r', 'LineWidth', 2); % Regression line
		hold off;
		%pos=find(~ContourToNodes(md.mesh.x,md.mesh.y,'/projects/domino/KNS/Exp/frontforcing/1920.exp',1) & md.mask.ice_levelset>0);
		%md.friction.coefficient(pos) = 0.1; %This is a test for the forced front. Depending on outcome, should consider a timeseries friction coefficient given 1946 and 1948 shallow sills are causing unrealistic stability in the advance.

		savemodel(org, md);

	end % }}}
		disp('  --spc nunatak');
	out=ContourToNodes(md.mesh.x,md.mesh.y,'/projects/domino/KNS/Exp/frontforcing/nunatak-ext.exp',1);
	in=ContourToNodes(md.mesh.x,md.mesh.y,'/projects/domino/KNS/Exp/frontforcing/nunatak-int.exp',1);
	md.stressbalance.spcvx(find(out))=md.initialization.vx(find(out)); md.stressbalance.spcvx(find(in))=NaN;
	md.stressbalance.spcvy(find(out))=md.initialization.vy(find(out)); md.stressbalance.spcvy(find(in))=NaN;
	md.stressbalance.spcvz(find(out))=md.initialization.vz(find(out)); md.stressbalance.spcvz(find(in))=NaN;
	plotmodel(md,'data','BC')
	expdisp('/projects/domino/KNS/Exp/frontforcing/nunatak-ext.exp')

	%How to deal with Nunataks
	md=loadresultsfromdisk(md,'frontforce.outbin');
	%%
	figure('WindowState', 'maximized')
	plotmodel(md,'data','transient_movie','transient_movie_field','Vel','log#all',10,'caxis#all',[1.5,6000])
	%%
	plotmodel(md,'data',md.results.TransientSolution(end).Vel,'title','Modelled velocity - first timestep', ...
		'caxis#all',([1.5,6000]),'log#all',10,'data',md.initialization.vel,'title','Initial velocity')
	%%
	plotmodel(md,'data',md.initialization.vel,'title','initial velocity (m/yr)', ...
		'caxis#1',([1.5,6000]),'log#1',10, 'data', md.geometry.thickness,'title','initial thickness (m)')
	%%
	figure('WindowState', 'maximized')
	plotmodel(md,'data','transient_movie','transient_movie_field','Thickness')
	%%
	plotmodel(md,'data',md.results.TransientSolution(end).Thickness);


cal = []; vel = [];
for i = 1:md.results.TransientSolution(end).step
	c = max(md.results.TransientSolution(i).CalvingCalvingrate);
	v = max(md.results.TransientSolution(i).Vel);
	cal = [cal c];
	vel = [vel v];
end

time = [];
for i = 1:length(cal)
	t =  md.results.TransientSolution(i).time;
	time = [time t];
end

scatter(time, cal, 5,'r','filled');
hold on;
scatter(time, vel, 5, 'b','filled');


	%DEBUG ====
	md.settings.waitonlock = 1;
	md.cluster.interactive = 1; %only needed if you are using the generic cluster
	md.cluster.np=30;
	md.verbose = verbose(0);
	md.verbose.solution = 1;
	md=solve(md,'Transient');
	error('STOP');
	%DEBUG ====

	pos=find(md.geometry.bed<0 & ContourToNodes(md.mesh.x,md.mesh.y,'/projects/domino/KNS/Exp/no-ice-mask_0603.exp',1));
	md.geometry.thickness(pos)=0;
	loc=find(md.geometry.thickness<=10 & ~(ContourToNodes(md.mesh.x,md.mesh.y,'/projects/domino/KNS/Exp/no-ice-mask_0603.exp',1)));
	pos=find(md.geometry.bed>0 & ContourToNodes(md.mesh.x,md.mesh.y,'/projects/domino/KNS/Exp/no-ice-mask_0603.exp',1));
	md.geometry.thickness(loc)=10;
	md.geometry.thickness(pos)=10;
	md.geometry.surface=md.geometry.thickness+md.geometry.base; %Surface
	
