function [output]   = DUSR(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by An Qi Zhang in MATLAB R2020b
% [output] = DUSR() => output = generates figure showing Fig1 oscillation
% at the nominal period
% [output] = DUSR('initial') => output = initial conditions in column vector
% [output] = DUSR('states') => output = state names in cell-array
% [output] = DUSR('parameters') => output = parameter names and values
% [output] = DUSR(time,vars) => output = time derivatives in column vector
%
% [output] = DUSR(D2AR, DAT, V0) => output = DAex, F value computed from the state variables
%
% [output]= DUSR(circReg,tranSim); => output = simulated ODEs
%   in which,
%          circReg [shape, amplitude]
%          tranSim [starttime, duration, amplitude]
%          output [struct(DAex, D2AR, TDA, V0, F)];
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameters = struct('alpha', 0.09, ...      uM
                 'kVmax',    2.63*3600, ... uM/h
                 'Km',       0.2, ...       uM
                 'beta',     144, ...       1/h
                 'D2tot',    0.1, ...       uM
                 'k',        10.46, ...     1/h/uM 
                 'a',        1.7, ...       1/h
                 'c',        3.62, ...      1/h
                 'b',        0.012, ...     mV
                 'kV',       2.73*3600, ... mV/h/uM
                 'Fmax',     15*3600, ...   1/h
                 'theta',    25, ...        mV
                 'sigma',    18, ...        mV
                 'deltaT',   1.8, ...       1
                 'D0',       0.04, ...      uM
                 'kT',       87.5, ...      1/uM
                 'tauT',     0.15 ...       h
        ); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HANDLE VARIABLE INPUT ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 0
    [t,x] = ode45( @(t,x) DUSR_ODE(t,x,parameters,0,0,0),...
        [0 24*3],DUSR('initial'));
    D2AR = x(:,1); TDA = x(:,2); V0 = x(:,3);
    DAF = DUSR(D2AR, TDA, V0);
    DAex = DAF(:,1); F = DAF(:,2); 
    output = figure(); hold on;
    
    yyaxis left
    pDA = plot(t,DAex,'-k','LineWidth', 2);
    pD2 = plot(t,D2AR,'-k');
    ylabel('$$D2_{\mathrm{AR}}$$, $$DA^{\mathrm{ex}}$$ ($\mu$M)','interpreter','latex')
    ylim([0 0.15])
    yticks([0 0.05 0.1 0.15]); yticklabels({'0','0.05','0.10','0.15'})
    
    yyaxis right
    F = F./parameters.Fmax;
    pDAT = plot(t,TDA,':k');
    pF = plot(t,F,'--k');
    ylabel('$$T_{\mathrm{DA}}$$, $$\frac{F}{F_{\mathrm{max}}}$$','interpreter','latex')
    xlim([0 12]); ylim([0 2]); yticks([0 1 2])
    
    legend([pDA,pD2,pDAT,pF],{'$$DA^{\mathrm{ex}}$$','$$D2_{\mathrm{AR}}$$',...
        '$$T_{\mathrm{DA}}$$','$$\frac{F}{F_{\mathrm{max}}}$$'},'Location','northoutside','Orientation','horizontal','interpreter','latex')
    ax1 = gca; ax1.YAxis(1).Color = 'k'; ax1.YAxis(2).Color = 'k';
    set(gca,'XTick',0:4:24) 
    
    xlabel('time (h)')
    return
elseif nargin == 1
	if strcmp(varargin{1},'initial')
		% Return state names in cell-array
		x0 = [0.1 1.1 0];
        output = x0;
    elseif strcmp(varargin{1},'stateVars')
		% Return state variable names as a cell array
		output = {'D2AR', 'TDA', 'V0'};
	elseif strcmp(varargin{1},'parameters')
		% Return parameter names and values as a struct
		output = parameters;
	else
		error('Invalid input argument. Please read the help text.');
    end
	return
elseif nargin == 2
    circReg = varargin{1};
    tranSim = varargin{2};
elseif nargin == 3
	D2AR = varargin{1}; D2AR = D2AR(:);
    TDA = varargin{2}; TDA = TDA(:);
    V0 = varargin{3}; V0 = V0(:);
    p = parameters;
    try 
        F = p.Fmax./(1+exp((p.theta-V0)/p.sigma));
        DAex = (p.alpha*F-p.beta*p.Km-p.kVmax*TDA+...
            sqrt((p.alpha*F-p.beta*p.Km-p.kVmax*TDA).^2 + 4*p.beta*p.alpha*F*p.Km))...
            /(2*p.beta);
    catch
        warning("Invalid input argument.")
        DAex = NaN; F = NaN;
    end
    output = [DAex F];
    return
else
	error('Invalid number of input arguments. Please read the help text.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIMULATION OF ODES BASED ON CONDITION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(circReg) && isempty(tranSim)
    % Simulate DUSR without circadian regulation or transient stimulus
    output = struct();
    [t,x] = ode45( @(t,x) DUSR_ODE(t,x,parameters,0,0,0),...
        [0 24*3],DUSR('initial'));
elseif isempty(circReg) && ~isempty(tranSim)
    % Simulate DUSR with circadian regulation but no transient stimulus
    tranSim_starttime = tranSim(1);
    tranSim_duration = tranSim(2);
    tranSim_amplitude = tranSim(3);
    tranSim_endtime = tranSim_starttime + tranSim_duration;
    
    % Restart simulation at start & end of tranSim
    [t1,x1] = ode45( @(t,x) DUSR_ODE(t,x,parameters,0,0,0),...
        [0 tranSim_starttime],DUSR('initial')); x0_2 = x1(end,:);
    [t2,x2] = ode45( @(t,x) DUSR_ODE(t,x,parameters,0,0,tranSim_amplitude),...
        [0 tranSim_duration],x0_2); x0_3 = x2(end,:);
    [t3,x3] = ode45( @(t,x) DUSR_ODE(t,x,parameters,0,0,0),...
        [0 24*3],x0_3);
    t = [t1;t2 + tranSim_starttime;t3 + tranSim_endtime];
    x = [x1;x2;x3];
elseif ~isempty(circReg) && isempty(tranSim)
    % Simulate DUSR without circadian regulation but with transient stimulus
    circReg_shape = circReg(1); % 1-block; 2-sine; 3-tilted sine
    circReg_Amp = circReg(2);
    % Restart simulation at start & end of tranSim
    if circReg_shape == 1
    elseif circReg_shape == 2 || circReg_shape == 3
        [t,x] = ode45( @(t,x) DUSR_ODE(t,x,parameters,circReg_shape,circReg_Amp,0),...
        [0 24*10],DUSR('initial'));
    else
        error("Invalid circadian regulation waveform.")
    end
elseif ~isempty(circReg) && ~isempty(tranSim)
    % Simulate DUSR with both circadian regulation and transient stimulus
    circReg_shape = circReg(1); % 1-block; 2-sine; 3-tilted sine
    circReg_Amp = circReg(2);
    tranSim_starttime = tranSim(1);
    tranSim_duration = tranSim(2);
    tranSim_amplitude = tranSim(3);
    tranSim_endtime = tranSim_starttime + tranSim_duration;
    
    if circReg_shape == 1
    elseif circReg_shape == 2 || circReg_shape == 3
        [t,x] = ode45( @(t,x) DUSR_ODE(t,x,parameters,circReg_shape,circReg_Amp,0),...
        [0 24*10],DUSR('initial'));
    else
        error("Invalid circadian regulation waveform.")
    end
else
    error('Invalid input combination.')
end
% record results
    D2AR = x(:,1); TDA = x(:,2); V0 = x(:,3);
    DAF = DUSR(D2AR, TDA, V0);
    output.t = t;
    output.DAex = DAF(:,1);
    output.D2AR = D2AR;
    output.TDA = TDA;
    output.V0 = V0;
    output.F = DAF(:,2);
end


function dxdt = DUSR_ODE(t, x, parameters, circShape, circAmp, tranSimAmp)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PARAMETERS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    alpha = parameters.alpha; kVmax = parameters.kVmax; Km = parameters.Km;
    beta = parameters.beta; D2tot = parameters.D2tot; k = parameters.k;
    a = parameters.a; c = parameters.c; b = parameters.b;
    kV = parameters.kV; Fmax = parameters.Fmax; theta = parameters.theta;
    sigma = parameters.sigma; deltaT = parameters.deltaT; 
    D0 = parameters.D0; kT = parameters.kT; tauT = parameters.tauT;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STATES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    StateVars = x;
    D2AR = StateVars(1);
    TDA = StateVars(2);
    V0 = StateVars(3);
    DAF = DUSR(D2AR, TDA, V0); DAex = DAF(1); F = DAF(2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CIRCADIAN REGULATION STRENGTH
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if circShape == 1 % block wave
        if rem(t,24) < 12
            circRegStrenth = circAmp;
        else
            circRegStrength = 0;
        end    
    elseif circShape == 2 % sine wave
            circRegStrength = circAmp/2*sin(2*pi*t/24)+circAmp/2;
    elseif circShape == 3 % tilted-sine wave
            xsin = 2*pi*t/24;
            circRegStrength = circAmp/2*sin(xsin-sin(xsin)/2)+circAmp/2;
    elseif circShape == 0
        circRegStrength = 0;
    else
        error("Invalid Circadian Regulation Waveform");
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DIFFERENTIAL EQUATIONS 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dD2AR_dt = k*(D2tot-D2AR)*DAex - a*D2AR;
    dTDA_dt = (1+(deltaT-1)/(1+exp(-kT*(D2AR-D0)))-TDA)/tauT;
    dV0_dt = -c*V0 + b*F - kV*D2AR - circRegStrength + tranSimAmp;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % RETURN VALUES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STATE ODEs
    dxdt(1) = dD2AR_dt;
    dxdt(2) = dTDA_dt;
    dxdt(3) = dV0_dt;
    dxdt = dxdt(:);
end
