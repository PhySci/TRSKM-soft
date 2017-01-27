function res= FitWave()
    global x1 Signal1
    
    fName = '160402-IMG-e.h5';
    img = Kerr_img2();
    img.fName = fName;
    img.open();
    %img.plotImg();
    
    
    imgInd = 1;
    y1 = 1;
    
    x0 = linspace(img.params{imgInd}.xMin,img.params{imgInd}.xMax,img.params{imgInd}.xSteps+1);
    Signal0 = mean(squeeze(img.Kerr1(imgInd,16:25,:)));
    
    x1 = x0(y1:end);
    Signal1 = Signal0(y1:end);
   
    
    A = [6.15521e-5 6.1e-5 6.2-5]; % amplitude 
    l = [19.336 18 21]; %wavelength
    k = [17.962 15 20]; % damping length
    
    % background and linear component
    C = [-3.789e-6 -3.9e-6 -3.7e-6];
    lin = [0 -1e-4 1e-4];
    
    % params of error function
    xErf = [13 10 16];
    cErf = [2 0 10];
    
    phi = [5.2 0 2*pi]; % phase
    
    options = saoptimset('PlotFcns',{@saplotbestx,...
                @saplotbestf,@saplotx,@saplotf,@saplottemperature},...
                'TolFun',1e-6,'InitialTemperature',10000); %,... % tolerance
                %'HybridFcn',@patternsearch,... %Hybrid function
                %'HybridInterval',10);
    
    parArr = [A; l; phi; k; C*1e3; xErf; cErf; lin*1e4];
    coeff = parArr(:,1);

    coeff = simulannealbnd(@err,coeff,parArr(:,2),parArr(:,3),options);
    
    h1 = figure();
    
    
    % function of linear dependence and background
    y(:,1) = coeff(5)/1e3+x0*coeff(8)*1e-4;
    
    % error function    
    y(:,2) = 0.5*coeff(1)*(1+erf((x0-coeff(6))/coeff(7)))+coeff(5)*1e-3+x0*coeff(8)*1e-4;
    
    % exponential decay
    y(:,3) = coeff(1)*exp(-(x0-x0(1))/coeff(4))+coeff(5)*1e-3+x0*coeff(8)*1e-4;
    y(:,4) = -coeff(1)*exp(-(x0-x0(1))/coeff(4))+coeff(5)*1e-3+x0*coeff(8)*1e-4;
    
    % whole curve
    y(:,5) = 0.5*coeff(1)*sin(2*pi*x0./coeff(2)+coeff(3)).*exp(-(x0-x0(1))/coeff(4)).*(1+erf((x0-coeff(6))/coeff(7)))+coeff(5)*1e-3+x0*coeff(8)*1e-4;
    
    plot(x0,Signal0,'ro',x0,y.');
    
    coeff(5) = coeff(5)*1e-3;
    coeff(8) = coeff(8)*1e-4; 
    
    text(50,0,num2str(coeff));
    legend('Experiment','Linear','Error function','Exponent','Exponent','Location','South');
    disp(['Error is ' num2str(err(coeff))]);        
end


function res = err(coeff)
    global Signal1
    err = (Signal1 - func(coeff)).*(abs(Signal1-coeff(5)/1e3).^1)*1e6;
    res = sum(err.^2,2);
end

function res = func(coeff)
    global x1
    res = 0.5*coeff(1)*sin(2*pi*x1./coeff(2)+coeff(3)).*exp(-(x1-x1(1))/coeff(4)).*(1+erf((x1-coeff(6))/coeff(7)))+coeff(5)*1e-3+x1*coeff(8)*1e-4;
end