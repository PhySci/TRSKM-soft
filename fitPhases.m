% try to improve a bad experiment and exclude thermal drift of the sample
function fitPhases
    clf();
    % base layer
    baseLayer = 3;
    xCoord = 18;
    y1 = 42;
 
    freq = 10.16e9;
    % initial and final time moments
    t1 = 1; t2 = 6;

    img = Kerr_img2();
    img.open('fName','170214-img7.h5');
    yScale = linspace(img.params{1,1}.yMin,img.params{1,1}.yMax,...
        img.params{1,1}.ySteps+1);
     
    
    %fill up arrays of data
    Monitor1 = [];
    Monitor2 = [];
    laserInt = [];
    getSignalsArrs();
    
    % correct power fluctuations
    [reflectivity,reflectivity2] = powerCorrection();
    
    
    %% correstion of thermal drift
    for phaseInd = 1:size(reflectivity2)
        shift(phaseInd,:) = fftshift(abs(ifft(fftshift(fft(reflectivity2(phaseInd,:))).*...
            conj(fftshift(fft(reflectivity2(baseLayer,:)))))));
    end
    [~,shiftPos] = max(shift,[],2);
    shiftPos = shiftPos - shiftPos(baseLayer);
    for phaseInd = 1:size(reflectivity2)
        reflectivity3(phaseInd,:) = circshift(reflectivity2(phaseInd,:),...
            [1,-shiftPos(phaseInd,1)]);
        Kerr1_2(phaseInd,:) = circshift(Kerr1(phaseInd,:),...
            [1,-shiftPos(phaseInd,1)]);
        
    end    
    
    % recalculate Kerr signal to complex plot
    tSlice = Kerr1(:,y1); % slice along time (or phase)
    phaseArr = getAngles(img);
    a = tSlice.*abs(cos(phaseArr).');
    b = tSlice.*abs(sin(phaseArr).');
    cKerr = a+j*b;
    
    % PLOT IMAGES
    figure(1);
        subplot(311);
            plot(yScale,reflectivity);
            title('Original reflectivity signal'); 
            legend(arrayfun(@num2str,t1:t2,'UniformOutput',false));
        subplot(312);
            plot(yScale,reflectivity2);
            title('Power correction');
            legend(arrayfun(@num2str,t1:t2,'UniformOutput',false));
        subplot(313);
            plot(yScale,reflectivity3);
            title('Shift correction');
            legend(arrayfun(@num2str,t1:t2,'UniformOutput',false));     
            
            
            
    figure(2);
        plot(1:61,Kerr1_2);
        legend(arrayfun(@num2str,phaseArr,'UniformOutput',false));
        title('Kerr rotation signal');
        legend(arrayfun(@num2str,t1:t2,'UniformOutput',false));     
    
        
    figure(3);
    
        hold on
        for ind = 1:size(cKerr,1)
            cmap = lines;
            plot(real(cKerr(ind)),imag(cKerr(ind)),'o','MarkerEdgeColor',cmap(ind,:),...
                'MarkerFaceColor',cmap(ind,:));
        end    
        mVal = max(abs(cKerr));
        xlim([-mVal mVal]);
        ylim([-mVal mVal]);
        hold off
        legend(arrayfun(@num2str,phaseArr/(2*pi),'UniformOutput',false));
        axis xy square
        grid on
    
    
    %x = linspace(-0.1,2*pi+0.1);
    %f = 1.8*cos(x+0.1)*1e-4;    
    %figure(4);
    %    plot(phaseArr,tSlice,'ro',...
    %        x,f,'-b');
    %    xlim([-0.1 2*pi+0.1])
    
    % read angles and calculate phase difference of images
    % input:
    %   t1, t2 are initial and final time indexes
    %   img is Kerr_img2 object
    %   freq is frequency of CW excitation
    % output:
    %   phaseArr is 1D array of relative phases
    function res = getAngles(img)
        dists = [];
        for phaseInd = t1:t2
            dists = [dists img.params{1,phaseInd}.ods];
        end
        timeArr = 2*dists/3e11;
        res = 2*pi*(timeArr - min(timeArr))*freq;
    end


        % linear interpolation of laser intensity
    % correction for power fluctuation
    % input:
    %   Monitor1, Monitor2 are 2D arrays of reflectivity
    %   laserInt is 2D array of laser intensity
    % output:
    %   reflectivity2 is corrected for laser power 
    function [r1, r2] = powerCorrection() 
        r1 = Monitor1+Monitor2;
        for phaseInd = 1:size(laserInt,1)
            coeff = polyfit(yScale,laserInt(phaseInd,:),1); %polynomial fit
            laserBaseLine = coeff(1)*yScale+coeff(2);
            %laserBaseLine(phaseInd,:) = laserBaseLine(phaseInd,:)/laserBaseLine(phaseInd,1); 
            signal = r1(phaseInd,:);
            laserBaseLine = laserBaseLine*signal(1)./max(laserBaseLine); 
            r2(phaseInd,:) = cleanCurve(signal-laserBaseLine);
        end
    end


    function getSignalsArrs()
        Monitor1 = squeeze(img.Monitor1(t1:t2,:,xCoord));
        Monitor2 = squeeze(img.Monitor2(t1:t2,:,xCoord));
        laserInt = squeeze(img.specAmp(t1:t2,:,xCoord)).*...
            squeeze(img.specWidth(t1:t2,:,xCoord));
        Kerr1 = squeeze(img.Kerr1(t1:t2,:,xCoord));
    end
        
end

% clean curve
function res = cleanCurve(curve)
    res = curve-min(curve);
    res = res./max(res);
end


