function adjustImages(img,varargin)
    
    p = inputParser();
    p.addParamValue('saveImg',false,@islogical);
    p.addParamValue('saveMat',false,@islogical);
    p.addParamValue('supportImage',0,@isnumeric);
    p.parse(varargin{:});
    inParams = p.Results; 
    
    % fName = '170224-img2.h5';
     freq = 10.56e9;
     x1 = 1; x2 = 21;
     
     % range of final plot of averaged Kerr rotation
     t1 = 1; t2 = 6;

     

     clf();
     %img = Kerr_img2();
     %img.open('fName',fName);

     %power correction
     for phaseInd = 1:size(img.Monitor1,1)
         for sliceInd = 1:size(img.Monitor1,3)
             M1 = squeeze(img.Monitor1(phaseInd,:,sliceInd));
             M2 = squeeze(img.Monitor2(phaseInd,:,sliceInd));
             laserInt = squeeze(img.specAmp(phaseInd,:,sliceInd).*...
                 img.specWidth(phaseInd,:,sliceInd));
             [R1(phaseInd,:,sliceInd),R2(phaseInd,:,sliceInd)]=...
                 powerCorrection(M1,M2,laserInt);
         end
     end
     
     % correct spatial shift
     if inParams.supportImage ~=0
         for phaseInd = 1:size(img.Monitor1,1)
             for sliceInd = 1:size(img.Monitor1,3)
                 [~,shiftPos(phaseInd,sliceInd)]= max(fftshift(abs(ifft(fftshift(fft(R2(phaseInd,:,sliceInd))).*...
                conj(fftshift(fft(R2(supportImage,:,sliceInd))))))));
             end    
         end
         tmp = shiftPos(supportImage,:);
         % shift procedure
         for phaseInd = 1:size(img.Monitor1,1)
             shiftPos(phaseInd,:) = shiftPos(phaseInd,:) - tmp;
             for sliceInd = 1:size(R2,3)
                 % shift reflectivity signal
                 R3(phaseInd,:,sliceInd)=circshift(squeeze(R2(phaseInd,:,sliceInd)),...
                     [1,shiftPos(phaseInd,sliceInd)]);

                 % shift Kerr rotation signal
                 Kerr(phaseInd,:,sliceInd)=circshift(squeeze(img.Kerr1(phaseInd,:,sliceInd)),...
                     [1,shiftPos(phaseInd,sliceInd)]);
             end
         end
     else
         R3 = R2;
         Kerr = img.Kerr1;
     end
     
     
     % name of file (to save pictures)
     [~,f,~] = fileparts(img.fName); 
     
     imgInd = 1;
     % loop over phases
     for phaseInd = 1:size(img.Monitor1,1)
         % get params and create X,Y rules
         params = img.params{1,phaseInd};
         xLabel = linspace(params.xMin,params.xMax,params.xSteps+1);
         yLabel = linspace(params.yMin,params.yMax,params.ySteps+1);
         
         
         
         slices = squeeze(Kerr(phaseInd,:,x1:x2));
         meanSlices(phaseInd,:) = mean(slices,2).';
         
         
         figure(imgInd);
           imgInd = imgInd+1;
           subplot(221);
             imagesc(xLabel,yLabel,squeeze(R2(phaseInd,:,:)));
             set(gca,'FontSize',12,'FontName','Times');
             axis xy equal; colormap(copper);
             xlim([min(xLabel) max(xLabel)]);
             ylim([min(yLabel) max(yLabel)]);
             xlabel('x (\mum)','FontSize',14,'FontName','Times');
             ylabel('y (\mum)','FontSize',14,'FontName','Times');
             title('Raw reflectivity');
             colorbar();
           subplot(222);
             imagesc(xLabel,yLabel,squeeze(img.Kerr1(phaseInd,:,:)));
             set(gca,'FontSize',12,'FontName','Times');
             axis xy equal;
             xlim([min(xLabel) max(xLabel)]);
             ylim([min(yLabel) max(yLabel)]);
             xlabel('x (\mum)','FontSize',14,'FontName','Times');
             ylabel('y (\mum)','FontSize',14,'FontName','Times');
             title('Raw Kerr rotation');
             colorbar()
           subplot(223);
             imagesc(xLabel,yLabel,squeeze(R3(phaseInd,:,:)));
             set(gca,'FontSize',12,'FontName','Times');
             axis xy equal; colormap(copper);
             xlim([min(xLabel) max(xLabel)]);
             ylim([min(yLabel) max(yLabel)]);
             xlabel('x (\mum)','FontSize',14,'FontName','Times');
             ylabel('y (\mum)','FontSize',14,'FontName','Times');
             title('Corrected reflectivity');
             colorbar();
           subplot(224);
             imagesc(xLabel,yLabel,squeeze(Kerr(phaseInd,:,:)));
             set(gca,'FontSize',12,'FontName','Times');
             axis xy equal;
             xlim([min(xLabel) max(xLabel)]);
             ylim([min(yLabel) max(yLabel)]);
             xlabel('x (\mum)','FontSize',14,'FontName','Times');
             ylabel('y (\mum)','FontSize',14,'FontName','Times');
             title('Corrected Kerr rotation');
             colorbar()
          if inParams.saveImg
              print(gcf,'-r600','-dpng',[f,'-kerrImg1-',num2str(phaseInd),'.png'])
          end
          
          figure(imgInd);
            imgInd = imgInd+1;
            subplot(221);
             imagesc(xLabel,yLabel,squeeze(R3(phaseInd,:,:)));
             set(gca,'FontSize',12,'FontName','Times');
             axis xy equal; colormap(copper);
             xlim([min(xLabel) max(xLabel)]);
             ylim([min(yLabel) max(yLabel)]);
             xlabel('x (\mum)','FontSize',14,'FontName','Times');
             ylabel('y (\mum)','FontSize',14,'FontName','Times');
             title('Adjusted reflectivity');
             colorbar();
           subplot(223);
             imagesc(xLabel,yLabel,squeeze(Kerr(phaseInd,:,:)));
             set(gca,'FontSize',12,'FontName','Times');
             axis xy equal;
             xlim([min(xLabel) max(xLabel)]);
             ylim([min(yLabel) max(yLabel)]);
             xlabel('x (\mum)','FontSize',14,'FontName','Times');
             ylabel('y (\mum)','FontSize',14,'FontName','Times');
             title('Adjusted Kerr rotation');
             colorbar()
             
           subplot(2,2,[2,4]);
             hold on
             plot(yLabel,slices);
             plot(yLabel,meanSlices(phaseInd,:),'LineWidth',3);
             hold off
             set(gca,'FontSize',12,'FontName','Times');   
             xlabel('y (\mum)'); ylabel('Kerr rotation (arb. units)');
             xlim([min(yLabel) max(yLabel)]);
             title('Kerr rotation');
             set(gca,'yaxislocation','Right','box','on')
           
           if inParams.saveImg  
               print(gcf,'-r600','-dpng',[f,'-kerrImg2-',num2str(phaseInd),'.png']);
           end
     end
     
     figure(imgInd);
         plot(yLabel,meanSlices(t1:t2,:),'LineWidth',2);
         legend(arrayfun(@num2str,getAngles(img,t1,t2,freq)/(2*pi),'UniformOutput',false));
         xlabel('y (\mum)'); xlim([min(yLabel) max(yLabel)]);
         
         print(gcf,'-r600','-dpng',[f,'kerrSlices-SI-',num2str(inParams.supportImage),'.png']);
        
    if inParams.saveMat
        angles = getAngles(img,t1,t2,freq)/(2*pi);
        sl = meanSlices(t1:t2,:);

        save('slices.mat', 'sl','angles','yLabel')
    end    
end


    % linear interpolation of laser intensity
% correction for power fluctuation
% input:
%   Monitor1, Monitor2 are 2D arrays of reflectivity
%   laserInt is 2D array of laser intensity
% output:
%   reflectivity2 is corrected for laser power 
function [r1, r2] = powerCorrection(Monitor1,Monitor2,laserInt) 
    r1 = Monitor1+Monitor2;
    yScale = 1:size(laserInt,2);
    for phaseInd = 1:size(laserInt,1)
        coeff = polyfit(yScale,laserInt(phaseInd,:),1); %polynomial fit
        laserBaseLine = coeff(1)*yScale+coeff(2);
        %laserBaseLine(phaseInd,:) = laserBaseLine(phaseInd,:)/laserBaseLine(phaseInd,1); 
        signal = r1(phaseInd,:);
        laserBaseLine = laserBaseLine*signal(1)./max(laserBaseLine); 
        r2(phaseInd,:) = cleanCurve(signal-laserBaseLine);
    end
end

% clean curve
function res = cleanCurve(curve)
    res = curve-min(curve);
    res = res./max(res);
end

% read angles and calculate phase difference of images
% input:
%   t1, t2 are initial and final time indexes
%   img is Kerr_img2 object
%   freq is frequency of CW excitation
% output:
%   phaseArr is 1D array of relative phases
function res = getAngles(img,t1,t2,freq)
    dists = [];
    for phaseInd = t1:t2
        dists = [dists img.params{1,phaseInd}.ods];
    end
    timeArr = 2*dists/3e11;
    res = 2*pi*(timeArr - min(timeArr))*freq;
end