classdef Kerr_img2 < hgsetget
    properties
        fName = '' % name of file
        
        % monitors
        Monitor1 = [];
        Monitor2 = [];
        signalMonitor = [];
        
        % Kerr rotation signal
        Kerr1 = [];
        Kerr2 = [];

        % parameters of the laser
        specAmp = [];
        specWidth = [];
        specCenter = [];
        
        % voltage of Hall sensors
        Hall1 = [];
        Hall2 = [];
        
        % corrected images (complex array of Kerr rotation)
        cKerr = []
        
        params
    end
    
    methods
        
        % constructor
        function obj = Kerr_img2
            disp('KERR_img2 object was created');
        end
        
        % open h5 file
        % params
        % fName is name of file
        function open(obj,varargin)
            p = inputParser();
            p.addParamValue('fName','',@isstr);
            p.parse(varargin{:});
            params = p.Results;
            
            % if strcmp(obj.fName,'')
            if isempty(params.fName)
                [fName,fPath,~] = uigetfile({'*.h5','*.*'});
                if (fName == 0)
                    suc = false;
                    return
                end
                suc = true;
                obj.fName = fullfile(fPath,fName);
            else
                obj.fName = params.fName;
            end
            
            % end
            
            info = h5info(obj.fName);
            h5disp(obj.fName);% scan datasets
            datasets = info.Datasets;
            
            if size(datasets,1)>0
                for datasetInd = 1:size(datasets,1)
                    switch datasets(datasetInd).Name
                        case 'Focus measure'
                            FM = h5read(obj.fName,'/Focus measure');
                            figure(1);
                            plot(FM(find(FM(:,1)>0),1),FM(find(FM(:,1)>0),2),'-rx');
                            xlabel('Focus distance (\num)');
                            ylabel('Focus measure (arb. units)');
                            
                        otherwise
                            disp('Unkwon dataset was found')
                    end
                end
            end
            
            obj.Monitor1 = [];
            obj.Monitor2 = [];
            obj.Kerr1 = [];
            obj.Kerr2 = [];
            
            imgGroup = h5info(obj.fName,'/images');
            
            % read images
            if (size(imgGroup.Groups,1)>0)
                for imgInd = 1:size(imgGroup.Groups,1)
                    params = struct();
                    groupName = imgGroup.Groups(imgInd).Name;
                    groupName = ['/images/',num2str(imgInd)];
                    params.xMin = h5readatt(obj.fName,groupName,'Initial x');
                    params.xMax = h5readatt(obj.fName,groupName,'Final x');
                    params.xSteps = h5readatt(obj.fName,groupName,'x steps');
                    params.yMin = h5readatt(obj.fName,groupName,'Initial y');
                    params.yMax = h5readatt(obj.fName,groupName,'Final y');
                    params.ySteps = h5readatt(obj.fName,groupName,'y steps');
                    params.ods = h5readatt(obj.fName,groupName,'ODS (mm)');
                    params.xScale = linspace(params.xMin,params.xMax,params.xSteps+1);
                    params.yScale = linspace(params.yMin,params.yMax,params.ySteps+1);
                    %params.delayLinePos = h5readatt(obj.fName,groupName,'ODS (mm)');
                    
                    obj.params{imgInd} = params;
                    
                    info =  h5info(obj.fName,groupName);
                    Names = info.Datasets;
                    for datasetInd = 1:size(Names,1)
                        dataset = Names(datasetInd);
                        switch dataset.Name
                            case 'monitor1'
                                obj.Monitor1(imgInd,:,:) = h5read(obj.fName,strcat(groupName,'/monitor1'));
                            case 'monitor2'
                                obj.Monitor2(imgInd,:,:) = h5read(obj.fName,strcat(groupName,'/monitor2'));
                            case 'kerr1'
                                obj.Kerr1(imgInd,:,:) = h5read(obj.fName,strcat(groupName,'/kerr1'));
                            case 'kerr2'
                                obj.Kerr2(imgInd,:,:) = h5read(obj.fName,strcat(groupName,'/kerr2'));
                            case 'spec amp'
                                obj.specAmp(imgInd,:,:) = h5read(obj.fName,strcat(groupName,'/spec amp'));
                            case 'spec width'
                                obj.specWidth(imgInd,:,:) = h5read(obj.fName,strcat(groupName,'/spec width'));
                            case 'spec center'
                                obj.specCenter(imgInd,:,:) = h5read(obj.fName,strcat(groupName,'/spec center'));
                            case 'Hall1'
                                obj.Hall1(imgInd,:,:) = h5read(obj.fName,strcat(groupName,'/Hall1'));
                            case 'Hall2'
                                obj.Hall2(imgInd,:,:) = h5read(obj.fName,strcat(groupName,'/Hall2')); 
                            case 'signal monitor'
                                obj.signalMonitor(imgInd,:,:) = h5read(obj.fName,strcat(groupName,'/signal monitor'));     
                                
                                
                        end
                    end
                end
            else
                disp('No images were found.')
            end
        end
        
        function plotAll(obj,varargin)
            p = inputParser();
            p.addParamValue('save',false,@islogical);
            p.parse(varargin{:});
            inParams = p.Results;
            
            % normalize images
            monLim = [min(min(obj.Monitor1(:)),min(obj.Monitor2(:)))...
                max(max(obj.Monitor1(:)),max(obj.Monitor2(:)))];
            monLim = [0.21 0.25]
            kerr1Lim = [min(obj.Kerr1(:)), max(obj.Kerr1(:))];
            kerr2Lim = [min(obj.Kerr2(:)), max(obj.Kerr2(:))];
            
            
            laserIntensity = obj.specAmp.*obj.specWidth;
            
            for imgInd = 1:size(obj.Monitor1,1)
                %  if ~any(size(obj.Kerr1)) && ~any(size(obj.Kerr2))
                
                % calculate rules
                
                params = obj.params{imgInd};
                xScale = linspace(params.xMin,params.xMax,params.xSteps);
                yScale = linspace(params.yMin,params.yMax,params.ySteps);
                
                

                hF1 = figure(imgInd);
                set(gca,'FontSize',16,'FontName','Times');
                subplot(222);
                    imagesc(xScale,yScale,squeeze(obj.Monitor1(imgInd,:,:)),monLim);
                    axis xy equal;
                    xlabel('x (\mum)','FontSize',22,'FontName','Times');
                    ylabel('y (\mum)','FontSize',22,'FontName','Times');
                    xlim([min(xScale) max(xScale)]);
                    title('Reflectivity 1','FontSize',22,'FontName','Times');
                    ylim([min(yScale) max(yScale)]);
                    t = colorbar('peer',gca);
                    set(get(t,'ylabel'),'FontSize',12,'FontName','Times','String', 'Voltage');
                
                set(gca,'FontSize',16,'FontName','Times');
                subplot(221);
                    imagesc(xScale,yScale,squeeze(obj.Kerr1(imgInd,:,:)),kerr1Lim);
                    axis xy equal;
                    xlabel('x (\mum)','FontSize',22,'FontName','Times');
                    ylabel('y (\mum)','FontSize',22,'FontName','Times');
                    xlim([min(xScale) max(xScale)]);
                    ylim([min(yScale) max(yScale)]);
                    title('Kerr 1','FontSize',22,'FontName','Times');
                    t = colorbar('peer',gca);
                    set(get(t,'ylabel'),'String', 'Voltage');
                    
                set(gca,'FontSize',16,'FontName','Times');
                subplot(223);
                    imagesc(xScale,yScale,squeeze(obj.Kerr2(imgInd,:,:)),kerr2Lim);
                    axis xy equal;
                    xlabel('x (\mum)','FontSize',22,'FontName','Times');
                    ylabel('y (\mum)','FontSize',22,'FontName','Times');
                    xlim([min(xScale) max(xScale)]);
                    ylim([min(yScale) max(yScale)]);
                    title('Kerr 2','FontSize',22,'FontName','Times');
                    t = colorbar('peer',gca);
                    set(get(t,'ylabel'),'String', 'Voltage');
                    colormap(copper)
                
                set(gca,'FontSize',16,'FontName','Times');    
                subplot(224);
                    imagesc(xScale,yScale,squeeze(laserIntensity(imgInd,:,:)));
                    axis xy equal;
                    xlabel('x (\mum)','FontSize',22,'FontName','Times');
                    ylabel('y (\mum)','FontSize',22,'FontName','Times');
                    xlim([min(xScale) max(xScale)]);
                    ylim([min(yScale) max(yScale)]);
                    title('Laser intensity','FontSize',22,'FontName','Times');
                    t = colorbar('peer',gca);
                    set(get(t,'ylabel'),'String', 'Voltage');
                    colormap(copper)    
                %
                %save images
                if inParams.save
                    [~,name,~] = fileparts(obj.fName)
                    imgName = strcat(name,'-kerr-',num2str(imgInd));
                    print(gcf,'-dpng',[imgName,'.png']);
                    savefig(gcf,[imgName,'.fig']);
                end    
                
                
                %   else
                if false
                    figure();
                    figTitle = [fName, ' Monitor 1. Focus distance is ',num2str(FM(imgInd,1)),' \mum,',...
                        'focus measure is ',num2str(FM(imgInd,2))];
                    imagesc(xScale,yScale,squeeze(m1(imgInd,:,:)),monLim);
                    axis xy square;
                    xlabel('x (\mum)','FontSize',14,'FontName','Times');
                    ylabel('y (\mum)','FontSize',14,'FontName','Times');
                    xlim([min(xScale) max(xScale)]);
                    ylim([min(yScale) max(yScale)]);
                    title(figTitle,'FontSize',14,'FontName','Times');
                    colorbar('peer',gca);
                    print(gcf,'-dpng',strcat(fName,'-',num2str(imgInd),'.png'));
                    savefig(gcf,strcat(fName,'-',num2str(imgInd),'.fig'));
                end
            end
        end
        
        function imgFFT(obj)
            i = 2;
            signal = squeeze(obj.Kerr1(i,:,:));
            monitor =  squeeze(obj.Monitor1(i,:,:));
            dx = abs(obj.params{i}.xMin-obj.params{i}.xMax)/size(signal,2);
            waveScale = linspace(-0.5/dx,0.5/dx,size(signal,2));
            
            figure(1);
            subplot(311);
            imagesc(obj.params{i}.xScale,obj.params{i}.yScale,monitor);
            xlabel('X, \mum'); ylabel('Y, \mum');
            title('Reflectivity','FontSize',14,'FontName','Times','FontWeight','bold');
            
            subplot(312);
            imagesc(obj.params{i}.xScale,obj.params{i}.yScale,signal);
            xlabel('X, \mum'); ylabel('Y, \mum');
            title('Kerr rotation','FontSize',14,'FontName','Times','FontWeight','bold');
            subplot(313);
            imagesc(waveScale,obj.params{i}.yScale,log10(fftshift(abs(fft(signal,[],2)),2)));
            xlabel('k_X, \mum^-^1'); ylabel('Y, \mum');
            title('FFT of Kerr rotation','FontSize',14,'FontName','Times','FontWeight','bold');
            xlim([-0.5 0.5]);
        end
        
        function waveFit(obj,varargin)
            p = inputParser();
            p.addParamValue('slices',[1 2],@isnumeric);
            p.parse(varargin{:});
            params = p.Results;
            
            
            Amp = squeeze(obj.Kerr1(1,1,:));
            A = [(max(Amp)-min(Amp))/2 0 2e-4]; % amplitude
            l = [4 3 6]; %wavelength
            k = [5 5 30]; % dampling length
            
            C = [1.5e-5 0 2e-4]; % background
            phi = [0.2 0 2*pi]; % phase
            xErf = [8 15 20];
            cErf = [1 1 10];
            
            options = saoptimset('PlotFcns',{@saplotbestx,...
                @saplotbestf,@saplotx,@saplotf},...
                'TolFun',1e-6);
            
            parArr = [A; l; phi; k; C; xErf; cErf];
            coeff = parArr(:,1);
            
            figure(1)
            coeff = simulannealbnd(@(x) err(x,Amp),coeff,parArr(:,2),parArr(:,3),options)
            h1 = figure();
            er = coeff(1)*erf((x-coeff(6))*coeff(7))+coeff(5);
            ep = coeff(1)*exp(-(x-x(1))/coeff(4))+coeff(5);
            plot(x,Amp,'-rx',x,func(coeff),x,er,x,ep);
            disp(['Error is ' num2str(err(coeff))]);
            %res(ind,:) = [num2str(err(coeff)) coeff];
            save res.mat res;
            
            function res = err(coeff,Amp)
                err = (Amp.' - func(coeff)).*(abs(Amp-coeff(5)).^1).'*1e9;
                res = sum(err.^2,2);
            end
            
            function res = func(coeff,Amp,x)
                A = coeff(1);
                l = coeff(2);
                phi = coeff(3);
                k = coeff(4);
                C = coeff(5);
                xErf = coeff(6);
                cErf = coeff(7);
                
                res = A*sin(2*pi*x./l+phi).*exp(-(x-x(1))/k).*erf((x-xErf)*cErf)+C;
            end
        end
        
        % plot one desired image
        % params
        %  index - number of the desired image (start from 1)
        %  shiftToZero - subtract constant value from the axis and set
        %                     origin of frame to zero
        %  normalize - divide Kerr signal on reflectivity signal
        function plotImg(obj,varargin)
            
            % read and parse params
            p = inputParser();
            
            p.addParamValue('index',1,@isnumeric);
            p.addParamValue('shiftToZero',false,@islogical);
            p.addParamValue('normalize',false,@islogical);
            p.addParamValue('xLim',[],@isnumeric);
            p.addParamValue('yLim',[],@isnumeric);
            p.addParamValue('scale',[],@isnumeric);
            p.addParamValue('clear',false,@islogical);
            p.addParamValue('blueMap',false,@islogical);
            p.addParamValue('subtract',false,@islogical);
            p.addParamValue('saveAs','',@isstr);
            
            p.parse(varargin{:});
            params = p.Results;
            
            % load parameters of the experiment
            expParams = obj.params{params.index};
            
            % create axis scales
            xScale = linspace(expParams.xMin,expParams.xMax,expParams.xSteps);
            yScale = linspace(expParams.yMin,expParams.yMax,expParams.ySteps)-0.29;
            
            if params.shiftToZero
                xScale = xScale -  expParams.xMin;
                yScale = yScale -  expParams.yMin;
            end
            
            if ~any(size(params.xLim) == [1 2])
                params.xLim = minmax(xScale);
            end
            
            if ~any(size(params.yLim) == [1 2])
                params.yLim = minmax(yScale);
            end
            
            % normalize (if required)
            if params.normalize
                mask = squeeze(obj.Monitor1(params.index,:,:));
                mask(find(mask<0.1*max(mask(:)))) = 0;
                kerrArr = squeeze(obj.Kerr1(params.index,:,:)).*mask;
            else
                kerrArr =  squeeze(obj.Kerr1(params.index,:,:));
            end
            
            % subtract background
            if params.subtract
                kerrArr = kerrArr - mean(kerrArr(:));
            end
            
            figure(1)
            if false
                subplot(211);
                
                imagesc(xScale,yScale,squeeze(obj.Monitor1(params.index,:,:)));
                axis xy equal;
                xlabel('x (\mum)','FontSize',14,'FontName','Times');
                ylabel('y (\mum)','FontSize',14,'FontName','Times');
                xlim(params.xLim); ylim(params.yLim);
                title(num2str(expParams.delayLinePos));
                title('Reflectivity','FontSize',14,'FontName','Times','FontWeight','bold');
                t = colorbar('peer',gca);
                set(get(t,'ylabel'),'FontSize',12,'FontName','Times','String', 'Voltage');
                
                subplot(212);
            end
            
            imagesc(xScale,yScale,kerrArr);
            axis xy equal;
            xlim(params.xLim); ylim(params.yLim);
            colormap(b2r(min(kerrArr(:)),max(kerrArr(:))));
            set(gca,'LineWidth',3,'TickLength',[0.015 0.015]);
            
            set(gca,'YTick',[0 5 10]);
            %set(gca,'XTick',[0 10 20 30 40]);
            if params.clear
                set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
            else
                xlabel('x (\mum)','FontSize',20,'FontName','Times','FontWeight','bold');
                ylabel('y (\mum)','FontSize',20,'FontName','Times','FontWeight','bold');
                title('Kerr rotation','FontSize',14,'FontName','Times','FontWeight','bold');
                set(get(colorbar('peer',gca),'ylabel'),'FontSize',12,'FontName','Times','String', 'Voltage');
                set(gca,'FontName','Times','FontSize',14,'FontWeight','bold');
            end
            
            % save image
            if ~strcmp(params.saveAs,'')
                print(gcf,'-depsc2',strcat(params.saveAs,'.eps'));
                print(gcf,'-dpng',strcat(params.saveAs,'.png'))
            end
            
            if false
                figure()
                lim = max(abs(kerrArr(:)));
                imagesc(xScale,yScale,kerrArr,[-lim lim]);
                axis xy equal;
                xlabel('x (\mum)','FontSize',14,'FontName','Times');
                ylabel('y (\mum)','FontSize',14,'FontName','Times');
                colormap(jet);
                xlim([min(xScale) max(xScale)]);
                ylim([min(yScale) max(yScale)]);
                
                print(gcf,'-depsc2',strcat(obj.fName,'-',num2str(params.index),'.eps'));
            end
        end
        
        
        % plot images of TRSCM and slices along OX axis
        function plotSlices(obj,varargin)
            % read and parse params
            p = inputParser();
            p.addParamValue('saveImg',false,@islogical);
            p.parse(varargin{:});
            params = p.Results;
            
            % normalize reflectivity
            reflect = sqrt(obj.Monitor1.^2+obj.Monitor2.^2);
            power = obj.specAmp.*obj.specWidth;
            power = power./max(power(:));
            normReflect = reflect./power;
            
            % find max and min of reflectivity, normalized reflectivity
            % and Kerr rotation
            refLimit = [min(reflect(:)), max(max(reflect(:)))];
            normRefLimit = [min(normReflect(:)), max(normReflect(:))];
            kerrLimit = [min(obj.Kerr1(:)), max(obj.Kerr1(:))];
            
            clf();
            for imgInd = 1:size(reflect,1)
                hF = figure(imgInd);
                
%TODO                xScale = linspace(obj.params{1,imgInd},)
                subplot(231);
                   imagesc(squeeze(reflect(imgInd,:,:)),refLimit);
                   axis xy square; 
                   title('Reflectivity','FontSize',20,'FontName','Times');
                   ylabel('y (\mum)'); xlabel('x (\mum)');
                   colorbar();

                subplot(234);
                   plot(squeeze(reflect(imgInd,:,:)));
                   xlabel('x (\mum)');
                   xlim([0 100]); ylim(refLimit);

                subplot(232);
                   imagesc(squeeze(normReflect(imgInd,:,:)),normRefLimit);
                   title('Normalized reflectivity','FontSize',20,'FontName','Times');
                   axis xy square;
                   xlabel('x (\mum)');

                subplot(235);
                   plot(squeeze(normReflect(imgInd,:,:)));
                   xlabel('x (\mum)');
                   xlim([0 100]); ylim(normRefLimit);

                subplot(233);
                   imagesc(squeeze(obj.Kerr1(imgInd,:,:)),kerrLimit);
                   title('Kerr rotation','FontSize',20,'FontName','Times');
                   axis xy square; xlabel('x (\mum)');

                subplot(236);
                   plot(squeeze(obj.Kerr1(imgInd,:,:)));   
                   xlabel('x (\mum)');
                   xlim([0 100]); ylim(kerrLimit);
                   
               if params.saveImg
                   [~,name,~] = fileparts(obj.fName); 
                   print(gcf,'-r600','-dpng',[name,'-',num2str(imgInd),'.png']);
                   %print(gcf,'-depsc2',[num2str(imgInd),'.eps']);
               end    

            end    
        end
        
        
        function plotMonitors(obj)
            
            % normalize images
            if ~any(size(obj.Monitor1)) && ~any(size(obj.Monitor2))
                monLim = [min(min(m1(:)),min(m2(:))) max(max(m1(:)),max(m2(:)))];
                monLim = [0 0.4]
            end
            
            laserIntensity = obj.specAmp.*obj.specWidth;
            
            for imgInd = 1:size(obj.Monitor1,1)
                %  if ~any(size(obj.Kerr1)) && ~any(size(obj.Kerr2))
                
                % calculate rules
                
                params = obj.params{imgInd};
                xScale = linspace(params.xMin,params.xMax,params.xSteps);
                yScale = linspace(params.yMin,params.yMax,params.ySteps);
                
                

                hF1 = figure(imgInd);
                
                subplot(221);
                    imagesc(xScale,yScale,squeeze(obj.Monitor1(imgInd,:,:)));
                    axis xy equal;
                    xlabel('x (\mum)','FontSize',10,'FontName','Times');
                    ylabel('y (\mum)','FontSize',10,'FontName','Times');
                    xlim([min(xScale) max(xScale)]);
                    title('Reflectivity 1','FontSize',10,'FontName','Times');
                    ylim([min(yScale) max(yScale)]);
                    t = colorbar('peer',gca);
                    set(get(t,'ylabel'),'FontSize',10,'FontName','Times','String', 'Voltage');
                
                subplot(222);
                    imagesc(xScale,yScale,squeeze(obj.Monitor2(imgInd,:,:)));
                    axis xy equal;
                    xlabel('x (\mum)','FontSize',10,'FontName','Times');
                    ylabel('y (\mum)','FontSize',10,'FontName','Times');
                    xlim([min(xScale) max(xScale)]);
                    title('Reflectivity 2','FontSize',10,'FontName','Times');
                    ylim([min(yScale) max(yScale)]);
                    t = colorbar('peer',gca);
                    set(get(t,'ylabel'),'FontSize',10,'FontName','Times','String', 'Voltage');    
                
                     
                subplot(223);
                    imagesc(xScale,yScale,squeeze(laserIntensity(imgInd,:,:)));
                    axis xy equal;
                    xlabel('x (\mum)','FontSize',10,'FontName','Times');
                    ylabel('y (\mum)','FontSize',10,'FontName','Times');
                    xlim([min(xScale) max(xScale)]);
                    ylim([min(yScale) max(yScale)]);
                    title('Laser intensity','FontSize',10,'FontName','Times');
                    t = colorbar('peer',gca);
                    set(get(t,'ylabel'),'FontSize',10,'FontName','Times','String', 'Voltage');
                    colormap(copper)    
                %print(gcf,'-dpng',strcat(fName,'-',num2str(imgInd),'-whole.png'));
                
                subplot(224);
                    imagesc(xScale,yScale,...
                        squeeze(obj.Monitor1(imgInd,:,:)+obj.Monitor2(imgInd,:,:))./squeeze(laserIntensity(imgInd,:,:)));
                    axis xy equal;
                    xlabel('x (\mum)','FontSize',10,'FontName','Times');
                    ylabel('y (\mum)','FontSize',10,'FontName','Times');
                    xlim([min(xScale) max(xScale)]);
                    ylim([min(yScale) max(yScale)]);
                    title('Normalized reflectivity','FontSize',14,'FontName','Times');
                    t = colorbar('peer',gca);
                    set(get(t,'ylabel'),'FontSize',10,'FontName','Times','String', 'Voltage');
                    colormap(copper)    
                
               [~,name,~] = fileparts(obj.fName)  
               print(gcf,'-dpng',strcat(name,'-',num2str(imgInd),'.png'));
               savefig(gcf,strcat(name,'-',num2str(imgInd),'.fig'));
            end
        end
        
        
        function rotateImage(obj,varargin)
            % read and parse params
            p = inputParser();
            p.addParamValue('imgInd',1,@isnumeric);
            p.parse(varargin{:});
            params = p.Results;
            
            
            expParams = obj.params{params.imgInd};
            xScale = linspace(expParams.xMin,expParams.xMax,expParams.xSteps);
            yScale = linspace(expParams.yMin,expParams.yMax,expParams.ySteps);
            
            % complex form of experimental signal
            cKerr = squeeze(obj.Kerr1(params.imgInd,:,:)+j*obj.Kerr2(params.imgInd,:,:));
            
            % вращение комплексного изображения
            ph = -4;
            amp0 = abs(cKerr); 
            phase0 = angle(cKerr);
            % complex form of corrected signal
            cKerr2 = amp0.*(cos(phase0+ph*2*pi/360)+j*sin(phase0+ph*2*pi/360));
            obj.cKerr(params.imgInd,:,:) = cKerr2;
            
            fH1 = figure(1);    
                subplot(241);
                    imagesc(xScale,yScale,squeeze(obj.Kerr1(params.imgInd,:,:)));
                    axis xy equal; colormap(copper)
                    xlabel('x (\mum)','FontSize',10,'FontName','Times');
                    ylabel('y (\mum)','FontSize',10,'FontName','Times');
                    xlim([min(xScale) max(xScale)]);
                    title('Kerr 1','FontSize',10,'FontName','Times');
                    ylim([min(yScale) max(yScale)]);
                    t = colorbar('peer',gca);
                    set(get(t,'ylabel'),'FontSize',10,'FontName','Times','String', 'Voltage');
                    colormap(flipud(copper));
                    freezeColors; cbfreeze;                    
                
                subplot(242);
                    imagesc(xScale,yScale,squeeze(obj.Kerr2(params.imgInd,:,:)));
                    axis xy equal;
                    xlabel('x (\mum)','FontSize',10,'FontName','Times');
                    ylabel('y (\mum)','FontSize',10,'FontName','Times');
                    xlim([min(xScale) max(xScale)]);
                    title('Kerr 2','FontSize',10,'FontName','Times');
                    ylim([min(yScale) max(yScale)]);
                    t = colorbar('peer',gca);
                    set(get(t,'ylabel'),'FontSize',10,'FontName','Times','String', 'Voltage');
                    colormap(flipud(copper));
                    freezeColors; cbfreeze;     
                    
                subplot(243);
                    imagesc(xScale,yScale,amp0);
                    axis xy equal; colormap(copper)
                    xlabel('x (\mum)','FontSize',10,'FontName','Times');
                    ylabel('y (\mum)','FontSize',10,'FontName','Times');
                    xlim([min(xScale) max(xScale)]);
                    title('Amplitude','FontSize',10,'FontName','Times');
                    ylim([min(yScale) max(yScale)]);
                    t = colorbar('peer',gca);
                    set(get(t,'ylabel'),'FontSize',10,'FontName','Times','String', 'Voltage');
                    colormap(jet); freezeColors; cbfreeze; 
                
                subplot(244);
                    imagesc(xScale,yScale,phase0,[-pi, pi]);
                    axis xy equal;
                    xlabel('x (\mum)','FontSize',10,'FontName','Times');
                    ylabel('y (\mum)','FontSize',10,'FontName','Times');
                    xlim([min(xScale) max(xScale)]);
                    title('Phase','FontSize',10,'FontName','Times');
                    ylim([min(yScale) max(yScale)]);
                    t = colorbar('peer',gca);
                    set(get(t,'ylabel'),'FontSize',10,'FontName','Times','String', 'Voltage');
                    colormap(hsv); freezeColors; cbfreeze;
                    
                subplot(245);
                    imagesc(xScale,yScale,real(cKerr2));
                    axis xy equal; colormap(copper)
                    xlabel('x (\mum)','FontSize',10,'FontName','Times');
                    ylabel('y (\mum)','FontSize',10,'FontName','Times');
                    xlim([min(xScale) max(xScale)]);
                    title('Kerr 1','FontSize',10,'FontName','Times');
                    ylim([min(yScale) max(yScale)]);
                    t = colorbar('peer',gca);
                    set(get(t,'ylabel'),'FontSize',10,'FontName','Times','String', 'Voltage');
                    colormap(flipud(copper));
                    freezeColors; cbfreeze;                    
                
                subplot(246);
                    imagesc(xScale,yScale,imag(cKerr2));
                    axis xy equal;
                    xlabel('x (\mum)','FontSize',10,'FontName','Times');
                    ylabel('y (\mum)','FontSize',10,'FontName','Times');
                    xlim([min(xScale) max(xScale)]);
                    title('Kerr 2','FontSize',10,'FontName','Times');
                    ylim([min(yScale) max(yScale)]);
                    t = colorbar('peer',gca);
                    set(get(t,'ylabel'),'FontSize',10,'FontName','Times','String', 'Voltage');
                    colormap(flipud(copper));
                    freezeColors; cbfreeze;     
                    
                subplot(247);
                    imagesc(xScale,yScale,abs(cKerr2));
                    axis xy equal; colormap(copper)
                    xlabel('x (\mum)','FontSize',10,'FontName','Times');
                    ylabel('y (\mum)','FontSize',10,'FontName','Times');
                    xlim([min(xScale) max(xScale)]);
                    title('Amplitude','FontSize',10,'FontName','Times');
                    ylim([min(yScale) max(yScale)]);
                    t = colorbar('peer',gca);
                    set(get(t,'ylabel'),'FontSize',10,'FontName','Times','String', 'Voltage');
                    colormap(jet); freezeColors; cbfreeze; 
                
                subplot(248);
                    imagesc(xScale,yScale,angle(cKerr2),[-pi, pi]);
                    axis xy equal;
                    xlabel('x (\mum)','FontSize',10,'FontName','Times');
                    ylabel('y (\mum)','FontSize',10,'FontName','Times');
                    xlim([min(xScale) max(xScale)]);
                    title('Phase','FontSize',10,'FontName','Times');
                    ylim([min(yScale) max(yScale)]);
                    t = colorbar('peer',gca);
                    set(get(t,'ylabel'),'FontSize',10,'FontName','Times','String', 'Voltage');
                    colormap(hsv); freezeColors; cbfreeze;
                    
              fH2 = figure(2)
                  lim = max(abs(cKerr(:)));
                  plot(cKerr);
                  xlim([-lim +lim]);
                  ylim([-lim +lim]);   

              fH3 = figure(3)   
                  lim = max(abs(cKerr2(:)));
                  plot(cKerr2);
                  xlim([-lim +lim]);
                  ylim([-lim +lim]);
             
            [~,shortFName,~] = fileparts(obj.fName);
             print(fH1,'-dpng',[shortFName,'kerrRot-img',num2str(params.imgInd),'.png']);     
             print(fH2,'-dpng',[shortFName,'kerrRot-img',num2str(params.imgInd),'plane1.png']);
             print(fH3,'-dpng',[shortFName,'kerrRot-img',num2str(params.imgInd),'plane2.png']);
                
        end    
        
    end
    
    methods (Access = protected)
    end
    
end

