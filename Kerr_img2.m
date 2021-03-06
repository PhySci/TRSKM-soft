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
    
    properties  (Access = protected)
        channels = {'Monitor1','Monitor2','Kerr1','Kerr2','realKerr','imagKerr'} 
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
            p.addParamValue('verbose',false,@islogical);
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
            if params.verbose
                h5disp(obj.fName);% scan datasets
            end    
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
            p.addParamValue('rotate',false,@islogical);
            p.parse(varargin{:});
            inParams = p.Results;
            
            
            % normalize images
            monLim = [min(min(obj.Monitor1(:)),min(obj.Monitor2(:)))...
                max(max(obj.Monitor1(:)),max(obj.Monitor2(:)))];
            %monLim = [0.21 0.25]
            kerr1Lim = [min(obj.Kerr1(:)), max(obj.Kerr1(:))];
            kerr2Lim = [min(obj.Kerr2(:)), max(obj.Kerr2(:))];
            
            plotSize = 4;
            
            
            laserIntensity = obj.specAmp.*obj.specWidth;
            
            for imgInd = 1:size(obj.Monitor1,1)
                %  if ~any(size(obj.Kerr1)) && ~any(size(obj.Kerr2))
                
                % calculate rules
                
                params = obj.params{imgInd};
                xScale = obj.getXScale();
                yScale = obj.getYScale();
                
                

                hF1 = figure(imgInd);
                clf(hF1)
                hS = []; %subplot handlers
                
                
                %subplot(411);
                %    hS = [hS imagesc(xScale,yScale,squeeze(obj.Monitor1(imgInd,:,:)),monLim)];
                %    axis xy equal;
                %    xlabel('x (\mum)','FontSize',22,'FontName','Times');
                %    ylabel('y (\mum)','FontSize',22,'FontName','Times');
                %    xlim([min(xScale) max(xScale)]);
                %    title('Reflectivity 1','FontSize',22,'FontName','Times');
                %    ylim([min(yScale) max(yScale)]);
                %    t = colorbar('peer',gca);
                %    set(get(t,'ylabel'),'FontSize',12,'FontName','Times','String', 'Voltage');
                
                if (plotSize == 2) 
                    subplot(121);
                        set(gca,'FontSize',16,'FontName','Times');
                        hS = [hS imagesc(xScale,yScale,squeeze(obj.Kerr1(imgInd,:,:)),kerr1Lim)];
                        axis xy equal;
                        xlabel('x (\mum)','FontSize',22,'FontName','Times');
                        ylabel('y (\mum)','FontSize',22,'FontName','Times');
                        xlim([min(xScale) max(xScale)]);
                        ylim([min(yScale) max(yScale)]);
                        title('Kerr 1','FontSize',22,'FontName','Times');
                        t = colorbar('peer',gca);
                        set(get(t,'ylabel'),'String', 'Voltage');

                    subplot(122);
                        set(gca,'FontSize',16,'FontName','Times');
                        hS = [hS imagesc(xScale,yScale,squeeze(obj.Kerr2(imgInd,:,:)),kerr2Lim)];
                        axis xy equal;
                        xlabel('x (\mum)','FontSize',22,'FontName','Times');
                        ylabel('y (\mum)','FontSize',22,'FontName','Times');
                        xlim([min(xScale) max(xScale)]);
                        ylim([min(yScale) max(yScale)]);
                        title('Kerr 2','FontSize',22,'FontName','Times');
                        t = colorbar('peer',gca);
                        set(get(t,'ylabel'),'String', 'Voltage');
                        colormap(copper)
                        
                elseif (plotSize == 4)
                     subplot(141);
                        set(gca,'FontSize',16,'FontName','Times');
                        hS = [hS imagesc(xScale,yScale,squeeze(obj.Monitor1(imgInd,:,:)),monLim)];
                        axis xy equal;
                        xlabel('x (\mum)','FontSize',22,'FontName','Times');
                        ylabel('y (\mum)','FontSize',22,'FontName','Times');
                        xlim([min(xScale) max(xScale)]);
                        ylim([min(yScale) max(yScale)]);
                        title('Monitor 1','FontSize',22,'FontName','Times');
                        t = colorbar('peer',gca);
                        set(get(t,'ylabel'),'String', 'Voltage');

                    subplot(142);
                        set(gca,'FontSize',16,'FontName','Times');
                        hS = [hS imagesc(xScale,yScale,squeeze(obj.Monitor2(imgInd,:,:)),monLim)];
                        axis xy equal;
                        xlabel('x (\mum)','FontSize',22,'FontName','Times');
                        ylabel('y (\mum)','FontSize',22,'FontName','Times');
                        xlim([min(xScale) max(xScale)]);
                        ylim([min(yScale) max(yScale)]);
                        title('Monitor 2','FontSize',22,'FontName','Times');
                        t = colorbar('peer',gca);
                        set(get(t,'ylabel'),'String', 'Voltage');
                        colormap(copper)
                        
                    subplot(143);
                        set(gca,'FontSize',16,'FontName','Times');
                        hS = [hS imagesc(xScale,yScale,squeeze(obj.Kerr1(imgInd,:,:)),kerr1Lim)];
                        axis xy equal;
                        xlabel('x (\mum)','FontSize',22,'FontName','Times');
                        ylabel('y (\mum)','FontSize',22,'FontName','Times');
                        xlim([min(xScale) max(xScale)]);
                        ylim([min(yScale) max(yScale)]);
                        title('Kerr 1','FontSize',22,'FontName','Times');
                        t = colorbar('peer',gca);
                        set(get(t,'ylabel'),'String', 'Voltage');

                    subplot(144);
                        set(gca,'FontSize',16,'FontName','Times');
                        hS = [hS imagesc(xScale,yScale,squeeze(obj.Kerr2(imgInd,:,:)),kerr2Lim)];
                        axis xy equal;
                        xlabel('x (\mum)','FontSize',22,'FontName','Times');
                        ylabel('y (\mum)','FontSize',22,'FontName','Times');
                        xlim([min(xScale) max(xScale)]);
                        ylim([min(yScale) max(yScale)]);
                        title('Kerr 2','FontSize',22,'FontName','Times');
                        t = colorbar('peer',gca);
                        set(get(t,'ylabel'),'String', 'Voltage');
                        colormap(copper)
                    
                end
                
                %set(gca,'FontSize',16,'FontName','Times');    
                %subplot(412);
                %    hS = [hS imagesc(xScale,yScale,squeeze(laserIntensity(imgInd,:,:)))];
                %    axis xy equal;
                %    xlabel('x (\mum)','FontSize',22,'FontName','Times');
                %    ylabel('y (\mum)','FontSize',22,'FontName','Times');
                %    xlim([min(xScale) max(xScale)]);
                %    ylim([min(yScale) max(yScale)]);
                %    title('Laser intensity','FontSize',22,'FontName','Times');
                %    t = colorbar('peer',gca);
                %    set(get(t,'ylabel'),'String', 'Voltage');
                %    colormap(copper)    
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
                %subplot(231);
                %   imagesc(squeeze(reflect(imgInd,:,:)),refLimit);
                %   axis xy square; 
                %   title('Reflectivity','FontSize',20,'FontName','Times');
                %   ylabel('y (\mum)'); xlabel('x (\mum)');
                %   colorbar();

                % raw reflectivity plot
                %subplot(221);
                %   plot(squeeze(reflect(imgInd,:,:)));
                %   xlabel('x (\mum)');
                %   xlim([0 100]); ylim(refLimit);

                % normalized reflectivity map   
                subplot(221);
                   imagesc(squeeze(normReflect(imgInd,:,:)),normRefLimit);
                   title('Normalized reflectivity','FontSize',20,'FontName','Times');
                   axis xy square;
                   xlabel('x (\mum)');

                % normalized reflectivity plot    
                subplot(223);
                   plot(squeeze(normReflect(imgInd,:,:)));
                   xlabel('x (\mum)');
                   xlim([0 100]); ylim(normRefLimit);
                   
                % Kerr map
                subplot(222);
                   imagesc(squeeze(obj.Kerr1(imgInd,:,:)),kerrLimit);
                   title('Kerr rotation','FontSize',20,'FontName','Times');
                   axis xy square; xlabel('x (\mum)');
                
                % Kerr plot   
                subplot(224);
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
                
                subplot(121);
                    imagesc(xScale,yScale,squeeze(obj.Monitor1(imgInd,:,:)));
                    axis xy equal;
                    xlabel('x (\mum)','FontSize',10,'FontName','Times');
                    ylabel('y (\mum)','FontSize',10,'FontName','Times');
                    xlim([min(xScale) max(xScale)]);
                    title('Reflectivity 1','FontSize',10,'FontName','Times');
                    ylim([min(yScale) max(yScale)]);
                    t = colorbar('peer',gca);
                    set(get(t,'ylabel'),'FontSize',10,'FontName','Times','String', 'Voltage');
                
                subplot(122);
                    imagesc(xScale,yScale,squeeze(obj.Monitor2(imgInd,:,:)));
                    axis xy equal;
                    xlabel('x (\mum)','FontSize',10,'FontName','Times');
                    ylabel('y (\mum)','FontSize',10,'FontName','Times');
                    xlim([min(xScale) max(xScale)]);
                    title('Reflectivity 2','FontSize',10,'FontName','Times');
                    ylim([min(yScale) max(yScale)]);
                    t = colorbar('peer',gca);
                    set(get(t,'ylabel'),'FontSize',10,'FontName','Times','String', 'Voltage');    
                
                 if false    
                subplot(143);
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
                
                subplot(144);
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
                 end
               [~,name,~] = fileparts(obj.fName)  
               print(gcf,'-dpng',strcat(name,'-',num2str(imgInd),'.png'));
               savefig(gcf,strcat(name,'-',num2str(imgInd),'.fig'));
            end
        end
        
        function rotateImage(obj,varargin)
            % read and parse params
            p = inputParser();
            p.addParamValue('imgInd',1,@isnumeric);
            p.addParamValue('angle',10,@isnumeric);
            p.parse(varargin{:});
            params = p.Results;
            
            % check params
            if params.imgInd > size(obj.Monitor1,1)
                error('Wrong input parameter: image index exceeds amount of images')
            end    
            
            
            expParams = obj.params{params.imgInd};
            xScale = linspace(expParams.xMin,expParams.xMax,expParams.xSteps);
            yScale = linspace(expParams.yMin,expParams.yMax,expParams.ySteps);
            
            % complex form of experimental signal
            cKerr = squeeze(obj.Kerr1(params.imgInd,:,:)+j*obj.Kerr2(params.imgInd,:,:));
            
            % вращение комплексного изображения
            ph = params.angle;
            amp0 = abs(cKerr); 
            phase0 = angle(cKerr);
            % complex form of corrected signal
            cKerr2 = amp0.*(cos(phase0+ph*2*pi/360)+j*sin(phase0+ph*2*pi/360));
            obj.cKerr(params.imgInd,:,:) = cKerr2;
            
            fH1 = figure(1);
                %clf(gcf) 
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
                    
              fH2 = figure(2);
                  clf(fH2)
                  lim = max(abs(cKerr(:)));
                  hold on
                  % add grid lines
                  plot([-lim, +lim],[0 0],':k')
                  plot([0 0],[-lim, +lim],':k')
                  plot(cKerr);
                  xlim([-lim +lim]);
                  ylim([-lim +lim]);
                  hold off

              fH3 = figure(3);
                  clf(fH3)
                  lim = max(abs(cKerr2(:)));
                  hold on
                  % add grid lines
                  plot([-lim, +lim],[0 0],':k')
                  plot([0 0],[-lim, +lim],':k')
                  plot(cKerr2);
                  xlim([-lim +lim]);
                  ylim([-lim +lim]);
                  hold off
             
             if false     
                 [~,shortFName,~] = fileparts(obj.fName);
                 print(fH1,'-dpng',[shortFName,'kerrRot-img',num2str(params.imgInd),'.png']);
                 print(fH2,'-dpng',[shortFName,'kerrRot-img',num2str(params.imgInd),'plane1.png']);
                 print(fH3,'-dpng',[shortFName,'kerrRot-img',num2str(params.imgInd),'plane2.png']);
             end   
        end 
        
        % make movie of Kerr rotation signal
        function makeMovie(obj,varargin)
            % read and parse params
            p = inputParser();
            %p.addParamValue('imgInd',1,@isnumeric);
            p.addParamValue('timeFrames',100,@isnumeric);
            p.addParamValue('fName','movie',@isstr);
            p.parse(varargin{:});
            params = p.Results;
            
            xScale = obj.getXScale();
            yScale = obj.getYScale();
            zScale = [min(real(obj.cKerr(:))) max(real(obj.cKerr(:)))]
            
            freq = 10e9;
            
            videoFile = fullfile(pwd,strcat(params.fName,'.avi'));
            writerObj = VideoWriter(videoFile);
            writerObj.FrameRate = 2;
            open(writerObj);
            
            fig=figure(1);
            for imgInd = 1:16
                ods = obj.params{1,imgInd}.ods;
                time = 2*(obj.params{1,imgInd}.ods-obj.params{1,1}.ods)/(3e11);
                phase = (time*freq);
                imagesc(xScale,yScale,squeeze(real(obj.cKerr(imgInd,:,:))),zScale);
                colormap(jet); colorbar();
                axis xy equal
                xlim([xScale(1) xScale(end)]);
                ylim([yScale(1) yScale(end)]);
                xlabel('x (\mum)','FontName','Times','FontSize',14);
                ylabel('y (\mum)','FontName','Times','FontSize',14);
                title(['L = ',num2str(ods),' mm. Time is ',num2str(time),...
                    '. Phase is ',num2str(phase)]);
                writeVideo(writerObj,getframe(fig));
            end
            
            close(writerObj);
        end    
                
        function res = getXScale(obj)
            ind = 1;
            res = linspace(obj.params{1,ind}.xMin,obj.params{1,ind}.xMax,...
                obj.params{1,ind}.xSteps+1);
        end 
        
        function res = getYScale(obj)
            ind = 1;
            res = linspace(obj.params{1,ind}.yMin,obj.params{1,ind}.yMax,...
                obj.params{1,ind}.ySteps+1);

        end 
        
        function plotKerrImages(obj,varargin)
        % Plot one selected channel for all dataframe with diffrent 
        % position of ODS
        
            % parse input parameters
            p = inputParser();
            p.addParamValue('channel','realKerr',@(x) any(strcmp(x,obj.channels)));
            p.addParamValue('saveAs','',@isstr)
            
            p.parse(varargin{:});
            params = p.Results;
            
            % get data
            switch params.channel
                case 'realKerr'
                    data = real(obj.cKerr);
                otherwise
                    error('Kerr_img:inputParameters','Unknown physical channel')
            end        
            
            zLim = [min(data(:)),max(data(:))];
            
            len = size(obj.Monitor1,1);
            xScale = obj.getXScale();
            yScale = obj.getYScale();
            
            figure(1)
            hS = [];
            for imgInd = 1:len
                hS = [hS,subplot(1,len,imgInd)];
                imagesc(xScale,yScale,...
                    squeeze(data(imgInd,:,:)),zLim);
                axis xy equal
                ylim([yScale(1),yScale(end)]); xlim([xScale(1),xScale(end)])
                title([num2str(obj.params{1,imgInd}.ods),' mm'])
            end    
             
            set(get(hS(1),'YLabel'),'String',' y (\mum)')
            
            if ~strcmp(params.saveAs,'')
                print(gcf,'-dpng','-r300',[params.saveAs,'.png']);
                savefig(gcf,[params.saveAs,'.fig'])
            end    
             
        end
        
        function save(obj,fName)
            %save(fName, obj)
        end
        
        
    end

     %methods (Access = protected)

    
end

