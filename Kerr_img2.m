classdef Kerr_img2 < hgsetget
    properties
        fName = '' % name of file
        Monitor1 = [];
        Monitor2 = [];
        Kerr1 = [];
        Kerr2 = [];
        params
    end
    
    methods
        
        % constructor
        function obj = Kerr_img2
            disp('KERR_img2 object was created');
        end
        
        % open file
        function open(obj)
           % if strcmp(obj.fName,'')
                [fName,fPath,~] = uigetfile({'*.h5','*.*'});
                if (fName == 0)
                    suc = false;
                    return
                end
                suc = true;
                fullName = fullfile(fPath,fName);
                obj.fName = fullName;
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
                    params.xMin = h5readatt(obj.fName,groupName,'Initial x');
                    params.xMax = h5readatt(obj.fName,groupName,'Final x');
                    params.xSteps = h5readatt(obj.fName,groupName,'x steps');
                    params.yMin = h5readatt(obj.fName,groupName,'Initial y');
                    params.yMax = h5readatt(obj.fName,groupName,'Final y');
                    params.ySteps = h5readatt(obj.fName,groupName,'y steps');
                    params.xScale = linspace(params.xMin,params.xMax,params.xSteps+1);
                    params.yScale = linspace(params.yMin,params.yMax,params.ySteps+1);
                    params.delayLinePos = h5readatt(obj.fName,groupName,'ODS (mm)');
                    
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
                        end
                    end
                end
            end
        end
        
        
        
        function plotAll(obj)
            
            % normalize images
            if ~any(size(obj.Monitor1)) && ~any(size(obj.Monitor2))
                monLim = [min(min(m1(:)),min(m2(:))) max(max(m1(:)),max(m2(:)))];
                monLim = [0 0.4]
            end
            if ~any(size(obj.Kerr1)) && ~any(size(obj.Kerr2))
                kerrLim = [min(min(k1(:)),min(k2(:))) max(max(k1(:)),max(k2(:)))];
                kerrLim = [-6e-5 6e-5];
            end
            
            
            
            for imgInd = 1:size(obj.Monitor1,1)
              %  if ~any(size(obj.Kerr1)) && ~any(size(obj.Kerr2))
                    clf();
                    figure(1)
                    subplot(221);
                    imagesc(xScale,yScale,squeeze(obj.Monitor2(imgInd,:,:)),monLim);
                    axis xy equal;
                    xlabel('x (\mum)','FontSize',14,'FontName','Times');
                    ylabel('y (\mum)','FontSize',14,'FontName','Times');
                    xlim([min(xScale) max(xScale)]);
                    title('Reflectivity 1','FontSize',14,'FontName','Times');
                    ylim([min(yScale) max(yScale)]);
                    t = colorbar('peer',gca);
                    set(get(t,'ylabel'),'FontSize',12,'FontName','Times','String', 'Voltage');
                    
                    
                    subplot(222);
                    imagesc(xScale,yScale,squeeze(m2(imgInd,:,:)),monLim);
                    axis xy equal;
                    xlabel('x (\mum)','FontSize',14,'FontName','Times');
                    ylabel('y (\mum)','FontSize',14,'FontName','Times');
                    xlim([min(xScale) max(xScale)]);
                    ylim([min(yScale) max(yScale)]);
                    title('Reflectivity 2','FontSize',14,'FontName','Times');
                    t = colorbar('peer',gca);
                    set(get(t,'ylabel'),'FontSize',12,'FontName','Times','String', 'Voltage');
                    
                    subplot(223);
                    imagesc(xScale,yScale,squeeze(k1(imgInd,:,:)),kerrLim);
                    axis xy equal;
                    xlabel('x (\mum)','FontSize',14,'FontName','Times');
                    ylabel('y (\mum)','FontSize',14,'FontName','Times');
                    xlim([min(xScale) max(xScale)]);
                    ylim([min(yScale) max(yScale)]);
                    title('Kerr 1','FontSize',14,'FontName','Times');
                    t = colorbar('peer',gca);
                    set(get(t,'ylabel'),'FontSize',12,'FontName','Times','String', 'Voltage');
                    
                    subplot(224);
                    imagesc(xScale,yScale,squeeze(k2(imgInd,:,:)),kerrLim);
                    axis xy equal;
                    xlabel('x (\mum)','FontSize',14,'FontName','Times');
                    ylabel('y (\mum)','FontSize',14,'FontName','Times');
                    xlim([min(xScale) max(xScale)]);
                    ylim([min(yScale) max(yScale)]);
                    title('Kerr 2','FontSize',14,'FontName','Times');
                    t = colorbar('peer',gca);
                    set(get(t,'ylabel'),'FontSize',12,'FontName','Times','String', 'Voltage');
                    
                    colormap(jet)
                    %print(gcf,'-dpng',strcat(fName,'-',num2str(imgInd),'-whole.png'));
                    
                    figure(2)
                    subplot(121);
                    imagesc(xScale,yScale,squeeze(m1(imgInd,:,:)),monLim);
                    axis xy equal;
                    %xlabel('x (\mum)','FontSize',14,'FontName','Times');
                    %ylabel('y (\mum)','FontSize',14,'FontName','Times');
                    xlim([min(xScale) max(xScale)]);
                    ylim([min(yScale) max(yScale)]);
                    title(num2str( delayLinePos));
                    %title('Reflectivity 1','FontSize',14,'FontName','Times');
                    %t = colorbar('peer',gca);
                    %set(get(t,'ylabel'),'FontSize',12,'FontName','Times','String', 'Voltage');
                    
                    subplot(122);
                    imagesc(xScale,yScale,squeeze(k1(imgInd,:,:)),kerrLim);
                    axis xy equal;
                    %xlabel('x (\mum)','FontSize',14,'FontName','Times');
                    %ylabel('y (\mum)','FontSize',14,'FontName','Times');
                    xlim([min(xScale) max(xScale)]);
                    ylim([min(yScale) max(yScale)]);
                    %title('Kerr 1','FontSize',14,'FontName','Times');
                    %t = colorbar('peer',gca);
                    %set(get(t,'ylabel'),'FontSize',12,'FontName','Times','String', 'Voltage');
                    %print(gcf,'-dpng',strcat(fName,'-',num2str(imgInd),'-m1k1.png'));
                    
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
        
        %% plot one desired image
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
            
    end
    
    methods (Access = protected)
    end
    
end
