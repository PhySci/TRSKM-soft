classdef Kerr_ODS < hgsetget
    %Kerr_ODS - matlab class for processing of time-resolved (ODS) signals
    %   Detailed explanation goes here
    
    properties
        
        fName = '';
        
        % arrays of signal 
        Chanel1 = 0;
        Chanel2 = 0;
        Monitor1 = 0;
        Monitor2 = 0;
        FFTspec = 0; % Fourier spectrum
        
        % params
        dt = 0;
        freqScale = 0;
        lengthScale = 0;
        timeScale = 0;
        position = [0 0 0];
        waitTime = 0; 
    end
    
    properties (Access = protected)
    end    
    
    methods
        
        % create object
        function obj = Kerr_ODS
            disp('Kerr_ODS object was created');
        end    
        
        % open new file file
        function res = open(obj)
            obj.fName = '';
            if obj.load() 
                obj.makeFFT('chanel','Chanel1');
                obj.plotKerr();
                obj.plotMonitors();
            end
        end
        
        % load file
        function suc = load(obj)
           if isempty(obj.fName)   
              [fName,fPath,~] = uigetfile({'*.h5';'*.txt'});
              if (fName == 0)
                  suc = false;
                  return
              end
              suc = true;
              fullName = fullfile(fPath,fName);
              obj.fName = fullName;
           end
          
          [~,~,ext] = fileparts(fName); 

          if strcmp(ext,'.txt')
              res = obj.readDataFile().';
              %readParamsFile(fullName);
          elseif strcmp(ext,'.h5')
              h5disp(obj.fName);
              obj.position = h5readatt(obj.fName,'/','Position');
              obj.waitTime = h5readatt(obj.fName,'/','time_wait_(ms) ');
              res = h5read(obj.fName,'/Signal').';
          end
          
          if res(1,1)>res(end,1)
              res = flipud(res);
          end
          
          obj.lengthScale = res(:,1);
          obj.Chanel1 = res(:,4);
          obj.Chanel2 = res(:,5);
          obj.Monitor1 = res(:,2);
          obj.Monitor2 = res(:,3);
         
              
          
          obj.calcFreqScale();         
        end    
        
        % calculate FFT spectra of the signal
        function makeFFT(obj,varargin)
            p = inputParser();
            p.addParamValue('chanel','Complex',@(x) any(strcmp(x,{'Complex','Chanel1','Chanel2'})))
            p.parse(varargin{:});
            params = p.Results;
            
            obj.FFTspec = zeros(size(obj.Chanel1,1),1);
            windArr = rectwin(size(obj.Chanel1,1));
            
            switch params.chanel
                case 'Complex'
                    signal = windArr.*(obj.Chanel1+i*obj.Chanel2);
                case 'Chanel1'
                    signal = windArr.*(obj.Chanel1);
                case 'Chanel2'
                    signal = windArr.*(obj.Chanel2);
            end
            
            obj.FFTspec(:,1) = fftshift(abs(fft(signal)));
            num2clip(obj.FFTspec(:,1));
        end
        
        % scan folder, calculate FFT for every file and average FFT spectra
        function FFTspec = scanFolder(obj)
            
            fList = dir(fullfile(pwd,'*.txt'));
            if (size(fList)==0)
                disp('Files havent been found');
                return;
            end
            
            FFTspecArr = [];
            for fInd = 1:size(fList)
                obj.fName = fullfile(pwd,fList(fInd).name);
                obj.load();
                obj.makeFFT();
                obj.plot();
                FFTspecArr(:,fInd) = obj.FFTspec;
            end

            FFTspec = mean(FFTspecArr,2);
            
            h2 = figure(2);
                plot(obj.freqScale,FFTspec);
                xlabel('Frequency (GHz)','FontSize',14,'FontName','Times');
                ylabel('FFT intensity (arb. units)','FontSize',14,'FontName','Times');
                xlim([0 20]);
        end 
        
        % plot signal and FFT spectra
        function plotKerr(obj,varargin)
            
            p = inputParser();
            p.addParamValue('saveAs','',@isstring);
            p.parse(varargin{:});
            params = p.Results; 
            
            [~,fName,~] = fileparts(obj.fName);
            figTitle = [fName, ' Position:(',num2str(obj.position.x_position,'%10.2f'),',',...
                num2str(obj.position.y_position,'% 10.2f'),',',...
                num2str(obj.position.z_position,'% 10.2f'),')'];
            
            hf = figure();
            subplot(211);
                plot(obj.timeScale,obj.Chanel1,'-r',obj.timeScale,obj.Chanel2,'-b');
                xlim([min(obj.timeScale) max(obj.timeScale)]);
                title(figTitle,'FontName','Times');
                xlabel('Delay (ns)','FontSize',14,'FontName','Times');
                ylabel('Signal (V)','FontSize',14,'FontName','Times');
                legend('Chanel 1','Chanel 2');
                
            subplot(212);
                plot(obj.freqScale, obj.FFTspec);
                title('FFT','FontName','Times');
                legend('rectangular window','Hamming window')
                xlim([0 20]);
                xlabel('Frequency (GHz)','FontSize',14,'FontName','Times');
                ylabel('FFT intensity (arb. units)','FontSize',14,'FontName','Times');
                
            [pathstr,fName,ext] = fileparts(obj.fName);
            savefig(hf,strcat(fName,'-kerr.fig'));
            print(hf,'-dpng','-r600',strcat(fName,'-kerr.png'));
        end 
        
        % plot monitors values and their FFT spectra
        function plotMonitors(obj,varargin)
            p = inputParser();
            p.addParamValue('saveAs','',@isstring);
            p.parse(varargin{:});
            params = p.Results; 
            
            [~,fName,~] = fileparts(obj.fName);
            figTitle = [fName, ' Position:(',num2str(obj.position.x_position,'% 10.2f'),',',...
                num2str(obj.position.y_position,'% 10.2f'),',',...
                num2str(obj.position.z_position,'% 10.2f'),')'];
            
            hf = figure(2);
            subplot(211);
                plot(obj.timeScale,obj.Monitor1,'-r',obj.timeScale,obj.Monitor2,'-b');
                xlim([min(obj.timeScale) max(obj.timeScale)]);
                title(figTitle);
                xlabel('Delay (ns)','FontSize',14,'FontName','Times');
                ylabel('Signal (V)','FontSize',14,'FontName','Times');
                legend('Monitor 1','Monitor 2');
            
            M1FFT = fftshift(abs(fft(obj.Monitor1)));
            M2FFT = fftshift(abs(fft(obj.Monitor2)));
            
            subplot(212);
                semilogy(obj.freqScale, M1FFT,obj.freqScale, M2FFT); title('FFT');
                xlim([-15 15]);
                xlabel('Frequency (GHz)','FontSize',14,'FontName','Times');
                ylabel('FFT intensity (arb. units)','FontSize',14,'FontName','Times');
                
            [pathstr,fName,ext] = fileparts(obj.fName);
            savefig(hf,strcat(fName,'-monitor.fig'));
            print(hf,'-dpng','-r600',strcat(fName,'-monitor.png'));
        end    
        
        function sinFit(obj,varargin)
            sinFunc  = @(b,x) (b(1)*sin(b(4)+2*pi*x/b(3))+b(2));
            fitRes = nlinfit(obj.lengthScale,obj.Chanel1,sinFunc,[1, 1,80,0.1]); %,...
            amp1 = fitRes(1);
            bias1 = fitRes(2);
            period1 = fitRes(3);
            shift1 = fitRes(4);
            yFit = amp1*sin(shift1+2*pi*obj.lengthScale/period1)+bias1;
            plot(obj.lengthScale,obj.Chanel1,'xr',obj.lengthScale,yFit,'-g')
        end
        
        function  fft = plotAverageSpec(obj,varargin)
            
            p = inputParser();
            p.addParamValue('saveAs','',@isstr);
            p.addParamValue('specNum',2,@isnumeric);
            p.addParamValue('chanel','Complex',@(x) any(strcmp(x,{'Complex','Chanel1','Chanel2'})))
            p.parse(varargin{:});
            params = p.Results; 
            
            % read files
            for specInd = 1:params.specNum
                obj.fName = '';
                obj.load();
                obj.makeFFT('chanel',params.chanel);
                signalArr(:,specInd) = obj.Chanel1;
                fftArr(:,specInd) = obj. FFTspec(:,1);
                
                %legStr{specInd} = ['(' num2str(obj.position.x_position,'%5.2f') ','...
                %    num2str(obj.position.y_position,'%5.2f') ','...
                %    num2str(obj.position.z_position,'%5.2f') ')'];
            end
            legStr = '';
            % calculate average value
            fft = mean(fftArr,2);
            
            subplot(211);
                plot(obj.timeScale,signalArr);
                xlabel('Delay (ns)','FontSize',14,'FontName','Times');
                ylabel('Kerr rotation (arb. units)','FontSize',14,'FontName','Times');
                legend(legStr,'Location','Best')
            subplot(212);
                plot(obj.freqScale, fft); title('FFT');
                xlim([0 20]);
                xlabel('Frequency (GHz)','FontSize',14,'FontName','Times');
                ylabel('FFT intensity (arb. units)','FontSize',14,'FontName','Times');
                
            
            % save img
            if ~isempty(params.saveAs)
                savefig(gcf,strcat(params.saveAs,'.fig'));
                print(gcf,'-dpng','-r600',strcat(params.saveAs,'.png'));
            end
            
            num2clip(fft);
            
        end
        
        % plot Kerr chanels as complex value
        function plotComplex(obj)
            figure(1);
                plot(obj.Chanel1,obj.Chanel2);
                xlabel('Chanel1'); ylabel('Chanel2');
                axis equal
            figure(2);
                subplot(211);
                title('Phase of signal')
                plot(angle(obj.Chanel1+i*obj.Chanel2));
                subplot(212);
                title('Phase of signal')
                plot(fftshift(abs(fft(angle(obj.Chanel1+i*obj.Chanel2)))));
                
            figure(3);
                 
                subplot(311);
                title('Phase of monitors')
                plot(angle(obj.Monitor1+i*obj.Monitor2));
                subplot(312);
                semilogy(obj.freqScale,abs(fftshift(fft(angle(obj.Monitor1+i*obj.Monitor2)))));
        end    
        
    end
    
    methods (Access = protected)
        
        % Read file and return array of data 
        function res = readDataFile(obj)
          fid = fopen(obj.fName);
          errorMsg = ferror(fid);
          if (~strcmp(errorMsg,''))
              disp(errorMsg);
              return;
          end    
          res=[];
          while (~feof(fid))
            str = fgetl(fid);
            [tmp, ~, errmsg] = sscanf(str, '%e');
            if (~strcmp(errmsg,''))
              disp('A problem occured with file reading');
              disp(errmsg);
              return;
            end    
            res = cat(2,res,tmp);
          end
          fclose(fid);
        end
        
        
        % Read file of parameters and print all information on the screen
        function readParamsFile(obj,fPath)
          % add "-i" suffix to name of file
          [pathstr,name,ext] = fileparts(fPath);
          name = strcat(name,'-i',ext);
          fPath = fullfile(pathstr,name);

          % open and read file
          fid = fopen(fPath);
          errorMsg = ferror(fid);
          if (~strcmp(errorMsg,''))
              disp(errorMsg);
              return;
          end
          disp(' ');
          disp(name);
          while (~feof(fid))
            disp(fgetl(fid));  
          end
          fclose(fid);
        end
        
        function calcFreqScale(obj)
            obj.timeScale = 2*(obj.lengthScale - min(obj.lengthScale))/3e2; % m
            obj.dt = abs(mean(diff(obj.timeScale)));
            obj.freqScale = linspace(-0.5/obj.dt,0.5/obj.dt,size(obj.timeScale,1)).';
        end
        
    end  
    
end

