classdef RsCorrection
    
	% performs offline series resistance correction for recorded currents
	% Input:
	%
	% - data = recorded current trace, can be single vector or multiple recordings as array or cell
	% - Rs = uncompensated series resistance (Rs) during the recording in Ohm, i.e. if the Rs during the experiment was 10 MOhm and online compensated by 50% by
	% the amplifier, the remaining uncompensated Rs will be 5 MOhm (5e6 Ohm = 5 MOhm)
	% - Cm = Membrane capacitance during the recording in Farad (e.g. 10e-12 F = 10 pF)
	% - Vhold = holding potential during the recording in Volts (e.g. -0.06 V =  -60 mV)
	% - Vrev = reversal potential of the recorded current in Volts (e.g. 0.01V = 10 mV)
	% - SR = Sampling rate during the recordings (in Hz)
	% - [optional] fractionV = voltage error to be compensated [0-1] (e.g. 1 if voltage error should be fully compensated)
	% - [optional] fractionC = fraction capacitative filtering error to be compensated [0-1] (e.g. 1 capacitative filtering error should be fully compensated)
	% - [optional] fc = cutoff frequency for filter to smooth capacitative current correction (in Hz) (if omitted, fc is calculated from the sampling interval as fc = 1/(2 * pi * si))

	% Based on: "Traynelis SF (1998) Software-based correction of single
	% compartment series resistance errors. J Neurosci Methods 86:25–34."
	%
	% EXAMPLE: RsCorrection(data, Rs, Cm, Vhold, Vrev, SR, 'fractionV', 1, 'fractionC', 1, 'fc', 10e3)
	
     
    
    properties
        dataRaw
        dataCorrected
        options
        
    end
    
    
    methods
        function obj = RsCorrection(data, Rs, Cm, Vhold, Vrev, SR, varargin)
            
            % CHECK INPUTS
            checkData = @(n)validateattributes(n, {'numeric','DimensionedVariable','cell'},{'nonempty'});
            checkNumericPos = @(n)validateattributes(n, {'numeric','DimensionedVariable'},{'nonnegative','nonnan','nonempty'});
            checkVoltage  = @(n)validateattributes(n, {'numeric','DimensionedVariable'},{'nonnan','nonempty'});
            
            P = inputParser;
            % REQUIRED INPUTS:
            P.addRequired('data',checkData)
            P.addRequired('Rs',checkNumericPos)
            P.addRequired('Cm',checkNumericPos)       
            P.addRequired('Vhold',checkVoltage)
            P.addRequired('Vrev',checkVoltage)
            P.addRequired('SR',checkNumericPos)
            
            % OPTIONAL INPUT
            P.addOptional('fractionC',1,checkNumericPos)
            P.addOptional('fractionV',1,checkNumericPos)
            P.addOptional('fc',NaN,checkNumericPos)

            P.parse(data, Rs, Cm, Vhold, Vrev, SR, varargin{:});
            opt = P.Results;
            obj.options = opt;
            obj.dataRaw = data;
%%            
            % INTERNALLY ALL DATA ARE TREATED AS CELLS                    
            dataIsCell = iscell(data);
                 
            if ~dataIsCell
                [~,dataDim] = max(size(data));
                data = num2cell(data,dataDim);
            end
            
            si = 1/SR;
            
            if isnan(opt.fc)
            tauLag = si;
            opt.fc = 1/(2*pi*tauLag);
            end
            
            filterfactor = (1-exp(-2*pi*si*opt.fc));
            
            for iTrace = 1:numel(data)
                nPoints = length(data{iTrace});
                
                Vlast = Vhold-data{iTrace}(1)*Rs; % INITIALIZE DC VOLTAGE AT MEMBRANE
                denominator = Vlast-Vrev;
                
                if denominator ~= 0
                    fracCorr = opt.fractionV*(1-(Vhold-Vrev)/denominator);
                else
                    fracCorr = 0;
                end
                data{iTrace}(1) = data{iTrace}(1)*(1-fracCorr); % DC CORRECTION FOR FIRST DATA POINT
                
                
                for i=2:nPoints-1 % LOOP THROUGH ALL OTHER DATA POINTS
                    Vthis = Vhold - data{iTrace}(i)*Rs;
                    if Vthis ~= Vrev
                        fracCorr = opt.fractionV*(1-(Vhold-Vrev)/(Vthis-Vrev));
                    else
                        fracCorr = 0;
                    end
                    Icap = Cm*(Vthis-Vlast)/si;
                    Icap = Icap*filterfactor;
                    data{iTrace}(i-1) = data{iTrace}(i-1)-(opt.fractionC*Icap);
                    data{iTrace}(i-1) = data{iTrace}(i-1)*(1-fracCorr);
                    Vlast = Vthis;
                end
                
                
            end
            
            if ~dataIsCell
                obj.dataCorrected = cell2mat(data);
            else 
                obj.dataCorrected = data;
            end

            
        end
    end
    
end


