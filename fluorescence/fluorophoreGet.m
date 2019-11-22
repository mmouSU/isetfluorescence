function val = fluorophoreGet(fl,param,varargin)
%
% Syntax
%   val = fluorophoreGet(fl,param,...)
% 
% Description
%   The getter method of the fluorophore object. This function is used to
%   extract different properties of the fluorophore object.
%
% Inputs:
%   fl - a fluorophore structure
%   param - a string describing the parameter of the fluorophore to be
%      returned. Param can have the following values
%
%      'name'                     - fluorophore name
%      'solvent'                  - fluorophore solvent
%      'comment'                  - comment string
%      'type'                     - always 'fluorophore'
%
%      'emission'                 - fluorophore's emission spectrum
%      'normalized emission'      - fluorophore's emission normalized to
%                                   unit amplitude
%      'peak emission'            - peak emission wavelength in nm
%
%      'excitation'               - fluorophore's excitation spectrum
%      'peak excitation'          - peak excitation wavelength in nm
%      'normalized excitation'    - fluorophore's excitation spectrum
%                                   normalized to unit amplitude
%
%      'eem'                      - Excitation-Emission matrix; also
%                                   called the Donaldson matrix
%
%      'Stokes shift'             - wavelength shift between excitation and
%                                   emission peaks
%      'qe'                       - fluorophore's quantum efficiency
%
%      'wave'                     - spectral sampling vector
%      'delta wave'               - interval between concescutive spectral
%                                   samples in nm
%      'nwave'                    - number of spectral samples
%
%      'photons'                  - photons emitted by the fluorophore
%                                   under a particular illuminant. The 
%                                   illuminant is an ISET structure and
%                                   needs to be passed as a parameter (see
%                                   example).
%
% Outputs:
%    val - the value of the requested property.
%
% Examples:
%   em = fluorophoreGet(fl,'emission');
%   name = fluorophoreGet(fl,'name');
% 
%   ill = illuminantCreate('d65');
%   ph = fluorophoreGet(fl,'photons',ill);
%
% Copyright Henryk Blasinski, 2016
%
% See also  fluorophoreCreate, fluorophoreSet, fluorophoreRead
%

%% Parameter checking
if ~exist('fl','var') || isempty(fl), error('Fluorophore structure required'); end
if ~exist('param','var') || isempty(param), error('param required'); end

val = [];

%% Main switch statement
param = ieParamFormat(param);
% param = lower(param);
% param = strrep(param,' ','');

switch param
    case 'name'
        val = fl.name;

    case 'type'
        % Should always be 'fluorophore'
        val = fl.type;
       
    case {'emission','emissionphotons'}
        % Data are stored in photon units
        if ~checkfields(fl,'emission'), val = []; return; end
        val = fl.emission(:);

    case {'emissionenergy'}
        val = fluorophoreGet(fl,'emission');  % Photons
        wave = fluorophoreGet(fl,'wave');     % nm
        val = Energy2Quanta(wave,val);
        
    case {'normemission','normalizedemission'}
        if ~checkfields(fl,'emission'), val = []; return; end
        val = fl.emission(:)/max(fl.emission);
        
    case {'excitation','excitationphotons'}
        % Data are stored in photon units (not energy)
        if ~checkfields(fl,'excitation'), val = []; return; end
        val = fl.excitation(:);
        
    case {'excitationenergy'}
        val = fluorophoreGet(fl,'excitation');  % Photons
        wave = fluorophoreGet(fl,'wave');     % nm
        val = Energy2Quanta(wave,val);
        
    case {'normexcitation','normalizedexcitation'}
        if ~checkfields(fl,'excitation'), val = []; return; end
        val = fl.excitation(:)/max(fl.excitation);
        
    case {'peakexcitation'}
        if ~checkfields(fl,'excitation'), val = []; return; end
        [~, id] = max(fl.excitation);
        val = fl.spectrum.wave(id);    
        
    case {'peakemission'}
        if ~checkfields(fl,'emission'), val = []; return; end
        [~, id] = max(fl.emission);
        val = fl.spectrum.wave(id);
        
    case {'stokesshift'}
        % Difference between the peak emission and peak excitation
        val = fluorophoreGet(fl,'peakemission') - fluorophoreGet(fl,'peakexcitation');
        
    case 'wave'
        if isfield(fl,'spectrum'), val = fl.spectrum.wave; end
        if isvector(val), val = val(:); end
        
    case {'deltawave'}
        wave = fluorophoreGet(fl,'wave');
        val = wave(2) - wave(1);
     
    case {'eem','excitationemissionmatrix','donaldsonmatrix',}
        % This is the excitation emission matrix. It is structured so that
        %
        %     fluorescenceSpectrum = dMatrix * illuminantPhotons(:)
        %
        % 
        deltaL = fluorophoreGet(fl,'delta wave');

        if isfield(fl,'donaldsonMatrix')
            % If the fluorophore is defined in terms of the Donaldson matrix,
            % then return the matrix.
            val = fl.donaldsonMatrix*deltaL;
        elseif isfield(fl,'eem')
            % We want this name instead.  Converting for now, will
            % eliminate donaldson path over time.
            val = fl.eem*deltaL;
        else
            % Otherwise compute the Donaldson matrix from excitation and
            % emission vectors.
            ex = fluorophoreGet(fl,'excitation photons');  % Column
            em = fluorophoreGet(fl,'emission photons');    % Column
            
            % Logic of the EEM is that an excitation photon at some
            % wavelength produces an emission spectrum.  We make the
            % emission vector sum to one (every photon goes somewhere, but
            % no amplificiation). 
            em = em/sum(em(:));
            
            % We scale the excitation vector so that it has a maximum of
            % one, though we could scale it so that it sums to one.  Either
            % way is probably OK because, well, there are no units here.
            ex = ex/max(ex(:));
            
            % We take the outer product so that every column of the
            % Donaldson matrix is the emission, scaled by the relative
            % excitability
            eem = em(:)*ex(:)';
            
            % Finally, we apply the Stoke's constraint so that only
            % photons with  energy lower than the energy of the
            % excitation wavelength (longer wavelengths than the
            % excitation wavelength) are emitted.
            %
            % We leave out the main diagonal (-1 argument) which would
            % normally contain the reflectance function
            %
            % deltaL is the wavelength spacing (delta lambda)
            val = tril(eem,-1) * deltaL;
            
            % Sometime the excitation and emission have NaNs.
            % Set NaNs to 0
            val(isnan(val)) = 0;
            
            % Historical note:  Delete some day.
            %
            % qe = fluorophoreGet(fl,'qe');
            % The fluorophore only sees a fraction of the incident
            % photons from the light.  We call that the quantum
            % efficiency of the fluorophore (a scalar).
            %
            % We have disabled this because it depends a great deal on the
            % solvent and general imaging conditions.
            
        end
        
        % No NaNs on the return.  Make the NaNs 0
        val(isnan(val)) = 0;
        
    case {'photons'}
        illWave  = illuminantGet(varargin{1},'wave');
        illSpd = illuminantGet(varargin{1},'photons');
        
        fl = fluorophoreSet(fl,'wave',illWave);
        DM = fluorophoreGet(fl,'Donaldson matrix');

        val = DM*illSpd;
        
    case 'nwave'
        val = length(fluorophoreGet(fl,'wave'));
        
    case 'comment'
        val = fl.comment;
    
    case 'qe'
        if isfield(fl,'qe')
            val = fl.qe;
        end
        
    case 'solvent'
        val = fl.solvent;
        
    otherwise
        error('Unknown fluorescence parameter %s\n',param)
end

end
