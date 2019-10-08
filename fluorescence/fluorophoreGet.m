function val = fluorophoreGet(fl,param,varargin)

% val = fluorophoreGet(fl,param,...)
% 
% The getter method of the fluorophore object. This function is used to
% extract different properties of the fluorophore object.
%
% Examples:
%   em = fluorophoreGet(fl,'emission');
%   name = fluorophoreGet(fl,'name');
% 
%   ill = illuminantCreate('d65');
%   ph = fluorophoreGet(fl,'photons',ill);
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
%      'peak emission'            - peak emission wavelenght in nm
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
% Copyright Henryk Blasinski, 2016

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
        
        if ~checkfields(fl,'emission'), val = []; return; end
        val = fl.emission(:);

    case {'normemission','normalizedemission'}
        if ~checkfields(fl,'emission'), val = []; return; end
        val = fl.emission(:)/max(fl.emission);
        
    case {'excitation','excitationphotons'}
        
        if ~checkfields(fl,'excitation'), val = []; return; end
        val = fl.excitation(:);
        
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
        
    case {'deltawave','deltaWave'}
        wave = fluorophoreGet(fl,'wave');
        val = wave(2) - wave(1);
     
    case {'eem','excitationemissionmatrix','donaldsonmatrix',}
        % This is also the excitation emission matrix.
        % The matrix is expected to be structured so that
        %
        %     fluorescenceSpectrum = dMatrix * illuminantPhotons(:)
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
            % Otherwise compute the Donaldson matrix from excitation and emission vectors. 
            % The matrix is the outer product
            %
            %       excitation amplitude times the single emission spectrum, and then
            % forcing the excitation emission matrix to lower
            % diagonal; only emissions at longer wavelengths than the
            % excitation wavelength.
            ex = fluorophoreGet(fl,'excitation photons');  % Column
            em = fluorophoreGet(fl,'emission photons');    % Column
            % qe = fluorophoreGet(fl,'qe');
            
            % For the separable case, a photon at some wavelength has
            % an efficacy of producing the emission spectrum.  That
            % efficacy is encoded by the excitation vector.
            %
            % The emission spectrum is encoded by the emission vector.
            %
            % These are combined so that the incident light at each
            % wavelength will be multiplied by the excitation efficacy
            % at each wavelength, and this scales the amount of the
            % emission spectrum produced.  The same calculation takes
            % place at every wavelength.
            %
            % The fluorophore only sees a fraction of the incident
            % photons from the light.  We call that the quantum
            % efficiency of the fluorophore (a scalar).
            %
            % Finally, we apply the Stoke's constraint so that only
            % photons with  energy lower than the energy of the
            % excitation wavelength (longer wavelengths than the
            % excitation wavelength) are emitted.
            %
            % deltaL is the wavelength spacing (delta lambda)
            
            % Each column is the emission spectrum.  We leave out the main
            % diagonal (-1 argument), which would normally contain the
            % reflectance function. We scale by delta lambda so that
            % different wavelength sampling schemes will return a similar
            % value.  But the reality is we have very little knowledge
            % about the absolute levels in here.
            val = tril(em*ex',-1) * deltaL;
            
            % Set NaNs to 0
            val(isnan(val)) = 0;
        end
        
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
