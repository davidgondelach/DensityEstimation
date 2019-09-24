function [GM, Re, maxdeg, C_gravmod, S_gravmod, sF_gravmod]= initgravitysphericalharmonic( varargin )
%  GRAVITYSPHERICALHARMONIC Implement a spherical harmonic representation
%   of planetary gravity. 
%   [GX GY GZ] = GRAVITYSPHERICALHARMONIC( P ) implements the mathematical
%   representation of spherical harmonic planetary gravity based on
%   planetary gravitational potential. Using P, a M-by-3 array of
%   Planet-Centered Planet-Fixed coordinates, GX, GY and GZ, arrays of M 
%   gravity values in the x-axis, y-axis and z-axis of the Planet-Centered
%   Planet-Fixed coordinates are calculated for planet using 120th degree 
%   and order spherical coefficients for EGM2008 by default. 
%
%   Alternate formats for calling spherical harmonic gravity are:
%   [GX GY GZ] = GRAVITYSPHERICALHARMONIC( P, DEGREE )   
%   [GX GY GZ] = GRAVITYSPHERICALHARMONIC( P, MODEL )   
%   [GX GY GZ] = GRAVITYSPHERICALHARMONIC( P, MODEL, DEGREE )   
%   [GX GY GZ] = GRAVITYSPHERICALHARMONIC( P, MODEL, DEGREE, ACTION )   
%   [GX GY GZ] = GRAVITYSPHERICALHARMONIC( P, 'Custom', DEGREE, {DATAFILE DFREADER}, ACTION )   
%
%   Inputs for spherical harmonic gravity are:
%   P        :a M-by-3 array of Planet-Centered Planet-Fixed coordinates in
%            meters where the z-axis is positive towards the North Pole. For
%            Earth this would be ECEF coordinates.
%   MODEL    :a string specifying the planetary model:
%            'EGM2008' (Earth), 'EGM96' (Earth), 'JGM-3' (Earth),'LP100K' 
%            (Moon), 'LP165P' (Moon), 'GMM2B' (Mars), 'Custom', or 
%            'EIGENGL04C' (EARTH). The default is 'EGM2008'. 
%   DEGREE   :a scalar value specifying the degree and order of the      
%            harmonic gravity model. For 'EGM2008', the maximum degree and
%            order is 2159 and the default degree and order is 120.   For
%            'EGM96', the maximum degree and order is 360 and the default
%            degree and order is 70.  For 'JGM3' the maximum degree and 
%            order is 70 and the default degree and order is 30. For 
%            'LP100K', the maximum degree and order is 100 and the default 
%            degree and order is 60.  For 'LP165P', the maximum degree and 
%            order is 165 and the default degree and order is 60.  For 
%            'GMM2B', the maximum degree and order is 80 and the default 
%            degree and order is 60.  For 'Custom', the default degree and 
%            order is the maximum degree. For 'EIGENGL04C', the maximum 
%            degree and order is 360 and the default degree and order is 70.
%   DATAFILE :a file containing the planetary gravitational parameter,
%            planet equatorial radius, maximum degree, and normalized 
%            spherical harmonic coefficient matrices.
%   DFREADER :a function handle to an MATLAB(R) function which reads
%            DATAFILE.  The MATLAB function must output planetary
%            gravitational parameter in meters cubed per second squared,
%            planet equatorial radius in meters, maximum degree, and the
%            normalized spherical harmonic coefficient matrices, C and S.
%   ACTION   :a string to determine action for out of range input. Specify
%            if out of range input invokes a 'Warning', 'Error', or no
%            action ('None'). The default is 'Warning'.
%
%   Output calculated for the spherical harmonic gravity includes:
%   GX     :an array of M gravity values in the x-axis of the
%          Planet-Centered Planet-Fixed coordinates in meters per second
%          squared.
%   GY     :an array of M gravity values in the y-axis of the
%          Planet-Centered Planet-Fixed coordinates in meters per second
%          squared. 
%   GZ     :an array of M gravity values in the z-axis of the
%          Planet-Centered Planet-Fixed coordinates in meters per second
%          squared. 
%
%   Limitations:                                                           
%
%   This function has the limitations of excluding the centrifugal effects
%   of planetary rotation, and the effects of a precessing reference frame.
%
%   Spherical harmonic gravity model is valid for radial positions greater
%   than the planet equatorial radius.  Using it near or at the planetary
%   surface can probably be done with negligible error.  The spherical
%   harmonic gravity model is not valid for radial positions less than
%   planetary surface. 
%
%   Examples:                                                              
%
%   Calculate the gravity in the x-axis at the equator on the surface of
%   Earth, using the 120 degree model of EGM2008 with warning actions:
%       gx = gravitysphericalharmonic( [-6378.1363e3 0 0] ) 
%
%   Calculate the gravity at 25000 meters over the south pole of Earth using
%   the 70 degree model of EGM96 with error actions: 
%       [gx, gy, gz] = gravitysphericalharmonic( [0 0 -6381.751e3], 'EGM96', 'Error' )   
%
%   Calculate the gravity at 15000 meters over the equator and 11000 meters
%   over the north pole using a 30th order GMM2B Mars model with warning
%   actions:
%       p  = [2412.648e3 -2412.648e3 0; 0 0 3376.2e3]
%       [gx, gy, gz] = gravitysphericalharmonic( p, 'GMM2B', 30, 'Warning' )   
%
%   Calculate the gravity at 15000 meters over the equator and 11000 meters
%   over the north pole using a 60th degree custom planetary model with no
%   actions:  
%       p       = [2412.648e3 -2412.648e3 0; 0 0 3376e3]
%       [gx, gy, gz] = gravitysphericalharmonic( p, 'custom', 60, ...
%                       {'GMM2BC80_SHA.txt' @astReadSHAFile}, 'None' )
%
%   See also GRAVITYWGS84, GRAVITYCENTRIFUGAL, GRAVITYZONAL, GEOIDEGM96

%   Copyright 2009-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.9.2.1 $  $Date: 2011/01/13 20:01:11 $

%   References:  
%   [1] Vallado, D. A., "Fundamentals of Astrodynamics and Applications",
%       McGraw-Hill, New York, 1997.  
%   [2] NIMA TR8350.2: "Department of Defense World Geodetic System 1984,
%       Its Definition and Relationship with Local Geodetic Systems." 
%   [3] Konopliv, A. S., S. W. Asmar, E. Carranza, W. L. Sjogen, D. N.
%       Yuan., "Recent Gravity Models as a Result of the Lunar Prospector
%       Mission", Icarus, Vol. 150, no. 1, pp 1?18, 2001.                    
%   [4] Lemoine, F. G., D. E. Smith, D.D. Rowlands, M.T. Zuber, G. A.
%       Neumann, and D. S. Chinn, "An improved solution of the gravity
%       field of Mars (GMM-2B) from Mars Global Surveyor", J. Geophys. Res.,
%       Vol. 106, No. E10, pp 23359-23376, October 25, 2001.   
%   [5] Kenyon S., J. Factor, N. Pavlis, and S. Holmes, "Towards the Next
%       Earth Gravitational Model", Society of Exploration Geophysicists
%       77th Annual Meeting, San Antonio, Texas, September 23-28, 2007.
%   [6] Pavlis, N.K., S.A. Holmes, S.C. Kenyon, and J.K. Factor, "An Earth
%       Gravitational Model to Degree 2160: EGM2008", presented at the 2008
%       General Assembly of the European Geosciences Union, Vienna,
%       Austria, April 13-18, 2008. 
%   [7] Grueber, T., and A. Kohl, "Validation of the EGM2008 Gravity Field
%       with GPS-Leveling and Oceanographic Analyses", presented at the IAG
%       International Symposium on Gravity, Geoid & Earth Observation 2008,
%       Chania, Greece, June 23-27, 2008.
%   [8] F?rste, C., Flechtner, F., Schmidt, R., K?nig, R., Meyer, U.,
%       Stubenvoll, R., Rothacher, M., Barthelmes, F., Neumayer, H.,
%       Biancale, R., Bruinsma, S., Lemoine, J.M., Loyer, S., "A Mean
%       Global Gravity Field Model From the Combination of Satellite
%       Mission and Altimetry/Gravmetry Surface Data - EIGEN-GL04C",
%       Geophysical Research Abstracts, Vol. 8, 03462, 2006 
%       http://icgem.gfz-potsdam.de/ICGEM/

error(nargchk(1, 5, nargin,'struct'));

% set default values
model  = 'EGM2008';
action = 'warning';


switch nargin                                                              
    case 2
        % set degree or model
        if ~(isnumeric( varargin{1} ) && isreal( varargin{1} ))
            if ~ischar( varargin{1} )
                error(message('aero:gravitysphericalharmonic:inputTypeVar'));
            else
                if strcmpi( varargin{1}, 'custom')
                    error(nargchk(4, 5, nargin,'struct'));
                else
                    % assign model or action
                    modeloraction( varargin{1} )
                    maxdeg = varargin{2};
                end
            end
        else
            maxdeg = varargin{1};
        end
    case 3
        if ~ischar( varargin{2} )
            % not setting action
            if ~ischar( varargin{1} ) || ...
              ~(isnumeric( varargin{2} ) && isreal( varargin{2} ))
                error('aero:gravitysphericalharmonic:inputTypeNoAction3', ...
                     ['Model input must be a string and maximum degree '...
                     'input must be numeric and real.'] );
            end
            if strcmpi( varargin{1}, 'custom')
                error(nargchk(4, 5, nargin,'struct'));
            end
            % set model and degree
            checkmodel( varargin{1} );
            maxdeg = varargin{2};    
        else
            if ~(isnumeric( varargin{1} ) && isreal( varargin{1} ))
                if ~ischar( varargin{1} )
                    error('aero:gravitysphericalharmonic:inputTypeVar3', ...
                         ['Model input must be a string or maximum degree',...
                         ' input must be numeric and real.'] );
                end
                if strcmpi( varargin{1}, 'custom')
                    error(nargchk(4, 5, nargin,'struct'));
                end
                % set model and action
                checkmodel( varargin{1} );
                checkaction( varargin{2} );              
            else
            % set degree and action
            maxdeg = varargin{1};
            checkaction( varargin{2} );
            end
        end
    case 4
        if ~ischar( varargin{1} ) || ...
                ~(isnumeric( varargin{2} ) && isreal( varargin{2} ))             
            error('aero:gravitysphericalharmonic:inputTypeVar4',...
                ['Model input must be a string' ...
                ' and maximum degree input must be numeric and real.'] );
        else
            if ischar( varargin{3} )
                % set model, degree and action
                checkmodel( varargin{1} );
                maxdeg = varargin{2};
                checkaction( varargin{3} );
            elseif iscell( varargin{3} ) && ischar(varargin{3}{1}) && ...
                    isa(varargin{3}{2}, 'function_handle') && ...
                    strcmpi( varargin{1}, 'custom')
                % set model, degree, data file and data file reader
                checkmodel( varargin{1} );
                maxdeg = varargin{2};
                datafile = varargin{3}{1};                                 
                dfreader = varargin{3}{2};
            else
                if strcmpi( varargin{1}, 'custom')
                    error('aero:gravitysphericalharmonic:inputTypeVar4Data', ...
                        ['Data file string and reader function handle inputs',... 
                         ' must be contained in a cell array.']);
                else
                    error(message('aero:gravitysphericalharmonic:inputTypeVar4Action'));
                end
            end
        end
    case 5
        if ischar( varargin{1} ) && ~strcmpi( varargin{1}, 'custom')
            % This option is only for 'custom'
            error(message('aero:gravitysphericalharmonic:wrongModel'));
        end
        if ~ischar( varargin{1} ) || ...
           ~(isnumeric( varargin{2} ) && isreal( varargin{2} )) || ...
           ~(iscell( varargin{3} ) && ischar(varargin{3}{1}) && ...
             isa(varargin{3}{2}, 'function_handle')) || ...
           ~(ischar( varargin{4} ))
            error('aero:gravitysphericalharmonic:inputTypeVar5',['Model and action ', ...
                  'inputs must be strings. Degree input must be numeric',...
                  ' and real. Data file string and reader function handle', ...
                  ' inputs must be contained in a cell array.'] );
        end
        % set model, degree, data file, data file reader, and action
        checkmodel( varargin{1} );
        maxdeg   = varargin{2};
        datafile = varargin{3}{1};                                         
        dfreader = varargin{3}{2}; 
        checkaction( varargin{4} );     
end
           
switch lower( model )                                                      
    case 'custom'
        if ~exist( datafile, 'file' )
            error(message('aero:gravitysphericalharmonic:noDataFile', datafile))
        end
        try
            [GM, Re, degree, C, S] = dfreader( datafile );
        catch MECustomFile
            try
                % try loading a mat-file
                dfreader( datafile );  % [GM, Re, degree, C, S] 
            catch MECustomFileLoad
                throwAsCaller(MECustomFile)
            end                                                            
            % check for the existence of correct variables in mat-file
            if ~exist( 'GM', 'var' ) 
                error('aero:gravitysphericalharmonic:noGM', ['Value for',...
                      ' planetary gravitational parameter, GM, must be ',...
                      'contained in mat-file, %s.'], datafile)
            end
            if ~exist( 'Re', 'var' ) 
                 error('aero:gravitysphericalharmonic:noRe', ['Value for',...
                       ' planet equatorial radius, Re, must be contained ',...
                       'in mat-file, %s.'], datafile)
            end
            if ~exist( 'degree', 'var' ) 
                 error('aero:gravitysphericalharmonic:noDegree', ['Value for',...
                       ' maximum degree, degree, must be contained ',...
                       'in mat-file, %s.'], datafile)
            end
            if ~exist( 'C', 'var' ) 
                 error('aero:gravitysphericalharmonic:noC', ['Value for',...
                       ' normalized spherical coefficient C, must be ',...
                       'contained in mat-file, %s.'], datafile)
            end
            if ~exist( 'S', 'var' ) 
                 error('aero:gravitysphericalharmonic:noS', ['Value for',...
                       ' normalized spherical coefficient S, must be ',...
                       'contained in mat-file, %s.'], datafile)
            end
        end
        
        % check data type and sizes of custom data
        if ~isscalar( GM ) || ~isscalar( Re ) || ~isscalar( degree ) || ~isnumeric( GM ) || ~isnumeric( Re ) || ~isnumeric( degree ) 
            error('aero:gravitysphericalharmonic:notScalar',['Values for ',...
                'planetary gravitational parameter, GM, planet equatorial ',...
                'radius, Re, and maximum degree, degree must be scalar numeric values.'])
        end
        if  (~isnumeric(C) || ~isnumeric(S)) || (~all( size(C) == [ degree+1 degree+1 ] ) || ~all( size(S) == [ degree+1 degree+1 ] ))
            error('aero:gravitysphericalharmonic:wrongMatrix',['Values ',...
                'for normalized spherical coefficients C and S must be ',...
                'numeric matrices of size (maximum degree+1)-by-(maximum degree+1).'])
        end
        
        % The recommended default degree for custom is unknown, so use the
        % maximum degree of custom unless maxdeg is input
        checkmaxdeg( degree, degree )                                                   
    case 'egm2008'                                                         
        % Earth
        % Load normalized coefficients and planetary constants
        %load('aeroegm2008.mat') % [GM, Re, degree, C, S] 
        load('EGM2008.mat') % [GM, Re, degree, C, S] 
        checkmaxdeg( degree, 120 ) %#ok<NODEF>
    case 'egm96'
        % Earth
        % Load normalized coefficients and planetary constants
        load('aeroegm96.mat') % [GM, Re, degree, C, S] 
        checkmaxdeg( degree, 70 ) %#ok<NODEF>
    case 'jgm3'                                                         
        % Earth
        % Load normalized coefficients and planetary constants
        load('JGM3.mat') % [GM, Re, degree, C, S] 
        checkmaxdeg( degree, 30 ) %#ok<NODEF>
    case 'lp100k' 
        % Moon
        % Load normalized coefficients and planetary constants
        load('aerolp100k.mat') % [GM, Re, degree, C, S] 
        checkmaxdeg( degree, 60 ) %#ok<NODEF>
    case 'lp165p'
        % Moon
        % Load normalized coefficients and planetary constants
        load('aerolp165p.mat') % [GM, Re, degree, C, S] 
        checkmaxdeg( degree, 60 ) %#ok<NODEF>
    case 'gmm2b'
        % Mars
        % Load normalized coefficients and planetary constants
        load('aerogmm2b.mat') % [GM, Re, degree, C, S] 
        checkmaxdeg( degree, 60 ) %#ok<NODEF>
    case 'eigengl04c'
        % Earth
        % Load normalized coefficients and planetary constants
        load('aeroeigengl04c.mat') % [GM, Re, degree, C, S] 
        checkmaxdeg( degree, 70 ) %#ok<NODEF>     
end

C_gravmod = C(1:maxdeg+2,1:maxdeg+2);
S_gravmod = S(1:maxdeg+2,1:maxdeg+2);

sF_gravmod = loc_gravLegendre_scaleFactor(maxdeg);

    function [scaleFactor] = loc_gravLegendre_scaleFactor( maxdegree )
        % loc_GRAVLEGENDRE internal function computing normalized associated
        % legendre polynomials, P, via recursion relations for spherical harmonic
        % gravity
        
        scaleFactor = zeros(maxdegree+3, maxdegree+3);
        
        % Seeds for recursion formula
        scaleFactor(1,1) = 0;
        scaleFactor(2,1) = 1;
        scaleFactor(2,2) = 0;
        
        for n = 2:maxdegree+2
            k = n + 1;
            
            for m = 0:n
                p = m + 1;
                % Scale Factor needed for normalization of dUdphi partial derivative
                
                if (n == m)
                    scaleFactor(k,k) = 0;
                elseif (m == 0)
                    scaleFactor(k,p) = sqrt( (n+1)*(n)/2);
                else
                    scaleFactor(k,p) = sqrt( (n+m+1)*(n-m));
                end
            end
        end
    end

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function checkmaxdeg( upperlimit, default )
        if exist('maxdeg','var') && maxdeg > upperlimit
            switch action
                case 'none'
                    % no passing of maxdeg messages
                case 'warning'
                    warning('aero:gravitysphericalharmonic:exceedMaxDeg', ...
                        ['Degree must be less than or equal to '...
                        'maximum degree of planetary model. Setting maximum', ...
                        ' degree to the maximum degree of the planetary model.']);
                case 'error'
                    error('aero:gravitysphericalharmonic:exceedMaxDeg', ...
                        ['Degree must be less than or equal to '...
                        'maximum degree of planetary model.']);
                otherwise
                    error('aero:gravitysphericalharmonic:unknownActionMaxDeg',...
                        ['Action for out of range input must be None, '...
                        'Warning or Error']);
            end
            maxdeg = upperlimit;
        else
            if ~exist('maxdeg','var')
            % maxdeg was not set
            maxdeg = default;            
            end
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function checkmodel( str )
        switch lower( str )
            case { 'egm2008', 'egm96', 'jgm3', 'lp100k', 'lp165p', 'gmm2b', 'custom', 'eigengl04c' }
                model = lower( str );
            otherwise
                error(message('aero:gravitysphericalharmonic:unknownModel'));
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function checkaction( str )
        switch lower( str )
            case { 'error', 'warning', 'none' }
                action = lower( str );
            otherwise
                error(message('aero:gravitysphericalharmonic:unknownAction'));
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function modeloraction( str )
        switch lower( str )
            case { 'egm2008', 'egm96', 'jgm3', 'lp100k', 'lp165p', 'gmm2b', 'eigengl04c' }
                model = lower( str );
            case { 'error', 'warning', 'none' }
                action = lower( str );
            otherwise
                error(message('aero:gravitysphericalharmonic:unknownString'));
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
end
