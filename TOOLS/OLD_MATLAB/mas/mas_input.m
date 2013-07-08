%==========================================================================
% input
%==========================================================================

%- set paths --------------------------------------------------------------

path(path,'../');
path(path,'/home/fichtner/Matlab/TFA');
path(path,'/home/fichtner/mmap/m_map');

path_data='/home/fichtner/Data_Europa/DATA_Turkey/';                       % path of real data

path_synth='/nfs/stig/fichtner/SES3D_3.MER.T1/DATA/OUTPUT/T10.42/';        % path of new synthetics
output_path='/nfs/stig/fichtner/SES3D_3.MER.T1/ADJOINT/0.selection.tf/';        % path of new adjoint sources

path_synth_old='/nfs/stig/fichtner/SES3D_3.MER.T1/DATA/OUTPUT/T10.0/';    % path of old synthetics
path_logfile='/nfs/stig/fichtner/SES3D_3.MER.T1/ADJOINT/39.10s.tf/';           % path of old logfile

%- parameters for the computation of phase misfits ------------------------

    taper_width=7.0;            % Breite des tapers [s]
    width=10.0;                 % width of the Gaussian for the tf transform [s]
    dt_new=1.5;                 % new time increment for the tf transform [s]
    nu_max=1/(3*dt_new);        % maximum frequency, only for plotting purposes
    threshold=0.0;
    misfit_upper_bound=10.7;     % misfits above this bound are set to zero
    %misfit_upper_bound=1.7;     % misfits above this bound are set to zero
    
%- data specifications ----------------------------------------------------
                                          
period_tag='.10s';         % period extension of data directory


event_list=[1002 1005 1009 1010 1012 1027 1029 1030 1014 1017 1018 1019 1020 1021 1037 1039 1048 1054 1062 1065 1066 1068 1069 1071 1072 1073 1074 1075 1076];
event_list=[1062];

%- rotation parameters ----------------------------------------------------

rot_flag='yes';
rot_axis=[0 1 0]';
rot_angle=57.5;

%- misfit specification ---------------------------------------------------

misfit='tf';

%- corrections ------------------------------------------------------------

station_correction=1;
