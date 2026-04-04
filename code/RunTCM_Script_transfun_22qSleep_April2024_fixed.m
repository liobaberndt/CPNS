function RunTCM_Script_transfun_22qSleep_April2024_fixed
    % Run subject-level thalamo-cortical DCM using transfer functions.
    %
    % This script applies the Shaw et al. (2020) thalamo-cortical neural
    % mass model to M/EEG data using a linearised transfer-function
    % on sleep EEG from individuals with 22q11.2DS (n=28) and siblings (n=17).
    %
    % Workflow summary:
    %   1) Configure local paths and dataset list
    %   2) Build and prepare one DCM per dataset
    %   3) Estimate parameters with spm_nlsi_GN
    %   4) Save fitted DCM outputs
    %
    % Requires atcm (thalamo cortical modelling package) and aoptim
    % (optimisation package)
    %
    % atcm: https://github.com/alexandershaw4/atcm
    % aoptim: https://github.com/alexandershaw4/aoptim
    
    %==========================================================================
    % SECTION 1: PATHS AND INPUT LIST
    %==========================================================================
    % Add the directory where this script is located.
    scriptPath = fileparts(mfilename('fullpath'));
    cd(scriptPath);
    addpath(scriptPath);
    addpath(fullfile(scriptPath, 'helperfunctions'));
    
    base_path = '/Users/liobaberndt/Dropbox/UoE';
    addpath(genpath(fullfile(base_path, 'DCM_TCM', 'atcm')));
    addpath(genpath(fullfile(base_path, 'DCM_TCM', 'aoptim')));
    
    data_path = fullfile(base_path, 'sleepdetectives', 'data', 'elifeNick');
    save_path = fullfile(base_path, 'sleepdetectives', 'dcm', 'eLifeNick', 'DCM_v01');
    Data.Datasets = fullfile(base_path, 'sleepdetectives', 'data', 'filelists', 'filelist_2.txt');
    
    % Ensure save directory exists
    if ~exist(save_path, 'dir')
        mkdir(save_path);
        fprintf('Created save directory: %s\n', save_path);
    end
    
    % Print active paths.
    fprintf('Data path: %s\n', data_path);
    fprintf('Save path: %s\n', save_path);
    fprintf('Dataset list: %s\n', Data.Datasets);
    
    % Switch to data directory for loading.
    cd(data_path);
    
    %==========================================================================
    % SECTION 2: DESIGN SETTINGS
    %==========================================================================
    Data.Design.X     = [];                % design matrix
    Data.Design.name  = {'undefined'};     % condition names
    
    % 1 = Wake, 2 = N1, 3 = N2, 4 = N3, 5 = REM
    Data.Design.tCode = 1;               % condition codes in SPM
    Data.Design.Ic    = 58; %find(strcmp(D.chanlabels,'Cz'));               % channel indices
    
    Data.Design.Sname = {'Cz'};            % channel (node) names
    Data.Prefix       = 'Wake_TCM_';      % outputted DCM prefix
    
    %==========================================================================
    % SECTION 3: MODEL STRUCTURE
    %==========================================================================
    % Model space: T is ns x ns (1 = forward, 2 = backward).
    T = [... % this is a 1-node model; nothing to put here...
        0];
    F = (T==1);
    B = (T==2);
    C = [1]';          % input(s)
    L = sparse(1,1);
    
    [p] = fileparts(which('atcm.integrate_1'));
    p = strrep(p,'+atcm','');
    addpath(p);
    
    %==========================================================================
    % SECTION 4: SUBJECT LOOP
    %==========================================================================
    for i = 1:length(Data.Datasets)
            
            %----------------------------------------------------------------------
            % 4.1 Build DCM structure for this dataset
            %----------------------------------------------------------------------
            DCM          = [];
            [fp fn fe]   = fileparts(Data.Datasets{i});
            DCM.name     = [Data.Prefix fn fe];
            
            DCM.xY.Dfile = Data.Datasets{i};  % original spm datafile
            Ns           = length(F);         % number of regions / modes
            DCM.xU.X     = Data.Design.X;     % design matrix
            DCM.xU.name  = Data.Design.name;  % condition names
            tCode        = Data.Design.tCode; % condition index (in SPM)
            DCM.xY.Ic    = Data.Design.Ic;    % channel indices
            DCM.Sname    = Data.Design.Sname; % channel names
            
            if exist(DCM.name);
                fprintf('Skipping model %d/%d - already exists!\n( %s )\n',i,length(Data.Datasets),DCM.name);
                continue;
            end
            
            %----------------------------------------------------------------------
            % 4.2 Set connectivity priors
            %----------------------------------------------------------------------
            DCM.A{1} = F;
            DCM.A{2} = B;
            DCM.A{3} = L;
            DCM.B{1} = DCM.A{1} | DCM.A{2};
            DCM.B(2:length(DCM.xU.X)) = DCM.B;
            DCM.C    = C;
            
            %----------------------------------------------------------------------
            % 4.3 Assign model functions
            %----------------------------------------------------------------------
            DCM.M.f  = @atcm.tc_hilge2;               % model function handle
            DCM.M.IS = @atcm.fun.alex_tf;            % Alex integrator/transfer function
            DCM.options.SpecFun = @atcm.fun.Afft;    % fft function for IS
            
            %----------------------------------------------------------------------
            % 4.4 Configure analysis options
            %----------------------------------------------------------------------
            fprintf('Running Dataset %d / %d\n',i,length(Data.Datasets));
            cd(data_path)
            fq =  [1 30];
            
            DCM.M.U            = sparse(diag(ones(Ns,1)));  %... ignore [modes]
            DCM.options.trials = tCode;                     %... trial code [GroupDataLocs]
            DCM.options.Tdcm   = [0 30e3];                   %... peristimulus time
            DCM.options.Fdcm   = fq;                    %... frequency window
            DCM.options.D      = 1;                         %... downsample
            DCM.options.han    = 1;                         %... apply hanning window
            DCM.options.h      = 4;                         %... number of confounds (DCT)
            DCM.options.DoData = 1;                         %... leave on [custom]
            %DCM.options.baseTdcm   = [-200 0];             %... baseline times [new!]
            DCM.options.Fltdcm = fq;                    %... bp filter [new!]
            DCM.options.UseButterband = fq;
    
            DCM.options.analysis      = 'CSD';              %... analyse type
            DCM.xY.modality           = 'LFP';              %... ECD or LFP data? [LFP]
            DCM.options.spatial       = 'LFP';              %... spatial model [LFP]
            DCM.options.model         = 'tc6';              %... neural model
            DCM.options.Nmodes        = length(DCM.M.U);    %... number of modes
            
            DCM.options.UseWelch      = 1010;
            DCM.options.FFTSmooth     = 1;
            DCM.options.BeRobust      = 0;
            DCM.options.FrequencyStep = 1/4;
            
            DCM.xY.name = DCM.Sname;
            
            %----------------------------------------------------------------------
            % 4.5 Prepare spectral data
            %----------------------------------------------------------------------
            fprintf('Preparing data...\n');
            DCM = atcm.fun.prepcsd(DCM);
            DCM.options.DATA = 1;
            fprintf('Processing data...\n');
            DCM.xY.y{1} = abs(DCM.xY.y{1});
            DCM.xY.y{1} = agauss_smooth(abs(DCM.xY.y{1}), .6)';
            fprintf('Data processing complete.\n');
            
            %----------------------------------------------------------------------
            % 4.6 Initialize parameters and priors
            %----------------------------------------------------------------------
            DCM = atcm.parameters(DCM,Ns);
                
            % Additional model initialization.
            DCM.M.solvefixed=0;      % oscillations == no fixed point search
            DCM.M.x = zeros(1,8,7);  % init state space: ns x np x nstates
            DCM.M.x(:,:,1)=-70;      % init pop membrane pot [mV]
            
            load([p '/newpoints3.mat'],'pE','pC')
    
            pE = spm_unvec(spm_vec(pE)*0,pE);
    
            pC.ID = pC.ID * 0;
            pC.T  = pC.T *0;
            
            pE.J = pE.J-1000;    
            pE.J(1:8) = log([.6 .8 .4 .6 .4 .6 .4 .4]);
            %pC.ID = pC.ID + 1/8;
            pE.L = 0;
            pC.a = pC.a*0;
    
            pE.Gb = pE.H;
            pC.Gb = [1   0   0   0   0   0   0   0;
                     0   1   1   0   0   0   0   0;
                     0   0   1   0   0   0   0   0;
                     0   0   0   1   1   0   0   0;
                     0   0   0   0   1   0   0   0;
                     0   0   0   0   1   1   0   0;
                     0   0   0   0   0   0   0   0;
                     0   0   0   0   0   0   1   0]/64;
    
            pC.J(1:8)=1/8;
            pC.d(1) = 1/8;

            % Fixed synaptic connections (set prior variance to zero).
            % Population index mapping used below:
            %   1=ss, 2=sp, 3=si, 4=dp, 5=di, 6=tp, 7=rt, 8=rl

            % AMPA fixed
            pC.H(5,4) = 0;  % AMPA_dp_to_di
            pC.H(1,6) = 0;  % AMPA_tp_to_ss
            pC.H(8,6) = 0;  % AMPA_tp_to_rl

            % GABA-A fixed
            pC.H(1,3) = 0;  % GABAA_si_to_ss

            % NMDA fixed
            pC.Hn(8,1) = 0; % NMDA_ss_to_rl
            pC.Hn(3,3) = 0; % NMDA_si_to_si
            pC.Hn(5,4) = 0; % NMDA_dp_to_di
            pC.Hn(6,4) = 0; % NMDA_dp_to_tp
            pC.Hn(1,6) = 0; % NMDA_tp_to_ss
            pC.Hn(8,6) = 0; % NMDA_tp_to_rl

            % GABA-B fixed
            pC.Gb(1,3) = 0; % GABAB_si_to_ss
            pC.Gb(4,5) = 0; % GABAB_di_to_dp
            pC.Gb(5,5) = 0; % GABAB_di_to_di
            pC.Gb(6,5) = 0; % GABAB_di_to_tp
            
            % Apply edited prior structures to DCM.
            DCM.M.pE = pE;
            DCM.M.pC = pC;
    
            %----------------------------------------------------------------------
            % 4.7 Estimate model parameters
            %----------------------------------------------------------------------
            w   = DCM.xY.Hz;
            Y   = DCM.xY.y{:};
            DCM.M.y  = DCM.xY.y;
            DCM.M.Hz = DCM.xY.Hz;
    
            ppE = DCM.M.pE;
            ppC = DCM.M.pC;
    
            fprintf('--------------- STATE ESTIMATION ---------------\n');
            fprintf('Search for a stable fixed point\n');
    
            xx = load([p '/newx.mat']); DCM.M.x = spm_unvec(xx.x,DCM.M.x);
            load('init_14dec','x');
            DCM.M.x = spm_unvec(x,DCM.M.x);
    
            x = atcm.fun.alexfixed(DCM.M.pE,DCM.M,1e-10);
            DCM.M.x = spm_unvec(x,DCM.M.x);
    
            norm(DCM.M.f(DCM.M.x,0,DCM.M.pE,DCM.M))
    
            fprintf('Finished...\n');
              
            fprintf('--------------- PARAM ESTIMATION ---------------\n');
            %fprintf('iteration %d\n',j);
    
            % Alex's version of the Levenberg-Marquard routine
            %M = AODCM(DCM);
    
            DCM.xY.Q = spm_Q(1/2,length(w)) * diag(w);
    
            %DCM.M.nograph = 1;
            DCM.M.nograph     =0;
            [Qp,Cp,Eh,F] = spm_nlsi_GN(DCM.M,DCM.xU,DCM.xY);
    
            %M.alex_lm;
    
            %M.compute_free_energy(M.Ep);
    
            %[Qp,Cp,Eh,F] = spm_nlsi_GN(DCM.M,DCM.xU,DCM.xY);
    
            %----------------------------------------------------------------------
            % 4.8 Store results and save
            %----------------------------------------------------------------------
            DCM.M.pE = ppE;
            DCM.Ep = Qp;%spm_unvec(M.Ep,DCM.M.pE);
            DCM.Cp = Cp;
    
            DCM.M.sim.dt  = 1./600;
            DCM.M.sim.pst = 1000*((0:DCM.M.sim.dt:(2)-DCM.M.sim.dt)');
    
            [y,w,G,s] = feval(DCM.M.IS,DCM.Ep,DCM.M,DCM.xU);
    
            DCM.pred = y;
            DCM.w = w;
            DCM.G = G;
            DCM.series = s;
            
            %DCM.Cp = atcm.fun.reembedreducedcovariancematrix(DCM,M.CP);
            %DCM.Cp = makeposdef(DCM.Cp);
            DCM.F  = F;%M.FreeEnergyF;
            %DCM.Cp = M.CP;
            cd(save_path)
            save(DCM.name); close all; clear global;
    end
end 