function results = mosComputations(data, results)
    % Calcul de la Margin of Stability (MoS) - AP et ML
    % Basé sur la méthode de Hof et al. (2005)
    
    if (~isfield(data, 'Left') || isempty(data.Left.stamps)) && (~isfield(data, 'Right') || isempty(data.Right.stamps))
        results.Left.mos = false;
        results.Right.mos = false;
        return;
    end
    
    g = 9810; % mm/s^2
    
    % Paramètres détection HO (heel-off cinématique)
    zThresh_mm = 5;
    minSlope_mmps = 5;
    minSustain_fr = 5;
    fallbackFrac = 0.40;
    
    % === PIED GAUCHE ===
    if isfield(data, 'Left') && ~isempty(data.Left.stamps)
        rate_pts = data.Left.angleInfos.frequency;
        
        % COM (déjà calculé dans votre structure)
        COM = data.Left.markers.CentreOfMass;
        COMx = COM(:,1); COMy = COM(:,2); COMz = COM(:,3);
        
        % Vitesse COM
        vx = diff(COMx) * rate_pts;
        vy = diff(COMy) * rate_pts;
        
        % Longueur pendule inversé (cheville -> COM)
        Lz_L = abs(data.Left.markers.LANK(1:end-1,3) - COMz(1:end-1));
        w0_L = sqrt(g ./ max(Lz_L, eps));
        
        % xCOM (extrapolated COM)
        xCOMx_L = COMx(1:end-1) + vx ./ w0_L;
        xCOMy_L = COMy(1:end-1) + vy ./ w0_L;
        
        % Événements
        HS_L = data.Left.stamps.Left_Foot_Strike.frameStamp;
        TO_L = data.Left.stamps.Left_Foot_Off.frameStamp;
        nL = min(numel(HS_L), numel(TO_L));
        HS_L = HS_L(1:nL); TO_L = TO_L(1:nL);
        
        % Initialisation des vecteurs de résultats
        results.Left.MoS_AP_Mean = nan(nL, 1);
        results.Left.MoS_ML_Mean = nan(nL, 1);
        results.Left.MoS_HS_AP = nan(nL, 1);
        results.Left.MoS_HS_ML = nan(nL, 1);
        
        for k = 1:nL
            hs = max(1, HS_L(k));
            to = min(TO_L(k), length(xCOMy_L));
            if to <= hs, continue; end
            
            % Estimation HO (heel-off)
            HO = estimateHO_heelLift(data.Left.markers.LHEE(:,3), hs, to, rate_pts, ...
                                     zThresh_mm, minSlope_mmps, minSustain_fr, fallbackFrac);
            
            % Direction ML
            dirML = sign(data.Left.markers.LM5(hs,1) - data.Left.markers.RM5(hs,1));
            if dirML == 0, dirML = 1; end
            
            % Fenêtres temporelles
            A1 = hs : max(hs, HO-1);  % HS -> HO
            A2 = HO : max(HO, to-1);  % HO -> TO
            
            % MoS AP
            MoS_AP_heel = data.Left.markers.LHEE(A1,2) - xCOMy_L(A1);
            MoS_AP_toe = data.Left.markers.LTOE(A2,2) - xCOMy_L(A2);
            results.Left.MoS_AP_Mean(k) = mean([mean(MoS_AP_heel), mean(MoS_AP_toe)]);
            results.Left.MoS_HS_AP(k) = MoS_AP_heel(1);
            
            % MoS ML
            MoS_ML_ankle = (data.Left.markers.LANK(A1,1) - xCOMx_L(A1)) * dirML;
            MoS_ML_m5 = (data.Left.markers.LM5(A2,1) - xCOMx_L(A2)) * dirML;
            mos_ml_total = [MoS_ML_ankle; MoS_ML_m5];
            results.Left.MoS_ML_Mean(k) = mean(mos_ml_total);
            results.Left.MoS_HS_ML(k) = mos_ml_total(1);
        end
    end
    
    % === PIED DROIT ===
    if isfield(data, 'Right') && ~isempty(data.Right.stamps)
        rate_pts = data.Right.angleInfos.frequency;
        
        COM = data.Right.markers.CentreOfMass;
        COMx = COM(:,1); COMy = COM(:,2); COMz = COM(:,3);
        
        vx = diff(COMx) * rate_pts;
        vy = diff(COMy) * rate_pts;
        
        Lz_R = abs(data.Right.markers.RANK(1:end-1,3) - COMz(1:end-1));
        w0_R = sqrt(g ./ max(Lz_R, eps));
        
        xCOMx_R = COMx(1:end-1) + vx ./ w0_R;
        xCOMy_R = COMy(1:end-1) + vy ./ w0_R;
        
        HS_R = data.Right.stamps.Right_Foot_Strike.frameStamp;
        TO_R = data.Right.stamps.Right_Foot_Off.frameStamp;
        nR = min(numel(HS_R), numel(TO_R));
        HS_R = HS_R(1:nR); TO_R = TO_R(1:nR);
        
        results.Right.MoS_AP_Mean = nan(nR, 1);
        results.Right.MoS_ML_Mean = nan(nR, 1);
        results.Right.MoS_HS_AP = nan(nR, 1);
        results.Right.MoS_HS_ML = nan(nR, 1);
        
        for k = 1:nR
            hs = max(1, HS_R(k));
            to = min(TO_R(k), length(xCOMy_R));
            if to <= hs, continue; end
            
            HO = estimateHO_heelLift(data.Right.markers.RHEE(:,3), hs, to, rate_pts, ...
                                     zThresh_mm, minSlope_mmps, minSustain_fr, fallbackFrac);
            
            dirML = sign(data.Right.markers.RM5(hs,1) - data.Right.markers.LM5(hs,1));
            if dirML == 0, dirML = 1; end
            
            A1 = hs : max(hs, HO-1);
            A2 = HO : max(HO, to-1);
            
            MoS_AP_heel = data.Right.markers.RHEE(A1,2) - xCOMy_R(A1);
            MoS_AP_toe = data.Right.markers.RTOE(A2,2) - xCOMy_R(A2);
            results.Right.MoS_AP_Mean(k) = mean([mean(MoS_AP_heel), mean(MoS_AP_toe)]);
            results.Right.MoS_HS_AP(k) = MoS_AP_heel(1);
            
            MoS_ML_ankle = (data.Right.markers.RANK(A1,1) - xCOMx_R(A1)) * dirML;
            MoS_ML_m5 = (data.Right.markers.RM5(A2,1) - xCOMx_R(A2)) * dirML;
            mos_ml_total = [MoS_ML_ankle; MoS_ML_m5];
            results.Right.MoS_ML_Mean(k) = mean(mos_ml_total);
            results.Right.MoS_HS_ML(k) = mos_ml_total(1);
        end
    end
end

%% Fonction auxiliaire
function HO = estimateHO_heelLift(heelZ, HS, TO, fs, zThresh, minSlope, minSustain, fallbackFrac)
    HS = max(1, HS); TO = max(HS+1, TO);
    idx = HS:TO;
    z = heelZ(:);
    zSeg = z(idx);
    
    % Filtre Butterworth
    fc = 6; order = 4;
    [b, a] = butter(order, fc/(fs/2), 'low');
    zSm = filtfilt(b, a, zSeg);
    
    % Baseline
    n20 = max(1, round(0.02*fs));
    win = 1:min(n20, numel(zSm));
    base = median(zSm(win));
    
    % Dérivée
    dz = gradient(zSm) * fs;
    
    % Conditions
    elev = (zSm - base) > zThresh;
    slope = dz > minSlope;
    mask = elev & slope;
    
    HOrel = firstSustained(mask, minSustain);
    
    if isempty(HOrel)
        HO = HS + max(1, round(fallbackFrac*(TO-HS)));
    else
        HO = HS + HOrel - 1;
    end
    HO = min(max(HO, HS+1), TO-1);
end

function idx = firstSustained(mask, L)
    if ~any(mask)
        idx = [];
        return;
    end
    run = 0;
    for i = 1:numel(mask)
        if mask(i), run = run+1; else, run = 0; end
        if run >= L
            idx = i - L + 1;
            return;
        end
    end
    idx = [];
end