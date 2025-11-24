function results = spatiotempComputations(data, results)
    
    % Si lorsqu'aucune donnée est envoyée (notamment pour le statique), il
    % manque des fields dans writeExcel, il est possible de les ajouter ici
    if (~isfield(data, 'Left') || isempty(data.Left.stamps)) && (~isfield(data, 'Right') || isempty(data.Right.stamps))
        results.Left.spatio = false;
        results.Right.spatio = false;
        return;
    end

    sides = fieldnames(data);
    for iS = 1:length(sides)
        side = sides{iS};
        if strcmp(side, 'Left')
            oppositeSide = 'Right';
        else
            oppositeSide = 'Left';
        end
        s = side(1);
        os = oppositeSide(1);

        if ~isfield(data.(side), 'stamps') || isempty(data.(side).stamps)
            continue
        end

        % Moment (en %) du toe off
        results.(side).pctToeOff = data.(side).stamps.(['CycleMarche' side])(data.(side).stamps.([side '_Foot_Off']).frameStamp);
        
        % Moment (en %) du toe off opposÃ©
        results.(side).pctToeOffOppose = data.(side).stamps.(['CycleMarche' side])(data.(side).stamps.([oppositeSide '_Foot_Off']).frameStamp);
        
        % Moment (en %) du contact talon opposÃ©
        results.(side).pctContactTalOppose = data.(side).stamps.(['CycleMarche' side])(data.(side).stamps.([oppositeSide '_Foot_Strike']).frameStamp);
        
        % Temps (en %) du simple appuie
        results.(side).pctSimpleAppuie = results.(side).pctContactTalOppose - results.(side).pctToeOffOppose;
        
        % Grandeur (en m) d'un pas et d'une foulÃ©e
        results.(side).distPas = abs(data.(side).markers.([s 'HEE'])(data.(side).stamps.([side '_Foot_Strike']).frameStamp(1),2) -  ... 
                         data.(side).markers.([os 'HEE'])(data.(side).stamps.([oppositeSide '_Foot_Strike']).frameStamp,2))/1000;
        results.(side).distFoulee = abs(diff(data.(side).markers.([s 'HEE'])(data.(side).stamps.([side '_Foot_Strike']).frameStamp,2)))/1000;

        % Temps (en s) d'une foulÃ©e
        results.(side).tempsFoulee = data.(side).tempsCycle;
        
        % Vitesse (en m/s) d'une foulÃ©e
        results.(side).vitFoulee = results.(side).distFoulee' ./ results.(side).tempsFoulee;
        
        % Calcul de la cadence (pas/minute) de marche
        results.(side).vitCadencePasParMinute = 1./results.(side).tempsFoulee * 60; % Convertir en "pas/minutes" 

        % Largeur de pas (mm) : Distance ML entre talons au contact 
        framesContact = data.(side).stamps.([side '_Foot_Strike']).frameStamp;
        framesContactOpp = data.(side).stamps.([oppositeSide '_Foot_Strike']).frameStamp;
        
        nStrides = length(framesContact) - 1; % Nombre de foulées
        stepWidthHeel = zeros(nStrides, 1);
        
        for iStride = 1:nStrides
            % Frame de contact du pied considéré
            frameCurrentSide = framesContact(iStride);
            frameNextSide = framesContact(iStride + 1);
            
            % Trouver le contact du pied opposé ENTRE ces deux contacts
            idxOppBetween = find(framesContactOpp > frameCurrentSide & ...
                                framesContactOpp < frameNextSide);
            
            if ~isempty(idxOppBetween)
                % Prendre le premier contact opposé dans cet intervalle
                frameOppSide = framesContactOpp(idxOppBetween(1));
                
                % Position ML (axe X) des talons
                posML_CurrentSide = data.(side).markers.([s 'HEE'])(frameCurrentSide, 1); % en mm
                posML_OppSide = data.(side).markers.([os 'HEE'])(frameOppSide, 1); % en mm
                
                % Distance ML = largeur de pas (en mm)
                stepWidthHeel(iStride) = abs(posML_CurrentSide - posML_OppSide);
            else
                % Si pas de contact opposé trouvé, mettre NaN
                stepWidthHeel(iStride) = NaN;
            end
        end
        
        % Stocker dans results (en mm pour cohérence avec autres mesures)
        results.(side).stepWidthHeel = stepWidthHeel;
    end
end