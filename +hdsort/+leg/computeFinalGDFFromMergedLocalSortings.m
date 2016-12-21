%function [R.gdf_merged, R.T_merged, R.localSorting, R.localSortingID R.gdf_discarded] = computeFinalGDFFromMergedLocalSortings(G)
function [R P] = computeFinalGDFFromMergedLocalSortings(G)
    P = struct();
    
    R.gdf_merged = [];
    R.T_merged = [];
    nextt = 1;
    
    R.gdf_discarded = [];
    R.T_discarded = [];
    nextt_dis = 1;
    
    R.localSorting = [];
    R.localSortingID = [];
    for g=1:length(G);
        localValidNeuronsIdx = find(G(g).templates.maxInThisGroup & cellfun(@(x) isempty(x), G(g).templates.masterTemplate)); 
        localDiscardedNeuronsIdx = find(G(g).templates.maxInThisGroup & cellfun(@(x) ~isempty(x), G(g).templates.masterTemplate));
     
        if ~isempty(localValidNeuronsIdx)
            localValidNeuronsId = localValidNeuronsIdx-1;
            validGdf = G(g).gdf(ismember(G(g).gdf(:,1), localValidNeuronsId),:); % the -1 is necessary since gdf ids start at zero!
            
            % set templates without any spikes in gdf as invalid
            validUnits = unique(validGdf(:,1));
            if length(validUnits) < length(localValidNeuronsId)
                warning('There were templates without spikes in the corresponding GDF!! They are removed. This should not happen');
                g
                G(g)
                validUnits
                localValidNeurons
            end
            localValidNeuronsId = intersect(localValidNeuronsId, validUnits);
            localValidNeuronsIdx = localValidNeuronsId+1;
            validGdf(:,1) = validGdf(:,1) + 1000*g;
            R.gdf_merged = [R.gdf_merged; validGdf];
            nT = length(localValidNeuronsId);
            R.T_merged(:,:, nextt:nextt+nT-1) = G(g).templates.wfs(:,:,localValidNeuronsIdx);
            nextt = nextt+nT;
            R.localSorting = [R.localSorting; g*ones(nT,1)];
            R.localSortingID = [R.localSortingID; localValidNeuronsId(:)];
        end
        if ~isempty(localDiscardedNeuronsIdx)
            localDiscardedNeuronsId = localDiscardedNeuronsIdx-1;
            discardedGdf = G(g).gdf(ismember(G(g).gdf(:,1), localDiscardedNeuronsId),:); % the -1 is necessary since gdf ids start at zero!
            
            % set templates without any spikes in gdf as invalid
            discardedUnits = unique(discardedGdf(:,1));
            
            if length(discardedUnits) < length(localDiscardedNeuronsId)
                warning('There were templates without spikes in the corresponding GDF!! They are removed. This should not happen');
                g
                G(g)
                discardedUnits
                localDiscardedNeurons
            end
            localDiscardedNeuronsId = intersect(localDiscardedNeuronsId, discardedUnits);
            localDiscardedNeuronsIdx = localDiscardedNeuronsId+1;
            discardedGdf(:,1) = discardedGdf(:,1) + 1000*g;
            R.gdf_discarded = [R.gdf_discarded; discardedGdf];
            nT = length(localDiscardedNeuronsId);
            R.T_discarded(:,:, nextt_dis:nextt_dis+nT-1) = G(g).templates.wfs(:,:,localDiscardedNeuronsIdx);
            nextt_dis = nextt_dis+nT;
        end
        
    end
    units = unique(R.gdf_merged(:,1));
    nU = length(units);
    assert(length(R.localSorting) == nU, 'must be identical');
    assert(length(R.localSortingID) == nU, 'must be identical');
    assert(size(R.T_merged,3) == nU, 'must be identical');
    assert(nextt-1 == nU, 'must be identical');
    
    units_dis = unique(R.gdf_discarded(:,1));
    nU_dis = length(units_dis);
    assert(size(R.T_discarded,3) == nU_dis, 'must be identical');
    assert(nextt_dis-1 == nU_dis, 'must be identical');
    