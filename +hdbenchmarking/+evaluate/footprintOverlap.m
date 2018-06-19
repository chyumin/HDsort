function [FPO, fpoFile] = footprintOverlap(artificialUnits, SpikeSortingResult, anaFolder)

%% Look at the footprint overlaps
fpoFile = fullfile(anaFolder, ['fpo.mat']);
try
    load(fpoFile)
catch
    %%
    [fp_sorted, resultsCutLeft] = SpikeSortingResult.getFootprint();
    assert( ~any(resultsCutLeft-resultsCutLeft(1)), 'Not all waveforms had the same cutLeft value!')
    resultsCutLeft = resultsCutLeft(1);
    
    vfp_sorted = hdsort.waveforms.t2v(fp_sorted);
    FPO.ku = []; FPO.dot_product = [];
    
    for ku = 1:artificialUnits.ST.nCells
        
        fp_gt = artificialUnits.FP.footprints(:,:,ku);
        
        cut0 = artificialUnits.FP.cutLeft-resultsCutLeft;
        cut1 = cut0 + size(fp_sorted,1) - 1;
        
        %fp_cut = fp((artificialUnits.FP.cutLeft-20):(artificialUnits.FP.cutLeft+54), :, :);
        fp_gt_cut = fp_gt(cut0:cut1, :, :);
        
        vfp_gt = hdsort.waveforms.t2v(fp_gt_cut);
        FPO.ku(:, :, ku) = fp_gt_cut;
        
        %% Dot product:
        FPO.dot_product(ku, :) = vfp_gt * vfp_sorted';
    end
    FPO.sorted = fp_sorted;
    
    [FPO.shiftedDistances, FPO.shiftedDistancesP] = ...
        hdsort.waveforms.tShiftedDistances(FPO.sorted, FPO.ku, 'maxShift', 10);
    
    save(fpoFile, 'FPO')
end

end