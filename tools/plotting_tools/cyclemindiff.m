function [diff2423,diff2422,diff2421] = cyclemindiff(spotRec,spotNorm, ...
                                        spotScale,tsi,low21,low22,low23, ...
                                        low24)
%CYCLEMINDIFF Given sunspot record and scaling, input TSI (from BTSI or
%various reconstructions), and the indices corresponding to the time of
%lowest sunspot counts, return the difference between cycle 24 and the
%cycle in question.
%
% November 22, 2024

%For xAll when spot variability can be sampled
if size(spotScale,3) == size(tsi,2) 
    xCorr=tsi-((repmat(spotRec,[1 size(tsi,2)])-min(spotRec))./...
        spotNorm)./squeeze(spotScale)';
else
    xCorr=tsi-((repmat(spotRec,[1 size(tsi,2)])-min(spotRec))./...
      spotNorm)./squeeze(mean(spotScale))';
end
diff2423=mean(xCorr(low24,:),1)-mean(xCorr(low23,:),1);
diff2422=mean(xCorr(low24,:),1)-mean(xCorr(low22,:),1);
diff2421=mean(xCorr(low24,:),1)-mean(xCorr(low21,:),1);
end

