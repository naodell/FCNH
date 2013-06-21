#!/bin/csh

set count = 5

while ($count)
    #hadd -f -k combined_histos/{$1}_cut{$count}_{$2}_{$3}.root batchStage/{$3}/fcncHistograms_cut{$count}_*root > /dev/null &
    hadd -f -k combined_histos/fcnh_cut{$count}_2012_{$1}.root batchStage/{$1}/fcncHistograms_cut{$count}_*root > /dev/null &
    @ count--
end

