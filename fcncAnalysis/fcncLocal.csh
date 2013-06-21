#!/bin/csh

cp fcncAnalyzer_Template.C fcncAnalyzer.C

set count = `ls -1 histos/ | wc -l`

if (count == 0) then
    rm histos/fcnc*.root
endif

sed -i "s/SUFFIX/$1/g" fcncAnalyzer.C
sed -i "s/SELECTION/$2/g" fcncAnalyzer.C
sed -i "s/PERIOD/$3/g" fcncAnalyzer.C

root -l -b -q run.C
