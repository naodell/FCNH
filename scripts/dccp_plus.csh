# !/bin/csh

set SOURCE = $1
set TARGET = $2

set files = (`ls $SOURCE`)

set DIR = `echo $SOURCE | rev | cut -d/ -f1 | rev`
mkdir $TARGET/$DIR

foreach file ($files)
    echo 'Copying '$SOURCE/$file' to '$TARGET/$DIR
    dccp $SOURCE/$file $TARGET/$DIR 
end

