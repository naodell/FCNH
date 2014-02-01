# !/bin/csh

set SOURCE = $1
set TARGET = $2

set files = (`ls $SOURCE`)
mkdir $TARGET

foreach file ($files)
    echo 'Copying '$SOURCE/$file' to '$TARGET
    dccp $SOURCE/$file $TARGET/.
end

