#!/bin/bash

echo "Rename in all .m files " + $1 + " to " + $2
#find . -type f -name '*.m' -exec sed -i '' s/$1/$2/ {} +


perl -pi -w -e 's/$1/$2/g;' *.m

