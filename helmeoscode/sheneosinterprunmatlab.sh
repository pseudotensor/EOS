#!/bin/bash

cd $1

export DISPLAY=

#with no debugging:
matlab -nosplash -nodesktop -r shen_interp $@

# with debugging
#matlab -nosplash -nodesktop -Dgdb -r shen_interp $@
