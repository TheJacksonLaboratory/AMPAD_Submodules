#!/usr/bin/env bash

# Usage:
# 	./update_sdata.sh
# 	Make sure to run the script from the AMPAD_Submodules/ directory!

dirname=${PWD##*/}
if [ $dirname = "AMPAD_Submodules" ]; then

	rsync --update -raz --progress ./* helix.jax.org:/sdata/carter-lab/carter/AMPAD/AMPAD_Submodules/
	ssh helix.jax.org chgrp -R carter_secure /sdata/carter-lab/carter/AMPAD/AMPAD_Submodules/
	ssh helix.jax.org chmod -R g+rx /sdata/carter-lab/carter/AMPAD/AMPAD_Submodules/

else

	echo "Please run this script from the AMPAD_Submodules/ directory!"

fi
