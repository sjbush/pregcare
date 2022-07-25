#!/bin/bash

for i in $(find /t1-data/project/pregcare/sbush/MinION8/1.fast5 -mindepth 1 -maxdepth 1 -exec basename {} \;)

        do
	sbatch 3.run_guppy.sh $i

done
