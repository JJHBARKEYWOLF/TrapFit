#!/bin/bash
#$ -S /bin/bash
#$ -q all.q
#$ -N trapfit_wrapper
#$ -cwd
#$ -j Y
#$ -V
#$ -m be

# R CMD BATCH --no-save --no-restore "--path='/exports/radiologie-hpc/JurriaanBarkeyWolf/vascular_reactivity/script/' --TR=3 --onduration=20 --blockduration=48 --multiple=TRUE --CVthresh=2 --outputname='testthres2'" trapfit.R
# R CMD BATCH "--path='/exports/radiologie-hpc/JurriaanBarkeyWolf/vascular_reactivity/script/' --TR=3 --onduration=20 --blockduration=48 --multiple=TRUE --CVthresh=2 --outputname='/exports/radiologie-hpc/JurriaanBarkeyWolf/vascular_reactivity/script/testthres2'" trapfit.R

R CMD BATCH "--args --path=/exports/radiologie-hpc/JurriaanBarkeyWolf/vascular_reactivity/script/data/ --TR=3 --onduration=20 --blockduration=48 --multiple=TRUE --CVthresh=2 --outputname=/exports/radiologie-hpc/JurriaanBarkeyWolf/vascular_reactivity/script/testthres2" trapfit.R

#example:
# $ R CMD BATCH --no-save --no-restore '--args a=1 b=c(2,5,6)' test.R test.out &