#!/bin/bash -l
# SBATCH --job-name=triplot_launcher
# SBATCH --nodes=1
# SBATCH --ntasks=60
# SBATCH --cpus-per-task=1
# SBATCH --mem=80G
# SBATCH --output=CMAverse_launcher.log
# SBATCH --time=8-24:00:00



mod_sizes="14"
#mod_sizes="5 7 8 10 14"
#deepSplit_values="0"
#deepSplit_values="2 1 0"
#MEDissThress="0.4"
#MEDissThress="0.25 0.4"
cohorts=""

treatments ="quantile10vs90 quantile25vs75 minVSmax"

workdir=""


source activate mediation


for cohort in $cohorts; do

mkdir -p ${workdir}/${cohort}


WGCNA_params () {
Rscript --no-save --no-restore --verbose $workdir/CMAverse_script_WGCNA_parameters_01.R\ 
                --workdir=${workdir}/ \
                --cohort=$cohort \
                --outfolder=${workdir}/${cohort}/ \
                --max_cores=15
}
#WGCNA_params               

wait

CMAverse_analysis () {
#for deepSplit in $deepSplit_values; do
        for num in $mod_sizes; do
        	for treatment in treatments; do
                #for MEDissThres in $MEDissThress; do


        Rscript --no-save --no-restore --verbose $workdir/CMAverse_script_analysis_02.R\
                --workdir="${workdir}/" \
                --cohort="$cohort" \
                --outfolder="${workdir}/${cohort}/" \
                --mod_size="$num" \
                --treatment="$treatment" \
                --max_cores=8 &
                #--deepSplit=$deepSplit \
                #--MEDissThres=$MEDissThres \

wait
                #done   

done

wait

done

wait

}
CMAverse_analysis


wait

done

wait

exit
