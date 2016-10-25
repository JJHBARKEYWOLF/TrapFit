#!/bin/bash
Usage() {
cat << EOF
Exracts mean timeseries of nonzero voxels using a z-statistic map as a mask and specified percentile bounds.

example usage: extract_TS.sh -featdirs ./data/ -lowperc 99 -topperc 100

basenames of the input feat directories are used as the output filenames.
OPTIONS:
    -featdirs 		Relative or absolute path to a directory containing input data (.feat directories)
    -lowperc		Desired lowest percentile. lowperc=70 will discard all voxels BELOW the 70th percentile.
    -topperc  		Desired upper percentile. topperc=80 will discard all voxels ABOVE the 80th percentile (defaults to 100).

EOF
    exit 1
}

if [ $# -lt 1 ]; then
    Usage
fi

# get options
while [ $# -gt 0 ] ;
do
    case $1 in
		-featdirs)
			featdirectories=$2
			shift 2
			;;
		-lowperc)
			lowperc=$2
			shift 2
			;;
		-topperc)
			topperc=$2
			shift 2
			;;
		-*) echo "Wrong option: <$1>"
			echo ""
			Usage
			;;
		*)   break
		;;
    esac
done

TSdir="./extract_TS_output/"
mkdir -p ${TSdir}


for inputdir in `ls -d ${featdirectories}*`
do
	echo ${inputdir}
	subjcode="$(basename "$inputdir")"
	
	echo "extracting timeseries of $subjcode"

	#make working directory
	stagingdir=${inputdir}/timeseries/topperc_${topperc}_lowperc_${lowperc}
	mkdir -p ${stagingdir}

	#copy relevant data files to dir
	cp ${inputdir}/thresh_zstat1.nii.gz ${stagingdir}/

	echo "calculating appropriate thresholding values from z-statistic map.."
	uthresh="$(fslstats ${stagingdir}/thresh_zstat1.nii.gz -P ${topperc})"
	lthresh="$(fslstats ${stagingdir}/thresh_zstat1.nii.gz -P ${lowperc})"
	echo "Thresholding values:"
	echo "upper Z:"
	echo "		" ${uthresh}
	echo "lower Z:"
	echo "		" ${lthresh}

	echo "taking voxels above these values and making a binary mask..."
	fslmaths ${stagingdir}/thresh_zstat1.nii.gz -uthr ${uthresh} -thr ${lthresh} -bin ${stagingdir}/thresh_zstat_nonzero_betw_percs_${lowperc}_${topperc}.nii.gz
	
	echo "binary mask contains (voxels, mm3) "
	fslstats ${stagingdir}/thresh_zstat_nonzero_betw_percs_${lowperc}_${topperc}.nii.gz -V
	echo "nonzero voxels"

	#extract timeseries of all voxels within binary mask
	fslmeants -i ${inputdir}/filtered_func_data.nii.gz -o ${stagingdir}/TS_top${perc}perc.txt -m ${stagingdir}/thresh_zstat_nonzero_betw_percs_${lowperc}_${topperc}.nii.gz

	#copy timeseries to aggregation folder
	cp ${stagingdir}/TS_top${perc}perc.txt ${TSdir}/${subjcode}_betw_percs_${lowperc}_${topperc}.txt

	echo "done for ${subjcode}"
done