# TrapFit (beta)
Many brain pathologies manifest with some degree of vascular dysfunction. Blood-oxygen-level dependent (BOLD) functional magnetic resonance imaging (fMRI) is frequently used to investigate the brain's haemodynamic response to stimuli. TrapFit allows an investigator to derive several parameters relevant to vascular physiology from block-design stimulus experiments. It works by generating a crude, trapezoid-shaped approximation of the haemodynamic response function from a BOLD timeseries so that the time-to-peak, time-to-baseline and amplitude of the response can be calculated. It is based primarily on the method described by [Dumas et al.(1)](http://onlinelibrary.wiley.com/doi/10.1002/ana.23566/abstract "Functional magnetic resonance imaging detection of vascular reactivity in cerebral amyloid angiopathy"), with some modifications. Many of the methodological considerations discussed in that paper  are also applicable for this method. If you want to cite this method, please consider citing the Dumas paper but be sure to consider the differences between the approaches. The actual curve fitting relies on the excellent ```stats::constrOptim``` function that comes with R. 

# Procedure
See ```theory.pdf```.

# Example Usage
Make sure you have the latest version of R installed. Note that working in a high performance computing environment like the SHARK-cluster is not necessary (see resource usage). However, when working on the SHARK-cluster, make sure the latest R module is imported. At the time of writing you can import the latest version with this command:

```
module load R/3.3.0
```

Download the entire repository by either pulling it or downloading it as a [.zip-file](https://github.com/JJHBARKEYWOLF/TrapFit/archive/master.zip) (see 'Clone or download' button).
Make sure the files ```TrapFit.R``` and ```functions.R``` are in the same directory.


The script is executed through ```Rscript```, which should be callable from the command line when you have R installed correctly (or when you imported the SHARK R-module).

To see the example usage and parameters, run the script as follows:
```
Rscript TrapFit.R
```

To run the script on your data, have a look at the example usage command.
```
Rscript ./TrapFit.R --path='./timeseries/' --TR=3 --onduration=20 --blockduration=48 --CVthresh=3 --restfirst=FALSE --outputname='test'
```
The path to the directory containing the source data files may be relative or absolute. 
Rscript is picky about parameters fed through the command line and I have not implemented any significant checks for this, so make sure the parameters you feed are correct and correctly formatted.



# Input data
The input data consist of a single (ROI-averaged) timeseries per subject. 
Some example datafiles are included in ```/example_data/``` to give you an idea of how they should be formatted. Basically, this amounts to running ```fslmeants``` from the [FSLutils](http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Fslutils) on your data with default parameters.
See the ```extract_TS.sh``` shell script for an example of this procedure.

# Output files
* ```$experiment$_fitparameters.csv``` this table contains the fit parameters and the calculated TTP and TTB values.
* ```$experiment$_raw_inputs.pdf``` plots of the mean timeseries that were provided as input data.
* ```$experiment$_allfits.pdf``` the trapezoid curves that were fit for every subject to evaluate fit quality.
* ```$experiment$_trapfit_settings.pdf``` logs the parameters you fed the script.
* ```$experiment$_blockaverages.pdf``` shows the average timeseries over all blocks for all subjects. This plot can be used to identify problems with your data.
* ```$experiment$_outlierblocks.pdf``` contains the blocks that were tagged as outliers because they exceeded the specified CV threshold.

# Resource use
CPU and memory use of the script is pretty low. Running on 20 subjects should take about 30 seconds on an average desktop computer.


# Bugs/issues/suggestions
At the time of writing, I have tested the script on only a single dataset. I'll be adding features and fixing bugs as I go along.
However, if you find something that doesn't work as you expected, feel free to contact me or open an issue.
If you have any suggestions regarding the implementation (used algorithms, parameters, etc.) please tell me. Any potential improvements are welcome!

##

# References
[1] Andrew Dumas, Gregory A Dierksen, M Edip Gurol, Amy Halpin, Sergi
Martinez-Ramirez, Kristin Schwab, Jonathan Rosand, Anand Viswanathan,
David H Salat, Jonathan R Polimeni, et al. Functional magnetic resonance
imaging detection of vascular reactivity in cerebral amyloid angiopathy. An-
nals of neurology, 72(1):76{81, 2012.
