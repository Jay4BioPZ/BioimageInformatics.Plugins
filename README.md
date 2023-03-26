# Bioimage Informatics 2023
## Basic Information
The code is mainly used for BII2023 homeworks and group project. Please refer to [BII2023](https://edu.epfl.ch/coursebook/fr/bioimage-informatics-BIO-410) for more details.

This is an example Maven project implementing an ImageJ command with a Pom adapted for BioImageIformatics at EPFL. It was initially cloned from https://github.com/BIOP/ijp-template-ij2 (14th September 2022). The basic programming framework was kept. 

New codes were included for resolving bioimaging tasks (i.e., multidimensional data pixel-wise operation, segmentation...) presented in lectures.

## Available Plugins
- <b>BleachingCorrection</b>: The fluorescence intensity of a life-imaging sample may decay exponentially as time goes because of photobleaching. The bleaching process can be represented as $I(t) = A*\exp{\(t/\tau\)}+C$. Given an image stack, this plugin fits the parameters $A$, $\tau$ and $C$. It would perform a pixel-wise computation to make the mean of each frame in a comparable level to compensate such gradual drop. A certain level of noise can be observed in later time frames due to higher SNR.
- <b>MultipleChannelsQuantification</b>: The plugin computes the gene expression (in form of fluorescence) of the cytoplasm at the very close periphery of every nucleus in a sequence of images. The input should includes two stacks, one for nucleus and the other for cytoplasm. The plugin first generates a mask of nucleus and then overlays ROIs onto the cytoplasm image. The measured values would be plotted into a time vs. intensity scatter plot.


