# Full-Spectrum-Monte-Carlo
Fast Monte Carlo simulation of full-spectrum backscattered diffuse reflectance (F-BDR) and FD-OCT

![full-wavelength FD-OCT B-scan image](https://github.com/Jianing-Mao/fullwaveOCT/blob/master/example/Bscan.png)
## Get Started
The code is builded in MATLAB R2020a and ANSI C. Please install the [MatScat](https://ww2.mathworks.cn/matlabcentral/fileexchange/36831-matscat) in your MATLAB.
## Prepare scattering parameters
Run [TT_para_gen.m](https://github.com/Jianing-Mao/fullwaveOCT/blob/master/Paras_gen/TT_para_gen.m) to generate the wavelength-dependent parameters for TT scattering phase functions.

Run [mus_gen.m](https://github.com/Jianing-Mao/fullwaveOCT/blob/master/Paras_gen/mus_gen.m) to generate the wavelength-dependent scattering coefficients.

Make sure the parameters (wavelength, angles, diameter, RI, number density...) are the same in both files. Copy the outputs to the [code](https://github.com/Jianing-Mao/fullwaveOCT/tree/master/code) directory.
## Begin Monte Carlo simulation
### Step I: prepare the tissue

Run [GenTissue.m](https://github.com/Jianing-Mao/fullwaveOCT/blob/master/code/GenTissue.m) to generate the simulation setting at all wavelengths.

### Step II: run the Monte Carlo

For the A-line F-BDR simulation:

Compile the [fullwaveMC.c](https://github.com/Jianing-Mao/fullwaveOCT/blob/master/code/fullwave_FBDR.c)
```sh
gcc fullwaveMC.c -lm -o test
```

Run the generated file with two arguments: (1) arg1:1/0, 1 means continuous simulation up to the wavelength sampling point specified in arg2, 0 means a single simulation of a wavelength sampling point specified in arg2; (2) arg2: 1-samplePoints, wavelength sampling point. For example, continuous simulation from sampling points 1-1024:
```sh
./test 1 1024
```
For the B-scan simulation:

Compile the [fullwaveMC_Bscan.c](https://github.com/Jianing-Mao/fullwaveOCT/blob/master/code/fullwaveMC_Bscan.c)
```sh
gcc fullwaveMC_Bscan.c -lm -o test
```

Run the generated file with one arguments: (1) arg1: the index of the detectors:
```sh
./test 128
```
You can generate a Bscan by running the simulation of all of the detector in parallel:
```sh
./test 1
./test 2
...
```

### Step III: show the results

For the A-line F-BDR simulation:

Run [lookFBDR.m](https://github.com/Jianing-Mao/fullwaveOCT/blob/master/code/lookFBDR.m) to show the F-BDR result.

For the B-scan simulation:

Run [lookFBDR_Bscan.m](https://github.com/Jianing-Mao/fullwaveOCT/blob/master/code/lookFBDR_Bscan.m) to show the F-BDR result.

Run [lookImage_Bscan.m](https://github.com/Jianing-Mao/fullwaveOCT/blob/master/code/lookImage_Bscan.m) to show the B-scan image with full-wavelength information.

# Example

# To be implemented
* Integrate all parameters and procedures into a calling function (now we need to carefully check the parameters of the generation, running, and viewing codes)
* More realistic light source
* ...

# Citation
If you use this code in your own research, please cite the following papers:

https://opg.optica.org/boe/fulltext.cfm?uri=boe-14-9-4644&id=536404

https://opg.optica.org/boe/fulltext.cfm?uri=boe-13-12-6317&id=518908

# Acknowledgement
The codes are built on [MatScat](https://ww2.mathworks.cn/matlabcentral/fileexchange/36831-matscat), open-source codes in [OMLC](https://omlc.org/software/mc/) and [OCT_MC](https://github.com/RMTariant/OCT_MC). We sincerely appreciate the authors for sharing their codes.
