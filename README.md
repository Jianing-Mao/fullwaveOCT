# Full-Wavelength-Monte-Carlo
Fast Monte Carlo simulation of full-wavelength backscattered diffuse reflectance (F-BDR) and FD-OCT
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

Compile the [fullwave_FBDR.c](https://github.com/Jianing-Mao/fullwaveOCT/blob/master/code/fullwave_FBDR.c)
```sh
gcc fullwave_FBDR.c -lm -o test
```

Run the generated file with two arguments: (1) arg1:1/0, 1 means continuous simulation up to the wavelength sampling point specified in arg2, 0 means a single simulation of a wavelength sampling point specified in arg2; (2) arg2: 1-samplePoints, wavelength sampling point. For example, continuous simulation from sampling points 1-1024:
```sh
./test 1 1024
```
For the B-scan simulation:

### Step III: show the results

For the A-line F-BDR simulation:

Run [lookFBDR.m](https://github.com/Jianing-Mao/fullwaveOCT/blob/master/code/lookFBDR.m) to show the F-BDR result.

# Example

# Acknowledgement
The codes are built on the open-source code for MC in [OMLC](https://omlc.org/software/mc/) and [OCT_MC](https://github.com/RMTariant/OCT_MC). We sincerely appreciate the authors for sharing their codes.
