# fullwaveOCT
Fast Monte Carlo simulation of full-wavelength backscattered diffuse reflectance and FD-OCT
## Get Started
The code is builded in MATLAB R2020a and ANSI C. Please install the [MatScat](https://ww2.mathworks.cn/matlabcentral/fileexchange/36831-matscat) in your MATLAB.
### Prepare scattering parameters
Run [TT_para_gen.m](https://github.com/Jianing-Mao/fullwaveOCT/blob/master/Paras_gen/TT_para_gen.m) to generate the wavelength-dependent parameters for TT scattering phase functions.

Run [mus_gen.m](https://github.com/Jianing-Mao/fullwaveOCT/blob/master/Paras_gen/mus_gen.m) to generate the wavelength-dependent scattering coefficients.

Make sure the parameters (wavelength, angles, diameter, RI, number density...) are the same in both files. Copy the outputs to the [code](https://github.com/Jianing-Mao/fullwaveOCT/tree/master/code) directory.
### Begin Monte Carlo simulation
Step I: prepare the tissue

Run [GenTissue.m](https://github.com/Jianing-Mao/fullwaveOCT/blob/master/code/GenTissue.m)
