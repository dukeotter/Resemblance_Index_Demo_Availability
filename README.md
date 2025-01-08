# Resemblance_Index_Demo_Availability
This is the documentation of code to run inter-channel fiber resemblance index (RI) analysis based on multiple morpho-structural features through multi-photon microscopy (MPM) images of distinct fibrous structures (e.g. elastin and collagen fibers in this paper).
## 1. Software
The code runs in the MATLAB programming language. To install, download MATLAB from: [http://www.mathworks.com/products/matlab/](url)
## 2. Data format
Theoretically, any format of data that is compatible with MATLAB is suitable to de operated by the code. We include the example data in the folder "Examples" accompanied with the code to demonstrate the quantification in a concise way.
## 3. Run the analysis
Here are a total of 4 MATLAB files, where 'RI_main.m' is the main program and the others are functions that will be called during the running of the main program. Explanations for variables within the code, as well as the ideas in organizing each part of the code, have been noted in the main program. Additionally, possible parameters that should be modified accordingly to the data sets have been highlighted as 'modify x' (x refers to the numbering).
Generally, the main program can be divided into 7 parts:
### 1)	Images loading
Here the images were loaded up for consequent multi-parametric morpho-structural quantification (here parameters include but not limited to local coverage, orientation/azimuthal angle(_θ_), directional variance and waviness).
### 2)	Create the binary mask selecting the fiber-only regions
Here binary masks were created mainly based on the signal intensity. The fiber-only regions of the images were identified as 1 and background regions were identified as 0 by this mask. The masks would be used in consequent quantification and analysis operations.
### 3)	Acquire the voxel-wise local coverage
Here the voxel-wise 3D local coverage of the 3D stack is acquired. The method is described in our previous papers (_Optics Express_ **30**, 25718-25733 (2022)).
### 4)	Acquire the voxel-wise orientation
Here the voxel-wise 3D orientation of the images was acquired. The method was described in our previous papers (_Biomed. Opt. Express_ **6**, 2294–2310 (2015); _Biomaterials_, **116**, 34-47 (2017); _Biomaterials_, **179**, 96-108 (2018)).
### 5)	Calculate the directional variance and waviness 
Here the voxel-wise 3D directional variance and waviness of fiber-like structure were obtained based on the extracted orientation information. The directional variance and waviness calculation method was described in our previous papers (_Laser & Photonics Review_ **16** (2022), 2100576;_ Opt. Lett_. **47** (2022) 357–360; _Biomaterials_, **116**, 34-47 (2017); _Biomaterials_, **179**, 96-108 (2018)).
### 6)	Closest inter-channel fiber pixel searching and difference normalization
Here the code searched the closest inter-fiber pixel for the two images' masks, generated the minimal distance map and compared morpho-structural features' difference. Resulted difference maps were then normalized and summarized into the final RI map.
###7)	Perform post-processing
Here 'pretty' images of RI maps were generated. In detail, the raw intensity image was used to modulate degree of inter-fiber resemblance, with resemblance level indicated by gradient colors to enhance visual contrast. 
## 4. Example
Combined images of elastin and collagen fibers in different channels from _ex vivo_ human lung cancer as examples to test code were saved in the folder 'Examples_cancer'. 
