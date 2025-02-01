# Image Analysis and Computer Vision Toolbox

This repository contains Matlab code that implements a series of computer vision techniques for analyzing an image. The primary goals of the project are to:

1. **Feature Extraction**: Extract relevant features from the given image using various image processing methods.
2. **2D Reconstruction**: Rectify (2D reconstruct) the horizontal plane from the useful selected image lines and features.
3. **Calibration**: First extract a vertical vanishing point and then use it together with useful information extracted during the rectification step, in order to estimate the calibration matrix K containing the intrinsic
parameters of the camera.
4. **Localization**: Determine the relative pose (i.e. position and orientation) between the reference attached to the horizontal plane and the camera reference.
5. **Reconstruction**: Use the knowledge of K to rectify also a vertical facade.
