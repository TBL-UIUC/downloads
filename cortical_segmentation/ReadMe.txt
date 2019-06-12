Project Name: 
Cortical Bone Segmentation 

Description: 
This code scans a thickness map obtained in ImageJ/Fiji from the BoneJ plugin in order to segment cortical/subchondral bone from 
trabecular. This is useful when determining cortical or trabecular morphological characteristics. 

Installation: 
Download corticalSegmentation.m

Usage: 
This code requires an ancillary code to read in ImageJ thickness maps. This can be self-generated in MATLAB based on user data. 

If the user has another method of finding the bone border, sections 3, 6, and 7, can be removed from the script and line 151 can be 
uncommented. Otherwise, strelVal 1, 2, and 3 can be adjusted to optimal values for different datasets. 

Credits: 
Author: Ida Ang and Mariana Kersh
Contact: ia267@cornell.edu and mkersh@illinois.edu
