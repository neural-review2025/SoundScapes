# SoundScapes

These scripts perform several steps:

1. S01_hyperaling.m: load the subject-wise response matrices X (M × N), where M=18 stimuli (6 objects per category) and N=number of vertices, for each one of the 360 ROIs of the HCP parcellation (Glasser et al., 2016), and apply hyperalignment (Haxby et al., 2011) (without scaling) across subjects. Data need to be organized in a way that there is one matrix (data.mat) for each ROI and hemisphere.
   
2. S02_PCA_ED.m: reduces data dimensionality via Principal Component Analysis (PCA; (Panichello & Buschman, 2021), and aligns PC-scores  using procrustes alignment (Barbosa et al., 2025) (without scaling) at the subject level, to ensure comparability across experiments. To quantify differences in representational geometry, in this script we compute Euclidean distances (and Principal Angle) at the object-level between OA object representations points and representations of scenes where the relevant object either was attended or a distractor in 3OA (using the first 3 PCs). To test whether the attended representation is closer than the distractor, we test euclidean difference scores using permutation-based paired t-tests at ROI-level.

data_demo.mat represents an example of the X matrix extracted from the left A5.
