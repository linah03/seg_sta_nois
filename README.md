# seg_sta_nois

Tools for segmentation and feature extraction of Lossless derivative EEGIP London dataset.

Use EEGLAB with Matlab working directory set to:

eegip/london/derivatives/lossless

Add to the Matlab path:

addpath('code/dependencies/eeglab_asr_amica/');
addpath('derivatives/seg_sta_nois/code/scripts/');

Then start EEGLAB:

eeglab redraw


To run the segmentation use the batch_context EEGLAB extension to execute the "seg_sta_nois.htb" script (located in 'derivatives/seg_sta_nois/code/scripts') on the *.qcr.set files in the lossless derivatives directory.

Use the getERPfeatures script to extract measures for the resulting segmented files.
