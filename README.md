# This is a demo code for PMC msDTI
/data：please download the test data from Data.txt

|-- PMC_msDWI
    |-- Code
    |   |-- src
    |   |-- tools
    |   |-- utils
    |   |-- multi_shot_recon_demo.m
    |   |-- single_shot_recon_demo.m
    |-- Data
    |   |-- ssDWI
    |   |-- msDWI

/utils: src

single_shot_recon_demo:

To demonstrate the capability that PMC can effectively removes motion-induced image blurring and artifacts

multi_shot_recon_demo:

To demonstrate the capability that PMC can effectively removes motion-induced image blurring and artifacts, while maintaining an invariance to the trace of the diffusion tensor. The background phases accompanying the use of large diffusion gradients were well tolerated by the LLR regularized reconstruction. 
