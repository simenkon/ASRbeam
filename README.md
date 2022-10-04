# ASRbeam
## Introduction
This is a 2D beam non-linear finite element software for reinforced concrete (RC) beams. Only material non-linearity is included,
but not geometric non-linearity. The purpose of the software is to calculate load effects owing
to ASR in statically indeterminate RC beams/frames.

The beam elements have 3 nodes (one mid node) and 7 degrees of freedom.

## INPUT
The input to the program are the files in the **INPUT** folder. A description of each input-file is listed below. 
- _nodes.txt_ -- defines the nodes' coordinate, and if they are free or fixed, i.e., the essential/kinematic boundary conditions. There are three essential conditions: 1) translation in x, 2) translation in z, and 3)rotation about y.  
- _eles.txt_ -- defines the element connectivity, the cross-section flag (id) for each element. 

