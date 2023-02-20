At one point or another, the following Toolboxes were required.
1. Computer Vision Toolbox
2. GPU Coder
3. Image Processing Toolbox
4. Parallel Computing Toolbox
5. Signal procssing Toolbox
6. Statistics and Machine Learning Toolbox



This folder contains some code for executing SPIDER on various datasets.

The basic idea is that there are scripts in each folder (besides SPIDER_functions) called lib_*.m,
where * is a wildcard. Running this scripts will produce the following important objects.

G      - a feature matrix contiaining weak-form integrals of various library terms
labels - LaTeX strings associated with the terms in G 
scales - a vector containing the physics-informed scales of each term.

Once these are obtained, running a subsection of sparse_regression.m will produce some sparse relations 
from this data. 







For questions or comments, email Matthew Golden at mgolden30@gatech.edu

Disclaimer:
SPIDER is not a fully automated process (although Daniel Gurevich is working on an automated version in Python).

The basic paradigm of SPIDER is to
1. identifying symmetries
2. identifying irreducible representations
3. create a library in each of these representations
4. Weak-form integration
5. Sparse regression

Only steps 4 and 5 are automated. The first 3 require thinking by the modeler.
