ofxCorrectPerspective
=====================

Introduction
------------
openFrameworks addon that performs automatic 2d rectification of images. It’s based on work done in “Shape from Angle Regularity” by Zaheer et al., ECCV 2012.
http://cvlab.lums.edu.pk/zaheer2012shape/

Demo Video: https://vimeo.com/95204456

Project Details: http://goo.gl/FvkFOh

Other libraries used in this addon:
- LSD (C), http://www.ipol.im/pub/art/2012/gjmr-lsd/
- dlib (C++), http://dlib.net/

License
-------
The code in this repository is available under the MIT License.
<br> Copyright (c) 2014 M. Haris Usmani, www.harisusmani.com

Installation
------------
- Copy to your openFrameworks/addons folder
- Generate a new project using projectGenerator
- AFTER generating project, download and copy dlib to ofxCorrectPerspective/src

Dependencies
------------
- ofxCV
- ofxOpenCv

Compatibility
-------------
openFrameworks V 0.8.0 Windows-CodeBlocks

Version History
---------------
### v1.0 5/15/2014
- initial version
