/*
 * Copyright (c) 2014 M. Haris Usmani <me@harisusmani.com>
 *
 * The code in this repository is available under the MIT License.
 * For information on usage and redistribution, and for a DISCLAIMER OF ALL
 * WARRANTIES, see the file, "LICENSE.txt," in this distribution.
 *
 * See https://github.com/harisusmani/ofxCorrectPerspective for documentation
 *
 */
#pragma once

#include "ofMain.h"
#include "ofxCorrectPerspective.h"

class testApp : public ofBaseApp{

	public:
		void setup();
		void update();
		void draw();

		void keyPressed(int key);
		void keyReleased(int key);
		void mouseMoved(int x, int y );
		void mouseDragged(int x, int y, int button);
		void mousePressed(int x, int y, int button);
		void mouseReleased(int x, int y, int button);
		void windowResized(int w, int h);
		void dragEvent(ofDragInfo dragInfo);
		void gotMessage(ofMessage msg);

		ofImage my_image;
		ofImage output_img;
		string filename;
		float focal_length;
		float sensor_width;

		ofxCorrectPerspective cPerspective;
};
