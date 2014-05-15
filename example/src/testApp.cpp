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
#include "testApp.h"

//--------------------------------------------------------------
void testApp::setup(){
    filename="dc_1.jpg"; //Input Image
    //filename="floor_hunt.jpg";
    /*cout << "Filename: ";
	cin >> filename;*/

	my_image.loadImage(filename);
	filename.insert(filename.size()-4,"_ofxOUT");

	//Import EXIF Data//...
        //Note: This is yet to be Implemented. Please enter Focal
        //Length and Sensor Width manually in same unit (mm).
	focal_length=22; //Canon EOS M for dc_1.jpg
	sensor_width=22.3;

	//focal_length=4; //Nexus 5, for floor_hunt.jpg
	//sensor_width=4.59;

	//...EXIF Data Completed//
	double alpha, beta;
	cPerspective.Rectify_Image(my_image,focal_length,sensor_width,alpha,beta);
	cout << alpha << "  " << beta;
}

//--------------------------------------------------------------
void testApp::update(){

}

//--------------------------------------------------------------
void testApp::draw(){
    my_image.draw(0, 0);
}

//--------------------------------------------------------------
void testApp::keyPressed(int key){
    output_img.grabScreen(0,0,my_image.width,my_image.height);
    output_img.saveImage(filename);
}

//--------------------------------------------------------------
void testApp::keyReleased(int key){

}

//--------------------------------------------------------------
void testApp::mouseMoved(int x, int y ){

}

//--------------------------------------------------------------
void testApp::mouseDragged(int x, int y, int button){

}

//--------------------------------------------------------------
void testApp::mousePressed(int x, int y, int button){

}

//--------------------------------------------------------------
void testApp::mouseReleased(int x, int y, int button){

}

//--------------------------------------------------------------
void testApp::windowResized(int w, int h){

}

//--------------------------------------------------------------
void testApp::gotMessage(ofMessage msg){

}

//--------------------------------------------------------------
void testApp::dragEvent(ofDragInfo dragInfo){

}
