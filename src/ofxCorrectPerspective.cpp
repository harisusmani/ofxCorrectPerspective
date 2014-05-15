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
#include "ofxCorrectPerspective.h"

using namespace cv;

typedef std::pair<double,int> mypair;
bool comparator ( const mypair& l, const mypair& r){
    return l.first < r.first;
    }

ofVec2f ofxCorrectPerspective::solveLinearSys(double a11,double a12,double a21,double a22,double b1,double b2){
    ofVec2f out;
    double det=(a11*a22)-(a12*a21);
    out.set((a22*b1-a12*b2)/det,(-a21*b1+a11*b2)/det);
    return out;
}

ofxCorrectPerspective::column_vector ofxCorrectPerspective::fitFunc4Lines(std::vector<std::vector<double> > L_vec,unsigned int r1,unsigned int r2,unsigned int r3,unsigned int r4, float f){
    column_vector starting_point(2);
    starting_point = ofRandom(PI/4),ofRandom(PI/4); // OR 0, 0
        find_min_bobyqa(cost_function(L_vec, r1, r2, r3, r4, f),
                        starting_point,
                        5,    // number of interpolation points
                        dlib::uniform_matrix<double>(2,1, -PI/2),  // lower bound constraint
                        dlib::uniform_matrix<double>(2,1, PI/2),   // upper bound constraint
                        PI/10,    // initial trust region radius
                        1e-100,  // stopping trust region radius
                        2000    // max number of objective function evaluations
        );
    return starting_point;
}

void ofxCorrectPerspective::distFunc(std::vector<std::vector<double> >L_vec,ofxCorrectPerspective::column_vector modelX,float f,double thresh,std::vector<int> & arIn,std::vector<int> & acIn){
    //MAKE ROTATION MATRIX
    ofMatrix4x4 R=ofMatrix4x4::newRotationMatrix(modelX(0)*180.0/PI, ofVec3f(-1, 0, 0), modelX(1)*180.0/PI, ofVec3f(0, -1, 0), 0, ofVec3f(0, 0, -1));
    double m[3][3] = {{R(0,0), R(0,1), R(0,2)}, {R(1,0), R(1,1), R(1,2)}, {R(2,0), R(2,1), R(2,2)}};
    cv::Mat R_mat = cv::Mat(3, 3, CV_64F, m);

    cv::Mat K_mat = (cv::Mat_<double>(3,3)<< f,0.0,0.0,0.0,f,0.0,0.0,0.0,1.0);

    cv::Mat K_c=K_mat.clone();
    K_c=K_c.inv();
    R_mat=R_mat.t();
    cv::Mat Hinv=K_mat*R_mat*K_c;

    double L_vec_D[3][L_vec[0].size()];
    for (int i = 0; i < L_vec[0].size(); ++i){
        L_vec_D[0][i]=L_vec[0][i];
        L_vec_D[1][i]=L_vec[1][i];
        L_vec_D[2][i]=L_vec[2][i];
    }

    cv::Mat L_vec_M=cv::Mat(3, L_vec[0].size(), CV_64F, L_vec_D);

    Hinv=Hinv.t();
    cv::Mat Lp=Hinv*L_vec_M;
    Lp.resize(2);

    double mag;
    for (int i=0; i<L_vec[0].size(); i++)
    {
        mag=(Lp.at<double>(0,i)*Lp.at<double>(0,i))+(Lp.at<double>(1,i)*Lp.at<double>(1,i));
        mag=sqrt(mag);
        Lp.at<double>(0,i)=Lp.at<double>(0,i)/mag;
        Lp.at<double>(1,i)=Lp.at<double>(1,i)/mag;
    }

    cv::Mat Lp_T=Lp.clone();
    Lp_T=Lp_T.t();
    cv::Mat C=Lp_T*Lp;
    multiply(C, C, C, 1.0);

    for (int i=0; i<L_vec[0].size(); i++)
    {
        for (int j=0; j<L_vec[0].size(); j++)
        {
            if (C.at<double>(i,j)<=thresh)
            {
                arIn.push_back(i);
                acIn.push_back(j);
            }
        }
    }
}

void ofxCorrectPerspective::Rectify_Image(ofImage & my_image,double focal_length,double sensor_width, double & alpha, double & beta){
    resize_image=1; //To Enable or Disable Resize

	//Rescale Image to be Max in 1000px
	if (resize_image && max(my_image.width,my_image.height)>1000)
    {
        float s=float(1000)/max(my_image.width,my_image.height);
        my_image.resize(floor(my_image.width*s),floor(my_image.height*s));
    }

    float f=focal_length*(float)max(my_image.width,my_image.height)/sensor_width;
    K.set(f,0.0,0.0,0.0,f,0.0,0.0,0.0,1.0);
    center.set(double(my_image.width)/2.0,double(my_image.height)/2.0);
    if (talk) cout << K;

    my_img_gray=my_image;
    my_img_gray.setImageType(OF_IMAGE_GRAYSCALE);

    //Converting Image to Image Double//
    image_double dub_image;
    ntuple_list lsd_out;
    unsigned int w=my_img_gray.width;
    unsigned int h=my_img_gray.height;
    unsigned char * imgP=my_img_gray.getPixels();

    // LSD parameters start
    double scale = 0.8;       // Scale the image by Gaussian filter to 'scale'.
    double sigma_scale = 0.6; // Sigma for Gaussian filter is computed as sigma = sigma_scale/scale.
    double quant = 2.0;       // Bound to the quantization error on the gradient norm.
    double ang_th = 22.5;     // Gradient angle tolerance in degrees.
    double eps = 0.0;         // Detection threshold, -log10(NFA).
    double density_th = 0.7;  // Minimal density of region points in rectangle.
    int n_bins = 1024;        // Number of bins in pseudo-ordering of gradient modulus.
    double max_grad = 255.0;  // Gradient modulus in the highest bin. The default value corresponds to the highest
                              // gradient modulus on images with gray levels in [0,255].
    // LSD parameters end

    bool verbose=0;

    dub_image = new_image_double(w,h);
    double px=0;
    cout << "\n--------\nInput data being written to image buffer \n";
    for(int j=0;j<(w*h);j++){
        px=imgP[j];
        dub_image->data[j] = px;
        if (verbose){
            cout << " " << dub_image->data[j];
        }
    }
    // Call LSD //
    lsd_out = LineSegmentDetection( dub_image, scale, sigma_scale, quant, ang_th, eps,
                               density_th, n_bins, max_grad, NULL );
    cout << "LSD has done it's thing!\n";
    if (talk) cout << "Number of Lines: "<< lsd_out->size << "Number of Dimensions: " << lsd_out->dim << "\n";

    if (verbose)
    {
    cout << "LSD Values: " << lsd_out->values[0] << " " << lsd_out->values[1] <<" " <<  lsd_out->values[2] <<" " <<  lsd_out->values[3] <<" " <<  lsd_out->values[4] <<" " <<  lsd_out->values[5] << "\n";
    }

    //Sorting in (Value, Index) pairs
        //http://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes
    std::vector<mypair> line_lengths;
    double sqd_distance;

    mesh.setMode(OF_PRIMITIVE_LINES);
    mesh.enableColors();

    ofVec3f first(0.0,0.0,0.0);
    ofVec3f second(0.0,0.0,0.0);
    double x1,x2,y1,y2;

    for(int j=0;j<(lsd_out->size*lsd_out->dim);j=j+5){
        x1=lsd_out->values[j];
        y1=lsd_out->values[j+1];
        x2=lsd_out->values[j+2];
        y2=lsd_out->values[j+3];
        sqd_distance=(x2-x1)*(x2-x1) + (y2-y1)*(y2-y1);

        line_lengths.push_back(make_pair(sqd_distance,j));

        //To Draw as Primitive Lines//
        /*first.set(x2,y2,0.0);
        second.set(x1,y1,0.0);
        mesh.addVertex(first);
        mesh.addColor(ofFloatColor(1.0, 0.0, 0.0));
        mesh.addVertex(second);
        mesh.addColor(ofFloatColor(1.0, 0.0, 0.0));*/
        //Lines Added, will be drawn//
       }
    sort(line_lengths.begin(),line_lengths.end());
    reverse(line_lengths.begin(), line_lengths.end());

    if (talk) cout << line_lengths[0].first << " " << line_lengths[0].second << " " << line_lengths[1].first << " " << line_lengths[1].second << "\n";

    unsigned int maxlines=700;
    unsigned int no_of_lines=min(lsd_out->size,maxlines);
    //Store these Lines pairs in a Matrix, in descending order of Distance
    cout << "Number of Lines: " << no_of_lines << "\n";
    L.resize(4);
    for (int i = 0; i < 4; ++i){
        L[i].resize(no_of_lines);
    }

    for (int j=0; j<no_of_lines; j++){
        L[0][j] = lsd_out->values[line_lengths[j].second];
        L[1][j] = lsd_out->values[(line_lengths[j].second)+1];
        L[2][j] = lsd_out->values[(line_lengths[j].second)+2];
        L[3][j] = lsd_out->values[(line_lengths[j].second)+3];
    }

    //LINE EXTENSION
    double extension_fac=.15; //For one side of the line
    double line_length;
    double line_gradient;
    double rise_angle;
    double delta_x;
    double delta_y;
    for (int j=0; j<no_of_lines; j++){
        x1=L[0][j];
        y1=L[1][j];
        x2=L[2][j];
        y2=L[3][j];
        line_length=sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
        line_gradient=(y2-y1)/(x2-x1);
        rise_angle=atan(abs(line_gradient));
        delta_x=cos(rise_angle)*(line_length*extension_fac);
        delta_y=sin(rise_angle)*(line_length*extension_fac);
        if (line_gradient<0)
        {
            if (x1>x2)
            {
                L[0][j]+=delta_x;
                L[1][j]-=delta_y;
                L[2][j]-=delta_x;
                L[3][j]+=delta_y;
            }
            else
            {
                L[2][j]+=delta_x;
                L[3][j]-=delta_y;
                L[0][j]-=delta_x;
                L[1][j]+=delta_y;
            }
        }
        else
        {
            if (x1>x2)
            {
                L[0][j]+=delta_x;
                L[1][j]+=delta_y;
                L[2][j]-=delta_x;
                L[3][j]-=delta_y;
            }
            else
            {
                L[2][j]+=delta_x;
                L[3][j]+=delta_y;
                L[0][j]-=delta_x;
                L[1][j]-=delta_y;
            }
        }
        x1=L[0][j];
        y1=L[1][j];
        x2=L[2][j];
        y2=L[3][j];

        L[2][j]-=center.x; //Move Origin to the Principle Point
        L[3][j]-=center.y;
        L[0][j]-=center.x;
        L[1][j]-=center.y;
    }

    //Finding ADJACENT Lines//
    bool adjflag=1;
    std::vector<int> ar; //To hold Adjacent Row Values
    std::vector<int> ac; //To hold Adjacent Column Values

    if (adjflag)
    {
    double athreshadj=10;

    std::vector<std::vector<bool> > adj; //Line x Line Inf Matrix Initialization, adj
                                            //Not Used at the Moment.
    adj.resize(no_of_lines); //Height
    for (int i = 0; i < no_of_lines; ++i){
        adj[i].resize(no_of_lines);
        for (int j = 0; j < no_of_lines; ++j){
            //adj[i][j]=1.0/0.0;
            adj[i][j]=0;
        }
    }

    ofVec2f v1,v2,x;
    athreshadj=abs(cos((athreshadj*PI)/180.0));
    for (int i = 0; i < no_of_lines; ++i){
        for (int j = i+1; j < no_of_lines; ++j){ //Everyline in-front
            v1.set(L[0][i]-L[2][i],L[1][i]-L[3][i]);
            v2.set(L[0][j]-L[2][j],L[1][j]-L[3][j]);
            v1.normalize();
            v2.normalize();
            if (abs(v1.dot(v2))<athreshadj) //acos(v1.dot(v2)) //So Angle is greater!
            {
               x=solveLinearSys(L[0][i]-L[2][i],-L[0][j]+L[2][j],L[1][i]-L[3][i],-L[1][j]+L[3][j],-L[2][i]+L[2][j],-L[3][i]+L[3][j]);
               if (not isinf(x.x) and not isinf(x.y))
               {
                   adj[i][j]=(x.x>=-DBL_EPSILON) && (x.x<=1+DBL_EPSILON) && (x.y>=-DBL_EPSILON) && (x.y<=1+DBL_EPSILON);
                   adj[j][i]=adj[i][j] || adj[j][i];
                   if (adj[i][j])
                   {
                      ar.push_back(i);
                      ac.push_back(j);
                   }
               }
            }
        }
    }
    }
    else
    {
    for (int i = 0; i < no_of_lines; ++i){ // nC2 Pairs
        for (int j = i+1; j < no_of_lines; ++j){
            ar.push_back(i);
            ac.push_back(j);
            }
        }
    }
    cout << "No. of Pairs: "<< ac.size() << "\n";

    //Adjacent Matrix FOUND//

    //Convert Line Segments to vector format for Rectification//
    std::vector<std::vector<double> > L_vec; //Line's Vector Form
    L_vec.resize(3); //Height
    for (int i = 0; i < 3; ++i){
        L_vec[i].resize(no_of_lines);
    }

    Mat A(3,3,CV_32F);
    Mat s, u, vt;
    A.at<float>(0,2)=1;
    A.at<float>(1,2)=1;
    A.at<float>(2,0)=0;
    A.at<float>(2,1)=0;
    A.at<float>(2,2)=0;
    for (int i = 0; i < no_of_lines; ++i){
            A.at<float>(0,0)=L[0][i];
            A.at<float>(0,1)=L[1][i];
            A.at<float>(1,0)=L[2][i];
            A.at<float>(1,1)=L[3][i];
            SVD::compute(A, s, u, vt);  //YY=U*S*V'
            vt=vt.t();
            vt.col(0)=vt.col(0)*-1;

            L_vec[0][i]=vt.at<float>(0,2);
            L_vec[1][i]=vt.at<float>(1,2);
            L_vec[2][i]=vt.at<float>(2,2);
    }

    //RANSAC//
    //int PairsC2=(ac.size()-1)*ac.size()/2;
    unsigned int maxTrials=400; //min(PairsC2*5,200);
    cout << "RANSAC: No. of Trials = " << maxTrials << endl;
    unsigned int trialcount=0;
    unsigned int r1, r2, r3, r4;
    unsigned int r_ind1, r_ind2;
    column_vector modelX;
    std::vector<int> arIn; //To hold Adjacent Row Values
    std::vector<int> acIn; //To hold Adjacent Column Values
    double thresh=0.001;

    std::vector<int> Best_arIn; //Best Adjacent Row Values
    std::vector<int> Best_acIn; //Best Adjacent Column Values
    unsigned int Bestscore=0; //Number of Inliers
    unsigned int score=0;
    column_vector Best_modelX;

    while (trialcount<maxTrials)
    {
        arIn.resize(0);
        acIn.resize(0);

        r_ind1=floor(ofRandom(ac.size()));
        r_ind2=floor(ofRandom(ac.size()));
        r1=ar[r_ind1];
        r2=ar[r_ind2];
        r3=ac[r_ind1];
        r4=ac[r_ind2];

        modelX=fitFunc4Lines(L_vec, r1, r2, r3, r4, f); //Fitting Function
        distFunc(L_vec,modelX,f,thresh,arIn,acIn);
        score=arIn.size();
        //score=distFunc2(L_vec,modelX,f,thresh,arIn,acIn, ar, ac); //Distance Function

        if (Bestscore<score)
        {
            Bestscore=score;
            Best_arIn=arIn;
            Best_acIn=acIn;
            Best_modelX=modelX;
            if (talk) cout << "No. of Inliers: "<< Bestscore<<endl;
        }

        trialcount++;
    }
    cout << "RANSAC Done." << endl;
    cout << "No. of Inliers: "<< Bestscore<<endl;

    column_vector solution=Best_modelX;

    if (talk) cout << "cost_function solution:\n" << Best_modelX << endl;
    cout << "cost_function solution:\n" << solution << endl;

    alpha=solution(0);
    beta=solution(1);

    //APPLYING HOMOGRAPHY for this SOLUTION//
    ofMatrix4x4 R=ofMatrix4x4::newRotationMatrix(solution(0)*180.0/PI, ofVec3f(-1, 0, 0), solution(1)*180.0/PI, ofVec3f(0, -1, 0), 0, ofVec3f(0, 0, -1));
    double m[3][3] = {{R(0,0), R(0,1), R(0,2)}, {R(1,0), R(1,1), R(1,2)}, {R(2,0), R(2,1), R(2,2)}};
    cv::Mat R_mat = cv::Mat(3, 3, CV_64F, m);

    cv::Mat K_mat = (cv::Mat_<double>(3,3)<< f,0.0,0.0,0.0,f,0.0,0.0,0.0,1.0);
    cv::Mat K_c= K_mat.clone();
    K_c=K_c.inv();

    cv::Mat C = (cv::Mat_<double>(3,3)<< 1,0,-center.x,0,1,-center.y,0,0,1);
    cv::Mat H=K_mat*R_mat*K_c*C;

    //Calclating Resultant Translation and Scale
    std::vector<Point2f> Ref_c;
    std::vector<Point2f> Ref_c_out;
    Ref_c.resize(4);
    Ref_c_out.resize(4);
    Ref_c[0].x=0;
    Ref_c[0].y=0;
    Ref_c[1].x=double(my_image.width);
    Ref_c[1].y=0;
    Ref_c[2].x=double(my_image.width);
    Ref_c[2].y=double(my_image.height);
    Ref_c[3].x=0;
    Ref_c[3].y=double(my_image.height);

    perspectiveTransform(Ref_c, Ref_c_out, H);

    if (talk) cout << "Ref Out New: " << Ref_c_out << endl;

    //Scalling:
    double scale_fac=abs((max(Ref_c_out[1].x,Ref_c_out[2].x)-min(Ref_c_out[0].x,Ref_c_out[3].x))/my_image.width); //Based on Length
    if (talk) cout << "Scale Factor: " << scale_fac << endl;

    Ref_c_out[0].x=Ref_c_out[0].x/scale_fac;
    Ref_c_out[0].y=Ref_c_out[0].y/scale_fac;
    Ref_c_out[1].x=Ref_c_out[1].x/scale_fac;
    Ref_c_out[1].y=Ref_c_out[1].y/scale_fac;
    Ref_c_out[2].x=Ref_c_out[2].x/scale_fac;
    Ref_c_out[2].y=Ref_c_out[2].y/scale_fac;
    Ref_c_out[3].x=Ref_c_out[3].x/scale_fac;
    Ref_c_out[3].y=Ref_c_out[3].y/scale_fac;

    Ref_c_out[1].x=Ref_c_out[1].x-Ref_c_out[0].x;
    Ref_c_out[1].y=Ref_c_out[1].y-Ref_c_out[0].y;
    Ref_c_out[2].x=Ref_c_out[2].x-Ref_c_out[0].x;
    Ref_c_out[2].y=Ref_c_out[2].y-Ref_c_out[0].y;
    Ref_c_out[3].x=Ref_c_out[3].x-Ref_c_out[0].x;
    Ref_c_out[3].y=Ref_c_out[3].y-Ref_c_out[0].y;
    Ref_c_out[0].x=Ref_c_out[0].x-Ref_c_out[0].x;
    Ref_c_out[0].y=Ref_c_out[0].y-Ref_c_out[0].y;
    if (talk) (cout << "Ref Out New: " << Ref_c_out << endl);

    H = getPerspectiveTransform( Ref_c, Ref_c_out ); //For the Translated/Scalled Image

    //Applying Homography//
    Mat src_img(cv:: Size (my_image.width, my_image.height),CV_8UC3,my_image.getPixels()); //OF to OpenCV
    cv::Mat dst_img;
	dst_img.create(src_img.size(), src_img.type());

	cv::warpPerspective(src_img, dst_img, H, src_img.size(), cv::INTER_LINEAR);

    //OpenCV to OF
    my_image.setFromPixels((unsigned char *) IplImage(dst_img). imageData,dst_img.size().width, dst_img.size().height,OF_IMAGE_COLOR);

    cout << endl <<"Press any key to Save Output Image." << endl;
}
