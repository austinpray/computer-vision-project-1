#include "opencv2/highgui.hpp"
#include <iostream>
using namespace cv;
using namespace std;
#include <math.h>       /* pow */

#include <Eigen/Dense>


Eigen::Vector3d RGB_to_Luv(Eigen::Vector3d RGB) {
    Eigen::Matrix3d rgbmat;
    rgbmat << 0.412453, 0.357580, 0.180423,
       0.212671, 0.715160, 0.072169,
       0.019334, 0.119193, 0.950227;
      Eigen::Vector3d XYZ = rgbmat * RGB;

      double X = XYZ(0);
      double Y = XYZ(1);
      double Z = XYZ(2);
      //cout << "X " << X << "\nY " << Y << "\nZ " << Z << "\n";
      double L;
      if (Y > 0.008856) {
        L = pow(116*Y, (1.0/3.0));
      } else {
        L = 903.3*Y;
      }
      double uprime = (4.0*X)/(X + (15.0*Y) + (3.0*Z));
      double vprime = (9.0*Y)/(X + (15.0*Y) + (3.0*Z));
      double un = 0.19793943;
      double vn = 0.46831096;
      double u = 13.0*L*(uprime - 0.19793943);
      double v = 13.0*L*(vprime - 0.46831096);

      Eigen::Vector3d out;
      out(0) = L;
      out(1) = u;
      out(2) = v;
      return out;
}


double getYfromL(double L) {
    double epsilon = 0.008856;
    double k = 903.3;
    if (L > k*epsilon) {
        return pow((L+16)/116, 3);
    }

    return L/k;
}

Eigen::Vector3d Luv_to_XYZ(Eigen::Vector3d Luv) {

        double L = Luv(0);
        double u = Luv(1);
        double v = Luv(2);

        double Xr = 0.9642;
        double Yr = 1.0;
        double Zr = 0.8249;
        double u0 = (4.0*Xr)/(Xr+15.0*Yr+3.0*Zr);
        double v0 = (9.0*Yr)/(Xr+15.0*Yr+3.0*Zr);
        double Y = getYfromL(L);
        double a = (1.0/3.0) * (((52.0 * L)/(u + (13.0 * L * u0))) - 1.0);
        double b = -5.0*Y;
        double c = -1.0 * (1.0 / 3.0);
        double d = Y*(((39*L)/(v+(13*L*v0)))-5);
        double X = (d - b)/(a - c);
        double Z = (X*a) - b;

        return Eigen::Vector3d(X, Y, Z);
}

Eigen::Vector3d XYZ_to_RGB(Eigen::Vector3d XYZ) {
    //X, Y and Z input refer to a D65/2° standard illuminant.
    //sR, sG and sB (standard RGB) output range = 0 ÷ 255

    double var_X = XYZ(0) / 100.0;
    double var_Y = XYZ(1) / 100.0;
    double var_Z = XYZ(2) / 100.0;

    double var_R = var_X *  3.2406 + var_Y * -1.5372 + var_Z * -0.4986;
    double var_G = var_X * -0.9689 + var_Y *  1.8758 + var_Z *  0.0415;
    double var_B = var_X *  0.0557 + var_Y * -0.2040 + var_Z *  1.0570;

    if ( var_R > 0.0031308 ) {
        var_R = 1.055 * pow( var_R , ( 1 / 2.4 ) ) - 0.055;
    } else {
       var_R = 12.92 * var_R;
    }
    if ( var_G > 0.0031308 ) {
        var_G = 1.055 * pow( var_G , ( 1 / 2.4 ) ) - 0.055;
    } else {
        var_G = 12.92 * var_G;
    }
    if ( var_B > 0.0031308 ) {
        var_B = 1.055 * pow( var_B , ( 1 / 2.4 ) ) - 0.055;
    } else {
        var_B = 12.92 * var_B;
    }

    double sR = var_R * 255;
    double sG = var_G * 255;
    double sB = var_B * 255;
    return Eigen::Vector3d(sR, sG, sB);
}


double Linear_Scale(double x, double a, double b, double A, double B) {
    return (((x-a)*(B-A))/(b-a))+A;
}

void runOnWindow(int W1,int H1, int W2,int H2, Mat inputImage, char *outName) {
  int rows = inputImage.rows;
  int cols = inputImage.cols;

  vector<Mat> i_planes;
  split(inputImage, i_planes);
  Mat iB = i_planes[0];
  Mat iG = i_planes[1];
  Mat iR = i_planes[2];

  // dynamically allocate RGB arrays of size rows x cols
  int** R = new int*[rows];
  int** G = new int*[rows];
  int** B = new int*[rows];
  for(int i = 0 ; i < rows ; i++) {
    R[i] = new int[cols];
    G[i] = new int[cols];
    B[i] = new int[cols];
  }

  for(int i = 0 ; i < rows ; i++)
    for(int j = 0 ; j < cols ; j++) {
      R[i][j] = iR.at<uchar>(i,j);
      G[i][j] = iG.at<uchar>(i,j);
      B[i][j] = iB.at<uchar>(i,j);
    }


  //	   The transformation should be based on the
  //	   historgram of the pixels in the W1,W2,H1,H2 range.
  //	   The following code goes over these pixels
  double largestL = -1;
  double smallestL = -1;
  for(int i = H1 ; i <= H2 ; i++) {
    for(int j = W1 ; j <= W2 ; j++) {
      double r = R[i][j];
      double g = G[i][j];
      double b = B[i][j];
      Eigen::Vector3d RGB(r/255.0, g/255.0, b/255.0);

      Eigen::Vector3d Luv = RGB_to_Luv(RGB);

      int getL = 0;
      if (largestL == -1) {
        largestL = Luv(getL);
      }
      if (smallestL == -1) {
        smallestL = Luv(getL);
      }
      if (Luv(getL) > largestL) {
        largestL = Luv(getL);
      }
      if (Luv(getL) < smallestL) {
        smallestL = Luv(getL);
      }

    }
  }
  for(int i = H1 ; i <= H2 ; i++) {
    for(int j = W1 ; j <= W2 ; j++) {
      double r = R[i][j];
      double g = G[i][j];
      double b = B[i][j];
      Eigen::Vector3d RGB(r/255.0, g/255.0, b/255.0);

      Eigen::Vector3d Luv = RGB_to_Luv(RGB);

      Luv(0) = Linear_Scale(Luv(0), smallestL, largestL, 0, 100);
      Eigen::Vector3d XYZ = Luv_to_XYZ(Luv);
      Eigen::Vector3d RGBdelta = XYZ_to_RGB(XYZ);

      R[i][j] = RGBdelta(0);
      G[i][j] = RGBdelta(1);
      B[i][j] = RGBdelta(2);
    }
  }

  Mat oR(rows, cols, CV_8UC1);
  Mat oG(rows, cols, CV_8UC1);
  Mat oB(rows, cols, CV_8UC1);
  for(int i = 0 ; i < rows ; i++)
    for(int j = 0 ; j < cols ; j++) {
      oR.at<uchar>(i,j) = R[i][j];;
      oG.at<uchar>(i,j) = G[i][j];;
      oB.at<uchar>(i,j) = B[i][j];;
    }

  Mat o_planes[] = {oB, oG, oR};
  Mat outImage;
  merge(o_planes, 3, outImage);
  
  imwrite(outName, outImage);
}

int main(int argc, char** argv) {
  if(argc != 7) {
    cerr << argv[0] << ": "
	 << "got " << argc-1 
	 << " arguments. Expecting six: w1 h1 w2 h2 ImageIn ImageOut." 
	 << endl ;
    cerr << "Example: proj1b 0.2 0.1 0.8 0.5 fruits.jpg out.bmp" << endl;
    return(-1);
  }
  double w1 = atof(argv[1]);
  double h1 = atof(argv[2]);
  double w2 = atof(argv[3]);
  double h2 = atof(argv[4]);
  char *inputName = argv[5];
  char *outputName = argv[6];

  if(w1<0 || h1<0 || w2<=w1 || h2<=h1 || w2>1 || h2>1) {
    cerr << " arguments must satisfy 0 <= w1 < w2 <= 1"
	 << " ,  0 <= h1 < h2 <= 1" << endl;
    return(-1);
  }

  Mat inputImage = imread(inputName, CV_LOAD_IMAGE_UNCHANGED);
  if(inputImage.empty()) {
    cout <<  "Could not open or find the image " << inputName << endl;
    return(-1);
  }

  if(inputImage.type() != CV_8UC3) {
    cout <<  inputName << " is not a standard color image  " << endl;
    return(-1);
  }

  int rows = inputImage.rows;
  int cols = inputImage.cols;
  int W1 = (int) (w1*(cols-1));
  int H1 = (int) (h1*(rows-1));
  int W2 = (int) (w2*(cols-1));
  int H2 = (int) (h2*(rows-1));

  runOnWindow(W1, H1, W2, H2, inputImage, outputName);

  return 0;
}
