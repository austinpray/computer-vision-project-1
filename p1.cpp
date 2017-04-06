#include "opencv2/highgui.hpp"
#include "opencv2/opencv.hpp"
#include <iostream>
using namespace cv;
using namespace std;
#include <cmath>

double getYfromL(double L) {
    double epsilon = 0.008856;
    double k = 903.3;
    if (L > k*epsilon) {
        return pow((L+16)/116, 3);
    }

    return L/k;
}

int main(int argc, char** argv) {
  if(argc != 3) {
    cout << argv[0] << ": "
	 << "got " << argc-1 << " arguments. Expecting two: width height."
	 << endl ;
    return(-1);
  }

  int width = atoi(argv[1]);
  int height = atoi(argv[2]);
  int** RED1 = new int*[height];
  int** GREEN1 = new int*[height];
  int** BLUE1 = new int*[height];
  int** RED2 = new int*[height];
  int** GREEN2 = new int*[height];
  int** BLUE2 = new int*[height];

  for(int i = 0 ; i < height ; i++) {
    RED1[i] = new int[width];
    GREEN1[i] = new int[width];
    BLUE1[i] = new int[width];
    RED2[i] = new int[width];
    GREEN2[i] = new int[width];
    BLUE2[i] = new int[width];
  }

  for(int i = 0 ; i < height ; i++)
    for(int j = 0 ; j < width ; j++)
      {
	int r1, g1, b1;
	int r2, g2, b2;

	double x = (double)j/(double)width;
	double y = (double)i/(double)height;
	double Y = 1.0;

	double L = 90;
	double u = x * 512 - 255;
	double v = y * 512 - 255;


	// translate xyY to byte sRGB and Luv to byte sRGB
	double X = x * Y / y;
	double Z = (1 - x - y) * Y / y;
	// Y
	r1 = (int) (2.5623 * X + (-1.1661) * Y + (-0.3962) * Z) * 255;
    g1 = (int) ((-1.0215) * X + 1.9778 * Y + 0.0437 * Z) * 255;
    b1 = (int) (0.0752 * X + (-0.2562) * Y + 1.1810 * Z) * 255;

    double Xr = 0.9642;
    double Yr = 1.0;
    double Zr = 0.8249;
    double u0 = (4.0*Xr)/(Xr+15.0*Yr+3.0*Zr);
    double v0 = (9.0*Yr)/(Xr+15.0*Yr+3.0*Zr);
    double Y2 = getYfromL(L);
    double a = (1.0/3.0) * (((52.0 * L)/(u + (13.0 * L * u0))) - 1.0);
    double b = -5.0*Y2;
    double c = -1.0 * (1.0 / 3.0);
    double d = Y2*(((39*L)/(v+(13*L*v0)))-5);
    double X2 = (d - b)/(a - c);
    double Z2 = (X2*a) - b;
    // cout << "diags\n";
    // cout << "Xr: " << Xr << "\n";
    // cout << "Yr: " << Yr << "\n";
    // cout << "Zr: " << Zr << "\n";
    // cout << "u0: " << u0 << "\n";
    // cout << "v0: " << v0 << "\n";
    // cout << "Y2: " << Y2 << "\n";
    // cout << "a: " << a << "\n";
    // cout << "b: " << b << "\n";
    // cout << "c: " << c << "\n";
    // cout << "d: " << d << "\n";
    // cout << "X2: " << X2 << "\n";
    // cout << "Z2: " << Z2 << "\n";

	r2 = (int) (2.5623 * X2 + (-1.1661) * Y2 + (-0.3962) * Z2) * 255;
    g2 = (int) ((-1.0215) * X2 + 1.9778 * Y2 + 0.0437 * Z2) * 255;
    b2 = (int) (0.0752 * X2 + (-0.2562) * Y2 + 1.1810 * Z2) * 255;
    //cout << "working\n";
    //cout << r1 << "\n" << X << "\n";
    //cout << "not working \n";
    //cout << r2 << "\n" << X2 << "\n";
	// end

	RED1[i][j] = r1;
	GREEN1[i][j] = g1;
	BLUE1[i][j] = b1;
	RED2[i][j] = r2;
	GREEN2[i][j] = g2;
	BLUE2[i][j] = b2;
      }


  Mat R1(height, width, CV_8UC1);
  Mat G1(height, width, CV_8UC1);
  Mat B1(height, width, CV_8UC1);

  Mat R2(height, width, CV_8UC1);
  Mat G2(height, width, CV_8UC1);
  Mat B2(height, width, CV_8UC1);

  for(int i = 0 ; i < height ; i++)
    for(int j = 0 ; j < width ; j++) {

      R1.at<uchar>(i,j) = RED1[i][j];
      G1.at<uchar>(i,j) = GREEN1[i][j];
      B1.at<uchar>(i,j) = BLUE1[i][j];

      R2.at<uchar>(i,j) = RED2[i][j];
      G2.at<uchar>(i,j) = GREEN2[i][j];
      B2.at<uchar>(i,j) = BLUE2[i][j];
    }

  Mat xyY;
  Mat xyY_planes[] = {B1, G1, R1};
  merge(xyY_planes, 3, xyY);
  namedWindow("xyY",CV_WINDOW_AUTOSIZE);
  imwrite("xyY.png", xyY);

  Mat Luv;
  Mat Luv_planes[] = {B2, G2, R2};
  merge(Luv_planes, 3, Luv);
  namedWindow("Luv",CV_WINDOW_AUTOSIZE);
  imwrite("Luv.png", Luv);
  return(0);
}