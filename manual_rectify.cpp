#pragma comment(lib, "gear_cvutil")
#include <cvlib.h>
#include <gear_cvutil.h>

#include <opencv2/opencv.hpp>
using namespace cv;

#include <iostream>
#include <fstream>
using namespace std;

namespace {
	vector<Point2f> pts1, pts2;
	const bool read_pts = true;
	const bool write_pts = false;
	const char * ofname = "matching_pts.txt";
	const char *ifname = "matching_point2.txt";
	const char *img1_name = "uimg1.png";
	const char *img2_name = "uimg2.png";
}

static void on_mouse(int event, int x, int y, int, void*data);

int main(int argc, char ** argv)
{
	Mat img1 = imread(img1_name, CV_LOAD_IMAGE_GRAYSCALE);
	Mat img2 = imread(img2_name, CV_LOAD_IMAGE_GRAYSCALE);
	Mat cimg;

	if (write_pts){
		connect_imgs(img1, img2, cimg);
		resize(cimg, cimg, Size(640, 240));
		namedWindow("disp");
		setMouseCallback("disp", on_mouse);
		imshow("disp", cimg);
		waitKey();
	
		ofstream ofs(ofname);
		ofs << pts1.size() << '\n';
		for (int i = 0; i < pts1.size(); ++i){
			cout << i << "th points" << endl;
			pts1[i].x *= img1.cols/(float)(cimg.cols/2);
			pts1[i].y *= img1.rows / (float)cimg.rows;
			pts2[i].x -= cimg.cols / 2;
			pts2[i].x *= img2.cols / (float)(cimg.cols/2);
			pts2[i].y *= img2.rows / (float)cimg.rows;
			ofs << pts1[i].x << " " << pts1[i].y << " "
				<< pts2[i].x << " " << pts2[i].y << '\n';
			cout << pts1[i].x << " " << pts1[i].y << " "
				<< pts2[i].x << " " << pts2[i].y << '\n';
		}
		ofs.close();
	
		cout << "pts1" << endl;
		print_pts2D(pts1);
		cout << "pts2" << endl;
		print_pts2D(pts2);
	}

	if (read_pts){
		vector<Point2f> pts3, pts4;
		read_matches(ifname, pts3, pts4);
		cout << "done read_matches" << endl;

		draw_matches(img1, pts3, img2, pts4, cimg);
		cout << "done << draw_matches" << endl;
		imwrite("cimg.png", cimg);
		resize(cimg, cimg, Size(640, 240));
		imshow("disp", cimg);
		waitKey();
	}

	return 0;
}

static void on_mouse(int event, int x, int y, int, void * data)
{
	static int counter = 0;
	if (event == EVENT_LBUTTONDOWN){
		if (counter % 2 == 0){
			cout << (counter+2)/2 << "th points" << endl;
			pts1.push_back(Point(x, y));
		}
		else{
			pts2.push_back(Point(x, y));
		}
		cout << "x, y : " << x << ", " << y << endl;
		++counter;
	}
}