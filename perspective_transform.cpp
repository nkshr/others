#ifdef _DEBUG
#pragma comment(lib, "opencv_core2410d")
#pragma comment(lib, "opencv_highgui2410d")
#pragma comment(lib, "opencv_imgproc2410d")
#else
#pragma comment(lib, "opencv_core2410")
#pragma comment(lib, "opencv_highgui2410")
#pragma comment(lib, "opencv_imgproc2410")
#endif

#define UPPER_LEFT 0
#define BOTTOM_LEFT 1
#define UPPER_RIGHT 2
#define BOTTOM_RIGHT 3

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <iostream>

using namespace cv;
using namespace std;

int counter;
bool plot = true;
vector<Point2f> src_pts;

void plot_pt(int event, int x, int y, int flag, void* param)
{

	switch (event)
	{

	case EVENT_LBUTTONUP:
		src_pts.push_back(Point2f(x, y));
		counter++;
		cout << "success: " << "(" 
			<< x << ", " << y << ")" << endl;

		switch (counter)
		{

		case BOTTOM_LEFT:
			cout << "Plot bottom left source corner." << endl;
			break;

		case UPPER_RIGHT:
			cout << "Plot bottom right source corner." << endl;
			break;

		case BOTTOM_RIGHT:
			cout << "Plot upper right source corner." << endl;
			break;

		default:
			cout << "Press space key to proceed next process." << endl;
			break;

		}

		break;

	}

}

vector<Point2f> dst_pts;
bool set = false;
Rect rect;
bool draw_rect = false;

void set_chess(int event, int x, int y, int flag, void *param)
{

	switch (event)
	{

	case EVENT_LBUTTONDOWN:
		draw_rect = true;
		rect.x = x;
		rect.y = y;
		//dst_pts.push_back(Point2f(x, y));
		break;

	case EVENT_MOUSEMOVE:

		if (draw_rect)
		{
			rect.width = x - rect.x;
			rect.height = y - rect.y;
		}

		break;

	case EVENT_LBUTTONUP:
		draw_rect = false;
		dst_pts.push_back(Point2f(rect.x, rect.y));
		dst_pts.push_back(Point2f(rect.x, y));
		dst_pts.push_back(Point2f(x, y));
		dst_pts.push_back(Point2f(x, rect.y));
		cout << "Destination rectangle setted." << endl
			<< "Press space key to confirm source points and destination rectangle." << endl;
		break;

	}

	set = true;
}

int main(int argc, char **argv)
{
	Mat img = imread("undistorted.png");
	Mat temp = img.clone();

	namedWindow("disp");
	counter = 0;
	cout << "Plot upper left source corner." << endl;
	setMouseCallback("disp", plot_pt, (void *)&img);
	imshow("disp", img);
	waitKey();

	for (int i = 0; i < src_pts.size(); i++)
	{
		circle(img, src_pts.at(i), 3, Scalar(200, 0, 0), -1, CV_AA);
		circle(temp, src_pts.at(i), 3, Scalar(200, 0, 0), -1, CV_AA);
	}

	cout << "Set destination rectangle by dragging and dropping." << endl;
	setMouseCallback("disp", set_chess, (void *)&temp);
	imshow("disp", temp);
	waitKey();

	rectangle(temp, rect, Scalar(0, 0, 0));
	imshow("disp", temp);
	cout << "Press space key to proceed warping process." << endl;
	waitKey();

	Mat M = getPerspectiveTransform(src_pts, dst_pts);
	Mat dst_img;
	warpPerspective(img, dst_img, M, img.size());
	rectangle(dst_img, rect, Scalar(0, 0, 0));
	cout << "Image warped." << endl;
	imshow("disp", dst_img);
	waitKey();

	return 0;
}