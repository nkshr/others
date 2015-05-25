#include <opencv2/opencv.hpp>
#include <iostream>

using namespace cv;
using namespace std;

void convertPoints(vector<Point3d> srcPoints, Mat rvec, Mat tvec, Mat cameraMatrix, vector<Point3d> &dstPoints)
{
	Mat temp1, temp2, src_pt_mat(4, 1, CV_64F), dst_pt_mat, rmat;
	Rodrigues(rvec, rmat);
	temp1 = (Mat_<double>(3, 4) <<
		rmat.at<double>(0, 0), rmat.at<double>(0, 1), rmat.at<double>(0, 2), tvec.at<double>(0, 0),
		rmat.at<double>(1, 0), rmat.at<double>(1, 1), rmat.at<double>(1, 2), tvec.at<double>(1, 0),
		rmat.at<double>(2, 0), rmat.at<double>(2, 1), rmat.at<double>(2, 2), tvec.at<double>(2, 0));
	temp2 = cameraMatrix * temp1;
	dstPoints.clear();
	for (unsigned int i = 0; i < srcPoints.size(); ++i)
	{
		src_pt_mat.at<double>(0, 0) = srcPoints[i].x;
		src_pt_mat.at<double>(1, 0) = srcPoints[i].y;
		src_pt_mat.at<double>(2, 0) = srcPoints[i].z;
		src_pt_mat.at<double>(3, 0) = 1;
		dst_pt_mat = temp2 * src_pt_mat;
		dstPoints.push_back(Point3d(dst_pt_mat.at<double>(0, 0), dst_pt_mat.at<double>(1, 0), dst_pt_mat.at<double>(2, 0)));
	}
}