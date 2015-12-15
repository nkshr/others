#ifdef WIN32
#include "cvlib.h"
#define _CRT_SECURE_NO_WARNINGS
#endif

#pragma comment(lib, "glfw3.lib")
#pragma comment(lib, "glu32.lib")
#pragma comment(lib, "opengl32.lib")
#include <GLFW/glfw3.h>

#include "gear_cvutil.h"

#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include <iostream>
#include <fstream>
#include <ios>

#define IFMFORD(m) (m.type() == CV_32F || m.type() == CV_64F) 

using namespace std;

void connect_imgs(const Mat &img1, const Mat &img2, Mat &cimg)
{
  CV_Assert(img1.size() == img2.size() &&
    ((img1.type() == CV_8U && img2.type() == CV_8U) ||
    (img1.type() == CV_8UC3 && img2.type() == CV_8UC3)));
  cimg.create(Size(img1.cols << 1, img1.rows), img1.type());
	for (int i = 0; i < cimg.rows; ++i){
		const uchar * pimg1 = img1.ptr<uchar>(i);
		const uchar * pimg2 = img2.ptr<uchar>(i);
		uchar * pcimg = cimg.ptr<uchar>(i);
		for (int j = 0; j < img1.cols; ++j){
			int j1 = j*cimg.channels();
			int j2 = j1 + img1.cols * cimg.channels();
			for (int k = 0; k < img1.channels(); ++k){
				pcimg[j1 + k] = pimg1[j1 + k];
				pcimg[j2 + k] = pimg2[j2 + k];
			}
		}
	}
}

void calcF(const vector<Point2f> &pts1,
  const vector<Point2f> &pts2, Mat &F)
{
  CV_Assert(pts1.size() == pts2.size() && pts1.size() > 7);
  Mat A((int)pts1.size(), 9, CV_32F);
  for (unsigned int i = 0; i < pts1.size(); ++i){
    float * pA = A.ptr<float>(i);
    pA[0] = static_cast<float>(pts2[i].x * pts1[i].x);
    pA[1] = static_cast<float>(pts2[i].x * pts1[i].y);
    pA[2] = static_cast<float>(pts2[i].x);
    pA[3] = static_cast<float>(pts2[i].y * pts1[i].x);
    pA[4] = static_cast<float>(pts2[i].y * pts1[i].y);
    pA[5] = static_cast<float>(pts2[i].y);
    pA[6] = static_cast<float>(pts1[i].x);
    pA[7] = static_cast<float>(pts1[i].y);
    pA[8] = 1.f;
  }
  SVD svdA(A);
  Mat Vt = svdA.vt;
  F.create(3, 3, CV_32F);
  float * pF = reinterpret_cast<float *>(F.data);
  float * pVt = Vt.ptr<float>(Vt.rows - 1);
  for (int i = 0; i < 9; ++i)
    pF[i] = pVt[i];
}


template void calcE(const vector<Point3_<double>> &pts1,
	const vector<Point3_<double>> &pts2, Mat &E);
template void calcE(const vector<Point3_<float>> &pts1,
	const vector<Point3_<float>> &pts2, Mat &E);
template <typename T>
void calcE(const vector<Point3_<T>> &pts1,
	const vector<Point3_<T>> &pts2, Mat &E)
{
	CV_Assert(pts1.size() == pts2.size() && pts1.size() > 7);
	Mat A((int)pts1.size(), 9, CV_64F);
	for (unsigned int i = 0; i < pts1.size(); ++i){
		double * pA = A.ptr<double>(i);
		pA[0] = pts2[i].x * pts1[i].x;
		pA[1] = pts2[i].x * pts1[i].y;
		pA[2] = pts2[i].x * pts1[i].z;
		pA[3] = pts2[i].y * pts1[i].x;
		pA[4] = pts2[i].y * pts1[i].y;
		pA[5] = pts2[i].y * pts1[i].z;
		pA[6] = pts2[i].z * pts1[i].x;
		pA[7] = pts2[i].z * pts1[i].y;
		pA[8] = pts2[i].z * pts1[i].z;
	}
	SVD svdA(A);
	Mat Vt = svdA.vt;
	E.create(3, 3, CV_64F);
	//float * pE = reinterpret_cast<float *>(E.data);
	double *pE = E.ptr<double>(0);
	double * pVt = Vt.ptr<double>(Vt.rows - 1);
	for (int i = 0; i < 9; ++i)
		pE[i] = pVt[i];
}

void draw_epi_line(const Mat &F,
  const vector<Point2f> &pts, Mat &img)
{
  CV_Assert(F.rows == 3 && F.cols == 3 &&
    !img.empty() && pts.size() > 0 && F.type() == CV_32F);
  float * l = new float[3];
  for (int i = 0; i < pts.size(); ++i){
    for (int j = 0; j < F.rows; ++j){
      const float * pF = F.ptr<float>(j);
      l[j] = pF[0] * pts[i].x + pF[1] * pts[i].y + pF[2];
    }
    Point2f start, end;
    start.x = 0;
    start.y = l[2] / -l[1];
    end.x = static_cast<float>(img.cols);
    end.y = (end.x * l[0] + l[2]) / -l[1];
    line(img, start, end, Scalar(255));
  }
  delete l;
}

void draw_epi_line(const Point2f &pt, const Point3f &_line,
  Mat &img1, Mat &img2)
{
  CV_Assert(!img1.empty() && !img2.empty());
  if (img1.type() == CV_8U || img2.type() == CV_8U){
    cvtColor(img1, img1, CV_GRAY2BGR);
    cvtColor(img2, img2, CV_GRAY2BGR);
  }
  static RNG rng = theRNG();
  Scalar color(rng(256), rng(256), rng(256));
  Point2f start, end;
  start.x = 0;
  start.y = -_line.z / _line.y;
  end.x = static_cast<float>(img2.cols);
  end.y = -(_line.x * end.x + _line.z) / _line.y;
  cv::line(img2, start, end, color);
  circle(img1, pt, 10, color);
}

template void draw_epi_lines(const vector<Point_<double>> &pts,
	const Mat &F, const int whichimg, Mat &img1, Mat &img2);
template void draw_epi_lines(const vector<Point_<float>> &pts,
	const Mat &F, const int whichimg, Mat &img1, Mat &img2);
template <typename T>
void draw_epi_lines(const vector<Point_<T>> &pts,
	const Mat &F, const int whichimg, Mat &img1, Mat&img2)
{
	for (int i = 0; i < pts.size(); ++i)
		draw_epi_line(pts[i], F, whichimg, img1, img2);
}

template void draw_epi_line(const Point_<double> &pt,
	const Mat &F, const int whichimg, Mat &img1, Mat &img2);
template void draw_epi_line(const Point_<float> &pt,
	const Mat &F, const int whichimg, Mat &img1, Mat &img2);
template <typename T>
void draw_epi_line(const Point_<T> &pt, const Mat &F, 
const int whichimg, Mat &img1, Mat &img2)
{
	CV_Assert(F.rows == 3 && F.cols == 3 &&
	 (F.type() == CV_32F || F.type() == CV_64F) &&
	 !img1.empty() && !img2.empty() &&
	  (whichimg == 1 || whichimg == 2));
	if (img1.type() == CV_8U || img2.type() == CV_8U){
		cvtColor(img1, img1, CV_GRAY2BGR);
		cvtColor(img2, img2, CV_GRAY2BGR);
	}
	static RNG rng = theRNG();
	Scalar color(rng(256), rng(256), rng(256));
	float * l = new float[3];
	Point2f start, end;
	switch (whichimg){
	case 1:
		if (F.type() == CV_32F){
			for (int i = 0; i < 3; ++i){
				const float * pF = F.ptr<float>(i);
				l[i] = static_cast<float>(pF[0] * pt.x
					+ pF[1] * pt.y + pF[2]);
			}
		}
		else{
			for (int i = 0; i < 3; ++i){
				const double * pF = F.ptr<double>(i);
				l[i] = static_cast<float>(pF[0] * pt.x +
					pF[1] * pt.y + pF[2]);
			}
		}
		start.x = 0;
		start.y = l[2] / -l[1];
		end.x = static_cast<float>(img1.cols);
		end.y = (end.x * l[0] + l[2]) / -l[1];
		line(img2, start, end, color);
		circle(img1, pt, 10, color);
		break;
	case 2:
		Mat Ft = F.t();
		if (Ft.type() == CV_32F){
			for (int i = 0; i < 3; ++i){
				const float * pFt = Ft.ptr<float>(i);
				l[i] = static_cast<float>(pFt[0] * pt.x
					+ pFt[1] * pt.y + pFt[2]);
			}
		}
		else{
			for (int i = 0; i < 3; ++i){
				const double * pFt = Ft.ptr<double>(i);
				l[i] = static_cast<float>(pFt[0] * pt.x +
					pFt[1] * pt.y + pFt[2]);
			}
		}
		start.x = 0;
		start.y = l[2] / -l[1];
		end.x = static_cast<float>(img2.cols);
		end.y = (end.x * l[0] + l[2]) / -l[1];
		line(img1, start, end, color);
		circle(img2, pt, 10, color);
		break;
	}
	delete l;
}

void draw_pline(Mat &img)
{
	int offset = img.rows/20;
	for (int i = 0; i < 20; ++i){
		Point start(0, offset*i), end(img.cols, offset*i);
		line(img, start, end, Scalar(255));
	}
}

void extractRTfromE(const Mat &E, int flag, Mat &R, Mat &T)
{
	CV_Assert(IFMFORD(E));
	printf("extract : %p\n", T.data);


	if (E.type() == CV_32F){
		T.create(Size(1, 3), CV_32F);
		Mat W(3, 3, CV_32F, Scalar(0));
		W.at<float>(0, 1) = -1;
		W.at<float>(1, 0) = 1;
		W.at<float>(2, 2) = 1;
		SVD decomp = SVD(E);
		Mat U = decomp.u;
		Mat Vt = decomp.vt;
		float *pT = T.ptr<float>(0);
		float *pU = U.ptr<float>(0);
		switch (flag){
		case 0:
			R = U * W * Vt;
			pT[0] = pU[2];
			pT[1] = pU[5];
			pT[2] = pU[8];
			break;
		case 1:
			R = U * W * Vt;
			pT[0] = pU[2];
			pT[1] = pU[5];
			pT[2] = pU[8];
			break;
		case 2:
			R = U * W.t() * Vt;
			pT[0] = pU[2];
			pT[1] = pU[5];
			pT[2] = pU[8];
			break;
		case 3:
			R = U * W.t() * Vt;
			pT[0] = pU[2];
			pT[1] = pU[5];
			pT[2] = pU[8];
			break;
		default:
			exit(EXIT_FAILURE);
		}
	}
	else{
		T.create(Size(1, 3), CV_64F);
		Mat W(3, 3, CV_64F, Scalar(0));
		W.at<double>(0, 1) = -1;
		W.at<double>(1, 0) = 1;
		W.at<double>(2, 2) = 1;
		SVD decomp = SVD(E);
		Mat U = decomp.u;
		Mat Vt = decomp.vt;
		double *pT = T.ptr<double>(0);
		double *pU = U.ptr<double>(0);
		switch (flag){
		case 0:
			R = U * W * Vt;
			pT[0] = pU[2];
			pT[1] = pU[5];
			pT[2] = pU[8];
			break;
		case 1:
			R = U * W * Vt;
			pT[0] = pU[2];
			pT[1] = pU[5];
			pT[2] = pU[8];
			break;
		case 2:
			R = U * W.t() * Vt;
			pT[0] = pU[2];
			pT[1] = pU[5];
			pT[2] = pU[8];
			break;
		case 3:
			R = U * W.t() * Vt;
			pT[0] = pU[2];
			pT[1] = pU[5];
			pT[2] = pU[8];
			break;
		default:
			exit(EXIT_FAILURE);
		}
	}
	printf("extract : %p\n", T.data);
}

template void cnv_pix_to_pt3D(const Point_<double> &pix,
	const double snsr_sz, const double flen,
	const Mat &img, Point3_<double> &pt3D);
template void cnv_pix_to_pt3D(const Point_<float> &pix,
	const float snsr_sz, const float flen,
	const Mat &img, Point3_<float> &pt3D);
template <typename T>
void cnv_pix_to_pt3D(const Point_<T> &pix,
 const T snsr_sz, const T flen,
 const Mat &img, Point3_<T> &pt3D)
{
	pt3D.x = pix.x * snsr_sz;
	pt3D.y = pix.y * snsr_sz;
	pt3D.z = snsr_sz * flen;
}

void decomp_aff(const Mat &aff, Mat &R, Mat &T)
{
	CV_Assert((aff.type() == CV_32F || aff.type() == CV_64F) 
	&& aff.cols == 4 && aff.rows == 3);
	if (aff.type() == CV_32F){
		R.create(3, 3, CV_32F);
		T.create(3, 1, CV_32F);
		for (int i = 0; i < 3; ++i){
			const float * paff = aff.ptr<float>(i);
			float * pR = R.ptr<float>(i);
			for (int j = 0; j < 3; ++j)
				pR[j] = paff[i];

			float * pT = T.ptr<float>(i);
			pT[0] = paff[3];
		}
	}
	else{
		R.create(3, 3, CV_64F);
		T.create(3, 1, CV_64F);
		for (int i = 0; i < 3; ++i){
			const double * paff = aff.ptr<double>(i);
			double * pR = R.ptr<double>(i);
			for (int j = 0; j < 3; ++j)
				pR[j] = paff[i];

			double * pT = T.ptr<double>(i);
			pT[0] = paff[3];
		}
	}
}

double calc_pts_ydiff(vector<Point2f> pts1, vector<Point2f> pts2)
{
	CV_Assert(pts1.size() == pts2.size());
	double ydiff = 0;
	for (int i = 0; i < pts1.size(); ++i)
		ydiff += pow(pts1[i].y - pts2[i].y, 2.0);

	ydiff = sqrt(ydiff);
	ydiff /= pts1.size();
	return ydiff;
}

double calc_pts_diff(vector<Point2f> pts1, vector<Point2f> pts2)
{
	CV_Assert(pts1.size() == pts2.size());
	double diff = 0;
	for (int i = 0; i < pts1.size(); ++i){
		double xdiff = pow(pts1[i].x - pts2[i].x,2.0);
		double ydiff = pow(pts1[i].y - pts2[i].y, 2.0);
		diff += xdiff + ydiff;
	}
	diff /= pts1.size();
	diff = sqrt(diff);
	return diff;
}

void print_pts2D(vector<Point2f> pts)
{
	for (int i = 0; i < pts.size(); ++i)
		cout << "x, y : " <<
		 pts[i].x << ", " << pts[i].y << endl;
}

double checkF(const vector<Point2f> &pts1,
	const vector<Point2f> &pts2, const Mat &F)
{
	CV_Assert(pts1.size() == pts2.size() && 
	F.rows == 3 && F.cols == 3 &&
	 (F.type() == CV_32F || F.type() == CV_64F));
	double err = 0;
	double * Fx = new double[3];
	if (F.type() == CV_64F){
		const double * pF = F.ptr<double>(0);
		for (int i = 0; i < pts1.size(); ++i){
			Fx[0] = pF[0] * pts1[i].x + pF[1] * pts1[i].y + pF[2];
			Fx[1] = pF[3] * pts1[i].x + pF[4] * pts1[i].y + pF[5];
			Fx[2] = pF[6] * pts1[i].x + pF[7] * pts1[i].y + pF[8];
			err += pts2[i].x*Fx[0] + pts2[i].y*Fx[1] + Fx[2];
		}
	}
	else{
		const float * pF = F.ptr<float>(0);
		for (int i = 0; i < pts1.size(); ++i){
			Fx[0] = pF[0] * pts1[i].x + pF[1] * pts1[i].y + pF[2];
			Fx[1] = pF[3] * pts1[i].x + pF[4] * pts1[i].y + pF[5];
			Fx[2] = pF[6] * pts1[i].x + pF[7] * pts1[i].y + pF[8];
			err += pts2[i].x*Fx[0] + pts2[i].y*Fx[1] + Fx[2];
		}
	}
	err /= static_cast<double>(pts1.size());
	delete[] Fx;
	return err;
}


template double checkE(const vector<Point3_<double>> &pts1,
	const vector<Point3_<double>> &pts2, const Mat &E);
template double checkE(const vector<Point3_<float>> &pts1,
	const vector<Point3_<float>> & pts2, const Mat &E);
template <typename T>
double checkE(const vector<Point3_<T>> &pts1,
 const vector<Point3_<T>> &pts2, const Mat &E)
{
	CV_Assert(pts1.size() == pts2.size() && E.size() == Size(3, 3)
	&& (E.type() == CV_64F || E.type() == CV_32F));
	double err = 0;
	double *Ex = new double[3];
	if (E.type() == CV_64F){
		const double *pE = E.ptr<double>(0);
		for (int i = 0; i < pts1.size(); ++i){
			Ex[0] = pts1[i].x * pE[0] + pts1[i].y * pE[1] + pts1[i].z * pE[2];
			Ex[1] = pts1[i].x * pE[3] + pts1[i].y * pE[4] + pts1[i].z * pE[5];
			Ex[2] = pts1[i].x * pE[6] + pts1[i].y * pE[7] + pts1[i].z * pE[8];
			err += pow(pts2[i].x * Ex[0] + pts2[i].y * Ex[1] + pts2[i].z * Ex[2], 2.0);
		}
		delete[] Ex;
	}
	else{
		const float *pE = E.ptr<float>(0);
		for (int i = 0; i < pts1.size(); ++i){
			Ex[0] = pts1[i].x * pE[0] + pts1[i].y * pE[1] + pts1[i].z * pE[2];
			Ex[1] = pts1[i].x * pE[3] + pts1[i].y * pE[4] + pts1[i].z * pE[5];
			Ex[2] = pts1[i].x * pE[5] + pts1[i].y * pE[6] + pts1[i].z * pE[7];
			err += pts2[i].x * Ex[0] + pts2[i].y * Ex[1] + pts2[i].z * Ex[2];
		}
		delete[] Ex;
	}
	err /= static_cast<double>(pts1.size());
	err = sqrt(err);
	return err;
}


template void calcE(const vector<Point_<double>> &pts1,
	const vector<Point_<double>> &pts2,
	const Mat &Mcam, const double snsr_sz, 
	Mat &E);
template void calcE(const vector<Point_<float>> &pts1,
	const vector<Point_<float>> &pts2,
	const Mat &Mcam, const double snsr_sz,
	Mat &E);
template <typename T>
void calcE(const vector<Point_<T>> &pts1,
	const vector<Point_<T>> &pts2, 
	const Mat &Mcam, const double snsr_sz, 
	Mat &E)
{
	CV_Assert(pts1.size() == pts2.size());
	double flen = (Mcam.at<double>(0, 0) + Mcam.at<double>(1, 1))/2;
	double cx = Mcam.at<double>(0, 2);
	double cy = Mcam.at<double>(1, 2);
	vector<Point3d> pts3D1, pts3D2;
	pts3D1.reserve(pts1.size());
	pts3D2.reserve(pts2.size());
	for (int i = 0; i < pts1.size(); ++i){
		Point3d pt;
		cnv_pix_to_pt3D(pts1[i], snsr_sz, flen, cx, cy, pt);
		pts3D1.push_back(pt);
		cnv_pix_to_pt3D(pts2[i], snsr_sz, flen, cx, cy, pt);
		pts3D2.push_back(pt);
	}
	calcE(pts3D1, pts3D2, E);
	cout << "Eerr: " << 
	checkE(pts3D1, pts3D2, E) << endl;
}

void cnv_pix_to_pt3D(const Point_<double> &pix, 
	const double snsr_sz, const double flen,
	const double cx, const double cy,
	Point3_<double> &pt3D)
{
	pt3D.x = (pix.x - cx)* snsr_sz;
	pt3D.y = (pix.y - cy) * snsr_sz;
	pt3D.z = snsr_sz * flen;
}

//https://en.wikipedia.org/wiki/Essential_matrix
template double calc_depth(const Point_<float> &pt1,
	const Point_<float> &pt2, const Mat &R, const Mat &T);
template double calc_depth(const Point_<double> &pt1,
	const Point_<double> &pt2, const Mat &R, const Mat &T);
template <typename T_>
double calc_depth(const Point_<T_> &pt1, const Point_<T_> &pt2,
const Mat &R, const Mat &T)
{
	CV_Assert(R.type() == T.type());
	double depth;
	double *r1 = new double[3], *r3 = new double[3], *vec = new double[3];

	if (R.type() == CV_64F && T.type() == CV_64F){
		const double * pR = R.ptr<double>(0);
		for (int i = 0; i < R.cols; ++i)
			r1[i] = pR[i];
		pR = R.ptr<double>(2);
		for (int i = 0; i < R.cols; ++i)
			r3[i] = pR[i];	
		for (int i = 0; i < 3; ++i)
			vec[i] = r1[i] - pt2.x * r3[i];
		const double * pT = T.ptr<double>(0);
		double temp1 = 0;
		for (int i = 0; i < T.rows;)
			temp1 += vec[i] * pT[i];
		double temp2 = vec[0] * pt1.x + vec[1] * pt1.y + vec[2];
		depth = temp1/temp2;
	}
	else{
		const float *pR = R.ptr<float>(0);
		for (int i = 0; i < R.cols; ++i)
			r1[i] = pR[i];
		pR = R.ptr<float>(2);
		for (int i = 0; i < R.cols; ++i)
			r3[i] = pR[i];
		for (int i = 0; i < 3; ++i)
			vec[i] = r1[i] - pt2.x * r3[i];
		const float * pT = T.ptr<float>(0);
		double temp1 = 0;
		for (int i = 0; i < T.rows;)
			temp1 += vec[i] * pT[i];
		double temp2 = vec[0] * pt1.x + vec[1] * pt1.y + vec[2];
		depth = temp1 / temp2;
	}

	delete[] r1, r3, vec;
	return depth;
}


template void draw_matches(const Mat &img1,
	const vector<Point_<float>> &pts1, const Mat &img2,
  const vector<Point_<float>> &pts2, Mat &dimg);
template void draw_matches(const Mat &img1,
	const vector<Point_<double>> &pts1, const Mat &img2,
	const vector<Point_<double>> &pts2, Mat &dimg);
template <typename T>
void draw_matches(const Mat &img1, const vector<Point_<T>> &pts1,
	const Mat &img2, const vector<Point_<T>> &pts2, Mat &dimg)
{
	CV_Assert(pts2.size() == pts1.size() &&
		pts2.size() > 0 && !img1.empty() && !img2.empty());
	connect_imgs(img1, img2, dimg);
	if (dimg.type() == CV_8U){
		cvtColor(dimg, dimg, CV_GRAY2BGR);
	}
	RNG rng = theRNG();
	for (int i = 0; i < pts1.size(); ++i){
		Scalar color(rng(256), rng(256), rng(256));
		Point2d start, end;
		start.x = pts1[i].x;
		start.y = pts1[i].y;
		end.x = pts2[i].x + img1.cols;
		end.y = pts2[i].y;
		line(dimg, start, end, color);
		circle(dimg, start, 5, color);
		circle(dimg, end, 5, color);
	}
}

template void read_matches(const char *fname,
vector<Point_<double>> &pts1, vector<Point_<double>> &pts2);
template void read_matches(const char *fname,
	vector<Point_<float>> &pts1, vector<Point_<float>> &pts2);
template <typename T>
void read_matches(const char *fname,
	vector<Point_<T>> &pts1, vector<Point_<T>> &pts2)
{
	CV_Assert(fname != NULL);
	ifstream ifs(fname);
	if(!ifs.is_open())
		exit(EXIT_FAILURE);
	char *buf = new char[1024];
	int sz;
	ifs.getline(buf, 1024);
	sz = atoi(buf);
	pts1.reserve(sz);
	pts2.reserve(sz);
	while (ifs.getline(buf, 1024)){
		//cout << "buf : " << buf << endl;
		Point2d pt1, pt2;
		pt1.x = atof(strtok(buf, " "));
		pt1.y = atof(strtok(NULL, " "));
		pt2.x = atof(strtok(NULL, " "));
		pt2.y = atof(strtok(NULL, " "));
		pts1.push_back(pt1);
		pts2.push_back(pt2);
	}
	//print_state(ifs);
	delete[]buf;
	ifs.close();
}

void print_state(const std::ios& stream) {
	std::cout << " good()=" << stream.good();
	std::cout << " eof()=" << stream.eof();
	std::cout << " fail()=" << stream.fail();
	std::cout << " bad()=" << stream.bad();
}

Camera::Camera() :cap(NULL), Mcam(Mat()), discoeff(Mat()){};

Camera::Camera(VideoCapture &_cap,
 const Mat &_Mcam, const Mat &_discoeff) : cap(_cap),
  Mcam(_Mcam), discoeff(_discoeff){};

EMPTS::EMPTS()
{
	matcher = BFMatcher(NORM_L2, true);
}

void EMPTS::operator()(const Mat &img1, const Mat &img2,
	double threshold = 2.5)
{
	sift(img1, Mat(), keypts1, descriptor1);
	sift(img2, Mat(), keypts2, descriptor2);
	matcher.match(descriptor1, descriptor2, matches);
	double min_dist = DBL_MAX;
	for (int i = 0; i < matches.size(); ++i){
		if (min_dist > matches[i].distance)
			min_dist = matches[i].distance;
	}

	double _threshold = 2.5f * min_dist;
	for (int i = 0; i < matches.size(); ++i){
		if (_threshold > matches[i].distance)
			good_matches.push_back(matches[i]);
	}

	for (int i = 0; i < good_matches.size(); ++i){
		pts1.push_back(keypts1[good_matches[i].queryIdx].pt);
		pts2.push_back(keypts2[good_matches[i].trainIdx].pt);
	}
}

void EMPTS::get_matches(vector<Point2f> &_pts1,
	vector<Point2f> &_pts2) const
{
	_pts1 = pts1;
	_pts2 = pts2;
}

template void cross(const double *a,
	const double *b, double *c);
template void cross(const float *a, 
	const float *b, float *c);
template <typename T>
void cross(const T *a, const T *b,
T *c)
{
	c[0] = a[1] * b[2] - a[2] * b[1];
	c[1] = a[2] * b[0] - a[0] * b[2];
	c[2] = a[0] * b[1] - a[1] * b[0];
}

template void calc_epipoles(const Mat &F,
	Point_<double> &e1, Point_<double> &e2);
template void calc_epipoles(const Mat &F,
	Point_<float> &e1, Point_<float> &e2);
template <typename T>
void calc_epipoles(const Mat &F,
Point_<T> &e1, Point_<T> &e2)
{
	CV_Assert(F.size() == Size(3, 3) &&
		(F.type() == CV_32F || F.type() == CV_64F));
	Point2d pts[2];
	pts[0].x = 0;
	pts[0].y = 0;
	pts[1].x = 1;
	pts[1].y = 1;
	double elines[2][3];
	double _e1[3], _e2[3];
	if (F.type() == CV_32F){
		//float *pF = F.ptr<float>(0);
		float *pF = (float*)F.data;
		for (int i = 0; i < 2; ++i){
			for (int j = 0; j < 3; ++j){
				elines[i][j] = pF[j*3] * pts[i].x +
					pF[j * 3 + 1] * pts[i].y + pF[j*3+2];
			}
		}
		cross(elines[0], elines[1], _e2);
		_e2[0] /= _e2[2];
		_e2[1] /= _e2[2];
		e2.x = static_cast<float>(_e2[0]);
		e2.y = static_cast<float>(_e2[1]);

		Mat Ft = F.t();
		float *pFt = (float*)Ft.data;
		for (int i = 0; i < 2; ++i){
			for (int j = 0; j < 3; ++j){
				elines[i][j] = static_cast<float>(pF[j * 3] * pts[i].x +
					pF[j * 3 + 1] * pts[i].y + pF[j*3+2]);
			}
		}
		cross(elines[0], elines[1], _e1);
		_e1[0] /= _e1[2];
		_e1[1] /= _e1[2];
		e1.x = static_cast<float>(_e1[0]);
		e1.y = static_cast<float>(_e1[1]);
	}
	else{
		double *pF = (double*)F.data;
		for (int i = 0; i < 2; ++i){
			for (int j = 0; j < 3; ++j){
				elines[i][j] = pF[j * 3] * pts[i].x +
					pF[j * 3 + 1] * pts[i].y + pF[j * 3 + 2];
			}
		}
		cross(elines[0], elines[1], _e2);
		_e2[0] /= _e2[2];
		_e2[1] /= _e2[2];
		e2.x = static_cast<float>(_e2[0]);
		e2.y = static_cast<float>(_e2[1]);

		Mat Ft = F.t();
		double *pFt = (double*)Ft.data;
		for (int i = 0; i < 2; ++i){
			for (int j = 0; j < 3; ++j){
				elines[i][j] = pF[j * 3] * pts[i].x +
					pF[j * 3 + 1] * pts[i].y + pF[j * 3 + 2];
			}
		}
		cross(elines[0], elines[1], _e1);
		_e1[0] /= _e1[2];
		_e1[1] /= _e1[2];
		e1.x = static_cast<float>(_e1[0]);
		e1.y = static_cast<float>(_e1[1]);
	}
}

template void shift_pts(const int x, const int y,
	vector<Point_<double>> &pts);
template void shift_pts(const int x, const int y,
	vector<Point_<float>> &pts);
template <typename T>
void shift_pts(const int x, const int y,
	vector<Point_<T>> &pts)
{
	for (int i = 0; i < pts.size(); ++i){
		pts[i].x += x;
		pts[i].y += y;
	}
}

//flen and img_width have same unit.
//baseline belongs to world cordinate.
void get_stereo_properties(const double baseline,
	const double flen, const int img_width,
	double &max_depth, double &min_depth)
{
	max_depth = baseline * flen;
	min_depth = baseline * flen / (img_width - 1);
}

double get_min_obj_sz_per_pix(const double flen, const double depth)
{
	return depth / flen;
}

void makeE(const Mat &R, const Mat &T, Mat &E){
	CV_Assert(IFMFORD(T) && IFMFORD(R) 
		&& R.type() == T.type());
	if (R.type() == CV_64F){
		const double *pT = T.ptr<double>(0);
		Mat Tx = (Mat_<double>(3, 3) 
			<< 0, -pT[2], pT[1],
			pT[2], 0, -pT[0]
			- pT[1], pT[0], 0);
		E = Tx * R;
	}
	else{
		const float *pT = T.ptr<float>(0);
		Mat Tx = (Mat_<float>(3, 3)
			<< 0, -pT[2], pT[1],
			pT[2], 0, -pT[0]
			- pT[1], pT[0], 0);
		E = Tx * R;
	}
}

template void dotm33pt3(const Mat &m,
 const Point3_<double> &src, 	Point3_<double> &dst);
template <typename T>
void dotm33pt3(const Mat &m,
	const Point3_<T> &src, Point3_<T> &dst){
	if (m.type() == CV_64F){
		const double *pm = m.ptr<double>(0);
		dst.x = pm[0] * src.x + pm[1] * src.y + pm[2] * src.z;
		dst.y = pm[3] * src.x + pm[4] * src.y + pm[5] * src.z;
		dst.z = pm[6] * src.x + pm[7] * src.y + pm[8] * src.z;
	}
	else{

	}
}

//flen belongs to pixel cordinates.. 
void rot_trans_img(const Mat &R, const Mat &T, 
	const double sensor_sz, 	const double flen,
	const Mat &src, Mat &dst)
{
	cout << "T in rot_trans" << endl << T << endl;
	dst = Mat::zeros(src.size(), CV_8U);
	int half_w = src.cols/2, half_h = src.rows/2;
	if (	R.type() == CV_64F && T.type() == CV_64F){
		cout << "enteing if" << endl;
		const double _flen = flen * sensor_sz;
//		flen = Mcam.at<double>(0, 0) * sensor_sz;
		const double *pT = T.ptr<double>(0);
		printf("rot_trans : %p\n", pT);
		cout << pT[0] << ", " << pT[1]
			<< ", " << pT[2] << endl;
		for (int i = 0; i < src.rows; ++i){
			const uchar *psrc = src.ptr<uchar>(i);
			for (int j = 0; j < src.cols; ++j){
				Point3d dst_pt;
				Point3d src_pt((j-half_w)*sensor_sz,
				 (i - half_h)*sensor_sz, _flen);
				dotm33pt3(R, src_pt, dst_pt);
				dst_pt.x += pT[0];
				dst_pt.y += pT[1];
				dst_pt.z += pT[2];
				src_pt = dst_pt;
				dst_pt.x = src_pt.x * flen + src_pt.z * half_w;
				dst_pt.y = src_pt.y * flen + src_pt.z * half_h;
				dst_pt.x /= dst_pt.z;
				dst_pt.y /= dst_pt.z;
				/*cout << "dst : " << dst_pt.x << ", "
					<< dst_pt.y << endl;*/
				if (0 <= dst_pt.x && dst_pt.x < src.cols-0.5 
					&& 0 <= dst_pt.y && dst_pt.y < src.rows-0.5){
					dst_pt.x = round(dst_pt.x);
					dst_pt.y = round(dst_pt.y);
					int index = static_cast<int>(dst_pt.x + dst_pt.y*dst.cols);
					dst.data[index] = psrc[j];
				}
			}
		}
	}
}

Mat getRx(const double alpha)
{
	double _alpha = alpha * (CV_PI/180.0);
	Mat Rx = (Mat_<double>(3, 3) <<
		1, 0, 0,
		0, cos(_alpha), -sin(_alpha),
		0, sin(_alpha), cos(_alpha));
	return Rx;
}

Mat getRy(const double alpha)
{
	double _alpha = alpha * (CV_PI/180);
	Mat Ry = (Mat_<double>(3, 3) << 
		cos(_alpha), 0, -sin(_alpha),
		0,           1, 0,
		sin(_alpha), 0, cos(_alpha));
	return Ry;
}

Mat getRz(const double alpha)
{
	double _alpha = alpha * (CV_PI/180);
	cout << _alpha << endl;
	Mat Rz = (Mat_<double>(3, 3) << 
		cos(_alpha), -sin(_alpha), 0,
		sin(_alpha), cos(_alpha), 0,
		0, 0, 1);
	return Rz;
}

void print_mat_propaties(const Mat &m)
{
	CV_Assert(!m.empty());
	printf("address : %p\n", &m);
	printf("data address : %p\n", m.data);
	cout << "total : "<< m.total() << endl;
	Size wholeSize;
	Point ofs;
	cout << "wholeSize : " << wholeSize << endl;
	cout << "offset : " << ofs << endl;
	cout << "colums : " << m.cols << endl;
	cout << "rows : " << m.rows << endl;
	switch (m.type()){
	case CV_8S:
		cout << "type : CV_8S" << endl;
		break;
	case CV_8SC2:
		cout << "type : CV_8SC2" << endl;
		break;
	case CV_8SC3:
		cout << "type : CV_8SC3" << endl;
		break;
	case CV_8SC4:
		cout << "type : CV_8SC4" << endl;
		break;
	case CV_8U:
		cout << "type : CV_8U" << endl;
		break;
	case CV_8UC2:
		cout << "type : CV_8UC2" << endl;
		break;
	case CV_8UC3:
		cout << "type : CV_8UC3" << endl;
		break;
	case	 CV_8UC4:
		cout << "type : CV_8UC4" << endl;
		break;
	case CV_16S:
		cout << "type : CV_16S" << endl;
		break;
	case CV_16SC2:
		cout << "type : CV_16SC2" << endl;
		break;
	case CV_16SC3:
		cout << "type : CV_16SC3" << endl;
		break;
	case CV_16SC4:
		cout << "type : CV_16SC4" << endl;
		break;
	case CV_16U:
		cout << "type : CV_16UC" << endl;
		break;
	case CV_16UC2:
		cout << "type : CV_16UC2" << endl;
		break;
	case CV_16UC3:
		cout << "type : CV_16UC3" << endl;
		break;
	case CV_16UC4:
		cout << "type : CV_16UC4" << endl;
		break;
	case CV_32F:
		cout << "type : CV_32F" << endl;
		break;
	case CV_32FC2:
		cout << "type : CV_32FC" << endl;
		break;
	case CV_32FC3:
		cout << "type : CV_32FC3" << endl;
		break;
	case CV_32FC4:
		cout << "type : CV_32FC4" << endl;
		break;
	case CV_64F:
		cout << "type : CV_64F" << endl;
		break;
	case CV_64FC2:
		cout << "type : CV_64FC2" << endl;
		break;
	case CV_64FC3:
		cout << "type : CV_64FC3" << endl;
		break;
	case CV_64FC4:
		cout << "type : CV_64FC4" << endl;
		break;
	}
}

void set_texture(const GLsizei width,
	const GLsizei height, const GLvoid *data,
	GLuint &texture)
{
	glEnable(GL_TEXTURE_2D);
	glGenTextures(1, &texture);
	glBindTexture(GL_TEXTURE_2D, texture);
	glTexParameterf(GL_TEXTURE_2D,
		GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameterf(GL_TEXTURE_2D,
		GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D, 0, 3,
		width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, data);
	glBindTexture(GL_TEXTURE_2D, 0);
}

void gl_line2D(const GLdouble x1, const GLdouble y1,
	const GLdouble x2, const GLdouble y2, 
	const GLfloat width, const GLdouble *color)
{
	GLdouble cur_color[4];
	glGetDoublev(GL_CURRENT_COLOR, cur_color);
	for (int i =0; i < 4; ++i)
		cout << cur_color[i] << ", ";
	cout << endl;
	glColor4d(color[0], color[1],
	 color[2], color[3]);
	glLineWidth(width);
	glBegin(GL_LINES);
	glVertex2d(x1, y1);
	glVertex2d(x2, y2);
	glEnd();
	glColor4d(cur_color[0], cur_color[1],
		cur_color[2], cur_color[3]);
}

void gl_line3D(const GLdouble x1, const GLdouble y1,
	const GLdouble z1, const GLdouble x2,
	const GLdouble y2, const GLdouble z2,
	const float width)
{
	glLineWidth(width);
	glBegin(GL_LINES);
	glVertex3d(x1, y1, z1);
	glVertex3d(x2, y2, z2);
	glEnd();
}


template void gl_line2D(const Point_<double> &start_pt,
	const Point_<double> &end_pt,
	const float width, const GLdouble *color);
template void gl_line2D(const Point_<float> &start_pt,
	const Point_<float> &end_pt,
	const float width, const GLdouble *color);
template <typename T>
void gl_line2D(const Point_<T> &start_pt,
	const Point_<T> &end_pt,
	const GLfloat width, const GLdouble *color)
{
	gl_line2D(start_pt.x, start_pt.y,
		end_pt.x, end_pt.y, width, color);
}

template void print_data_by_mat(const int cols,
	const int rows, double *data);
template void print_data_by_mat(const int cols,
	const int rows, float *data);
template <typename T>
void print_data_by_mat(const int cols,
 const int rows,  T *data)
{
	for (int i = 0; i < rows; ++i){
		for (int j = 0; j < cols; ++j){
			if (j == cols - 1){
				cout << data[i*rows + j] << ";\n";
			}
			else{
				cout << data[i*rows + j] << ", ";
			}
		}
	}
}

void print_gl_properties()
{
	GLdouble color[4];
	glGetDoublev(GL_CURRENT_COLOR, color);
	cout << "color : ";
	for (int i = 0; i < 4; ++i)
		cout << color[i];
	cout << endl;
	GLboolean enable;
	glGetBooleanv(GL_DEPTH_TEST, &enable);
	cout << "GL_DEPTH_TEST : " 
		<< enable << endl;
	GLfloat mat[16];
	glGetFloatv(GL_MODELVIEW_MATRIX, mat);
	cout << "MODELVIEW MATRIX" << endl;
	print_data_by_mat(4, 4, mat);
	glGetFloatv(GL_PROJECTION_MATRIX, mat);
	cout << "PROJECTION MATRIX" << endl;
	print_data_by_mat(4, 4, mat);
}