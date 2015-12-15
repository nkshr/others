#ifndef GEAR_CVUTIL
#define GEAR_CVUTIL

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/nonfree/nonfree.hpp>
#include <opencv2/features2d/features2d.hpp>
using namespace cv;

void connect_imgs(const Mat &img1, const Mat &img2, Mat &cimg);
void calcF(const vector<Point2f> &pts1,
  const vector<Point2f> &pts2, Mat &F);
void draw_epi_line(const Mat &F, const vector<Point2f> &pts,
  Mat &img);
void draw_epi_line(const Point2f &pt, const Point3f &line,
  Mat &img1, Mat &img2);
template <typename T>
void draw_epi_line(const Point_<T> &pt,
	const Mat &F, const int whichimg,
	Mat &img1, Mat &img2);
template <typename T>
void draw_epi_lines(const vector<Point_<T>> &pts,
	const Mat &F, const int whichimg,
	Mat &img1, Mat &img2);
void draw_pline(Mat &img);
void extractRTfromE(const Mat &E, int flag, Mat &R, Mat &T);
template <typename T>
void cnv_pix_to_pt3D(const Point_<T> &pix, const T snsr_sz,
	const T flen, const Mat &img, Point3_<T> &pt3D);
void cnv_pix_to_pt3D(Point_<double> &pix,
	const double snsr_sz, const double flen,
	const double cx, const double cy,
	Point3_<double> &pt3D);
void decomp_aff(const Mat &src, Mat &R, Mat &T);
template <typename T>
void calcE(const vector<Point3_<T>> &pts1,
	const vector<Point3_<T>> &pts2, Mat &E);
template <typename T>
void calcE(const vector<Point_<T>> &pts1,
	const vector<Point_<T>> &pts2,
	const Mat &Mcam, const double snsr_sz, Mat &E);
double calc_pts_ydiff(vector<Point2f> pts1, vector<Point2f> pts2);
double calc_pts_diff(vector<Point2f> pts1, vector<Point2f> pts2);
void print_pts2D(vector<Point2f> pts);
double checkF(const vector<Point2f> &pts1,
	const vector<Point2f> &pts2, const Mat &F);
template<typename T_>
double calc_depth(const Point_<T_> &pt1, const Point_<T_> &pt2,
const Mat &R, const Mat &T);
template <typename T>
void draw_matches(const Mat &img1, const vector<Point_<T>> &pts1,
	const Mat &img2, const vector<Point_<T>> &pts2, Mat &dimg);
template <typename T>
double checkE(const vector<Point3_<T>> &pts1,
	const vector<Point3_<T>> &pts2, const Mat &E);
template <typename T>
void read_matches(const char *fname,
	vector<Point_<T>> &pts1, vector<Point_<T>> &pts2);
void print_state(const std::ios& stream);

class Camera
{
public:
	VideoCapture cap;
	Mat Mcam, discoeff;
	explicit Camera();
	Camera(VideoCapture &cap,
	 const Mat &Mcam, const Mat &discoeff);
};

class EMPTS
{
private:
	SIFT sift;
	vector<KeyPoint> keypts1, keypts2;
	Mat descriptor1, descriptor2;
	vector<DMatch> matches, good_matches;
	BFMatcher matcher;
	vector<Point2f> pts1, pts2;

public:
	EMPTS();
	void operator()(const Mat &img1,
	 const Mat&img2,
		double threshold);
	void get_matches(vector<Point2f> &pts1,
		vector<Point2f> &pts2) const;
};

template <typename T>
void cross(const T *a, const T *b, T *c);
template <typename T>
void calc_epipoles(const Mat &F,
	Point_<T> &e1, Point_<T> &e2);
template <typename T>
void shift_pts(const int x, const int y,
	vector<Point_<T>> &pts);
void makeE(const Mat &R, const Mat &T,
	Mat &E);
template <typename T>
void dotm33pt3(const Mat &m,
	const Point3_<T> &src, Point3_<T> &dst);
void rot_trans_img(const Mat &R, const Mat &T,
	const double sensor_sz, const double flen,
	const Mat &src, Mat &dst);
Mat getRx(const double alpha);
Mat getRy(const double alpha);
Mat getRz(const double alpha);
void print_mat_propaties(const Mat &m);

#include <GLFW/glfw3.h>
void set_texture(const GLsizei width,
	const GLsizei height, const GLvoid *data, 
	GLuint &texture);
void gl_line2D(const GLdouble x1, const GLdouble y1,
	const GLdouble x2, const GLdouble y2,
	const GLfloat width, const GLdouble *color);
template <typename T>
void gl_line2D(const Point_<T> &start_pt,
	const Point_<T> &end_pt,
	const float width, const GLdouble *color);
template <typename T>
void print_data_by_mat(const int cols, 
	const int rows, T *data);
void print_gl_properties();
#endif