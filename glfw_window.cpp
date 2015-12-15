#pragma comment(lib, "glfw3.lib")
#pragma comment(lib, "glu32.lib")
#pragma comment(lib, "opengl32.lib")
#pragma comment(lib, "gear_cvutil")

#include <GLFW/glfw3.h>

#include <cvlib.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
using namespace cv;

#include <iostream>
using namespace std;

#include <gear_cvutil.h>

namespace
{
	const char *img_name1 = "uimg1.png";
	const char *img_name2 = "uimg2.png";
	const unsigned int width = 640;
	const unsigned int height = 480;
	GLdouble left_val, right_val,
	 bottom_val, top_val, centerx, centery;
	vector<Point2d> pts1, pts2;
	bool found = false;
	int counter = 0;
}

static void init_win(const Mat &img1, const Mat &img2,
	GLuint &texture1, GLuint &texture2);
static void mb_cb(GLFWwindow* window, int button,
	int action, int mods);
static void scr_cb(GLFWwindow* window, double xoffset, double yoffset);
static void key_cb(GLFWwindow* window, int key, int scancode,
 int action, int mods);
static void cp_cb(GLFWwindow* window, double xpos, double ypos);

int main(int argc, char **argv)
{
	GLFWwindow* window;

	if (!glfwInit())
		return -1;

	window = glfwCreateWindow(width, height, "Hello World", NULL, NULL);
	if (!window)
	{
		glfwTerminate();
		return -1;
	}

	glfwMakeContextCurrent(window);
	Mat img1 = imread(img_name1);
	Mat img2 = imread(img_name2);
	if (img1.empty() || img2.empty()){
		cerr << "Error: one of images is empty." << endl;
		return -1;
	}
	Mat tmp;
	resize(img1, img1, Size(640, 480));
	resize(img2, img2, Size(640, 480));
	GLuint texture1, texture2;
	init_win(img1, img2, texture1, texture2);
	glfwSetMouseButtonCallback(window, mb_cb);
	glfwSetScrollCallback(window, scr_cb);
	glfwSetKeyCallback(window, key_cb);
	glfwSetCursorPosCallback(window, cp_cb);
	left_val = 0.0;
	right_val = width;
	bottom_val = 0.0;
	top_val = height;
	centerx = 0.0;
	centery = 0.0;

	while (!glfwWindowShouldClose(window))
	{
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glOrtho(left_val, right_val, bottom_val, top_val, -1.0, 1.0);
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, texture1);
		glEnable(GL_ALPHA_TEST);
		glBegin(GL_POLYGON);
		glTexCoord2f(0.0f, 0.0f); glVertex2d(0, img1.rows);
		glTexCoord2f(0.0f, 1.0f); glVertex2d(0, 0);
		glTexCoord2f(1.0f, 1.0f); glVertex2d(img1.cols, 0);
		glTexCoord2f(1.0f, 0.0f); glVertex2d(img1.cols, img1.rows);
		glEnd();
		glDisable(GL_ALPHA_TEST);
		glDisable(GL_TEXTURE_2D);
		double xpos, ypos;
		glfwGetCursorPos(window, &xpos, &ypos);
		if (found){
			GLdouble color[4] = {1.0, 0.0, 0.0, 1.0};
			Point2d start_pt(pts1[counter].x,
				height - pts1[counter].y), end_pt(xpos, height- ypos);
			gl_line2D(start_pt, end_pt, 10.f, color);
		}
		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	glfwTerminate();
	return 0;
}

static void init_win(const Mat &img1, const Mat &img2,
	GLuint &texture1, GLuint &texture2)
{
	glClearColor(1.0, 0.0, 0.0, 1.0);
	GLvoid *data = reinterpret_cast<void*>(img1.data);
	set_texture((GLsizei)img1.cols, (GLsizei)img1.rows,
		data, texture1);
	data = reinterpret_cast<GLvoid*>(img2.data);
	set_texture((GLsizei)img2.cols, (GLsizei)img2.rows,
		data, texture2);
}

static void mb_cb(GLFWwindow* window, int button,
	int action, int mods)
{
	static int pre_action = GLFW_RELEASE;
	double xpos, ypos;
	glfwGetCursorPos(window, &xpos, &ypos);
	if (button == GLFW_MOUSE_BUTTON_LEFT
		&& action == GLFW_RELEASE 
		&& pre_action == GLFW_PRESS){
		cout << "x, y : " << xpos 
			<< ", "  << ypos << endl;
		if (found == false){
			found = true;
			pts1.push_back(Point2d(xpos, ypos));
		}
		else{
			pts2.push_back(Point2d(xpos, ypos));
			found = false;
			counter++;
		}
	}
	else if(button == GLFW_MOUSE_BUTTON_LEFT
		&& action == GLFW_PRESS){
		pre_action = GLFW_PRESS;
	}
}

static void scr_cb(GLFWwindow* window,
	double xoffset, double yoffset)
{
	static double offset = 0;
	const double scale = 3.0;
	double centerx = right_val - left_val;
	double centery = bottom_val - top_val;
	
	yoffset *= scale;
	left_val += yoffset;
	right_val -= yoffset;
	bottom_val += yoffset;
	top_val -= yoffset;
	cout << "offset : " << offset << endl;
	cout << "left : " << left_val << endl;
	cout << "right : " << right_val << endl;
	cout << "bottom : "<< bottom_val << endl;
	cout << "top : " << top_val << endl;
}

static void key_cb(GLFWwindow* window,
	int key, int scancode, int action, int mods)
{
	const double scale = 3.0;
	if (action == GLFW_PRESS){
		switch (key){
		case GLFW_KEY_UP:
			bottom_val -= scale;
			top_val -= scale;
			break;
		case GLFW_KEY_RIGHT:
			left_val += scale;
			right_val += scale;
			break;
		case GLFW_KEY_DOWN:
			bottom_val += scale;
			top_val += scale;
			break;
		case GLFW_KEY_LEFT:
			left_val -= scale;
			right_val -= scale;
			break;
		}
	}
}


static void cp_cb(GLFWwindow* window, double xpos,
	double ypos)
{
	if (found){
		cout << pts1[counter].x << ", " 
			<< pts1[counter].y << endl;
	}
}
