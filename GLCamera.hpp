#include <GL/freeglut.h>
#include <quaternion.hpp>

class GLCamera
{
private:
	quat qeye, qcenter, qup;
	void update();
public:
	GLdouble eyex, eyey, eyez,
		centerx, centery, centerz,
		upx, upy, upz;
	GLCamera(GLdouble _eyex, GLdouble _eyey, GLdouble _eyez,
		GLdouble _centerx, GLdouble _centery, GLdouble _centerz, 
		GLdouble _upx, GLdouble _upy, GLdouble _upz);
	void rotate_arb(vec3d axisr, double angle);
};

GLCamera::GLCamera(GLdouble _eyex, GLdouble _eyey, GLdouble _eyez,
	GLdouble _centerx, GLdouble _centery, GLdouble _centerz,
	GLdouble _upx, GLdouble _upy, GLdouble _upz) : eyex(_eyex), eyey(_eyey), eyez(_eyez),
	centerx(_centerx), centery(_centery), centerz(_centerz),
	upx(_upx), upy(_upy), upz(_upz)
{
	qeye = quat(0, eyex, eyey, eyez);
	qcenter = quat(0, centerx, centery, centerz);
	qup = quat(0, upx, upy, upz);
}

void GLCamera::rotate_arb(vec3d axisr, double angle)
{
	quat temp, temp2;

	temp2.s = 0;
	temp2.u = qcenter.u - qeye.u;
	temp2.v = qcenter.v - qeye.v;
	temp2.w = qcenter.w - qeye.w;

	rotateq(temp2, &temp, axisr, angle);
	qcenter.s = temp.s;
	qcenter.u = temp.u + qeye.u;
	qcenter.v = temp.v + qeye.v;
	qcenter.w = temp.w + qeye.w;

	rotateq(qup, &temp, axisr, angle);
	qup.s = temp.s;
	qup.u = temp.u;
	qup.v = temp.v;
	qup.w = temp.w;
	update();
}

void GLCamera::update()
{
	eyex = qeye.u;
	eyey = qeye.v;
	eyez = qeye.w;

	centerx = qcenter.u;
	centery = qcenter.v;
	centerz = qcenter.w;

	upx = qup.u;
	upy = qup.v;
	upz = qup.w;
}