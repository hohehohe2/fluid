#include <stdlib.h>
#include <GL/glut.h>
#include "ExampleSphApp.h"

using namespace hohehohe2;

ExampleSphApp esa;

float  camera_x = 0.0;
float  camera_y = 0.0;
float  camera_z = -30.0;
float  camera_pitch_x = 0.0;
float  camera_pitch_y = 0.0;

int drag_mouse_l = 0;
int zoom_mouse_l = 0;
int move_mouse_l = 0;
int last_mouse_x;
int last_mouse_y;

void  display( void )
{
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity();
	glTranslatef(camera_x, camera_y, camera_z);
	glRotatef(-camera_pitch_y, 1.0, 0.0, 0.0);
	glRotatef(-camera_pitch_x, 0.0, 1.0, 0.0);

	//float  light0_position[] = { 10.0, 10.0, 10.0, 1.0 };
	//glLightfv( GL_LIGHT0, GL_POSITION, light0_position );

	glEnable(GL_DEPTH_TEST);
	esa.draw();
	glDisable(GL_DEPTH_TEST);

	glutSwapBuffers();
}


void  reshape( int w, int h )
{
	glViewport(0, 0, w, h);

	glMatrixMode( GL_PROJECTION );
	glLoadIdentity();
	gluPerspective( 45, (double)w/h, 1, 500 );
}


void  mouse( int button, int state, int mx, int my )
{
	if ( ( button == GLUT_LEFT_BUTTON ) && ( state == GLUT_DOWN ) )
		drag_mouse_l = 1;
	else if ( ( button == GLUT_LEFT_BUTTON ) && ( state == GLUT_UP ) )
		drag_mouse_l = 0;

	if ( ( button == GLUT_MIDDLE_BUTTON ) && ( state == GLUT_DOWN ) )
		move_mouse_l = 1;
	else if ( ( button == GLUT_MIDDLE_BUTTON ) && ( state == GLUT_UP ) )
		move_mouse_l = 0;

	if ( ( button == GLUT_RIGHT_BUTTON ) && ( state == GLUT_DOWN ) )
		zoom_mouse_l = 1;
	else if ( ( button == GLUT_RIGHT_BUTTON ) && ( state == GLUT_UP ) )
		zoom_mouse_l = 0;

	last_mouse_x = mx;
	last_mouse_y = my;
}

void keyboard(unsigned char key, int x, int y)
{
	esa.onKey(key);
	display();
}

void  motion( int mx, int my )
{
	if ( drag_mouse_l == 1 )
	{
		camera_pitch_y -= ( my - last_mouse_y );
		if ( camera_pitch_y < -90.0 )
			camera_pitch_y = -90.0;

		camera_pitch_x -= ( mx - last_mouse_x );
		if ( camera_pitch_x < -90.0 )
			camera_pitch_x = -90.0;
	}

	if ( move_mouse_l == 1 )
	{
		camera_y -= ( my - last_mouse_y ) * 0.1f;
		camera_x += ( mx - last_mouse_x ) * 0.1f;
	}

	if ( zoom_mouse_l == 1 )
	{
		camera_z += my - last_mouse_y;
	}

	last_mouse_x = mx;
	last_mouse_y = my;

	glutPostRedisplay();
}


void  idle( void )
{
}


//void  initEnvironment( void )
//{
//	float  light0_position[] = { 10.0, 10.0, 10.0, 1.0 };
//	float  light0_diffuse[] = { 0.8, 0.8, 0.8, 1.0 };
//	float  light0_specular[] = { 1.0, 1.0, 1.0, 1.0 };
//	float  light0_ambient[] = { 0.1, 0.1, 0.1, 1.0 };
//	glLightfv( GL_LIGHT0, GL_POSITION, light0_position );
//	glLightfv( GL_LIGHT0, GL_DIFFUSE, light0_diffuse );
//	glLightfv( GL_LIGHT0, GL_SPECULAR, light0_specular );
//	glLightfv( GL_LIGHT0, GL_AMBIENT, light0_ambient );
//	glEnable( GL_LIGHT0 );
//
//	glEnable( GL_LIGHTING );
//
//	glEnable( GL_COLOR_MATERIAL );
//
//	glEnable( GL_DEPTH_TEST );
//
//	glClearColor( 1, 1, 1, 0.0 );
//}


int main( int argc, char ** argv )
{
	glutInit( &argc, argv );
	glutInitDisplayMode( GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH );
	glutInitWindowSize( 640, 640 );
	glutInitWindowPosition( 0, 0 );
	glutCreateWindow( "OpenGL & GLUT sample program" );

	glutDisplayFunc( display );
	glutReshapeFunc( reshape );
	glutMouseFunc( mouse );
	glutKeyboardFunc( keyboard );
	glutMotionFunc( motion );
	glutIdleFunc( idle );

	//initEnvironment();

	esa.reset();

	glutMainLoop();

	return 0;
}
