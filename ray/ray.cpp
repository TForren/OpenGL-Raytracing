#include <stdio.h>
#include <stdlib.h>

#include "udray.h"
#include "glm.h"
#include <GL/glut.h>
#include "QFHeader.h"

#include <iostream>
using namespace std;
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

extern Camera *ray_cam;       // camera info
extern int image_i, image_j;  // current pixel being shaded
extern bool wrote_image;      // has the last pixel been shaded?

// reflection/refraction recursion control

extern int maxlevel;          // maximum depth of ray recursion 
extern double minweight;      // minimum fractional contribution to color

// these describe the scene

extern vector < GLMmodel * > model_list;
extern vector < Sphere * > sphere_list;
extern vector < Light * > light_list;


//=========================VALUES TO CHANGE======================================
bool enableShadows = true;

bool enableAmbiant = true;
bool enableDiffuse = true;
bool enableSpecular = true;

//Depth visualization
bool depthCue = false;
//==============================================================================

double maxdepth = 1;
double mindepth = -1;


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

// intersect a ray with the entire scene (.obj models + spheres)

// x, y are in pixel coordinates with (0, 0) the upper-left hand corner of the image.
// color variable is result of this function--it carries back info on how to draw the pixel

void trace_ray(int level, double weight, Ray *ray, Vect color)
{
	Intersection *nearest_inter = NULL;
	Intersection *inter = NULL;
	int i;

	// test for intersection with all .obj models

	for (i = 0; i < model_list.size(); i++) {
		inter = intersect_ray_glm_object(ray, model_list[i]);
		update_nearest_intersection(&inter, &nearest_inter);
	}

	// test for intersection with all spheres

	for (i = 0; i < sphere_list.size(); i++) {
		inter = intersect_ray_sphere(ray, sphere_list[i]);
		update_nearest_intersection(&inter, &nearest_inter);
	}

	// "color" the ray according to intersecting surface properties

	// choose one of the simpler options below to debug or preview your scene more quickly.
	// another way to render faster is to decrease the image size.

	if (nearest_inter) {
		//shade_ray_false_color_normal(nearest_inter, color);
		//    shade_ray_intersection_mask(color);  
		shade_ray_diffuse(ray, nearest_inter, color);
		shade_ray_recursive(level, weight, ray, nearest_inter, color);
	}

	// color the ray using a default

	else
		shade_ray_background(ray, color); 
}

//----------------------------------------------------------------------------

// test for ray-sphere intersection; return details of intersection if true

Intersection *intersect_ray_sphere(Ray *ray, Sphere *S)
{

	Intersection *sphereIntersec = make_intersection();
	double discriminant;
	double eps = 0.0000001;

	// a = |E|^2
	// b = 2E * (D - C) 
	// c = |D - C|^2 - R^2

	//a = ray->dir^2
	//b = 2 * (ray->dir dot eye) - 2 * (ray->dir dot sphere->P)
	//c = eye^2 + sphere->P^2 - sphere->radius^2 - 2 * (eye dot sphere->P)

	Vect t4e;
	Vect t4dc;
	double a = VectDotProd(ray->dir, ray->dir);

	double b = 2 * VectDotProd(ray->dir, ray_cam->eye) - (2 * VectDotProd(ray->dir, S->P));

	double c = VectDotProd(ray_cam->eye, ray_cam->eye) + VectDotProd(S->P, S->P) - pow(S->radius, 2) - (2 * VectDotProd(ray_cam->eye, S->P));

	discriminant = (pow(b,2) - 4 * a * c);

	double x1 = (-b + sqrt(b * b - 4 * a * c)) / (2 * a);
	double x2 = (-b - sqrt(b * b - 4 * a * c)) / (2 * a);


	if (discriminant < 0.0){ // no results so return NULL
		return NULL;
	}
	else if (discriminant > -eps && discriminant < eps){ // one result

		sphereIntersec->t = x1; //can use just x1 because if discriminant == 0, then x1 and x2 are the same
	}
	else if (discriminant > 0){
		if (x1 < eps && x2 < eps){ // both values are too close to 0 to care about
			return NULL;
		}
		else if (x1 < eps && x2 > eps) { //x1 is too small to care about but x2 isn't.
			sphereIntersec->t = x2;
		}
		else if (x2 < eps && x1 > eps) { //x2 is too small to care about but x1 isn't.
			sphereIntersec->t = x1;
		}
		else if (x1 > eps && x2 > eps){ // both values are greater than eps 

			sphereIntersec->t = ((x1 < x2) ? x1 : x2);
		}
	}
	Vect p;
	VectAddS(sphereIntersec->t, ray->dir, ray_cam->eye, p);

	Vect nom; 
	VectSub(p, S->P, nom);
	
	VectUnit(nom);

	VectCopy(sphereIntersec->P, p);
	VectCopy(sphereIntersec->N, nom);

	sphereIntersec->surf = S->surf;

	//depth cue stuff
	if (depthCue){

		
	}
	return sphereIntersec;
}

//----------------------------------------------------------------------------

// only local, ambient + diffuse lighting (no specular, shadows, reflections, or refractions)

void shade_ray_diffuse(Ray *ray, Intersection *inter, Vect color)
{
	Vect L;
	double diff_factor;
	Vect t4dl;
	double t4dnl;

	// iterate over lights

	for (int i = 0; i < light_list.size(); i++) {

		if (depthCue){ // in depth view
			double z;
			z = inter->P[2];

			if (z < mindepth){ //clamping
				z = mindepth;
			}
			else if (z > maxdepth){ //clamping
				z = maxdepth;
			}

			color[R] = z;
			color[G] = z;
			color[B] = z;

			//cout << "P: " << inter->P[2] << endl;
		}
		else{ //not in depth visualization

			// AMBIENT
			if (enableAmbiant){
			color[R] += inter->surf->amb[R] * light_list[i]->amb[R];
			color[G] += inter->surf->amb[G] * light_list[i]->amb[G];
			color[B] += inter->surf->amb[B] * light_list[i]->amb[B];
			}
			// DIFFUSE 
			//N = normal
			//Kd = k_diffuse
			//Kd (N dot L) Id

			//calculating L by subtracting our intersection point from our current light pos
			VectSub(light_list[i]->P, inter->P, t4dl);
			//normalize
			VectUnit(t4dl);

			//dot product of N and L 
			t4dnl = VectDotProd(inter->N, t4dl);
		
			//allow depth simulation with light distances (further = darker. Closer = lighter)
			//t4dnl = t4dnl / (VectMag(t4dl) * VectMag(inter->N));
			//t4dnl = t4dnl / 4;

			//clamping 0-1
			t4dnl = (t4dnl > 1.0) ? 1.0 : t4dnl;
			t4dnl = (t4dnl < 0.0) ? 0.0 : t4dnl;

			if (enableDiffuse){
						// Kd						//Id				//(N*L)
			color[R] += inter->surf->diff[R] * light_list[i]->diff[R] * t4dnl;
			color[G] += inter->surf->diff[G] * light_list[i]->diff[G] * t4dnl;
			color[B] += inter->surf->diff[B] * light_list[i]->diff[B] * t4dnl;
			}
		}
	
	}

	// clamp color to [0, 1]
	VectClamp(color, 0, 1);
}
//---------------------------------------------------------------------------

//Checking if ray from surface to light is occluded
bool intersect_ray_shadow(Ray *ray, Sphere *S){

	double discriminant;
	double eps = 1.0e-7;
	double t4diror;

	bool result;

	Vect t4e;
	Vect t4dc;
	
	double a = VectDotProd(ray->dir, ray->dir);

	double b = 2 * VectDotProd(ray->dir, ray->orig) - 2 * VectDotProd(ray->dir, S->P);

	double c = VectDotProd(ray->orig, ray->orig) + VectDotProd(S->P, S->P) - pow(S->radius, 2) - 2 * VectDotProd(ray->orig, S->P);

	discriminant = (pow(b,2) - (4 * a * c));
	//cout << "Discriminant: " << discriminant << endl;

	double x;
	double x1;
	double x2;


	if (discriminant < 0.0){ // no results
		//cout << "Hit Nothing" << endl;
		result = false; //not Occluded
	}
	else if (discriminant > -eps && discriminant < eps){ //one result
		cout << "Hit Edge" << endl;
		x = (-b / (2 * a));
		if (x > 0.0){
			result = true; //Occluded
		}
		else {
			result = false;
		}
	}
	else if (discriminant > 0.0){ // two results
		x1 = (-b + sqrt(b * b - 4 * a * c)) / (2 * a);
		x2 = (-b - sqrt(b * b - 4 * a * c)) / (2 * a);

		if (x1 > 0.0 || x2 > 0.0){
			//cout << "Hit Sphere Center" << endl;
			result = true; //Occluded
		}
		else{
			result = false;
		}
	}
	return result;
}
//----------------------------------------------------------------------------

// same as shade_ray_diffuse(), but add specular lighting + shadow rays (i.e., full Phong illumination model)

void shade_ray_local(Ray *ray, Intersection *inter, Vect color)
	//Ks = inner surface specular
{	//V = unit vector in direction of viewer
	//R = mirior reflectance dir
	//R calculated from incoming light dir and surf norm
	// I_Specular = k_s * I_light ( V dot R )^n_shiney

	Vect t4sl;
	double t4snl;
	double t4vrs;
	double t4eyer;
	double scalar;
	Vect t4n;
	Vect t4sr;
	Vect t4eye;
	Vect t4shadowl;
	
	//shadows stff

	//iterating over all lights
	for (int i = 0; i < light_list.size(); i++) {

		

		//calculating L
		VectSub(light_list[i]->P, inter->P, t4sl);
		VectUnit(t4sl);
		

		//SHADOWS
		//For shadows, tracing rays from our intersection point to our lights. Seeing if occluded.
		Ray *shadowRay = make_ray();
		//setting ray origin to our intersection point
		VectCopy(shadowRay->orig, inter->P);
		//adding eps so it does not intersec with itself
		VectAddS(0.00001, t4sl, inter->P, shadowRay->orig);
		//direction from intersection point to light
		VectSub(light_list[i]->P, shadowRay->orig, shadowRay->dir);
		VectUnit(shadowRay->dir);

		//SPEC
		
		//making copy of ray cam eye
		VectCopy(t4eye, ray_cam->eye);
		VectUnit(t4eye);

		//calculating N dot L * 2
		t4snl = (2 * (VectDotProd(inter->N, t4sl)));

		//calculating 2(N dot L) * N - L putting result in t4sr
		VectNegate(t4sl, t4sl);
		VectAddS(t4snl, inter->N, t4sl, t4sr); 
		VectUnit(t4sr);

		//calculating (V dot R)^shiney
		t4eyer = VectDotProd(t4eye, t4sr);
		t4vrs = pow(t4eyer, inter->surf->spec_exp);

		// Kd						//Id							//(V dot R) ^shiney
		if ((t4eyer > 0) && enableSpecular){
			color[R] += inter->surf->spec[R] * light_list[i]->spec[R] * t4vrs;
			color[G] += inter->surf->spec[G] * light_list[i]->spec[G] * t4vrs;
			color[B] += inter->surf->spec[B] * light_list[i]->spec[B] * t4vrs;
		}
		else {
			//do nothing or else the backsides will be specular highlighted
		}



		// test for intersection with all spheres

	for (i = 0; i < sphere_list.size(); i++) {
			if (intersect_ray_shadow(shadowRay, sphere_list[i]) && enableShadows){// we hit something along the way to the light
				color[R] = 0;
				color[G] = 0;
				color[B] = 0;
				//cout << "Hit something" << endl;
			}
			else{
				//don't shadow point
			}
	
		}

	// test for intersection with all .obj models

	for (int m = 0; m < model_list.size(); m++) {
		if ((intersect_ray_glm_object(shadowRay, model_list[m]) != NULL) && enableShadows){ //we hit something along the way to the light
			color[R] = 0;
			color[G] = 0;
			color[B] = 0;
		}
	}

	}
	VectClamp(color, 0, 1);
}


//----------------------------------------------------------------------------

// full shading model: ambient/diffuse/specular lighting, shadow rays, recursion for reflection, refraction

// level = recursion level (only used for reflection/refraction)

void shade_ray_recursive(int level, double weight, Ray *ray, Intersection *inter, Vect color)
{
	Surface *surf;
	int i;

	// initialize color to Phong reflectance model

	shade_ray_local(ray, inter, color);

	// if not too deep, recurse

	if (level + 1 < maxlevel) {

		// add reflection component to color

		if (surf->reflectivity * weight > minweight) {

			// FILL IN CODE

		}

		// add refraction component to color

		if (surf->transparency * weight > minweight) {

			// GRAD STUDENTS -- FILL IN CODE

		}
	}
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

// ray trace another pixel if the image isn't finished yet

void idle()
{
	if (image_j < ray_cam->im->h) {

		raytrace_one_pixel(image_i, image_j);

		image_i++;

		if (image_i == ray_cam->im->w) {
			image_i = 0;
			image_j++;
		}    
	}

	// write rendered image to file when done

	else if (!wrote_image) {

		write_PPM("output.ppm", ray_cam->im);

		wrote_image = true;
	}

	glutPostRedisplay();
}

//----------------------------------------------------------------------------

// show the image so far

void display(void)
{
	// draw it!

	glPixelZoom(1, -1);
	glRasterPos2i(0, ray_cam->im->h);

	glDrawPixels(ray_cam->im->w, ray_cam->im->h, GL_RGBA, GL_FLOAT, ray_cam->im->data);

	glFlush ();
}

//----------------------------------------------------------------------------

void init()
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0, ray_cam->im->w, 0.0, ray_cam->im->h);
}

//----------------------------------------------------------------------------

int main(int argc, char** argv)
{
	//QuadFormula qf;
	//qf.quadFormula(9, 12, 4);
	//cout << "x: " << qf.x << " x1: " << qf.x1 << " x2: " << qf.x2 << endl;
	
	glutInit(&argc, argv);

	// initialize scene (must be done before scene file is parsed)

	init_raytracing();

	if (argc == 2)
		parse_scene_file(argv[1], ray_cam);
	else {
		printf("missing .scene file\n");
		exit(1);
	}

	// opengl business
	
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(ray_cam->im->w, ray_cam->im->h);
	glutInitWindowPosition(500, 300);
	glutCreateWindow("hw3");
	init();

	glutDisplayFunc(display); 
	glutIdleFunc(idle);

	glutMainLoop();

	return 0; 
}

//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
