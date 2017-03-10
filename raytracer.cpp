#include "raytracer.h"
#include "scene_types.h"
#include "ray.h"
#include "image.h"
#include "kdtree.h"
#include <stdio.h>
float pi=M_PI;
float IOR=1;
#define ANTIA 2.f
#define TMAX 100000

/// acne_eps is a small constant used to prevent acne when computing intersection
//  or boucing (add this amount to the position before casting a new ray !
const float acne_eps = 1e-4;

bool intersectPlane(Ray *ray, Intersection *intersection, Object *obj) {
	//! \todo : compute intersection of the ray and the plane object
	if(dot(obj->geom.plane.normal,ray->dir)==0)
		return false;
	float t= -(dot(obj->geom.plane.normal,ray->orig)+obj->geom.plane.dist)/dot(obj->geom.plane.normal,ray->dir);
	if(t<=ray->tmin || t>ray->tmax)
		return false;
	ray->tmax=t;
	point3 p=rayAt(*ray,t);
	intersection->normal=(dot(obj->geom.plane.normal,ray->orig-p)>0)?obj->geom.plane.normal:-obj->geom.plane.normal;//mettre le vecteur normal du bon coté
	intersection->mat=&(obj->mat);
	intersection->position=p;
	return true;
}

bool intersectSphere(Ray *ray, Intersection *intersection, Object *obj) {
	//! \todo : compute intersection of the ray and the sphere object
	vec3 co= ray->orig-obj->geom.sphere.center;

	float r=obj->geom.sphere.radius, b=2.0*dot(ray->dir,co), c=dot(co,co)-r*r,
		delta=b*b-4.0*c;
	float t1,t2;
	if(delta<=0)
		return false;
	float t;
	t1=(-b+sqrt(delta))/2.f;
	t2=(-b-sqrt(delta))/2.f;
	if(t1*t2<0) t= t1<t2? t2:t1;
	else t= t1<t2? t1:t2;
	if(t<=ray->tmin || t>ray->tmax)
		return false;
	ray->tmax=t;
	point3 P=rayAt(*ray,t);
	intersection->normal=normalize(P-obj->geom.sphere.center);
	if(t1*t2<0) intersection->normal= -intersection->normal;//flip the normal vector if the camera is inside the sphere
	intersection->position=P+acne_eps*intersection->normal;
	intersection->mat=&(obj->mat);
	return true;
}

bool intersectCylinder(Ray *ray, Intersection *intersection, Object *obj){
	vec3 D = obj->geom.cylinder.direction;
	vec3 C = obj->geom.cylinder.center;
	float r = obj->geom.cylinder.radius;
	vec3 d= ray->dir-dot(ray->dir,D)*D;
	point3 dp= ray->orig-C;
	vec3 cc=dp-dot(dp,D)*D;
	float a= dot(d,d),
		b=2.f*dot(d,cc),
		c= dot(cc,cc)-r*r;
	float delta = b*b -4.f*a*c,
		t1,t2;
	if(delta<0) return false;
	t1=(-b+sqrt(delta))/(2.0*a);
	t2=(-b-sqrt(delta))/(2.0*a);
	float t;
	if(t1*t2<0) t= t1<t2? t2:t1;
	else t= t1<t2? t1:t2;
	if(t<=ray->tmin || t>ray->tmax)
		return false;
	point3 P=rayAt(*ray,t);
	point3 pc=C+(dot(P-C,D)/dot(D,D))*D;
	float l=obj->geom.cylinder.length/2;
	if(length(pc-C)>l){
		Ray rTest; Intersection iTest;
		rayInit(&rTest,ray->orig,ray->dir);
		bool b= intersectPlane(&rTest, &iTest, initPlane(D,-dot(D,C+l*D),obj->mat));
		b= intersectPlane(&rTest, &iTest, initPlane(-D,-dot(-D,C-l*D),obj->mat))||b;
		if(length((C+l*D)-iTest.position)<=r||length((C-l*D)-iTest.position)<=r){
			ray->tmax=length(iTest.position-ray->orig)/length(ray->dir);
			intersection->position=iTest.position;
			intersection->normal=iTest.normal;
			intersection->mat=&(obj->mat);
			return true;
		}
		return false;
	}
	intersection->normal=normalize(P-pc);
	if(t1*t2<0) intersection->normal= -intersection->normal;//flip the normal vector if the camera is inside the cylinder
	intersection->position=P+acne_eps*intersection->normal;
	ray->tmax=t;
	intersection->mat=&(obj->mat);
	return true;
}

float cos2(float f){
	return .5+.5*cos(2.f*f);
}
float sin2(float f){
	return .5-.5*cos(2.f*f);
}

bool intersectCone(Ray *ray, Intersection *intersection, Object *obj){
	float alpha=obj->geom.cone.alpha;
	point3 dp=ray->orig-obj->geom.cone.top;
	vec3 va=normalize(obj->geom.cone.base-obj->geom.cone.top);
	vec3 v1=ray->dir-dot(ray->dir,va)*va,
		 v2=dp-dot(dp,va)*va;
	float a=cos2(alpha)*dot(v1,v1)-sin2(alpha)*dot(ray->dir,va)*dot(ray->dir,va),
		  b=2.f*cos2(alpha)*dot(v1,v2)-2.f*sin2(alpha)*dot(ray->dir,va)*dot(dp,va),
		  c=cos2(alpha)*dot(v2,v2)-sin2(alpha)*dot(dp,va)*dot(dp,va);
	float delta=b*b-4.f*a*c, t,t1,t2;
	if(delta<0) return false;
	t1=(-b+sqrt(delta))/(2.0*a);
	t2=(-b-sqrt(delta))/(2.0*a);
	if(t1*t2<0) t= t1<t2? t2:t1;
	else t= t1<t2? t1:t2;
	if(t<=ray->tmin || t>ray->tmax)
		return false;
	point3 P=rayAt(*ray,t);
	point3 pc=obj->geom.cone.top+(dot(P-obj->geom.cone.top,va)/dot(va,va))*va;
	if(dot(va,pc-obj->geom.cone.top)<=0)
		return false;
	float l = length(obj->geom.cone.top-obj->geom.cone.base);
	if(length(pc-obj->geom.cone.top)>l){
		Ray rTest; Intersection iTest;
		rayInit(&rTest,ray->orig,ray->dir);
		if(intersectPlane(&rTest, &iTest, initPlane(va,-dot(va,obj->geom.cone.base),obj->mat))&&length(obj->geom.cone.base-iTest.position)<=tan(alpha)*l){
			ray->tmax=length(iTest.position-ray->orig)/length(ray->dir);
			intersection->position=iTest.position;
			intersection->normal=iTest.normal;
			intersection->mat=&(obj->mat);
			return true;
		}
		return false;
	}
	intersection->normal=normalize(P-pc);
	intersection->position=P+acne_eps*intersection->normal;
	intersection->mat=&(obj->mat);
	ray->tmax=t;
	return true;
}

bool intersectScene(const Scene *scene, Ray *ray, Intersection *intersection) {
	bool hasIntersection = false;
	int objectCount = scene->objects.size();
	//!\todo loop on each object of the scene to compute intersection
	for(int i=0;i<objectCount;i++){
		switch(scene->objects[i]->geom.type){
			case PLANE:
				hasIntersection|=intersectPlane(ray,intersection,scene->objects[i]);
			break;
			case SPHERE:
				hasIntersection|=intersectSphere(ray,intersection,scene->objects[i]);
			break;
			case CYLINDER:
				hasIntersection|=intersectCylinder(ray,intersection,scene->objects[i]);
			break;
			case CONE:
				hasIntersection|=intersectCone(ray,intersection,scene->objects[i]);
			break;
			default:break;
		}
	}
	return hasIntersection;
}

/* --------------------------------------------------------------------------- */
/*
 *	The following functions are coded from Cook-Torrance bsdf model description and are suitable only
 *  for rough dielectrics material (RDM. Code has been validated with Mitsuba renderer)
 */

float RDM_chiplus(float c){
	return (c > 0.f) ? 1.f : 0.f;
}

/** Normal Distribution Function : Beckmann
 * NdotH : Norm . Half
 */
float RDM_Beckmann(float NdotH, float alpha){
	//! \todo compute Beckmann normal distribution
	float cos2x=NdotH*NdotH;
	float tan2x=(1-cos2x)/cos2x;
	return exp(-tan2x/(alpha*alpha))/(pi*alpha*alpha*cos2x*cos2x);
}

// Fresnel term computation. Implantation of the exact computation. we can use the Schlick approximation
// LdotH : Light . Half
float RDM_Fresnel(float LdotH, float extIOR, float intIOR){
	//! \todo compute Fresnel term
	float f=1-((extIOR/intIOR)*(extIOR/intIOR)*(1-LdotH*LdotH));
	if(f<0)
		return 1.f;
	float t=sqrt(f);
	float f1=extIOR*LdotH-intIOR*t,
	f2=extIOR*LdotH+intIOR*t,
	f3=extIOR*t-intIOR*LdotH,
	f4=extIOR*t+intIOR*LdotH,
	Rs=(f1*f1)/(f2*f2),
	Rp=(f3*f3)/(f4*f4);
	return .5*(Rs+Rp);
}

float RDM_G1(float DdotH, float DdotN, float alpha) {
	//!\todo compute G1 term of the Smith fonction
	float b=1.f/(alpha*(sqrt(1-DdotN*DdotN)/DdotN));
	float k=DdotH/DdotN;
	if(k<=0)
		return 0.f;
	if(b<1.6)
		return 1.f;
	return (3.535*b+ 2.181*b*b)/(1+2.276*b+2.577*b*b);
}

float RDM_Smith(float LdotH, float LdotN, float VdotH, float VdotN, float alpha){
	//!\todo the Smith fonction
	return RDM_G1(LdotH,LdotN,alpha)*RDM_G1(VdotH,VdotN,alpha);
}

// Specular term of the Cook-torrance bsdf
// LdotH : Light . Half
// NdotH : Norm . Half
// VdotH : View . Half
// LdotN : Light . Norm
// VdotN : View . Norm
color3 RDM_bsdf_s(float LdotH, float NdotH, float VdotH, float LdotN, float VdotN, Material *m){
	//!\todo specular term of the bsdf, using D = RDB_Beckmann, F = RDM_Fresnel, G = RDM_Smith
	return m->specularColor*(
		RDM_Beckmann(NdotH,m->roughness)*
		RDM_Fresnel(LdotH,IOR,m->IOR)*
		RDM_Smith(LdotH,LdotN,VdotH,VdotN, m->roughness)
		)/(4.f*LdotN*VdotN);
}
// diffuse term of the cook torrance bsdf
color3 RDM_bsdf_d(Material *m) {
	//!\todo compute diffuse component of the bsdf
	return m->diffuseColor/pi;
}

// The full evaluation of bsdf(wi, wo) * cos (thetai)
// LdotH : Light . Half
// NdotH : Norm . Half
// VdotH : View . Half
// LdotN : Light . Norm
// VdtoN : View . Norm
// compute bsdf * cos(Oi)
color3 RDM_bsdf(float LdotH, float NdotH, float VdotH, float LdotN, float VdotN, Material *m){
	//! \todo compute bsdf diffuse and specular term
	return clamp(RDM_bsdf_d(m)+RDM_bsdf_s(LdotH,NdotH,VdotH,LdotN,VdotN,m),0.f,1.f);
}

/* --------------------------------------------------------------------------- */

color3 shade(vec3 n, vec3 v, vec3 l, color3 lc, Material *mat){
	//! \todo compute bsdf, return the shaded color taking into account the lightcolor
	if(dot(l,n)<0)
		return color3(0.f);
	vec3 h= normalize((v+l)/length(v+l));
	float LdotH=dot(l,h), NdotH=dot(n,h), VdotH=dot(v,h), LdotN=dot(l,n), VdotN=dot(v,n);
	return clamp(lc*RDM_bsdf(LdotH,NdotH,VdotH,LdotN,VdotN,mat)*LdotN,0.f,1.f);
}

//! if tree is not null, use intersectKdTree to compute the intersection instead of intersect scene
color3 trace_ray(Scene * scene, Ray *ray, KdTree *tree) {
	color3 ret = scene->skyColor;
	Intersection intersection;
	if(intersectScene(scene, ray, &intersection)){
		if(intersection.mat->IOR==5040)
			return intersection.normal;
		ret=color3(0.f);
		point3 P=intersection.position;
		for(unsigned int i=0;i<scene->lights.size();i++){
			point3 L= scene->lights[i]->position;
			vec3 l =normalize(L-P);
			Ray rOmbre; Intersection dummyInter;
			rayInit(&rOmbre,P+acne_eps*l,l,0,length(P-L));
			if(!intersectScene(scene,&rOmbre,&dummyInter))
				ret+=shade(intersection.normal,-ray->dir, l, scene->lights[i]->color,intersection.mat);
		}
		if(ray->depth>0){
			vec3 ref=normalize(reflect(ray->dir,intersection.normal));
			rayInit(ray,P+acne_eps*ref,ref,ray->tmin,TMAX,ray->depth-1);
			ret+=RDM_Fresnel(dot(ref,intersection.normal),IOR,intersection.mat->IOR)*trace_ray(scene,ray,tree);
		}
	}
	return clamp(ret,0.f,1.f);
}

float plusOuMoins(){
	float r= (float)rand()/RAND_MAX;
	return r/ANTIA-1.f/(ANTIA*2.f);
}

void renderImage(Image *img, Scene *scene) {
	//! This function is already operational, you might modify it for antialiasing and kdtree initializaion
	float aspect = 1.f/scene->cam.aspect;

	KdTree *tree =  NULL;

	//! \todo initialize KdTree

	float delta_y = 1.f / (img->height * 0.5f); //! one pixel size
	vec3 dy = delta_y * aspect * scene->cam.ydir; //! one pixel step
	vec3 ray_delta_y = (0.5f - img->height * 0.5f) / (img->height * 0.5f) * aspect * scene->cam.ydir;

	float delta_x = 1.f / (img->width * 0.5f);
	vec3 dx = delta_x * scene->cam.xdir;
	vec3 ray_delta_x = (0.5f - img->width * 0.5f) / (img->width * 0.5f) *scene->cam.xdir;

	for(size_t j=0; j<img->height; j++){
		if(j!=0) printf("\033[A\r");
		float progress = (float)j/img->height*100.f;
		printf("progress\t[");
		int cpt = 0;
		for(cpt = 0; cpt<progress; cpt+=5) printf("-");
		for(       ; cpt<100; cpt+=5) printf(" ");
		printf("]\n");
		#pragma omp parallel for
		for(size_t i=0; i<img->width; i++){
			Ray rx;
			color3 c=color3(0);
			for(int x=0;x<ANTIA*ANTIA;x++){
				vec3 ray_dir = scene->cam.center + ray_delta_x + ray_delta_y + float(i+plusOuMoins()+((x%int(ANTIA))-1)/ANTIA)*dx + float(j+plusOuMoins()+((x/int(ANTIA))-1)/ANTIA)*dy;
				rayInit(&rx, scene->cam.position, normalize(ray_dir),0,TMAX,3);
				c+= trace_ray(scene,&rx,tree);
			}
			color3 *ptr = getPixelPtr(img, i,j);
			*ptr = c/(ANTIA*ANTIA);
		}
	}
}
