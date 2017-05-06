#define _USE_MATH_DEFINES

#include <cmath>
#include <iostream>
#include <algorithm>
#include <math.h>
#include <glm/glm.hpp>
#include <X11/Xlib.h>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "../bin/TestModel2.h"

using namespace std;
using glm::vec3;
using glm::mat3;

/* ----------------------------------------------------------------------------*/
/* GLOBAL VARIABLES AND STRUCTS                                                        */

const int SCREEN_WIDTH = 500;
const int SCREEN_HEIGHT = 500;
SDL_Surface* screen;
int t;

vector<Triangle> triangles;

struct Intersection
{
	vec3 position;
	float distance;
	int triangleIndex;
};

vec3 cameraPos( 0.f, 0.f, -4.f );
vec3 ceilingLight( 0, -0.5, -0.8 );
vec3 ceilingLColor = 9.f * vec3( 1, 0.5, 0 );
// vec3 ceilingLColor = 10.f * vec3( 0.3, 0.4, 0.9 );
vec3 secondLight( -0.58, 0.1, -0.40);
// vec3 secondLColor = 9.f * vec3( 1, 1, 1);
vec3 secondLColor = 10.f * vec3( 0.1, 0.3, 0.8);
vec3 thirdLight( 0.58, 0, -0.40);
// vec3 thirdLColor = 9.f * vec3( 1, 1, 1);
vec3 thirdLColor = 10.f * vec3( 0.1, 0.3, 0.8);


// vec3 indirectLight = 0.5f*vec3( 1, 1, 1 );

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

void Update();
void Draw();
void Rotate( mat3& R, float yaw );
void Interpolate( vec3 a, vec3 b, vector<vec3>& result );
vec3 DirectLight( const Intersection& i, vec3 lightPos, vec3 lightColor, float diff, float spec, float m );
vec3 ReflectedLight( vec3& rColor, const Intersection& i, int depth );
bool ClosestIntersection( vec3 start, vec3 dir, Intersection& closestIntersection );
float Fresnel( vec3 L, vec3 normal, float eta1, float eta2 );
float TorranceCook( vec3 normal, vec3 V, vec3 L, float m );

// bool distanceCompare( const Intersection &x, const Intersection &y );

int main( int argc, char* argv[] )
{
	// mat3 R;
	// float yaw = 0.0;
	XInitThreads();

	screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT );
	t = SDL_GetTicks();	// Set start value for timer.

	LoadTestModel( triangles );

	while( NoQuitMessageSDL() )
	{
		Update();
		Draw();
	}

	SDL_SaveBMP( screen, "raytracerTS7.bmp" );
	return 0;
}

void Update()
{
	mat3 R;
	float yaw = 0.0;
	// Compute frame time:
	int t2 = SDL_GetTicks();
	float dt = float(t2-t);
	t = t2;
	// cout << "Render time: " << dt << " ms." << endl;

	Uint8* keystate = SDL_GetKeyState( 0 );
	if( keystate[SDLK_UP] )
	{
		cameraPos.z += 0.3;
	}
	if( keystate[SDLK_DOWN] )
	{
		cameraPos.z -= 0.3;
	}
	if( keystate[SDLK_LEFT] )
	{
		yaw += M_PI/15.f;
		Rotate( R, yaw );
		// cout << "x: " << cameraPos.x << "\ny: " << cameraPos.y<< "\nz: " << cameraPos.z << endl;
	}
	if( keystate[SDLK_RIGHT] )
	{
		yaw -= M_PI/15.f;
		Rotate( R, yaw );
		// cout << "x: " << cameraPos.x << "\ny: " << cameraPos.y<< "\nz: " << cameraPos.z << endl;
	}

}

void Rotate( mat3& R, float yaw )
{
	R[0].x = cos( yaw );
	R[0].z = - sin( yaw );
	R[2].x = sin( yaw );
	R[2].z = cos( yaw );

	cameraPos = R * cameraPos;
}

void Draw()
{
	SDL_FillRect( screen, 0, 0 );

	if( SDL_MUSTLOCK(screen) )
		SDL_LockSurface(screen);

	// Intersection closestIntersection;
	// float focalLength,centerx,centery;
	// vec3 d, color;
	// vec3 tempcol = vec3(0,0,0);
	// bool hasIntercept;
	// focalLength = SCREEN_HEIGHT*1.5;
	// centerx = SCREEN_WIDTH/2;
	// centery = SCREEN_HEIGHT/2;

	#pragma omp parallel for schedule(dynamic)
	for( int n=0; n<SCREEN_HEIGHT*SCREEN_WIDTH; ++n )
	{
		int x = n/SCREEN_WIDTH;
		int y = n%SCREEN_WIDTH;
		vec3 avColor( 0, 0, 0 );
		Intersection closestIntersection;
		float focalLength,centerx,centery;
		focalLength = SCREEN_HEIGHT*1.5;
		centerx = SCREEN_WIDTH/2;
		centery = SCREEN_HEIGHT/2;

		// shoot multiple rays for anti-alias
		for( int s = 0; s < 3; s++ )
		{
			for( int t = 0; t < 3; t++ )
			{
				vec3 d = vec3( x-centerx+0.5*s, y-centery+0.5*t, focalLength );
				// hasIntercept = ClosestIntersection( cameraPos, d, triangles, closestIntersection);
				if( ClosestIntersection( cameraPos, d, closestIntersection ) )
				{
					float diff_val = triangles[closestIntersection.triangleIndex].diffuse;
					float spec_val = triangles[closestIntersection.triangleIndex].specular;
					float m_val = triangles[closestIntersection.triangleIndex].m;
					// vec3 triColor = triangles[closestIntersection.triangleIndex].color;
					vec3 color = triangles[closestIntersection.triangleIndex].color *
						( DirectLight( closestIntersection, ceilingLight, ceilingLColor, diff_val, spec_val, m_val ) +
							DirectLight( closestIntersection, secondLight, secondLColor, diff_val, spec_val, m_val ) +
							DirectLight( closestIntersection, thirdLight, thirdLColor, diff_val, spec_val, m_val )
						);
					avColor += color;
				}
				else
					PutPixelSDL( screen, x, y, vec3(0,0,0) );
			}
		}
		avColor.x = avColor.x/9.f;
		avColor.y = avColor.y/9.f;
		avColor.z = avColor.z/9.f;
		PutPixelSDL( screen, x, y, avColor );
	}


	if( SDL_MUSTLOCK(screen) )
		SDL_UnlockSurface(screen);

	SDL_UpdateRect( screen, 0, 0, 0, 0 );
}

vec3 DirectLight( const Intersection& i, vec3 lightPos, vec3 lightColor, float diff, float spec, float m )
{
	vec3 r, l_pos, l_dir, normal, B, color, diff_color, spec_color;
	float l_dist;
	Intersection nearestObject;

	for( int s = 0; s < 5; s++ )
	{
		for( int t = 0; t < 5; t++ )
		{

			l_pos = vec3( lightPos.x + 0.00005*s, lightPos.y, lightPos.z + 0.00005*t );
			r = l_pos - i.position;
			l_dist = glm::length( r );
			l_dir = glm::normalize( r );
			normal = triangles[i.triangleIndex].normal;

			ClosestIntersection( i.position+l_dir*0.0001f, l_dir, nearestObject );
			if( nearestObject.distance < l_dist )
			{
				color = vec3( 0.0, 0.0, 0.0 );
				return color;
			}

			B = (lightColor / (float)(4*M_PI*pow(l_dist, 2)));
			diff_color += B * std::max( glm::dot( normal, l_dir ), 0.f );

			vec3 V = cameraPos - i.position;
			vec3 L = i.position - l_pos;
			L = glm::normalize( L );
			V = glm::normalize( V );

			// Phong reflection
			// vec3 R = L - 2 * glm::dot( L, normal ) * normal;
			// spec_color += float(pow(std::max( glm::dot( V, R ), 0.f ), 50)) * B;

			spec_color += TorranceCook( normal, V, L, m );
		}
	}
	color = ((diff * diff_color) + (spec * spec_color/float(M_PI)))/25.f;

	return color;
}

float Fresnel( vec3 L, vec3 normal, float eta1, float eta2 )
{
	vec3 normal_n;
	float theta_i, sin_theta_t, cos_theta_t, cos_theta_i, rs, rp, rf;
	theta_i = glm::acos(glm::dot(L,normal_n));
	cos_theta_i = glm::cos(theta_i);
	sin_theta_t = pow(eta1/eta2, 2) * (1 -  pow(cos_theta_i, 2));
	cos_theta_t = sqrt( 1 - sin_theta_t );

	rs = pow( ( ( eta1*cos_theta_i - eta2*cos_theta_t ) / ( eta1*cos_theta_i + eta2*cos_theta_t ) ) , 2 );
	rp = pow( ( ( eta2*cos_theta_i - eta1*cos_theta_t ) / ( eta2*cos_theta_i + eta1*cos_theta_t ) ) , 2 );
	rf = rs + rp / 2;
	// tf = 1 - rf;

	return rf;
}

float TorranceCook( vec3 normal, vec3 V, vec3 L, float m )
{
	// Torrance-Cook BRDF
	normal = glm::normalize(normal);
	float NdotV = glm::dot( normal, V );
	float F = Fresnel(L, normal, 1.f, 1.49f);
	vec3 H = glm::normalize( V + L );
	float NdotH = glm::dot( normal, H );
	float NdotL = glm::dot( normal, L );
	float VdotH = glm::dot( V, H );

	float G = std::min(std::min(1.f, ( 2 * NdotH * NdotV / VdotH )), ( 2 * NdotH * NdotL / VdotH ));

	// float m = 0.4;
	float alpha = glm::acos(NdotH);
	float D = 0.8*exp(-(alpha*alpha)/(m*m));

	return ( F * D * G ) / ( M_PI * NdotL * NdotV );

}

bool ClosestIntersection( vec3 start, vec3 dir, Intersection& closestIntersection )
{
	closestIntersection.distance = numeric_limits<float>::max();
	bool hasIntercept = false;

	for( int s=0; s<triangles.size(); ++s)
	{

		vec3 e1 = triangles[s].v1 - triangles[s].v0;
		vec3 e2 = triangles[s].v2 - triangles[s].v0;
		vec3 b = start - triangles[s].v0;
		mat3 A( -dir, e1, e2 );
		vec3 x = inverse( A ) * b;
		// cout << s << ", ";
		if( (x.x >= 0) && (x.y >= 0) && (x.z >= 0) && ((x.y + x.z) <= 1) )	// x.x = t, x.y = u, x.z = v
		{
			if( closestIntersection.distance > x.x )
			{
				closestIntersection.position = start + x.x * dir;
				closestIntersection.distance = x.x;
				closestIntersection.triangleIndex = s;
				hasIntercept = true;
			}
		}
	}

	return hasIntercept;
}

void Interpolate( vec3 a, vec3 b, vector<vec3>& result)
{
	float size,intervalx,intervaly,intervalz;
	size = result.size()-1;
	intervalx = (b.x - a.x)/size;
	intervaly = (b.y - a.y)/size;
	intervalz = (b.z - a.z)/size;

	for ( int i=0; i<result.size(); ++i)
	{
		result[i].x = a.x + intervalx*i;
		result[i].y = a.y + intervaly*i;
		result[i].z = a.z + intervalz*i;
	}
}
