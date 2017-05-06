#include <iostream>
#include <vector>
#include <functional>
#include <algorithm>
#include <glm/glm.hpp>
#include <math.h>
#include <SDL.h>
#include "limits.h"
#include "SDLauxiliary.h"
#include "TestModel.h"

using namespace std;
using namespace glm;
using glm::vec3;
using glm::mat3;

/* ----------------------------------------------------------------------------*/
/* GLOBAL VARIABLES                                                            */

const int SCREEN_WIDTH = 500;
const int SCREEN_HEIGHT = 500;
SDL_Surface* screen;
int t;
float focalLength = SCREEN_HEIGHT*1.5;
vec3 cameraPos( 0, 0, -4.001f);
vec3 black( 0, 0, 0 );
vec3 white( 1, 1, 1 );
vec3 lightPos( 0, -0.5, -0.7);
vec3 lightColor  = 14.f * vec3 ( 1, 1, 1 );
vector<Triangle> triangles;
vec3 indirectLight = 0.5f * vec3( 1, 1, 1 );

struct Intersection
{
	vec3 position;
	float distance;
	int triangleIndex;
};

// Intersection closestIntersection;

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

void Update();
void Draw(float yaw);
void Interpolate(vec3 a, vec3 b, vector<vec3>& result);
void LoadTestModel(vector<Triangle>& triangles);
bool ClosestIntersection(vec3 start, vec3 dir, const vector<Triangle>& triangles, Intersection& closestIntersection);
vec3 DirectLight( const Intersection& closestIntersection );

int main( int argc, char* argv[] )
{
	screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT );
	t = SDL_GetTicks();	// Set start value for timer.

	LoadTestModel(triangles);

	float yaw = 0;

	while( NoQuitMessageSDL() )
	{
		Uint8* keystate = SDL_GetKeyState( 0 );

		if (keystate[SDLK_j])
			yaw += -0.1;
		if (keystate[SDLK_l])
			yaw += 0.1;

		Update();
		Draw(yaw);
	}

	SDL_SaveBMP( screen, "screenshot.bmp" );

	return 0;
}

bool ClosestIntersection(vec3 start, vec3 dir, const vector<Triangle>& triangles, Intersection& closestIntersection)
{
	closestIntersection.distance = std::numeric_limits<float>::max();
	bool foundIntersect = false;
	for( size_t i=0; i<triangles.size(); i++ )
	{
		vec3 v0 = triangles[i].v0;
		vec3 v1 = triangles[i].v1;
		vec3 v2 = triangles[i].v2;
		vec3 e1 = v1 - v0;
	 	vec3 e2 = v2 - v0;
		vec3 b = start - v0;

		mat3 A(-dir, e1, e2);
		vec3 x = glm::inverse(A) * b;

		if (x.y >= 0 && x.z >= 0 && (x.y + x.z) <= 1 && x.x >= 0)
		{
			if (x.x < closestIntersection.distance)
			{
				closestIntersection.distance = x.x;
				closestIntersection.position = start + x.x * dir;
				closestIntersection.triangleIndex = i;
				foundIntersect = true;
			}
		}
	}

	return foundIntersect;
}

void Update()
{
	// Compute frame time:
	int t2 = SDL_GetTicks();
	float dt = float(t2-t);
	t = t2;
	//cout << "Render time: " << dt << " ms." << endl;
	//rotation
	Uint8* keystate = SDL_GetKeyState( 0 );

	if( keystate[SDLK_UP] )
	{
		vec3 forward(0, 0, 0.1);
		cameraPos = cameraPos + forward;
		//Move camera forward
	}
	if( keystate[SDLK_DOWN] )
	{
		vec3 backwards(0, 0, -0.1);
		cameraPos = cameraPos + backwards;
		//Move camera backward
	}
	if( keystate[SDLK_LEFT] )
	{
		vec3 left(-0.1, 0, 0);
		cameraPos += left;
	}
	if( keystate[SDLK_RIGHT] )
	{
		vec3 right(0.1, 0, 0);
		cameraPos = cameraPos + right;
		//Move camera to the right
	}
	if (keystate[SDLK_i])
	{
		vec3 up(0, -0.1, 0);
		cameraPos += up;
	}
	if (keystate[SDLK_k])
	{
		vec3 down(0, 0.1, 0);
		cameraPos +=down;
	}
	// move light source ------------------------------------
	if( keystate[SDLK_w] )
	{
		vec3 forward(0, 0, 0.1);
		lightPos = lightPos + forward;
		//Move light forward
	}
	if( keystate[SDLK_s] )
	{
		vec3 backwards(0, 0, -0.1);
		lightPos = lightPos + backwards;
		//Move light backward
	}
	if( keystate[SDLK_a] )
	{
		vec3 left(-0.1, 0, 0);
		lightPos = lightPos + left;
		//Move light to the left
	}
	if( keystate[SDLK_d] )
	{
		vec3 right(0.1, 0, 0);
		lightPos = lightPos + right;
		//Move light to the right
	}
}

vec3 DirectLight(const Intersection& closestIntersection)
{
	vec3 directIllum(0, 0, 0);

	vec3 r = lightPos - closestIntersection.position;
	vec3 n = triangles[closestIntersection.triangleIndex].normal;

	float lightDistance = glm::length(r);
	vec3 lightDirection = glm::normalize(r);

	float Area = 4 * 3.14159265359 * lightDistance * lightDistance;
	vec3 powerPerArea = lightColor / Area;

	float i = (lightDirection.x * n.x) + (lightDirection.y * n.y) + (lightDirection.z * n.z);

	if (i > 0)
		directIllum = powerPerArea * i;

	return directIllum;

}

void Draw(float yaw)
{
	if( SDL_MUSTLOCK(screen) )
		SDL_LockSurface(screen);

	Intersection closestIntersection;
	Intersection lightIntersection;

	for( int y = 0; y < SCREEN_HEIGHT; y++ )
	{
		for( int x=0; x< SCREEN_WIDTH; x++ )
		{
			vec3 d( x - (SCREEN_WIDTH/2), y - (SCREEN_WIDTH/2), focalLength );
			d = d / sqrt(d.x * d.x + d.y * d.y + d.z * d.z);

			mat3 R(cos(yaw), 0, sin(yaw), 0, 1, 0, -sin(yaw), 0, cos(yaw));
			d = d * R;

			if (ClosestIntersection(cameraPos, d, triangles, closestIntersection))
			{
				vec3 position = closestIntersection.position;
				vec3 info = DirectLight(closestIntersection);
				float l_dist = glm::length(lightPos - position);
				vec3 l_dir = glm::normalize(lightPos - position);
				vec3 newColor = triangles[closestIntersection.triangleIndex].color * (info + indirectLight);

				if (ClosestIntersection(position+l_dir * 0.0001f, l_dir, triangles, lightIntersection)) {
					if (lightIntersection.distance < l_dist)
						newColor = black * triangles[closestIntersection.triangleIndex].color;
				}
				PutPixelSDL(screen, x, y, newColor);
			}
			else
			{
				PutPixelSDL(screen, x, y, black);
			}
		}
	}

	if( SDL_MUSTLOCK(screen) )
		SDL_UnlockSurface(screen);

	SDL_UpdateRect( screen, 0, 0, 0, 0 );
}
