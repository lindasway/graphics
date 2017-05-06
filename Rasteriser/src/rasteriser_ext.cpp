#include <iostream>
#include <algorithm>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"

using namespace std;
using glm::vec2;
using glm::vec3;
using glm::mat3;
using glm::ivec2;


/* ----------------------------------------------------------------------------*/
/* GLOBAL VARIABLES AND STRUCTS                                                        */

const int SCREEN_WIDTH = 500;
const int SCREEN_HEIGHT = 500;
SDL_Surface* screen;
int t;

vector<Triangle> triangles;

struct Intersection
{
	glm::vec3 position;
	float distance;
	int triangleIndex;
};

struct Pixel
{
	int x;
	int y;
	float zinv;
	vec3 position3d;
};

struct Vertex
{
	vec3 position;
};

float focal = SCREEN_HEIGHT;
float centerx = SCREEN_WIDTH/2;
float centery = SCREEN_HEIGHT/2;
vec3 cameraPos( 0.f, 0.f, -3.001 );
// mat3 R;
// float yaw = 0;

float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];

vec3 lightPos(0,-0.5,-0.7);
vec3 lightPower = 10.f*vec3( 1, 1, 1 );
vec3 indirectLightPowerPerArea = 0.5f*vec3( 1, 1, 1 );

// vec3 currentNormal;
// vec3 currentReflectance;

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

void Update();
void Draw( float yaw );
// void Rotate( float angle );
void VertexShader( const Vertex& v, Pixel& p );
void PixelShader( const Pixel& p, vec3 currentNormal, vec3 currentReflectance );
void DrawLineSDL( SDL_Surface* surface, ivec2 a, ivec2 b, vec3 currentNormal, vec3 currentReflectance );
// void DrawPolygonEdges( const vector<Vertex>& vertices );
void ComputePolygonRows( const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels );
void DrawRows( const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels, vec3 currentNormal, vec3 currentReflectance );
void DrawPolygon( const vector<Vertex>& vertices, vec3 currentNormal, vec3 currentReflectance );
void Interpolate( Pixel a, Pixel b, vector<Pixel>& result );
// bool ClosestIntersection( vec3 start, vec3 dir, const vector<Triangle>& triangles, Intersection& closestIntersection );
// bool distanceCompare( const Intersection &x, const Intersection &y );

int main( int argc, char* argv[] )
{

	screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT );
	t = SDL_GetTicks();	// Set start value for timer.

	LoadTestModel( triangles );

	// // -- Test ComputePolygonRows --
	// vector<Pixel> vertexPixels(3);
	// vertexPixels[0].x = 10;
	// vertexPixels[0].y = 5;
	// vertexPixels[0].zinv = 0.2;
	// vertexPixels[1].x = 5;
	// vertexPixels[1].y = 10;
	// vertexPixels[1].zinv = 0.3;
	// vertexPixels[2].x = 15;
	// vertexPixels[2].y = 15;
	// vertexPixels[2].zinv = 0.4;
	// vector<Pixel> leftPixels;
	// vector<Pixel> rightPixels;
	// ComputePolygonRows( vertexPixels, leftPixels, rightPixels );
	// for( int row=0; row<leftPixels.size(); ++row )
	// {
	// 	cout << "Start: ("
	// 	<< leftPixels[row].x << ","
	// 	<< leftPixels[row].y << "). "
	// 	<< "End: ("
	// 	<< rightPixels[row].x << ","
	// 	<< rightPixels[row].y << "). " << endl;
	// }
	// //------------------------------

	Uint8* keystate = SDL_GetKeyState( 0 );
	float yaw = 0;

	while( NoQuitMessageSDL() )
	{

		// if ( keystate[SDLK_LEFT] )
		// 	yaw += 0.01f;
		// if (keystate[SDLK_RIGHT])
		// 	yaw -= 0.01f;
		// Update();
		Draw(yaw);
	}

	SDL_SaveBMP( screen, "screenshot.bmp" );
	return 0;
}

void Update()
{
	// Compute frame time:
	int t2 = SDL_GetTicks();
	float dt = float(t2-t);
	t = t2;
	cout << "Render time: " << dt << " ms." << endl;

	Uint8* keystate = SDL_GetKeyState( 0 );
	// if( keystate[SDLK_UP] )
	// {
	// 	vec3 forward(0, 0, 0.01);
	// 	cameraPos = cameraPos + forward;
	// 	//Move camera forward
	// }
	// if( keystate[SDLK_DOWN] )
	// {
	// 	vec3 backwards(0, 0, -0.01);
	// 	cameraPos += backwards;
	// 	//Move camera backward
	// }
	if( keystate[SDLK_w] )
	{
		lightPos.z += 0.2f;
	}
	if( keystate[SDLK_s] )
	{
		lightPos.z -= 0.2f;
	}
	if( keystate[SDLK_a] )
	{
		lightPos.x -= 0.2f;
	}
	if( keystate[SDLK_d] )
	{
		lightPos.x += 0.2f;
	}


}

// void Rotate( float angle )
// {
// 	yaw += angle;
// 	R[0].x = cos( angle );
// 	R[0].z = - sin( angle );
// 	R[2].x = sin( angle );
// 	R[2].z = cos( angle );
//
// 	cameraPos = R * cameraPos;
// }

void Draw( float yaw )
{
	SDL_FillRect( screen, 0, 0 );

	if( SDL_MUSTLOCK(screen) )
		SDL_LockSurface(screen);

	// vec3 color ( 1.0, 1.0, 1.0 );

	// vec3 color( 1, 1, 1 );
	// DrawLineSDL( screen, ivec2( 100, 300 ), ivec2( 350, 100 ), color );
	//
	// vec3 color2( 1, 0.1, 0.3 );
	// DrawLineSDL( screen, ivec2( 400, 250 ), ivec2( 450, 200 ), color2 );


	#pragma omp parallel for
	for( int y=0; y<SCREEN_HEIGHT; ++y )
		for( int x=0; x<SCREEN_WIDTH; ++x )
			depthBuffer[y][x] = 0;

	#pragma omp parallel for
	for( int s=0; s<triangles.size(); s++ )
	{
		vector<Vertex> vertices(3);
		mat3 R(cos(yaw), 0, sin(yaw), 0, 1, 0, -sin(yaw), 0, cos(yaw));

		// vec3 currentColor = triangles[s].color;

		vertices[0].position = triangles[s].v0 * R;
		vertices[1].position = triangles[s].v1 * R;
		vertices[2].position = triangles[s].v2 * R;
		vec3 currentNormal = triangles[s].normal;
		vec3 currentReflectance = triangles[s].color;


		DrawPolygon( vertices, currentNormal, currentReflectance );

		// for( int v=0; v<3; v++ )
		// {
		// 	glm::ivec2 projPos;
		// 	VertexShader( vertices[v], projPos );
		// 	vec3 color( 1, 1, 1 );
		// 	PutPixelSDL( screen, projPos.x, projPos.y, color );
		// }
	}

	if( SDL_MUSTLOCK(screen) )
		SDL_UnlockSurface(screen);

	SDL_UpdateRect( screen, 0, 0, 0, 0 );
}

void VertexShader( const Vertex& v, Pixel& p )
{

	vec3 pos = v.position - cameraPos;
	// cout << "Pos.z is: " << pos.z << endl;
	p.zinv = 1 / pos.z;
	p.x = int(focal*(pos.x/pos.z)) + centerx;
	p.y = int(focal*(pos.y/pos.z)) + centery;
	p.position3d = v.position;
	// cout << "P.position3d is: " << p.position3d.x << ", " << p.position3d.y << ", " << p.position3d.z << endl;

	// cout << "pixel  illumination: " << (p.illumination).x << ", " << (p.illumination).y << ", " << (p.illumination).z << endl;

}

void PixelShader( const Pixel& p, vec3 currentNormal, vec3 currentReflectance )
{
	int x = p.x;
	int y = p.y;
	vec3 diff,spec;
	if( p.zinv > depthBuffer[y][x] )
	{
		depthBuffer[y][x] = p.zinv;

		vec3 r = lightPos - p.position3d;
		float l_dist = glm::length( r );
		vec3 l_dir = glm::normalize( r );
		vec3 B = lightPower / (float)(4*M_PI*pow(l_dist, 2));
		vec3 D = B * std::max( glm::dot( currentNormal, l_dir ), 0.f );
		diff += ( D + indirectLightPowerPerArea ) * currentReflectance;

		// Phong reflection
		// TODO: add Torrance-Cook brdf stuff
		vec3 V = cameraPos - p.position3d;
		vec3 L = p.position3d - lightPos;
		L = glm::normalize( L );
		V = glm::normalize( V );
		vec3 R = L - 2 * glm::dot( L, currentNormal ) * currentNormal;
		spec += float(pow(std::max( glm::dot( V, R ), 0.f ), 10)) * B;

		vec3 color = (0.8f * diff) + (0.2f * spec);

		PutPixelSDL( screen, x, y, color );
	}
}

void DrawLineSDL( SDL_Surface* surface, Pixel a, Pixel b, vec3 currentNormal, vec3 currentReflectance )
{
	ivec2 delta;
	delta.x = glm::abs( a.x - b.x );
	delta.y = glm::abs( a.y - b.y );
	int pixels = glm::max( delta.x, delta.y ) + 1;

	// cout << "pixel a illumination: " << (a.illumination).x << ", " << (a.illumination).y << ", " << (a.illumination).z << endl;

	vector<Pixel> line( pixels );
	Interpolate( a, b, line );

	for( int i=0; i<line.size(); i++ )
	{
		PixelShader( line[i], currentNormal, currentReflectance );
	}
}

// void DrawPolygonEdges( const vector<Vertex>& vertices )
// {
// 	int V = vertices.size();
//
// 	// Transform each vertex from 3D world position to 2D image position:
// 	vector<Pixel> projectedVertices( V );
// 	for( int i=0; i<V; ++i )
// 	{
// 		VertexShader( vertices[i], projectedVertices[i] );
// 	}
//
// 	// Loop over all vertices and draw the edge from it to the next vertex:
// 	for( int i=0; i<V; ++i )
// 	{
// 		int j = (i+1)%V; // The next vertex
// 		vec3 color( 1, 1, 1 );
// 		DrawLineSDL( screen, projectedVertices[i], projectedVertices[j] );
// 	}
//
// }

void ComputePolygonRows( const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels )
{

	// 1. Find max and min y-value of the polygon
	// and compute the number of rows it occupies.
	vector<int> yvals;
	for( int s=0; s<vertexPixels.size(); s++)
	{
		yvals.push_back(vertexPixels[s].y);
	}

	int ymin = *min_element( yvals.begin(), yvals.end());
	int ymax = *max_element( yvals.begin(), yvals.end());

	int nrows = ymax - ymin + 1;

	// 2. Resize leftPixels and rightPixels
	// so that they have an element for each row.

	leftPixels.resize(nrows);
	rightPixels.resize(nrows);

	// 3. Initialize the x-coordinates in leftPixels
	// to some really large value and the x-coordinates
	// in rightPixels to some really small value.

	for( int s=0; s<nrows; s++ )
	{
		leftPixels[s].x = numeric_limits<int>::max();
		rightPixels[s].x = -numeric_limits<int>::max();
	}

	// 4. Loop through all edges of the polygon and use
	// linear interpolation to find the x-coordinate for
	// each row it occupies. Update the corresponding
	// values in rightPixels and leftPixels.

	vector< vector<Pixel> > edgeVals(vertexPixels.size());

	// Interpolate edge values
	int index = 0;
	for( int s=0; s<vertexPixels.size()-1; s++ )
	{
		for( int t=s+1; t<vertexPixels.size(); t++ )
		{
			// cout << "zinv before interpolation: " << vertexPixels[t].zinv << endl;
			edgeVals[index].resize( max( abs( vertexPixels[s].x - vertexPixels[t].x ), abs( vertexPixels[s].y - vertexPixels[t].y ) ) + 1 );
			Interpolate( vertexPixels[s], vertexPixels[t], edgeVals[index] );
			index++;
		}
	}

	// Get values of left and right pixels
	for( int s=0; s<edgeVals.size(); s++ )
	{
		for( int t=0; t<edgeVals[s].size(); t++ )
		{
			// cout << "in loop iter" << endl;
			int i = edgeVals[s][t].y - ymin;
			// cout << "Edge no. " << s << ", Edge val: " << edgeVals[s][t].x << ", " << edgeVals[s][t].y << endl;
			// cout << "i is " << i << endl;
			if( edgeVals[s][t].x < leftPixels[i].x )
			{
				leftPixels[i].x = edgeVals[s][t].x;
				leftPixels[i].y = edgeVals[s][t].y;
				leftPixels[i].zinv = edgeVals[s][t].zinv;
				leftPixels[i].position3d = edgeVals[s][t].position3d;
			}

			if( edgeVals[s][t].x > rightPixels[i].x )
			{
				rightPixels[i].x = edgeVals[s][t].x;
				rightPixels[i].y = edgeVals[s][t].y;
				rightPixels[i].zinv = edgeVals[s][t].zinv;
				rightPixels[i].position3d = edgeVals[s][t].position3d;
			}
		}
	}
}

void DrawRows( const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels, vec3 currentNormal, vec3 currentReflectance )
{
	#pragma omp parallel for
	for( int s=0; s<leftPixels.size(); s++)
	{
		DrawLineSDL( screen, leftPixels[s], rightPixels[s], currentNormal, currentReflectance );
	}

}

void DrawPolygon( const vector<Vertex>& vertices, vec3 currentNormal, vec3 currentReflectance )
{

	int V = vertices.size();
	vector<Pixel> vertexPixels( V );

	#pragma omp parallel for
	for( int i=0; i<V; ++i )
		VertexShader( vertices[i], vertexPixels[i] );

	vector<Pixel> leftPixels;
	vector<Pixel> rightPixels;
	ComputePolygonRows( vertexPixels, leftPixels, rightPixels );
	DrawRows( leftPixels, rightPixels, currentNormal, currentReflectance );

}

void Interpolate( Pixel a, Pixel b, vector<Pixel>& result )
{

	// cout << "a.zinv is: " << a.zinv << endl;
	int N = result.size();
	float step_x = (b.x-a.x) / float(max(N-1,1));
	float step_y = (b.y-a.y) / float(max(N-1,1));
	float step_z = (b.zinv - a.zinv) / float(max(N-1,1));
	// vec3 step_3d = (b.position3d - a.position3d) / float(max(N-1,1));
	float step_3dx = (b.position3d.x*b.zinv - a.position3d.x*a.zinv) / float(max(N-1,1));
	float step_3dy = (b.position3d.y*b.zinv - a.position3d.y*a.zinv) / float(max(N-1,1));
	float step_3dz = (b.position3d.z*b.zinv - a.position3d.z*a.zinv) / float(max(N-1,1));

	float current_x = a.x;
	float current_y = a.y;
	float current_z = a.zinv;
	// vec3 current_3d = a.position3d;
	float current_3dx = a.position3d.x*a.zinv;
	float current_3dy = a.position3d.y*a.zinv;
	float current_3dz = a.position3d.z*a.zinv;


	for( int i=0; i<N; ++i )
	{
		result[i].x = (int) current_x;
		result[i].y = (int) current_y;
		result[i].zinv = current_z;
		result[i].position3d = vec3(current_3dx/current_z, current_3dy/current_z, current_3dz/current_z);

 		// current += step;
		current_x += step_x;
		current_y += step_y;
		current_z += step_z;
		current_3dx += step_3dx;
		current_3dy += step_3dy;
		current_3dz += step_3dz;

		// cout << "result.zinv is: " << result[i].zinv << endl;
	}
}
