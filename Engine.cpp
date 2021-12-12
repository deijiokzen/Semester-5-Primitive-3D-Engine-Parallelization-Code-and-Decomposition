#include "olcConsoleGameEngine.h"
#include <omp.h>
#include <fstream>
#include<string>
using namespace std;


struct vec3d
{
	float coordinates[3];
};

struct triangle
{
	vec3d p[3];
};

struct mesh
{
	vector<triangle> tris;
};

struct matrixarray
{
	float m[4][4] = { 0 };
};

class Engine3D : public olcConsoleGameEngine
{
public:
	Engine3D()
	{
		m_sAppName = L"Saud/Umer/Zaid 3D ENGINE";
	}


private:
	mesh CubeMesh;
	matrixarray ProjectionMatrix;

	float fTheta, w;

	void MultiplyMatrixVector(vec3d& b, vec3d& o, matrixarray& m)
	{
		std::ofstream out("MultiplicationDecomposition.txt", ios::out | ios::app);
		string input;
		float keep = 0;

		for (int i = 0; i < 4; i++)
		{
			keep = 0;
			float vec3dvalue = 0;
#pragma omp parallel for reduction(+:keep) firstprivate(i) num_threads(3) 
			for (int j = 0; j < 3; j++)
			{
				//input = "\nMultiplication Thread #"+ to_string(omp_get_thread_num()) +"\nVector Value = ";
				if (j == 0)
				{
					vec3dvalue = b.coordinates[0] * m.m[j][i];
					//input += to_string(vec3dvalue) + "\n";
				}
				if (j == 1)
				{
					vec3dvalue = b.coordinates[1] * m.m[j][i];
					//input += to_string(vec3dvalue) + "\n";
				}
				if (j == 2)
				{
					vec3dvalue = b.coordinates[2] * m.m[j][i];
					//input += to_string(vec3dvalue) + "\n";
				}
				keep += vec3dvalue;
				//input += "Individual local multiplication value: " + to_string(keep) + "\n";

			}
			keep += m.m[3][i];
			//input += "Final multiplication value: " + to_string(keep) + "\n";
			//out << input;
			if (i == 0)
			{
				o.coordinates[0] = keep;
			}
			if (i == 1)
			{
				o.coordinates[1] = keep;
			}
			if (i == 2)
			{
				o.coordinates[2] = keep;
			}
			if (i == 3)
			{
				w = keep;
			}
		}/*
		o.x = b.x * m.m[0][0] + b.y * m.m[1][0] + b.z * m.m[2][0] + m.m[3][0];
		o.y = b.x * m.m[0][1] + b.y * m.m[1][1] + b.z * m.m[2][1] + m.m[3][1];
		o.z = b.x * m.m[0][2] + b.y * m.m[1][2] + b.z * m.m[2][2] + m.m[3][2];
		w = b.x * m.m[0][3] + b.y * m.m[1][3] + b.z * m.m[2][3] + m.m[3][3];
*/
		if (w != 0.0f)
		{
			o.coordinates[0] /= w; o.coordinates[1] /= w; o.coordinates[2] /= w;
		}
	}

	void ProjectionMatrix_(float& fNear, float& fFar, float& fFov, float& fAspectRatio, float& fFovRad)
	{
		fNear = 0.1f;
		fFar = 1000.0f;
		fFov = 90.0f;
		fAspectRatio = (float)ScreenHeight() / (float)ScreenWidth();
		fFovRad = 0.5f / tanf(fFov * 0.5f / 180.0f * 3.14159f);
	}
	void ZRotate(matrixarray& m)
	{
		m.m[0][0] = cosf(fTheta);
		m.m[0][1] = sinf(fTheta);
		m.m[1][0] = -sinf(fTheta);
		m.m[1][1] = cosf(fTheta);
		m.m[2][2] = 1;
		m.m[3][3] = 1;
	}
	void YRotate(matrixarray& m)
	{
		m.m[0][0] = cosf(fTheta * 0.5f);
		m.m[1][1] = 1;
		m.m[0][2] = sinf(fTheta * 0.5f);
		m.m[2][0] = -sinf(fTheta * 0.5f);
		m.m[2][2] = cosf(fTheta * 0.5f);
		m.m[3][3] = 1;
	}

	void XRotate(matrixarray& m)
	{
		m.m[0][0] = 1;
		m.m[1][1] = cosf(fTheta * 0.5f);
		m.m[1][2] = sinf(fTheta * 0.5f);
		m.m[2][1] = -sinf(fTheta * 0.5f);
		m.m[2][2] = cosf(fTheta * 0.5f);
		m.m[3][3] = 1;
	}

public:
	bool OnUserCreate() override
	{
		//string all_faces[6] = { "SOUTH", "EAST", "NORTH", "WEST", "TOP", "BOTTOM" };
		std::fstream myfile("CubeCoordinates.txt", std::ios_base::in);
		int i = 0, j = 0, k = 0;
		float a;
		triangle b;
		while (myfile >> a)
		{
			k++; //if (i == 1 && k <72) <36 <52 for open rectangle calculations
			if (i == 1 && k < 72)
			{
				b.p[j].coordinates[i] = a * 2;
			}
			else/* FOR RECTANGLE */
				b.p[j].coordinates[i] = a;

			i++;

			if (i == 3)
			{
				i = 0;
				j++;
			}

			if (j == 3)
			{
				i = 0, j = 0;
				CubeMesh.tris.push_back(b);
			}
		}


		float fNear, fFar, fFov, fAspectRatio, fFovRad;
		// Projection Matrix
		ProjectionMatrix_(fNear, fFar, fFov, fAspectRatio, fFovRad);
		/*
		static void my_PerspectiveFOV(double fov, double aspect, double near, double far, double* mret)
		{
			double D2R = M_PI / 180.0;
			double yScale = 1.0 / tan(D2R * fov / 2);
			double xScale = yScale / aspect;
			double nearmfar = near - far;
			double m[] = {
			xScale, 0, 0, 0,
			0, yScale, 0, 0,
			0, 0, (far + near) / nearmfar, -1,
			0, 0, 2*far*near / nearmfar, 0
			};
			memcpy(mret, m, sizeof(double)*16);
		}
		*/

		ProjectionMatrix.m[0][0] = fAspectRatio * fFovRad;
		ProjectionMatrix.m[1][1] = fFovRad;
		ProjectionMatrix.m[2][2] = fFar / (fFar - fNear);
		ProjectionMatrix.m[3][2] = (-fFar * fNear) / (fFar - fNear);
		ProjectionMatrix.m[2][3] = 1.0f;
		ProjectionMatrix.m[3][3] = 0.0f;

		return true;
	}

	bool OnUserUpdate(float fElapsedTime) override
	{
		// Clear Screen
		Fill(0, 0, ScreenWidth(), ScreenHeight(), PIXEL_SOLID, FG_DARK_CYAN);

		// Set up rotation matrices
		matrixarray matrixRotationZ, matrixRotationX;
		fTheta += 1.0f * fElapsedTime * 3;

		// Rotation Z
		ZRotate(matrixRotationZ);


		// Rotation X
		YRotate(matrixRotationX);

		std::ofstream out("IntermediateDecomposition.txt");


#pragma omp parallel for num_threads(CubeMesh.tris.size())
		for (int i = 0; i < CubeMesh.tris.size(); i++)
		{
			triangle Triangle_Proj, Triangle_Translated, Triangle_RotatedZ, Triangle_RotatedZX;


			// Rotate in Z-Axis
			MultiplyMatrixVector(CubeMesh.tris[i].p[0], Triangle_RotatedZ.p[0], matrixRotationZ);
			MultiplyMatrixVector(CubeMesh.tris[i].p[1], Triangle_RotatedZ.p[1], matrixRotationZ);
			MultiplyMatrixVector(CubeMesh.tris[i].p[2], Triangle_RotatedZ.p[2], matrixRotationZ);
			string input = "Mesh Thread#" + to_string(omp_get_thread_num()) +"\nRotated Point X: "+ to_string(CubeMesh.tris[i].p[0].coordinates[0])+ 
				", Point Y: " + to_string(CubeMesh.tris[i].p[0].coordinates[1]) +  
				", Point Z: " + to_string(CubeMesh.tris[i].p[0].coordinates[2]) + 
				"--> Point X: " + to_string(Triangle_RotatedZ.p[0].coordinates[0]) +
                ",Point Y: " + to_string(Triangle_RotatedZ.p[0].coordinates[1]) +  
			    ",Point Z: " + to_string(Triangle_RotatedZ.p[0].coordinates[1]) +" About the Z Axis.\n";
			out << input;
			// Rotate in X-Axis
			MultiplyMatrixVector(Triangle_RotatedZ.p[0], Triangle_RotatedZX.p[0], matrixRotationX);
			MultiplyMatrixVector(Triangle_RotatedZ.p[1], Triangle_RotatedZX.p[1], matrixRotationX);
			MultiplyMatrixVector(Triangle_RotatedZ.p[2], Triangle_RotatedZX.p[2], matrixRotationX);
			input = "Mesh Thread#" + to_string(omp_get_thread_num()) + "\nRotated Point X: " + to_string(Triangle_RotatedZ.p[0].coordinates[0]) +
				", Point Y: " + to_string(Triangle_RotatedZ.p[0].coordinates[1]) +
				", Point Z: " + to_string(Triangle_RotatedZ.p[0].coordinates[2]) +
				"--> Point X: " + to_string(Triangle_RotatedZX.p[0].coordinates[0]) +
				",Point Y: " + to_string(Triangle_RotatedZX.p[0].coordinates[1]) +
				",Point Z: " + to_string(Triangle_RotatedZX.p[0].coordinates[1]) + " About the ZX Axis.\n";
			out << input;

			// Offset into the screen
			Triangle_Translated = Triangle_RotatedZX;
			Triangle_Translated.p[0].coordinates[2] = Triangle_RotatedZX.p[0].coordinates[2] + 3.0f;
			Triangle_Translated.p[1].coordinates[2] = Triangle_RotatedZX.p[1].coordinates[2] + 3.0f;
			Triangle_Translated.p[2].coordinates[2] = Triangle_RotatedZX.p[2].coordinates[2] + 3.0f;

			// Project triangles from 3D --> 2D
			MultiplyMatrixVector(Triangle_Translated.p[0], Triangle_Proj.p[0], ProjectionMatrix);
			MultiplyMatrixVector(Triangle_Translated.p[1], Triangle_Proj.p[1], ProjectionMatrix);
			MultiplyMatrixVector(Triangle_Translated.p[2], Triangle_Proj.p[2], ProjectionMatrix);
			input = "Mesh Thread#" + to_string(omp_get_thread_num()) + "\nTranslated Point X: " + to_string(Triangle_Translated.p[0].coordinates[0]) +
				", Point Y: " + to_string(Triangle_Translated.p[0].coordinates[1]) +
				", Point Z: " + to_string(Triangle_Translated.p[0].coordinates[2]) +
				"--> Point X: " + to_string(Triangle_Proj.p[0].coordinates[0]) +
				",Point Y: " + to_string(Triangle_Proj.p[0].coordinates[1]) +
				",Point Z: " + to_string(Triangle_Proj.p[0].coordinates[1]) +"\n";
			// Scale into view
			out << input;

			Triangle_Proj.p[0].coordinates[0] += 1.0f; Triangle_Proj.p[0].coordinates[1] += 1.0f;
			Triangle_Proj.p[1].coordinates[0] += 1.0f; Triangle_Proj.p[1].coordinates[1] += 1.0f;
			Triangle_Proj.p[2].coordinates[0] += 1.0f; Triangle_Proj.p[2].coordinates[1] += 1.0f;
			Triangle_Proj.p[0].coordinates[0] *= 0.5f * (float)ScreenWidth();
			Triangle_Proj.p[0].coordinates[1] *= 0.5f * (float)ScreenHeight();
			Triangle_Proj.p[1].coordinates[0] *= 0.5f * (float)ScreenWidth();
			Triangle_Proj.p[1].coordinates[1] *= 0.5f * (float)ScreenHeight();
			Triangle_Proj.p[2].coordinates[0] *= 0.5f * (float)ScreenWidth();
			Triangle_Proj.p[2].coordinates[1] *= 0.5f * (float)ScreenHeight();

			// Rasterize triangle
			DrawTriangle(Triangle_Proj.p[0].coordinates[0], 
				Triangle_Proj.p[0].coordinates[1],
				Triangle_Proj.p[1].coordinates[0], 
				Triangle_Proj.p[1].coordinates[1],
				Triangle_Proj.p[2].coordinates[0], 
				Triangle_Proj.p[2].coordinates[1],
				PIXEL_SOLID, FG_WHITE);
		}


		return true;
	}

};




int main()
{
	Engine3D Saud3D;
	if (Saud3D.ConstructConsole(256, 240, 4, 4))
		Saud3D.Start();
	return 0;
}

