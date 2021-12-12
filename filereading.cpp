#include <iostream>
#include <fstream>
#include<vector>
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
int main()
{
    mesh meshCube;
    std::fstream myfile("RectangleCoordinates.txt", std::ios_base::in);
    int i = 0,j=0;
    float a;
    triangle b;
    while (myfile >> a)
    {
        b.p[j].coordinates[i] = a;
        printf("%f ", a);
        i++;
        if (i == 3)
        {
            i = 0;
            j++;
        }
        if (j == 3)
        {
            i = 0, j = 0;
            meshCube.tris.push_back(b);
        }
    }



    return 0;
}