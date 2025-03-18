// ConsoleGameEngine.cpp : This file contains the 'main' function. Program execution begins and ends there.
#include <iostream>
#define OLC_PGE_APPLICATION
#include "olcPixelGameEngine.h"
struct vec3D
{
    float x, y, z;
};
struct triangle
{
    vec3D p[3];
};
struct mesh
{
    std::vector<triangle> tris;
};
struct mat4x4
{
    float m[4][4] = { 0 };
};
class Engine3D:public olc::PixelGameEngine
{
public:
    Engine3D()
    {

        sAppName = "3D Demo";
    }
private:
    mesh meshCube;
    mat4x4 matProj;
    float fTheta = 0.0f;
    void MultiplyMatrix(vec3D& i, vec3D& output, mat4x4& m)
    {
        output.x = i.x * m.m[0][0] + i.y * m.m[1][0] + i.z * m.m[2][0] + m.m[3][0];
        output.y = i.x * m.m[0][1] + i.y * m.m[1][1] + i.z * m.m[2][1] + m.m[3][1];
        output.z = i.x * m.m[0][2] + i.y * m.m[1][2] + i.z * m.m[2][2] + m.m[3][2];
        float w =  i.x * m.m[0][3] + i.y * m.m[1][3] + i.z * m.m[2][3] + m.m[3][3];
        if (w != 0.0f)
        {
            output.x /= w;
            output.y /= w;
            output.z /= w;
        }
    }
public:
    bool OnUserCreate()override
    {
        meshCube.tris = {
            //SOUTH
            {0.0f,0.0f,0.0f, 0.0f,1.0f,0.0f, 1.0f,1.0f,0.0f},
            {0.0f,0.0f,0.0f, 1.0f,1.0f,0.0f, 1.0f,0.0f,0.0f},

            //EAST
            {1.0f,0.0f,0.0f, 1.0f,1.0f,0.0f, 1.0f,1.0f,1.0f},
            {1.0f,0.0f,0.0f, 1.0f,1.0f,1.0f, 1.0f,0.0f,1.0f},

            //NORTH
            { 1.0f,0.0f,1.0f, 1.0f,1.0f,1.0f, 0.0f,1.0f,1.0f },
            { 1.0f,0.0f,1.0f, 0.0f,1.0f,1.0f, 0.0f,0.0f,1.0f },

            //WEST
            {0.0f,0.0f,1.0f, 0.0f,1.0f,1.0f, 0.0f,1.0f,0.0f},
            {0.0f,0.0f,1.0f, 0.0f,1.0f,0.0f, 0.0f,0.0f,0.0f},

            //TOP
            {0.0f,1.0f,0.0f, 0.0f,1.0f,1.0f, 1.0f,1.0f,1.0f},
            {0.0f,1.0f,0.0f, 1.0f,1.0f,1.0f, 1.0f,1.0f,0.0f},

            //BOTTOM
            {1.0f,0.0f,1.0f, 0.0f,0.0f,1.0f, 0.0f,0.0f, 0.0f },
            {1.0f,0.0f,1.0f, 0.0f,0.0f,0.0f, 1.0f,0.0f, 0.0f },
        };

        //populate the perspective matrix
        float fNear = 0.1f;
        float fFar = 1000.0f;
        float fFov = 90.0f;
        float fAspect = (float)ScreenHeight() / (float)ScreenWidth();
        float fFovRad = 1.0 / tanf(fFov * 0.5f / 180.0f * 3.14159265358979323846f);//TODO: (fFov * 0.5f/ 180.0f * 3.14159f);
        matProj.m[0][0] = fAspect * fFovRad;
        matProj.m[1][1] = fFovRad;
        matProj.m[2][2] = fFar / (fFar - fNear);
        matProj.m[3][2] = (- fFar * fNear )/ (fFar - fNear);
        matProj.m[2][3] = 1.0f;
        matProj.m[3][3] = 0.0f;
        return true;
    }
    bool OnUserUpdate(float fElapsedTime)override
    {
        for (int x = 0; x < ScreenWidth(); x++)
        {
            for (int y = 0; y < ScreenHeight(); y++)
            {
                Draw(x, y,olc::BLACK);
            }
        }

        //making rotation matrix
        mat4x4 matRotZ, matRotX;
        fTheta += 1.0f * fElapsedTime;
        //Rotation Z 
        matRotZ.m[0][0] = cos(fTheta);
        matRotZ.m[0][1] = -sin(fTheta);
        matRotZ.m[1][0] = sin(fTheta);
        matRotZ.m[1][1] = cos(fTheta);
        matRotZ.m[2][2] = 1.0f;
        matRotZ.m[3][3] = 1.0f;

        //Rotation X
        matRotX.m[0][0] = 1.0f;
        matRotX.m[1][1] = cos(fTheta * 0.5f);
        matRotX.m[1][2] = -sin(fTheta * 0.5f);
        matRotX.m[2][1] = sin(fTheta * 0.5f);
        matRotX.m[2][2] = cos(fTheta * 0.5f);
        matRotX.m[3][3] = 1.0f;

        //Draw Triangle
        for (auto &tri : meshCube.tris)
        {
            triangle triProjected,triTranslated,triRotatedZ,triRotatedZX;
            //Rotating the triangle by z then x axis
            MultiplyMatrix(tri.p[0], triRotatedZ.p[0], matRotZ);
            MultiplyMatrix(tri.p[1], triRotatedZ.p[1], matRotZ);
            MultiplyMatrix(tri.p[2], triRotatedZ.p[2], matRotZ);

            MultiplyMatrix(triRotatedZ.p[0], triRotatedZX.p[0], matRotX);
            MultiplyMatrix(triRotatedZ.p[1], triRotatedZX.p[1], matRotX);
            MultiplyMatrix(triRotatedZ.p[2], triRotatedZX.p[2], matRotX);
            //here we offset the triangle into the z axis
            triTranslated = triRotatedZX;
            triTranslated.p[0].z = triRotatedZX.p[0].z + 3.0f;
            triTranslated.p[1].z = triRotatedZX.p[1].z + 3.0f;
            triTranslated.p[2].z = triRotatedZX.p[2].z + 3.0f;
            MultiplyMatrix(triTranslated.p[0], triProjected.p[0], matProj);
            MultiplyMatrix(triTranslated.p[1], triProjected.p[1], matProj);
            MultiplyMatrix(triTranslated.p[2], triProjected.p[2], matProj);
            //scalling the triangle
            triProjected.p[0].x += 1.0f; triProjected.p[0].y += 1.0f;
            triProjected.p[1].x += 1.0f; triProjected.p[1].y += 1.0f;
            triProjected.p[2].x += 1.0f; triProjected.p[2].y += 1.0f;

            triProjected.p[0].x *= 0.5f * (float)ScreenWidth();
            triProjected.p[0].y *= 0.5f * (float)ScreenHeight();
            triProjected.p[1].x *= 0.5f * (float)ScreenWidth();
            triProjected.p[1].y *= 0.5f * (float)ScreenHeight();
            triProjected.p[2].x *= 0.5f * (float)ScreenWidth();
            triProjected.p[2].y *= 0.5f * (float)ScreenHeight();

            DrawTriangle(triProjected.p[0].x, triProjected.p[0].y,
                triProjected.p[1].x, triProjected.p[1].y,
                triProjected.p[2].x, triProjected.p[2].y,olc::RED);
        }
        return true;
    }
};
int main()
{
    Engine3D demo;
    if (demo.Construct(256, 240, 4, 4))
        demo.Start();
    std::cout << "Hello World!\n";
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
