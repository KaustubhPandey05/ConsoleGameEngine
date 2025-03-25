// ConsoleGameEngine.cpp : This file contains the 'main' function. Program execution begins and ends there.
#include <iostream>
#define OLC_PGE_APPLICATION
#include "olcPixelGameEngine.h"
#include <fstream>
#include <sstream>
#include <algorithm>
struct vec3D
{
    float x = 0;
    float y = 0;
    float z = 0;
    float w = 1;// Need a 4th term to perform sensible matrix vector multiplication

    vec3D operator+(vec3D& v2)
    {
        return { x + v2.x, y + v2.y, z + v2.z };
    }
    vec3D operator-(vec3D& v2)
    {
        return { x - v2.x, y - v2.y, z - v2.z };
    }
    vec3D operator*(float k)
    {
        return { x * k, y * k, z * k };
    }
    vec3D operator/(float k)
    {
        return { x / k, y / k, z / k };
    }
};
struct triangle
{
    vec3D p[3];
    olc::Pixel col;//this should be stored per vertex instead
};
struct mesh
{
    std::vector<triangle> tris;
    bool loadFromObjFile(std::string fileName)
    {
        std::ifstream f(fileName);
        if (!f.is_open())
            return false;
        std::vector<vec3D> verts;
        while (!f.eof())
        {
            char line[128];
            f.getline(line,128);
            std::stringstream s;
            //std::strstream s;
            s << line;
            char junk;
            if (line[0] == 'v')
            {
                vec3D v;
                s >> junk >> v.x >> v.y >> v.z;
                verts.push_back(v);
            }
            if (line[0] == 'f')
            {
                int f[3];
                s >> junk >> f[0] >> f[1] >> f[2];
                tris.push_back({ verts[f[0] - 1],verts[f[1] - 1],verts[f[2] - 1] });
            }
        }
        return true;
    }
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

    vec3D vCamera;
    vec3D vLookDir;

    float fYaw;// FPS Camera rotation in XZ plane
    float fTheta = 0.0f;
    vec3D Matrix_MultiplyVector(mat4x4& m, vec3D& i)
    {
        vec3D v;
        v.x = i.x * m.m[0][0] + i.y * m.m[1][0] + i.z * m.m[2][0] + i.w * m.m[3][0];
        v.y = i.x * m.m[0][1] + i.y * m.m[1][1] + i.z * m.m[2][1] + i.w * m.m[3][1];
        v.z = i.x * m.m[0][2] + i.y * m.m[1][2] + i.z * m.m[2][2] + i.w * m.m[3][2];
        v.w = i.x * m.m[0][3] + i.y * m.m[1][3] + i.z * m.m[2][3] + i.w * m.m[3][3];
        return v;
    }
    mat4x4 Matrix_MakeIdentity()
    {
        mat4x4 matrix;
        matrix.m[0][0] = 1.0f;
        matrix.m[1][1] = 1.0f;
        matrix.m[2][2] = 1.0f;
        matrix.m[3][3] = 1.0f;
        return matrix;
    }
    mat4x4 Matrix_MakeRotationX(float fAngleRad)
    {
        mat4x4 matrix;
        matrix.m[0][0] = 1.0f;
        matrix.m[1][1] = cosf(fAngleRad);
        matrix.m[1][2] = sinf(fAngleRad);
        matrix.m[2][1] = -sinf(fAngleRad);
        matrix.m[2][2] = cosf(fAngleRad);
        matrix.m[3][3] = 1.0f;
        return matrix;
    }
    mat4x4 Matrix_MakeRotationY(float fAngleRad)
    {
        mat4x4 matrix;
        matrix.m[0][0] = cosf(fAngleRad);
        matrix.m[0][2] = sinf(fAngleRad);
        matrix.m[2][0] = -sinf(fAngleRad);
        matrix.m[1][1] = 1.0f;
        matrix.m[2][2] = cosf(fAngleRad);
        matrix.m[3][3] = 1.0f;
        return matrix;
    }

    mat4x4 Matrix_MakeRotationZ(float fAngleRad)
    {
        mat4x4 matrix;
        matrix.m[0][0] = cosf(fAngleRad);
        matrix.m[0][1] = sinf(fAngleRad);
        matrix.m[1][0] = -sinf(fAngleRad);
        matrix.m[1][1] = cosf(fAngleRad);
        matrix.m[2][2] = 1.0f;
        matrix.m[3][3] = 1.0f;
        return matrix;
    }

    mat4x4 Matrix_MakeTranslation(float x, float y, float z)
    {
        mat4x4 matrix;
        matrix.m[0][0] = 1.0f;
        matrix.m[1][1] = 1.0f;
        matrix.m[2][2] = 1.0f;
        matrix.m[3][3] = 1.0f;
        matrix.m[3][0] = x;
        matrix.m[3][1] = y;
        matrix.m[3][2] = z;
        return matrix;
    }
    mat4x4 Matrix_MakeProjection(float fFovDegrees, float fAspectRatio, float fNear, float fFar)
    {
        float fFovRad = 1.0f / tanf(fFovDegrees * 0.5f / 180.0f * 3.14159f);
        mat4x4 matrix;
        matrix.m[0][0] = fAspectRatio * fFovRad;
        matrix.m[1][1] = fFovRad;
        matrix.m[2][2] = fFar / (fFar - fNear);
        matrix.m[3][2] = (-fFar * fNear) / (fFar - fNear);
        matrix.m[2][3] = 1.0f;
        matrix.m[3][3] = 0.0f;
        return matrix;
    }
    mat4x4 Matrix_MultiplyMatrix(mat4x4& m1, mat4x4& m2)
    {
        mat4x4 matrix;
        for (int c = 0; c < 4; c++)
            for (int r = 0; r < 4; r++)
                matrix.m[r][c] = m1.m[r][0] * m2.m[0][c] + m1.m[r][1] * m2.m[1][c] + m1.m[r][2] * m2.m[2][c] + m1.m[r][3] * m2.m[3][c];
        return matrix;
    }

    mat4x4 Matrix_PointAt(vec3D& pos, vec3D& target, vec3D& up)
    {
        // Calculate new forward direction
        vec3D newForward = target-pos;
        newForward = Vector_Normalise(newForward);
        // Calculate new Up direction
        vec3D a = newForward * Vector_DotProduct(up, newForward);
        vec3D newUp = up - a;
        newUp = Vector_Normalise(newUp);

        // New Right direction is easy, its just cross product
        vec3D newRight = Vector_CrossProduct(newUp, newForward);
      
        // Construct Dimensioning and Translation Matrix	
        mat4x4 matrix;
        matrix.m[0][0] = newRight.x;	matrix.m[0][1] = newRight.y;	matrix.m[0][2] = newRight.z;	matrix.m[0][3] = 0.0f;
        matrix.m[1][0] = newUp.x;		matrix.m[1][1] = newUp.y;		matrix.m[1][2] = newUp.z;		matrix.m[1][3] = 0.0f;
        matrix.m[2][0] = newForward.x;	matrix.m[2][1] = newForward.y;	matrix.m[2][2] = newForward.z;	matrix.m[2][3] = 0.0f;
        matrix.m[3][0] = pos.x;			matrix.m[3][1] = pos.y;			matrix.m[3][2] = pos.z;			matrix.m[3][3] = 1.0f;
        return matrix;

    }
    mat4x4 Matrix_QuickInverse(mat4x4& m) // Only for Rotation/Translation Matrices
    {
        mat4x4 matrix;
        matrix.m[0][0] = m.m[0][0]; matrix.m[0][1] = m.m[1][0]; matrix.m[0][2] = m.m[2][0]; matrix.m[0][3] = 0.0f;
        matrix.m[1][0] = m.m[0][1]; matrix.m[1][1] = m.m[1][1]; matrix.m[1][2] = m.m[2][1]; matrix.m[1][3] = 0.0f;
        matrix.m[2][0] = m.m[0][2]; matrix.m[2][1] = m.m[1][2]; matrix.m[2][2] = m.m[2][2]; matrix.m[2][3] = 0.0f;
        matrix.m[3][0] = -(m.m[3][0] * matrix.m[0][0] + m.m[3][1] * matrix.m[1][0] + m.m[3][2] * matrix.m[2][0]);
        matrix.m[3][1] = -(m.m[3][0] * matrix.m[0][1] + m.m[3][1] * matrix.m[1][1] + m.m[3][2] * matrix.m[2][1]);
        matrix.m[3][2] = -(m.m[3][0] * matrix.m[0][2] + m.m[3][1] * matrix.m[1][2] + m.m[3][2] * matrix.m[2][2]);
        matrix.m[3][3] = 1.0f;
        return matrix;
    }

    float Vector_DotProduct(vec3D& v1, vec3D& v2)
    {
        return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
    }

    float Vector_Length(vec3D& v)
    {
        return sqrtf(Vector_DotProduct(v, v));
    }

    vec3D Vector_Normalise(vec3D& v)
    {
        float l = Vector_Length(v);
        return { v.x / l, v.y / l, v.z / l };
    }

    vec3D Vector_CrossProduct(vec3D& v1, vec3D& v2)
    {
        vec3D v;
        v.x = v1.y * v2.z - v1.z * v2.y;
        v.y = v1.z * v2.x - v1.x * v2.z;
        v.z = v1.x * v2.y - v1.y * v2.x;
        return v;
    }

    vec3D Vector_IntersectPlane(vec3D& plane_p, vec3D& plane_n, vec3D& lineStart, vec3D& lineEnd)
    {
        plane_n = Vector_Normalise(plane_n);
        float plane_d = -Vector_DotProduct(plane_n, plane_p);
        float ad = Vector_DotProduct(lineStart, plane_n);
        float bd = Vector_DotProduct(lineEnd, plane_n);
        float t = -(plane_d - ad) / (bd - ad);
        vec3D lineStartToEnd = lineEnd - lineStart;
        vec3D lineToIntersect = lineStartToEnd * t;
        return lineStart + lineStartToEnd;
    }

    int Triangle_ClippingPlane(vec3D plane_p, vec3D plane_n,triangle &in_tri,triangle &out_tri1,triangle &out_tri2)//for a given clipping we are giving out two output traingles atmost
    {
        // Make sure plane normal is indeed normal
        plane_n = Vector_Normalise(plane_n);

        // Return signed shortest distance from point to plane, plane normal must be normalised
        auto dist = [&](vec3D& p)
            {
                vec3D n = Vector_Normalise(p);
                return (plane_n.x * p.x + plane_n.y * p.y + plane_n.z * p.z - Vector_DotProduct(plane_n, plane_p));
            };
        // Create two temporary storage arrays to classify points either side of plane
        // If distance sign is positive, point lies on "inside" of plane
        vec3D* inside_points[3];  int nInsidePointsCount = 0;
        vec3D* outside_points[3]; int nOutsidePointsCount = 0;

        float d0 = dist(in_tri.p[0]);
        float d1 = dist(in_tri.p[1]);
        float d2 = dist(in_tri.p[2]);

        if (d0 >= 0) { inside_points[nInsidePointsCount++] = &in_tri.p[0]; }
        else { outside_points[nOutsidePointsCount++] = &in_tri.p[0]; }
        if (d1 >= 0) { inside_points[nInsidePointsCount++] = &in_tri.p[1]; }
        else { outside_points[nOutsidePointsCount++] = &in_tri.p[1]; }
        if(d2>=0){ inside_points[nInsidePointsCount++] = &in_tri.p[2]; }
        else { outside_points[nOutsidePointsCount++] = &in_tri.p[2]; }

        // Now classify triangle points, and break the input triangle into 
        // smaller output triangles if required. There are four possible
        // outcomes...

        if (nInsidePointsCount == 0)
        {
            // All points lie on the inside of plane, so do nothing
            // and allow the triangle to simply pass through
            out_tri1 = in_tri;

            return 1; // Just the one returned original triangle is valid
        }
        if (nInsidePointsCount == 1 && nOutsidePointsCount == 2)
        {
            // Triangle should be clipped. As two points lie outside
            // the plane, the triangle simply becomes a smaller triangle

            // Copy appearance info to new triangle

            out_tri1.col = in_tri.col;

            // The inside point is valid, so keep that...
            out_tri1.p[0] = *inside_points[0];

            // but the two new points are at the locations where the 
            // original sides of the triangle (lines) intersect with the plane
            out_tri1.p[1] = Vector_IntersectPlane(plane_p, plane_n, *outside_points[0], *outside_points[0]);
            out_tri1.p[2] = Vector_IntersectPlane(plane_p, plane_n, *outside_points[0], *outside_points[1]);
            return 1;//Return the newly formed single triangle
        }
        if (nInsidePointsCount == 2 && nOutsidePointsCount == 1)
        {
            // Triangle should be clipped. As two points lie inside the plane,
            // the clipped triangle becomes a "quad". Fortunately, we can
            // represent a quad with two new triangles

            // Copy appearance info to new triangles
            out_tri1.col = in_tri.col;
            out_tri2.col = in_tri.col;

            // The first triangle consists of the two inside points and a new
            // point determined by the location where one side of the triangle
            // intersects with the plane
            out_tri1.p[0] = *inside_points[0];
            out_tri1.p[1] = *inside_points[1];
            out_tri1.p[2] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[0]);

            // The second triangle is composed of one of he inside points, a
            // new point determined by the intersection of the other side of the 
            // triangle and the plane, and the newly created point above
            out_tri2.p[0] = *inside_points[1];
            out_tri2.p[1] = out_tri1.p[2];
            out_tri2.p[2] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[1], *outside_points[0]);
            return 2; // Return two newly formed triangles which form a quad
        }
    }
public:
    olc::Pixel GetColour(float lum)
    {
        int i = (int)(lum * 255.0f);
        i = std::max(0, std::min(255, i));
        return olc::Pixel( i,i,i );
    }
    bool OnUserCreate()override
    {
        meshCube.loadFromObjFile("axis.obj");
        //populate the perspective matrix
        matProj = Matrix_MakeProjection(90.0f, (float)ScreenHeight() / (float)ScreenWidth(), 0.1f, 1000.0f);
        return true;
    }
    bool OnUserUpdate(float fElapsedTime)override
    {
        if (GetKey(olc::Key::UP).bHeld)
            vCamera.y += 8.0f * fElapsedTime;	// Travel Upwards

        if (GetKey(olc::Key::DOWN).bHeld)
            vCamera.y -= 8.0f * fElapsedTime;	// Travel Downwards


        // Dont use these two in FPS mode, it is confusing :P
        if (GetKey(olc::Key::LEFT).bHeld)
            vCamera.x -= 8.0f * fElapsedTime;	// Travel Along X-Axis

        if (GetKey(olc::Key::RIGHT).bHeld)
            vCamera.x += 8.0f * fElapsedTime;	// Travel Along X-Axis

        vec3D vForward = vLookDir * (8.0f * fElapsedTime);
        // Standard FPS Control scheme, but turn instead of strafe
        if (GetKey(olc::Key::W).bHeld)
            vCamera = (vCamera + vForward);

        if (GetKey(olc::Key::S).bHeld)
            vCamera = (vCamera - vForward);

        if (GetKey(olc::Key::A).bHeld)
            fYaw -= 2.0f * fElapsedTime;

        if (GetKey(olc::Key::D).bHeld)
            fYaw += 2.0f * fElapsedTime;
        
        for (int x = 0; x < ScreenWidth(); x++)
        {
            for (int y = 0; y < ScreenHeight(); y++)
            {
                Draw(x, y,olc::BLACK);
            }
        }

        //making rotation matrix
        mat4x4 matRotZ, matRotX;
        //fTheta += 1.0f * fElapsedTime;// Uncomment to spin me right round baby right round
        //Rotation Z 
        matRotZ = Matrix_MakeRotationZ(fTheta * 0.5f);

        //Rotation X
        matRotX = Matrix_MakeRotationX(fTheta);
        mat4x4 matTrans;
        matTrans = Matrix_MakeTranslation(0.0f, 0.0f, 5.0f);

        mat4x4 matWorld;
        matWorld = Matrix_MakeIdentity();
        matWorld = Matrix_MultiplyMatrix(matRotZ, matRotX);
        matWorld = Matrix_MultiplyMatrix(matWorld, matTrans);


        // Create "Point At" Matrix for camera
        vec3D vUp = { 0,1,0 };
        vec3D vTarget = { 0,0,1 };
        mat4x4 matCameraRot = Matrix_MakeRotationY(fYaw);
        vLookDir = Matrix_MultiplyVector(matCameraRot, vTarget);
        vTarget = vCamera + vLookDir;
        mat4x4 matCamera = Matrix_PointAt(vCamera, vTarget, vUp);

        // Make view matrix from camera
        mat4x4 matView = Matrix_QuickInverse(matCamera);

        // Store triagles for rastering later
        std::vector<triangle> vecTrianglesToRaster;
        //Draw Triangle
        for (auto &tri : meshCube.tris)
        {
            triangle triProjected, triTransformed, triViewed;
            // World Matrix 
            triTransformed.p[0] = Matrix_MultiplyVector(matWorld, tri.p[0]);
            triTransformed.p[1] = Matrix_MultiplyVector(matWorld, tri.p[1]);
            triTransformed.p[2] = Matrix_MultiplyVector(matWorld, tri.p[2]);

            vec3D normal, line1, line2;
            // Get lines either side of triangle
            line1 = triTransformed.p[1] - triTransformed.p[0];
            line2 = triTransformed.p[2] - triTransformed.p[0];
            
            // Take cross product of lines to get normal to triangle surface
            normal = Vector_CrossProduct(line1, line2);
           
            // You normally need to normalise a normal!
            normal = Vector_Normalise(normal);

            // Get Ray from triangle to camera
            vec3D vCameraRay = triTransformed.p[0] - vCamera;

            if (Vector_DotProduct(normal, vCameraRay) < 0.0f)
            {
                //Illumination
                vec3D light_direction = { 0.0f,0.0f,-1.0f };

                //normalise the light vector
                light_direction = Vector_Normalise(light_direction);

                // How "aligned" are light direction and triangle surface normal?
                float dp = std::max(0.1f, Vector_DotProduct(light_direction, normal));

                olc::Pixel c = GetColour(dp);//get the gray shade by looking at the dot product value
                triTransformed.col = c;//this is very console specific stuff

                // Convert World Space --> View Space
                triViewed.p[0] = Matrix_MultiplyVector(matView, triTransformed.p[0]);
                triViewed.p[1] = Matrix_MultiplyVector(matView, triTransformed.p[1]);
                triViewed.p[2] = Matrix_MultiplyVector(matView, triTransformed.p[2]);
                triViewed.col = triTransformed.col;

                // Project triangles from 3D --> 2D
                triProjected.p[0] = Matrix_MultiplyVector(matProj, triViewed.p[0]);
                triProjected.p[1] = Matrix_MultiplyVector(matProj, triViewed.p[1]);
                triProjected.p[2] = Matrix_MultiplyVector(matProj, triViewed.p[2]);
                triProjected.col = triViewed.col;

                // Scale into view, we moved the normalising into cartesian space
                // out of the matrix.vector function from the previous videos, so
                // do this manually
                triProjected.p[0] = (triProjected.p[0]/ triProjected.p[0].w);
                triProjected.p[1] = (triProjected.p[1]/ triProjected.p[1].w);
                triProjected.p[2] = (triProjected.p[2]/ triProjected.p[2].w);

                // Offset verts into visible normalised space
                vec3D vOffsetView = { 1,1,0 };
                triProjected.p[0] = (triProjected.p[0] + vOffsetView);
                triProjected.p[1] = (triProjected.p[1] + vOffsetView);
                triProjected.p[2] = (triProjected.p[2] + vOffsetView);
                triProjected.p[0].x *= 0.5f * (float)ScreenWidth();
                triProjected.p[0].y *= 0.5f * (float)ScreenHeight();
                triProjected.p[1].x *= 0.5f * (float)ScreenWidth();
                triProjected.p[1].y *= 0.5f * (float)ScreenHeight();
                triProjected.p[2].x *= 0.5f * (float)ScreenWidth();
                triProjected.p[2].y *= 0.5f * (float)ScreenHeight();
                vecTrianglesToRaster.push_back(triProjected);
            }
        }   

        //sort the triangles to draw them first from the furthest away
        std::sort(vecTrianglesToRaster.begin(), vecTrianglesToRaster.end(),[](triangle t1, triangle t2)
        {
                float z1 = (t1.p[0].z + t1.p[1].z + t1.p[2].z) / 3.0f;
                float z2 = (t2.p[0].z + t2.p[1].z + t2.p[2].z) / 3.0f;
                return z1 > z2;
        });
        for (auto& triProjected : vecTrianglesToRaster)
        {
            FillTriangle(triProjected.p[0].x, triProjected.p[0].y,
                triProjected.p[1].x, triProjected.p[1].y,
                triProjected.p[2].x, triProjected.p[2].y, triProjected.col);

            /*DrawTriangle(triProjected.p[0].x, triProjected.p[0].y,
                triProjected.p[1].x, triProjected.p[1].y,
                triProjected.p[2].x, triProjected.p[2].y, olc::BLACK);*/
        }
        return true;
    }
};
int main()
{
    Engine3D demo;
    if (demo.Construct(256, 240, 4, 4))
        demo.Start();
}
