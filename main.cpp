#define NOMINMAX

#include <Novice.h>
#define _USE_MATH_DEFINES
#include <algorithm>
#include <assert.h>
#include <cmath>
#include <cstdint>
#include <imgui.h>
#include <math.h>

struct Vector3 {
    float x;
    float y;
    float z;
};
struct Matrix4x4 {
    float m[4][4];
};
struct Sphere {
    Vector3 center;
    float radius;
    int color;
};

static const int krowheight = 20;
static const int kcolumnwidth = 60;

void MatrixScreenPrintf(int x, int y, const Matrix4x4& matrix, const char* num)
{

    Novice::ScreenPrintf(x, y, "%s", num);
    for (int row = 0; row < 4; ++row) {
        for (int column = 0; column < 4; ++column) {
            Novice::ScreenPrintf(x + column * kcolumnwidth, (y + row * krowheight) + 20, "%6.03f", matrix.m[row][column]);
        }
    }
}

Vector3 add(const Vector3& v1, const Vector3& v2)
{
    Vector3 num;
    num.x = v1.x + v2.x;
    num.y = v1.y + v2.y;
    num.z = v1.z + v2.z;
    return num;
}
Vector3 Subtract(const Vector3& v1, const Vector3& v2)
{
    Vector3 num;
    num.x = v1.x - v2.x;
    num.y = v1.y - v2.y;
    num.z = v1.z - v2.z;
    return num;
}
float Dot(const Vector3& v1, const Vector3& v2)
{
    float num;
    return num = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}
float Length(const Vector3& v)
{
    float num;
    return num = sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
}
Matrix4x4 Multiply(const Matrix4x4& m1, const Matrix4x4& m2)
{
    Matrix4x4 num;
    num.m[0][0] = m1.m[0][0] * m2.m[0][0] + m1.m[0][1] * m2.m[1][0] + m1.m[0][2] * m2.m[2][0] + m1.m[0][3] * m2.m[3][0];
    num.m[0][1] = m1.m[0][0] * m2.m[0][1] + m1.m[0][1] * m2.m[1][1] + m1.m[0][2] * m2.m[2][1] + m1.m[0][3] * m2.m[3][1];
    num.m[0][2] = m1.m[0][0] * m2.m[0][2] + m1.m[0][1] * m2.m[1][2] + m1.m[0][2] * m2.m[2][2] + m1.m[0][3] * m2.m[3][2];
    num.m[0][3] = m1.m[0][0] * m2.m[0][3] + m1.m[0][1] * m2.m[1][3] + m1.m[0][2] * m2.m[2][3] + m1.m[0][3] * m2.m[3][3];

    num.m[1][0] = m1.m[1][0] * m2.m[0][0] + m1.m[1][1] * m2.m[1][0] + m1.m[1][2] * m2.m[2][0] + m1.m[1][3] * m2.m[3][0];
    num.m[1][1] = m1.m[1][0] * m2.m[0][1] + m1.m[1][1] * m2.m[1][1] + m1.m[1][2] * m2.m[2][1] + m1.m[1][3] * m2.m[3][1];
    num.m[1][2] = m1.m[1][0] * m2.m[0][2] + m1.m[1][1] * m2.m[1][2] + m1.m[1][2] * m2.m[2][2] + m1.m[1][3] * m2.m[3][2];
    num.m[1][3] = m1.m[1][0] * m2.m[0][3] + m1.m[1][1] * m2.m[1][3] + m1.m[1][2] * m2.m[2][3] + m1.m[1][3] * m2.m[3][3];

    num.m[2][0] = m1.m[2][0] * m2.m[0][0] + m1.m[2][1] * m2.m[1][0] + m1.m[2][2] * m2.m[2][0] + m1.m[2][3] * m2.m[3][0];
    num.m[2][1] = m1.m[2][0] * m2.m[0][1] + m1.m[2][1] * m2.m[1][1] + m1.m[2][2] * m2.m[2][1] + m1.m[2][3] * m2.m[3][1];
    num.m[2][2] = m1.m[2][0] * m2.m[0][2] + m1.m[2][1] * m2.m[1][2] + m1.m[2][2] * m2.m[2][2] + m1.m[2][3] * m2.m[3][2];
    num.m[2][3] = m1.m[2][0] * m2.m[0][3] + m1.m[2][1] * m2.m[1][3] + m1.m[2][2] * m2.m[2][3] + m1.m[2][3] * m2.m[3][3];

    num.m[3][0] = m1.m[3][0] * m2.m[0][0] + m1.m[3][1] * m2.m[1][0] + m1.m[3][2] * m2.m[2][0] + m1.m[3][3] * m2.m[3][0];
    num.m[3][1] = m1.m[3][0] * m2.m[0][1] + m1.m[3][1] * m2.m[1][1] + m1.m[3][2] * m2.m[2][1] + m1.m[3][3] * m2.m[3][1];
    num.m[3][2] = m1.m[3][0] * m2.m[0][2] + m1.m[3][1] * m2.m[1][2] + m1.m[3][2] * m2.m[2][2] + m1.m[3][3] * m2.m[3][2];
    num.m[3][3] = m1.m[3][0] * m2.m[0][3] + m1.m[3][1] * m2.m[1][3] + m1.m[3][2] * m2.m[2][3] + m1.m[3][3] * m2.m[3][3];

    return num;
}
Vector3 Multiply(const float& m2, const Vector3& m1)
{
    Vector3 num;
    num.x = m1.x * m2;
    num.y = m1.y * m2;
    num.z = m1.z * m2;

    return num;
}
Vector3 Multiply(const Vector3& m2, const Vector3& m1)
{
    Vector3 num;
    num.x = m1.x * m2.x;
    num.y = m1.y * m2.y;
    num.z = m1.z * m2.z;

    return num;
}

Matrix4x4 Inverse(const Matrix4x4& m)
{
    float determinant;
    Matrix4x4 num;

    determinant = m.m[0][0] * m.m[1][1] * m.m[2][2] * m.m[3][3] + m.m[0][0] * m.m[1][2] * m.m[2][3] * m.m[3][1] + m.m[0][0] * m.m[1][3] * m.m[2][1] * m.m[3][2]
        - m.m[0][0] * m.m[1][3] * m.m[2][2] * m.m[3][1] - m.m[0][0] * m.m[1][2] * m.m[2][1] * m.m[3][3] - m.m[0][0] * m.m[1][1] * m.m[2][3] * m.m[3][2]
        - m.m[0][1] * m.m[1][0] * m.m[2][2] * m.m[3][3] - m.m[0][2] * m.m[1][0] * m.m[2][3] * m.m[3][1] - m.m[0][3] * m.m[1][0] * m.m[2][1] * m.m[3][2]
        + m.m[0][3] * m.m[1][0] * m.m[2][2] * m.m[3][1] + m.m[0][2] * m.m[1][0] * m.m[2][1] * m.m[3][3] + m.m[0][1] * m.m[1][0] * m.m[2][3] * m.m[3][2]
        + m.m[0][1] * m.m[1][2] * m.m[2][0] * m.m[3][3] + m.m[0][2] * m.m[1][3] * m.m[2][0] * m.m[3][1] + m.m[0][3] * m.m[1][1] * m.m[2][0] * m.m[3][2]
        - m.m[0][3] * m.m[1][2] * m.m[2][0] * m.m[3][1] - m.m[0][2] * m.m[1][1] * m.m[2][0] * m.m[3][3] - m.m[0][1] * m.m[1][3] * m.m[2][0] * m.m[3][2]
        - m.m[0][1] * m.m[1][2] * m.m[2][3] * m.m[3][0] - m.m[0][2] * m.m[1][3] * m.m[2][1] * m.m[3][0] - m.m[0][3] * m.m[1][1] * m.m[2][2] * m.m[3][0]
        + m.m[0][3] * m.m[1][2] * m.m[2][1] * m.m[3][0] + m.m[0][2] * m.m[1][1] * m.m[2][3] * m.m[3][0] + m.m[0][1] * m.m[1][3] * m.m[2][2] * m.m[3][0];

    if (determinant == 0.0f) {
        return m;
    };

    num.m[0][0] = (m.m[1][1] * m.m[2][2] * m.m[3][3] + m.m[1][2] * m.m[2][3] * m.m[3][1] + m.m[1][3] * m.m[2][1] * m.m[3][2] - m.m[1][3] * m.m[2][2] * m.m[3][1] - m.m[1][2] * m.m[2][1] * m.m[3][3] - m.m[1][1] * m.m[2][3] * m.m[3][2]) / determinant;
    num.m[0][1] = (-m.m[0][1] * m.m[2][2] * m.m[3][3] - m.m[0][2] * m.m[2][3] * m.m[3][1] - m.m[0][3] * m.m[2][1] * m.m[3][2] + m.m[0][3] * m.m[2][2] * m.m[3][1] + m.m[0][2] * m.m[2][1] * m.m[3][3] + m.m[0][1] * m.m[2][3] * m.m[3][2]) / determinant;
    num.m[0][2] = (m.m[0][1] * m.m[1][2] * m.m[3][3] + m.m[0][2] * m.m[1][3] * m.m[3][1] + m.m[0][3] * m.m[1][1] * m.m[3][2] - m.m[0][3] * m.m[1][2] * m.m[3][1] - m.m[0][2] * m.m[1][1] * m.m[3][3] - m.m[0][1] * m.m[1][3] * m.m[3][2]) / determinant;
    num.m[0][3] = (-m.m[0][1] * m.m[1][2] * m.m[2][3] - m.m[0][2] * m.m[1][3] * m.m[2][1] - m.m[0][3] * m.m[1][1] * m.m[2][2] + m.m[0][3] * m.m[1][2] * m.m[2][1] + m.m[0][2] * m.m[1][1] * m.m[2][3] + m.m[0][1] * m.m[1][3] * m.m[2][2]) / determinant;

    num.m[1][0] = (-m.m[1][0] * m.m[2][2] * m.m[3][3] - m.m[1][2] * m.m[2][3] * m.m[3][0] - m.m[1][3] * m.m[2][0] * m.m[3][2] + m.m[1][3] * m.m[2][2] * m.m[3][0] + m.m[1][2] * m.m[2][0] * m.m[3][3] + m.m[1][0] * m.m[2][3] * m.m[3][2]) / determinant;
    num.m[1][1] = (m.m[0][0] * m.m[2][2] * m.m[3][3] + m.m[0][2] * m.m[2][3] * m.m[3][0] + m.m[0][3] * m.m[2][0] * m.m[3][2] - m.m[0][3] * m.m[2][2] * m.m[3][0] - m.m[0][2] * m.m[2][0] * m.m[3][3] - m.m[0][0] * m.m[2][3] * m.m[3][2]) / determinant;
    num.m[1][2] = (-m.m[0][0] * m.m[1][2] * m.m[3][3] - m.m[0][2] * m.m[1][3] * m.m[3][0] - m.m[0][3] * m.m[1][0] * m.m[3][2] + m.m[0][3] * m.m[1][2] * m.m[3][0] + m.m[0][2] * m.m[1][0] * m.m[3][3] + m.m[0][0] * m.m[1][3] * m.m[3][2]) / determinant;
    num.m[1][3] = (m.m[0][0] * m.m[1][2] * m.m[2][3] + m.m[0][2] * m.m[1][3] * m.m[2][0] + m.m[0][3] * m.m[1][0] * m.m[2][2] - m.m[0][3] * m.m[1][2] * m.m[2][0] - m.m[0][2] * m.m[1][0] * m.m[2][3] - m.m[0][0] * m.m[1][3] * m.m[2][2]) / determinant;

    num.m[2][0] = (m.m[1][0] * m.m[2][1] * m.m[3][3] + m.m[1][1] * m.m[2][3] * m.m[3][0] + m.m[1][3] * m.m[2][0] * m.m[3][1] - m.m[1][3] * m.m[2][1] * m.m[3][0] - m.m[1][1] * m.m[2][0] * m.m[3][3] - m.m[1][0] * m.m[2][3] * m.m[3][1]) / determinant;
    num.m[2][1] = (-m.m[0][0] * m.m[2][1] * m.m[3][3] - m.m[0][1] * m.m[2][3] * m.m[3][0] - m.m[0][3] * m.m[2][0] * m.m[3][1] + m.m[0][3] * m.m[2][1] * m.m[3][0] + m.m[0][1] * m.m[2][0] * m.m[3][3] + m.m[0][0] * m.m[2][3] * m.m[3][1]) / determinant;
    num.m[2][2] = (m.m[0][0] * m.m[1][1] * m.m[3][3] + m.m[0][1] * m.m[1][3] * m.m[3][0] + m.m[0][3] * m.m[1][0] * m.m[3][1] - m.m[0][3] * m.m[1][1] * m.m[3][0] - m.m[0][1] * m.m[1][0] * m.m[3][3] - m.m[0][0] * m.m[1][3] * m.m[3][1]) / determinant;
    num.m[2][3] = (-m.m[0][0] * m.m[1][1] * m.m[2][3] - m.m[0][1] * m.m[1][3] * m.m[2][0] - m.m[0][3] * m.m[1][0] * m.m[2][1] + m.m[0][3] * m.m[1][1] * m.m[2][0] + m.m[0][1] * m.m[1][0] * m.m[2][3] + m.m[0][0] * m.m[1][3] * m.m[2][1]) / determinant;

    num.m[3][0] = (-m.m[1][0] * m.m[2][1] * m.m[3][2] - m.m[1][1] * m.m[2][2] * m.m[3][0] - m.m[1][2] * m.m[2][0] * m.m[3][1] + m.m[1][2] * m.m[2][1] * m.m[3][0] + m.m[1][1] * m.m[2][0] * m.m[3][2] + m.m[1][0] * m.m[2][2] * m.m[3][1]) / determinant;
    num.m[3][1] = (m.m[0][0] * m.m[2][1] * m.m[3][2] + m.m[0][1] * m.m[2][2] * m.m[3][0] + m.m[0][2] * m.m[2][0] * m.m[3][1] - m.m[0][2] * m.m[2][1] * m.m[3][0] - m.m[0][1] * m.m[2][0] * m.m[3][2] - m.m[0][0] * m.m[2][2] * m.m[3][1]) / determinant;
    num.m[3][2] = (-m.m[0][0] * m.m[1][1] * m.m[3][2] - m.m[0][1] * m.m[1][2] * m.m[3][0] - m.m[0][2] * m.m[1][0] * m.m[3][1] + m.m[0][2] * m.m[1][1] * m.m[3][0] + m.m[0][1] * m.m[1][0] * m.m[3][2] + m.m[0][0] * m.m[1][2] * m.m[3][1]) / determinant;
    num.m[3][3] = (m.m[0][0] * m.m[1][1] * m.m[2][2] + m.m[0][1] * m.m[1][2] * m.m[2][0] + m.m[0][2] * m.m[1][0] * m.m[2][1] - m.m[0][2] * m.m[1][1] * m.m[2][0] - m.m[0][1] * m.m[1][0] * m.m[2][2] - m.m[0][0] * m.m[1][2] * m.m[2][1]) / determinant;

    return num;
}
Vector3 Normalize(const Vector3& v)
{
    float Normalize;
    Vector3 num;
    Normalize = sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
    num.x = v.x / Normalize;
    num.y = v.y / Normalize;
    num.z = v.z / Normalize;
    return num;
}
Vector3 Cross(const Vector3& v1, const Vector3& v2)
{
    Vector3 num;
    num = { v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x };
    return num;
}

Matrix4x4 MakeRotateXMatrix(float radian)
{
    Matrix4x4 num;
    num = { 1, 0, 0, 0,
        0, std::cos(radian), std::sin(radian), 0,
        0, std::sin(-radian), std::cos(radian), 0,
        0, 0, 0, 1 };
    return num;
}
Matrix4x4 MakeRotateYMatrix(float radian)
{
    Matrix4x4 num;
    num = { std::cos(radian), 0, std::sin(-radian), 0,
        0, 1, 0, 0,
        std::sin(radian), 0, std::cos(radian), 0,
        0, 0, 0, 1 };
    return num;
}
Matrix4x4 MakeRotateZMatrix(float radian)
{
    Matrix4x4 num;
    num = { std::cos(radian), std::sin(radian), 0, 0,
        std::sin(-radian), std::cos(radian), 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1 };
    return num;
}

Matrix4x4 MakeAffineMatrix(const Vector3& scale, const Vector3& rotate, const Vector3& translate)
{
    Matrix4x4 rotateX = MakeRotateXMatrix(rotate.x);
    Matrix4x4 rotateY = MakeRotateYMatrix(rotate.y);
    Matrix4x4 rotateZ = MakeRotateZMatrix(rotate.z);
    Matrix4x4 rotateXYZ = Multiply(rotateX, Multiply(rotateY, rotateZ));

    Matrix4x4 num;
    num.m[0][0] = scale.x * rotateXYZ.m[0][0];
    num.m[0][1] = scale.x * rotateXYZ.m[0][1];
    num.m[0][2] = scale.x * rotateXYZ.m[0][2];
    num.m[0][3] = 0.0f * 0.0f * 0.0f * 0.0f;
    num.m[1][0] = scale.y * rotateXYZ.m[1][0];
    num.m[1][1] = scale.y * rotateXYZ.m[1][1];
    num.m[1][2] = scale.y * rotateXYZ.m[1][2];
    num.m[1][3] = 0.0f * 0.0f * 0.0f * 0.0f;
    num.m[2][0] = scale.z * rotateXYZ.m[2][0];
    num.m[2][1] = scale.z * rotateXYZ.m[2][1];
    num.m[2][2] = scale.z * rotateXYZ.m[2][2];
    num.m[2][3] = 0.0f * 0.0f * 0.0f * 0.0f;
    num.m[3][0] = translate.x;
    num.m[3][1] = translate.y;
    num.m[3][2] = translate.z;
    num.m[3][3] = 1.0f;
    return num;
}
Matrix4x4 MakePrespectiveFovMatrix(float fovY, float aspectRatio, float nearClip, float farClip)
{
    Matrix4x4 num;
    num = { (1 / aspectRatio) * (1 / tanf(fovY / 2)), 0, 0, 0, 0, (1 / tanf(fovY / 2)), 0, 0, 0, 0, farClip / (farClip - nearClip), 1, 0, 0, (-nearClip * farClip) / (farClip - nearClip) };
    return num;
}
Matrix4x4 MakeOrthographicMatrix(float left, float top, float right, float bottom, float nearClip, float farClip)
{
    Matrix4x4 num;
    num = { 2 / (right - left), 0, 0, 0, 0, 2 / (top - bottom), 0, 0, 0, 0, 1 / (farClip - nearClip), 0, (left + right) / (left - right),
        (top + bottom) / (bottom - top),
        nearClip / (nearClip - farClip), 1 };
    return num;
}
Matrix4x4 MakeViewportMatrix(float left, float top, float width, float height, float minDepth, float maxDepth)
{
    Matrix4x4 num;
    num = { width / 2, 0, 0, 0, 0, -(height / 2), 0, 0, 0, 0, maxDepth - minDepth, 0, left + (width / 2), top + (height / 2), minDepth, 1 };
    return num;
}

Vector3 operator-(const Vector3& v) { return Vector3(-v.x, -v.y, -v.z); }

Vector3 operator+(const Vector3& v1, const Vector3& v2) { return add(v1, v2); }
Vector3 operator-(const Vector3& v1, const Vector3& v2) { return Subtract(v1, v2); }
Vector3 operator*(float v1, const Vector3& v2) { return Multiply(v1, v2); }
Vector3 operator*(const Vector3& v, float s) { return s * v; }
Vector3 operator/(const Vector3& v, float s) { return Multiply(1.0f / s, v); }
Matrix4x4 operator*(const Matrix4x4& m1, const Matrix4x4& m2) { return Multiply(m1, m2); }
Vector3 operator*(const Vector3& m1, const Vector3& m2) { return Multiply(m1, m2); }

Matrix4x4 MakeRotateAxisAngle(const Vector3& axis, float angle)
{
    Vector3 n = Normalize(axis);
    float cos = cosf(angle);
    float sin = sinf(angle);
    float i = 1.0f - cos;

    Matrix4x4 m {};

    m.m[0][0] = n.x * n.x * i + cos;
    m.m[0][1] = n.x * n.y * i + n.z * sin;
    m.m[0][2] = n.x * n.z * i - n.y * sin;
    m.m[0][3] = 0.0f;

    m.m[1][0] = n.y * n.x * i - n.z * sin;
    m.m[1][1] = n.y * n.y * i + cos;
    m.m[1][2] = n.y * n.z * i + n.x * sin;
    m.m[1][3] = 0.0f;

    m.m[2][0] = n.z * n.x * i + n.y * sin;
    m.m[2][1] = n.z * n.y * i - n.x * sin;
    m.m[2][2] = n.z * n.z * i + cos;
    m.m[2][3] = 0.0f;

    m.m[3][0] = m.m[3][1] = m.m[3][2] = 0.0f;
    m.m[3][3] = 1.0f;

    return m;
}

Matrix4x4 DirectionToDirection(const Vector3& from, const Vector3& to)
{
    // 1 c*v
    //  fromとtoのクロス積の結果を入れたvector3の変数を作る

    // 2 cos
    //  floatの変数に、uとvの内積を入れる

    // 3 sin
    //  floatの変数に、uとvのクロス積の長さを入れる

    // 4
    //  if文を作成
    //  反転行列なら(u.x ≠ 0 || u.y ≠0の場合)と(u.x ≠ 0 || u.z ≠0の場合)でif文を書く
    //  どちらも当てはならないならassretではじく
    //  反転行列じゃないなら
    // 1の求めたクロス積を正規化下内容がnになる。

    Vector3 u = from;
    Vector3 v = to;

    Matrix4x4 m {};

    Vector3 cross = Cross(u, v);
    float cos = Dot(u, v);
    float sin = Length(Cross(u, v));
    float i = 1.0f - cos;

    if (cos == -1.0f) {
        if (u.x != 0 || u.y != 0) {
            Vector3 axis;
            axis = { u.y, -u.x, 0.0f };

            axis = Normalize(axis);

            Vector3 n = axis;

            m.m[0][0] = (n.x * n.x) * i + cos;
            m.m[0][1] = (n.x * n.y) * i + n.z * sin;
            m.m[0][2] = (n.x * n.z) * i - n.y * sin;
            m.m[0][3] = 0.0f;

            m.m[1][0] = n.y * n.x * i - n.z * sin;
            m.m[1][1] = n.y * n.y * i + cos;
            m.m[1][2] = n.y * n.z * i + n.x * sin;
            m.m[1][3] = 0.0f;

            m.m[2][0] = n.z * n.x * i + n.y * sin;
            m.m[2][1] = n.z * n.y * i - n.x * sin;
            m.m[2][2] = n.z * n.z * i + cos;
            m.m[2][3] = 0.0f;

            m.m[3][0] = m.m[3][1] = m.m[3][2] = 0.0f;
            m.m[3][3] = 1.0f;

            return m;
        } else if (u.z != 0 || u.x != 0) {
            Vector3 axis;
            axis = { u.z, 0.0f, -u.x };

            axis = Normalize(axis);

            Vector3 n = axis;

            m.m[0][0] = (n.x * n.x) * i + cos;
            m.m[0][1] = (n.x * n.y) * i + n.z * sin;
            m.m[0][2] = (n.x * n.z) * i - n.y * sin;
            m.m[0][3] = 0.0f;

            m.m[1][0] = n.y * n.x * i - n.z * sin;
            m.m[1][1] = n.y * n.y * i + cos;
            m.m[1][2] = n.y * n.z * i + n.x * sin;
            m.m[1][3] = 0.0f;

            m.m[2][0] = n.z * n.x * i + n.y * sin;
            m.m[2][1] = n.z * n.y * i - n.x * sin;
            m.m[2][2] = n.z * n.z * i + cos;
            m.m[2][3] = 0.0f;

            m.m[3][0] = m.m[3][1] = m.m[3][2] = 0.0f;
            m.m[3][3] = 1.0f;
            return m;
        }
    }

    Vector3 n = Normalize(cross);

    m.m[0][0] = (n.x * n.x) * i + cos;
    m.m[0][1] = (n.x * n.y) * i + n.z * sin;
    m.m[0][2] = (n.x * n.z) * i - n.y * sin;
    m.m[0][3] = 0.0f;

    m.m[1][0] = n.y * n.x * i - n.z * sin;
    m.m[1][1] = n.y * n.y * i + cos;
    m.m[1][2] = n.y * n.z * i + n.x * sin;
    m.m[1][3] = 0.0f;

    m.m[2][0] = n.z * n.x * i + n.y * sin;
    m.m[2][1] = n.z * n.y * i - n.x * sin;
    m.m[2][2] = n.z * n.z * i + cos;
    m.m[2][3] = 0.0f;

    m.m[3][0] = m.m[3][1] = m.m[3][2] = 0.0f;
    m.m[3][3] = 1.0f;

    return m;
}

const char kWindowTitle[] = "LE2C_03_イケダ_タクミ";

// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int)
{

    // ライブラリの初期化
    Novice::Initialize(kWindowTitle, 1280, 720);

    // キー入力結果を受け取る箱
    char keys[256] = { 0 };
    char preKeys[256] = { 0 };

    Vector3 from0 = Normalize(Vector3 { 1.0f, 0.7f, 0.5f });
    Vector3 to0 = -from0;
    Vector3 from1 = Normalize(Vector3 { -0.6f, 0.9f, 0.2f });
    Vector3 to1 = Normalize(Vector3 { 0.4f, 0.7f - 0.5f });
    Matrix4x4 rotateMatrix0 = DirectionToDirection(Normalize(Vector3 { 1.0f, 0.0f, 0.0f }), Normalize(Vector3 { -1.0f, 0.0f, 0.0f }));
    Matrix4x4 rotateMatrix1 = DirectionToDirection(from0, to0);
    Matrix4x4 rotateMatrix2 = DirectionToDirection(from1, to1);

    // ウィンドウの×ボタンが押されるまでループ
    while (Novice::ProcessMessage() == 0) {
        // フレームの開始
        Novice::BeginFrame();

        // キー入力を受け取る
        memcpy(preKeys, keys, 256);
        Novice::GetHitKeyStateAll(keys);

        ///
        /// ↓更新処理ここから
        ///

        ///
        /// ↑更新処理ここまで
        ///

        ///
        /// ↓描画処理ここから
        ///

        MatrixScreenPrintf(0, 0, rotateMatrix0, "rotateMatrix0");
        MatrixScreenPrintf(0, krowheight * 5, rotateMatrix1, "rotateMatrix1");
        MatrixScreenPrintf(0, krowheight * 10, rotateMatrix2, "rotateMatrix2");

        ///
        /// ↑描画処理ここまで
        ///

        // フレームの終了
        Novice::EndFrame();

        // ESCキーが押されたらループを抜ける
        if (preKeys[DIK_ESCAPE] == 0 && keys[DIK_ESCAPE] != 0) {
            break;
        }
    }

    // ライブラリの終了
    Novice::Finalize();
    return 0;
}
