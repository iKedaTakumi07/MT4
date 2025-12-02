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
struct Quaternion {
    float x;
    float y;
    float z;
    float w;
};

static const int krowheight = 20;
static const int kcolumnwidth = 60;

static const int kcolumnwith = 60;
void VectorScreenPrintf(int x, int y, const Vector3& vector, const char* label)
{
    Novice::ScreenPrintf(x, y, "%0.2f", vector.x);
    Novice::ScreenPrintf(x + kcolumnwith, y, "%0.2f", vector.y);
    Novice::ScreenPrintf(x + kcolumnwith * 2, y, "%0.2f", vector.z);

    Novice::ScreenPrintf(x + kcolumnwith * 4, y, "%s", label);
}
void QuaternionScreenPrintf(int x, int y, const Quaternion& vector, const char* label)
{
    Novice::ScreenPrintf(x, y, "%0.2f", vector.x);
    Novice::ScreenPrintf(x + kcolumnwith, y, "%0.2f", vector.y);
    Novice::ScreenPrintf(x + kcolumnwith * 2, y, "%0.2f", vector.z);
    Novice::ScreenPrintf(x + kcolumnwith * 3, y, "%0.2f", vector.w);
    Novice::ScreenPrintf(x + kcolumnwith * 4, y, "%s", label);
}
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
Quaternion add(const Quaternion& v1, const Quaternion& v2)
{
    Quaternion num;
    num.x = v1.x + v2.x;
    num.y = v1.y + v2.y;
    num.z = v1.z + v2.z;
    num.w = v1.w + v1.w;
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
float Dot(const Quaternion& a, const Quaternion& b)
{
    float num;
    return num = a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
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

Quaternion IdentityQuaternion()
{
    return { 0.0f, 0.0f, 0.0f, 1.0f };
}
Quaternion Conjugation(const Quaternion& v)
{
    return { -v.x, -v.y, -v.z, v.w };
}
float Norm(const Quaternion& v)
{
    float norm;
    norm = sqrtf(v.w * v.w + v.x * v.x + v.y * v.y + v.z * v.z);

    return norm;
}
Quaternion Normalize(const Quaternion& v)
{
    float norm;
    Quaternion num;
    norm = sqrtf(v.x * v.x + v.y * v.y + v.z * v.z + v.w * v.w);
    num.x = v.x / norm;
    num.y = v.y / norm;
    num.z = v.z / norm;
    num.w = v.w / norm;

    return num;
}
Quaternion Inverse(const Quaternion& v)
{
    float normSq = v.x * v.x + v.y * v.y + v.z * v.z + v.w * v.w;
    if (normSq == 0.0f) {
        // アサ―と
        assert("0.0f!!!!!");
    }

    Quaternion result;
    result.x = -v.x / normSq;
    result.y = -v.y / normSq;
    result.z = -v.z / normSq;
    result.w = v.w / normSq;
    return result;
}
Quaternion Multiply(const Quaternion& v1, const Quaternion& v2)
{
    Quaternion result;
    result.x = v1.w * v2.x + v1.x * v2.w + v1.y * v2.z - v1.z * v2.y;
    result.y = v1.w * v2.y - v1.x * v2.z + v1.y * v2.w + v1.z * v2.x;
    result.z = v1.w * v2.z + v1.x * v2.y - v1.y * v2.x + v1.z * v2.w;
    result.w = v1.w * v2.w - v1.x * v2.x - v1.y * v2.y - v1.z * v2.z;
    return result;
}

Quaternion operator-(const Quaternion& v) { return Quaternion(-v.x, -v.y, -v.z, -v.z); }
Quaternion operator*(const Quaternion& m1, const Quaternion& m2) { return Multiply(m1, m2); }
Quaternion operator+(const Quaternion& v1, const Quaternion& v2) { return add(v1, v2); }
Vector3 operator-(const Vector3& v) { return Vector3(-v.x, -v.y, -v.z); }
Vector3 operator+(const Vector3& v1, const Vector3& v2) { return add(v1, v2); }
Vector3 operator-(const Vector3& v1, const Vector3& v2) { return Subtract(v1, v2); }
Vector3 operator*(float v1, const Vector3& v2) { return Multiply(v1, v2); }
Vector3 operator*(const Vector3& v, float s) { return s * v; }
Vector3 operator/(const Vector3& v, float s) { return Multiply(1.0f / s, v); }
Matrix4x4 operator*(const Matrix4x4& m1, const Matrix4x4& m2) { return Multiply(m1, m2); }
Vector3 operator*(const Vector3& m1, const Vector3& m2) { return Multiply(m1, m2); }

// 任意軸回転を表すQuaternionの生成
Quaternion MakeRotateAxisAngleQuaternion(const Vector3& axis, float angle)
{
    float sin = std::sinf(angle / 2.0f);
    float cos = std::cosf(angle / 2.0f);
    Quaternion num;

    num = { axis.x * sin, axis.y * sin, axis.z * sin, cos };

    return num;
}
// ベクトルをQuaternionで回転させた結果のベクトルを求める
Vector3 RotateVector(const Vector3& vector, const Quaternion& quaternion)
{
    Quaternion qv = { vector.x, vector.y, vector.z, 0 };
    Quaternion qc = Conjugation(quaternion);

    Quaternion t;
    t.x = quaternion.w * qv.x + quaternion.x * qv.w + quaternion.y * qv.z - quaternion.z * qv.y;
    t.y = quaternion.w * qv.y - quaternion.x * qv.z + quaternion.y * qv.w + quaternion.z * qv.x;
    t.z = quaternion.w * qv.z + quaternion.x * qv.y - quaternion.y * qv.x + quaternion.z * qv.w;
    t.w = quaternion.w * qv.w - quaternion.x * qv.x - quaternion.y * qv.y - quaternion.z * qv.z;

    Quaternion result;
    result.w = t.w * qc.w - t.x * qc.x - t.y * qc.y - t.z * qc.z;
    result.x = t.w * qc.x + t.x * qc.w + t.y * qc.z - t.z * qc.y;
    result.y = t.w * qc.y - t.x * qc.z + t.y * qc.w + t.z * qc.x;
    result.z = t.w * qc.z + t.x * qc.y - t.y * qc.x + t.z * qc.w;

    return { result.x, result.y, result.z };
}
// Quaternionから回転行列を求める
Matrix4x4 MakeRotateMatrix(const Quaternion& quaternion)
{
    Matrix4x4 num;
    num.m[0][0] = quaternion.w * quaternion.w + quaternion.x * quaternion.x - quaternion.y * quaternion.y - quaternion.z * quaternion.z;
    num.m[0][1] = (quaternion.x * quaternion.y + quaternion.w * quaternion.z) * 2.0f;
    num.m[0][2] = (quaternion.x * quaternion.z - quaternion.w * quaternion.y) * 2.0f;
    num.m[0][3] = 0.0f;

    num.m[1][0] = (quaternion.x * quaternion.y - quaternion.w * quaternion.z) * 2.0f;
    num.m[1][1] = quaternion.w * quaternion.w - quaternion.x * quaternion.x + quaternion.y * quaternion.y - quaternion.z * quaternion.z;
    num.m[1][2] = (quaternion.y * quaternion.z + quaternion.w * quaternion.x) * 2.0f;
    num.m[1][3] = 0.0f;

    num.m[2][0] = (quaternion.x * quaternion.z + quaternion.w * quaternion.y) * 2.0f;
    num.m[2][1] = (quaternion.x * quaternion.z - quaternion.w * quaternion.x) * 2.0f;
    num.m[2][2] = quaternion.w * quaternion.w - quaternion.x * quaternion.x - quaternion.y * quaternion.y + quaternion.z * quaternion.z;
    num.m[2][3] = 0.0f;

    num.m[3][0] = 0.0f;
    num.m[3][1] = 0.0f;
    num.m[3][2] = 0.0f;
    num.m[3][3] = 1.0f;

    return num;
}

// 球面線形補間
Quaternion Slerp(const Quaternion& q0, const Quaternion& q1, float t)
{
    float dot = Dot(q0, q1); // Q0とq1の内積
    Quaternion q;

    q = q0;

    if (dot < 0) {
        q = -q0; // もう片方の回転を利用する
        dot = -dot;
    }

    // なす角を求める
    float theta = std::acos(dot);
    float sinTheta = std::sin(theta);

    // thetaとsinを使って線形保管Scale0,scale1を求める

    // ここで求めている係数をそのまま計算
    float s0 = std::sin((1.0f - t) * theta) / sinTheta; // sin((1-t)θ)/sinθ
    float s1 = std::sin(t * theta) / sinTheta; // sin(tθ)/sinθ

    // それぞれの補間係数を利用して補間後のQuaternionを求める
    Quaternion num; 

    num = { s0 * q.x + s1 * q1.x,
        s0 * q.y + s1 * q1.y,
        s0 * q.z + s1 * q1.z,
        s0 * q.w + s1 * q1.w };

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
Vector3 Transform(const Vector3& vector, const Matrix4x4& matrix)
{
    Vector3 result;
    result.x = vector.x * matrix.m[0][0] + vector.y * matrix.m[1][0] + vector.z * matrix.m[2][0] + 1.0f * matrix.m[3][0];
    result.y = vector.x * matrix.m[0][1] + vector.y * matrix.m[1][1] + vector.z * matrix.m[2][1] + 1.0f * matrix.m[3][1];
    result.z = vector.x * matrix.m[0][2] + vector.y * matrix.m[1][2] + vector.z * matrix.m[2][2] + 1.0f * matrix.m[3][2];
    float w = vector.x * matrix.m[0][3] + vector.y * matrix.m[1][3] + vector.z * matrix.m[2][3] + 1.0f * matrix.m[3][3];
    assert(w != 0.0f);
    result.x /= w;
    result.y /= w;
    result.z /= w;

    return result;
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

    Vector3 u = Normalize(from);
    Vector3 v = Normalize(to);

    Matrix4x4 m {};

    Vector3 cross = Normalize(Cross(u, v));
    float cos = Dot(u, v);
    float sin = Length(Cross(u, v));
    float i = 1.0f - cos;

    if (cos == -1.0f) {
        if (u.x != 0.0f || u.y != 0.0f) {
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
        } else {
            assert(cos);
        }
    }

    m.m[0][0] = (cross.x * cross.x) * i + cos;
    m.m[0][1] = (cross.x * cross.y) * i + cross.z * sin;
    m.m[0][2] = (cross.x * cross.z) * i - cross.y * sin;
    m.m[0][3] = 0.0f;

    m.m[1][0] = (cross.y * cross.x) * i - cross.z * sin;
    m.m[1][1] = (cross.y * cross.y) * i + cos;
    m.m[1][2] = (cross.y * cross.z) * i + cross.x * sin;
    m.m[1][3] = 0.0f;

    m.m[2][0] = (cross.z * cross.x) * i + cross.y * sin;
    m.m[2][1] = (cross.z * cross.y) * i - cross.x * sin;
    m.m[2][2] = (cross.z * cross.z) * i + cos;
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

    Quaternion rotation0 = MakeRotateAxisAngleQuaternion({ 0.71f, 0.71f, 0.0f }, 0.3f);
    Quaternion rotation1 = MakeRotateAxisAngleQuaternion({ 0.71f, 0.0f, 0.71f }, 3.141692f);

    Quaternion interpolate0 = Slerp(rotation0, rotation1, 0.0f);
    Quaternion interpolate1 = Slerp(rotation0, rotation1, 0.3f);
    Quaternion interpolate2 = Slerp(rotation0, rotation1, 0.5f);
    Quaternion interpolate3 = Slerp(rotation0, rotation1, 0.7f);
    Quaternion interpolate4 = Slerp(rotation0, rotation1, 1.0f);

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

        QuaternionScreenPrintf(0, krowheight * 0, interpolate0, ": interpolate0, Slerp(rotation0, rotation1, 0.0f);");
        QuaternionScreenPrintf(0, krowheight * 1, interpolate1, ": interpolate0, Slerp(rotation0, rotation1, 0.3f);");
        QuaternionScreenPrintf(0, krowheight * 2, interpolate2, ": interpolate0, Slerp(rotation0, rotation1, 0.5f);");
        QuaternionScreenPrintf(0, krowheight * 3, interpolate3, ": interpolate0, Slerp(rotation0, rotation1, 0.7f);");
        QuaternionScreenPrintf(0, krowheight * 4, interpolate4, ": interpolate0, Slerp(rotation0, rotation1, 1.0f);");
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
