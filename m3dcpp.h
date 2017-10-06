#ifndef M3DCPP_LIBRARY_H
#define M3DCPP_LIBRARY_H

#include <string>
#include <sstream>
#include <cmath>
#include <array>
#include <experimental/array>

namespace M3D {
    template<int comp, typename T>
    struct Vector {

        Vector<comp,T>(){}
       // Vector<comp,T>(std::array<T,comp> arr) : m{arr}{};
        //Vector<comp,T>(const T (&arr)[comp]): m{std::experimental::to_array(arr)}{};


        Vector<4,float>(float x,float y, float z, float w) {
            m[0] = x;
            m[1] = y;
            m[2] = z;
            m[3] = w;
        };
        Vector<3,float>(float x,float y, float z) {
            m[0] = x;
            m[1] = y;
            m[2] = z;
        };
        Vector<2,float>(float x,float y) {
            m[0] = x;
            m[1] = y;
        };
/*
        Vector<4,int>(int x,int y, int z, int w) {
            m[0] = x;
            m[1] = y;
            m[2] = z;
            m[3] = w;
        };
        Vector<3,int>(int x,int y, int z) {
            m[0] = x;
            m[1] = y;
            m[2] = z;
        };
        Vector<2,int>(int x,int y) {
            m[0] = x;
            m[1] = y;
        };
*/
        Vector<comp,T>(T *arr){
            bool eoa = false;
            for(int i = 0; i < comp; i++){
                if(eoa){
                    m[i] = 0;
                }else {
                    eoa = (arr == nullptr);
                    m[i] = *arr++;
                }
            }
        };


        std::array<T,comp> m = {0};

        Vector<comp, T> operator+(const Vector<comp, T> &other) {
            Vector<comp, T> res;
            for (int i = 0; i < comp; ++i) {
                res.m[i] = this->m[i] + other.m[i];
            }
            return res;
        };
        Vector<comp, T> operator-(const Vector<comp, T> &other) {
            Vector<comp, T> res;
            for (int i = 0; i < comp; ++i) {
                res.m[i] = this->m[i] - other.m[i];
            }
            return res;
        };
        Vector<comp, T> operator*(const Vector<comp, T> &other) {
            Vector<comp, T> res;
            for (int i = 0; i < comp; ++i) {
                res.m[i] = this->m[i] * other.m[i];
            }
            return res;
        };
        Vector<comp, T> operator/(const Vector<comp, T> &other) {
            Vector<comp, T> res;
            for (int i = 0; i < comp; ++i) {
                if (other[i] == 0) return Vector<comp, T>();
                res.m[i] = (*this)[i] / other[i];
            }
            return res;
        };
        Vector<comp, T> operator+(float other) {
            Vector<comp, T> res;
            for (int i = 0; i < comp; ++i) {
                res.m[i] = this->m[i] + other;
            }
            return res;
        };
        Vector<comp, T> operator-(float other) {
            Vector<comp, T> res;
            for (int i = 0; i < comp; ++i) {
                res.m[i] = this->m[i] - other;
            }
            return res;
        };
        Vector<comp, T> operator*(float other) {
            Vector<comp, T> res;
            for (int i = 0; i < comp; ++i) {
                res.m[i] = this->m[i] * other;
            }
            return res;
        };
        Vector<comp, T> operator/(float other) {
            Vector<comp, T> res;
            if (other == 0) return res;
            for (int i = 0; i < comp; ++i) {
                res.m[i] = this->m[i] / other;
            }
            return res;
        };


        T &operator[](int at) {
            return m[at];
        };
        T operator[](int at) const {
            return m[at];
        }
        T &operator[](char at) {
            switch (at) {
                case 'x':
                case 'r':
                case 'u':
                    return m[0];
                case 'y':
                case 'g':
                case 'v':
                    return m[1];
                case 'z':
                case 'b':
                    return m[2];
                case 'w':
                case 'a':
                    return m[3];
                default:
                    return m[0];
            }
        }
        T operator[](char at) const {
            switch (at) {
                case 'x':
                case 'r':
                case 'u':
                    return m[0];
                case 'y':
                case 'g':
                case 'v':
                    return m[1];
                case 'z':
                case 'b':
                    return m[2];
                case 'w':
                case 'a':
                    return m[3];
                default:
                    return m[0];
            }
        }

        std::string ToString() const {
            std::stringstream ss;
            ss << "(";
            for (int i = 0; i < comp; i++) {
                if (i == comp - 1) ss << m[i] << ")\n";
                else
                    ss << m[i] << ",";
            }
            return ss.str();
        }

        template<int n, typename tT>
        friend std::ostream &operator<<(std::ostream &os, const Vector<n, tT> &vec);

        template<int n, typename tT>
        float Dot(const Vector<n, tT> &other) const {
            if (n != comp)
            {
                return 0;
            }else{
                float sum = 0;
                for (int i = 0; i < n; ++i) {
                    sum += m[i] * other[i];
                }
                return sum;
            }
        }
        Vector<3, T> Cross(const Vector<3, T> &other) const {
            Vector<3, T> res;
            res[0] = (*this)[1] * other[2] - (*this)[2] * other[1];
            res[1] = (*this)[2] * other[0] - (*this)[0] * other[2];
            res[2] = (*this)[0] * other[1] - (*this)[1] * other[0];
            return res;
        };

        float Length() const {
            return sqrtf((*this).Dot(*this));
        }

        void Normalize() {
            float len = Length();
            if (len == 0) return;
            *this = *this / len;

        }
        Vector<comp, T> Normalize_Copy() const {
            Vector<comp, T> res;
            float len = Length();
            if (len == 0) return res;
            for (int i = 0; i < comp; i++) {
                res[i] = m[i] / len;
            }
            return res;
        }

        constexpr int Size() const {
            return comp;
        }


    };
    template<int comp, typename T>
    struct Matrix {
        T m[comp][comp] = {0};

        Vector<comp, T> GetColumn(int col) {
            Vector<comp, T> vec;
            for (int i = 0; i < comp; i++) {
                vec[i] = m[i][col];
            }
            return vec;
        };

        Vector<comp, T> GetRow(int row) {
            Vector<comp, T> vec;
            for (int i = 0; i < comp; i++) {
                vec[i] = m[row][i];
            }
            return vec;
        };

        void MakeIdentity(T diagonal) {
            for (int i = 0; i < comp; i++) {
                m[i][i] = diagonal;
            }
        };

        void MakeTranslation(Vector<3, float> offset) {
            Identity(1);
            m[0][3] = offset[0];
            m[1][3] = offset[1];
            m[2][3] = offset[2];
        };

        void MakeRotation(float angle_in_rad, Vector<3, float> axis) {
            Matrix<4, float> res = Matrix<4, float>::Rotation(angle_in_rad, axis);
            *this = res;
        };

        void MakeScaling(Vector<3, float> scale) {
            Identity(1);
            for (int i = 0; i < 3; ++i) {
                m[i][i] = scale[i];

            }
        };

        void MakePerspective(float fovy, float aspect, float near, float far) {
            Matrix<4, float> res = Matrix<4, float>::Perspective(fovy, aspect, near, far);
            *this = res;
        };

        void MakeLookAt(Vector<3, float> from, Vector<3, float> to, Vector<3, float> up) {
            Matrix<4, float> res = Matrix<4, float>::LookAt(from, to, up);
            *this = res;
        };


        static Matrix<comp, T> Identity(T diagonal) {
            Matrix<comp, T> matrix;
            for (int i = 0; i < comp; i++) {
                matrix.m[i][i] = diagonal;
            }
            return matrix;
        };

        static Matrix<4, float> Translation(Vector<3, float> offset) {
            Matrix<4, float> res = {
                    1, 0, 0, offset[0],
                    0, 1, 0, offset[1],
                    0, 0, 1, offset[2],
                    0, 0, 0, 1
            };
            return res;
        };

        static Matrix<4, float> Rotation(float angle_in_rad, Vector<3, float> axis) {
            Vector<3, float> normalized_axis = axis.Normalize_Copy();
            float x = normalized_axis[0], y = normalized_axis[1], z = normalized_axis[2];
            float c = cosf(angle_in_rad), s = sinf(angle_in_rad);

            return Matrix<4, float>{
                    c + x * x * (1 - c), x * y * (1 - c) - z * s, x * z * (1 - c) + y * s, 0,
                    y * x * (1 - c) + z * s, c + y * y * (1 - c), y * z * (1 - c) - x * s, 0,
                    z * x * (1 - c) - y * s, z * y * (1 - c) + x * s, c + z * z * (1 - c), 0,
                    0, 0, 0, 1
            };

        };

        static Matrix<4, float> Scaling(Vector<3, float> scale) {
            Matrix<4, float> res = Matrix::Identity(1);

            for (int i = 0; i < 3; ++i) {
                res.m[i][i] = scale[i];

            }

            return res;
        };

        static Matrix<4, float> Perspective(float fovy, float aspect, float near, float far) {
            float f = 1.0f / tanf(fovy / 2.0f);
            float ar = aspect;
            float nd = near, fd = far;

            return Matrix<4, float>{
                    f / ar, 0, 0, 0,
                    0, f, 0, 0,
                    0, 0, (fd + nd) / (nd - fd), (2 * fd * nd) / (nd - fd),
                    0, 0, -1, 0
            };
        };

        static Matrix<4, float> Ortho(float top, float bottom, float left, float right, float near, float far) {
            float l = left, r = right, b = bottom, t = top, n = near, f = far;
            float tx = -(r + l) / (r - l);
            float ty = -(t + b) / (t - b);
            float tz = -(f + n) / (f - n);
            return Matrix<4, float>{
                    2 / (r - l), 0, 0, tx,
                    0, 2 / (t - b), 0, ty,
                    0, 0, 2 / (f - n), tz,
                    0, 0, 0, 1
            };
        };

        static Matrix<4, float> LookAt(Vector<3, float> from, Vector<3, float> to, Vector<3, float> up) {
            Vector<3, float> zz = (to - from);
            Vector<3, float> z = zz.Normalize_Copy() * -1;
            Vector<3, float> x = up.Cross(z).Normalize_Copy();
            Vector<3, float> y = z.Cross(x);

            return Matrix<4, float>{
                    x[0], x[1], x[2], -from.Dot(x),
                    y[0], y[1], y[2], -from.Dot(y),
                    z[0], z[1], z[2], -from.Dot(z),
                    0, 0, 0, 1
            };
        };

        template<int n, typename tT>
        friend std::ostream &operator<<(std::ostream &os, const Matrix<n, tT> &mat);


        Matrix<comp, T> operator+(const Matrix<comp, T> &other) {
            Matrix<comp, T> res;
            for (int i = 0; i < comp; ++i) {
                for (int j = 0; j < comp; ++j) {
                    res.m[i][j] = m[i][j] + other.m[i][j];
                }
            }
            return res;
        };

        Matrix<comp, T> operator-(const Matrix<comp, T> &other) {
            Matrix<comp, T> res;
            for (int i = 0; i < comp; ++i) {
                for (int j = 0; j < comp; ++j) {
                    res.m[i][j] = m[i][j] - other.m[i][j];
                }
            }
            return res;
        };

        Matrix<comp, T> operator*(const Matrix<comp, T> &other) {

            Matrix<comp, T> res;
            for (int i = 0; i < comp; i++) {
                for (int j = 0; j < comp; j++) {
                    float sum = 0;
                    for (int k = 0; k < comp; k++) {
                        sum += m[k][j] * other.m[i][k];
                    }
                    res.m[i][j] = sum;
                }
            }
            return res;

        };

        Vector<4, float> operator*(const Vector<4, float> &vec) {
            Vector<4, float> res;
            res = GetColumn(0) * vec[0] + GetColumn(1) * vec[1] + GetColumn(2) * vec[2] + GetColumn(3) * vec[3];
            return res;
        };

        Vector<3, float> operator*(const Vector<3, float> &vec) {
            Vector<3, float> res;
            res = GetColumn(0) * vec[0] + GetColumn(1) * vec[1] + GetColumn(2) * vec[2];
            return res;
        };

        Matrix<comp, T> operator*(T val) {
            Matrix<comp, T> res;
            for (int i = 0; i < comp; ++i) {
                for (int j = 0; j < comp; ++j) {
                    res.m[i][j] = m[i][j] * val;
                }
            }
        };

        Matrix<comp, T> operator+(T val) {
            Matrix<comp, T> res;
            for (int i = 0; i < comp; ++i) {
                for (int j = 0; j < comp; ++j) {
                    res.m[i][j] = m[i][j] + val;
                }
            }
        };

        Matrix<comp, T> operator-(T val) {
            Matrix<comp, T> res;
            for (int i = 0; i < comp; ++i) {
                for (int j = 0; j < comp; ++j) {
                    res.m[i][j] = m[i][j] - val;
                }
            }
        };

        Matrix<comp, T> operator/(T val) {
            Matrix<comp, T> res;
            for (int i = 0; i < comp; ++i) {
                for (int j = 0; j < comp; ++j) {
                    res.m[i][j] = m[i][j] / val;
                }
            }
        };

    };


    typedef Matrix<4, float> Mat4f;
    typedef Matrix<3, float> Mat3f;
    typedef Vector<4, float> Vec4f;
    typedef Vector<3, float> Vec3f;
    typedef Vector<2, float> Vec2f;
    typedef Vector<4, int> Vec4i;
    typedef Vector<3, int> Vec3i;
    typedef Vector<2, int> Vec2i;

    template<int n, typename tT>
    std::ostream &operator<<(std::ostream &os, const M3D::Vector<n, tT> &vec) {
        os << vec.ToString();
        return os;
    };

    template<int comp, typename T>
    std::ostream &operator<<(std::ostream &os, const M3D::Matrix<comp, T> &mat) {
        for (int i = 0; i < comp; i++) {
            os << "|";
            for (int j = 0; j < comp; j++) {
                if (j == comp - 1)
                    os << mat.m[i][j];
                else
                    os << mat.m[i][j] << ",";
            }
            os << "|\n";
        }
        return os;
    };


    bool InRange(float value, float min, float max);
    float Map(float value, float in_min, float in_max, float out_min, float out_max);
    float Constrain(float value, float min, float max);
    float InterpolateLinear(float start, float end, float factor);

}



#endif











