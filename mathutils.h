/**
 *  \brief  Mathematical utility functions.
 *
 *  \file   mathutils.h
 *
 *  \author Michal Butterweck
 *  \date   2015
 *
 *  \copyright Copyright (c) INNTEC.PL Sp. z o.o. All rights reserved.
 *             Prawa autorskie zarezerwowane - INNTEC.PL Sp. z o.o.
 *
 *  \note INNTEC.PL Sp. z o.o. Company Confidential Information.
 *        Informacje Poufne firmy INNTEC.PL Sp. z o.o.
 *
 */

#ifndef __common_mathutils_h
#define __common_mathutils_h

#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <float.h>
#include <stdio.h>
#include <string.h>
#include "error.h"

// types

/** \struct float2_t
 *
 *  \brief Basic 2-vector data structure.
 */
typedef struct float2_t
{
    union {
        struct {
            float x;    ///< x-component
            float y;    ///< y-component
        };
        struct {
            float s;    ///< s-component
            float t;    ///< t-component
        };
    };
} float2_t;

/** \struct float3_t
 *
 *  \brief Basic 3-vector data structure.
 */
typedef struct float3_t
{
    float x;    ///< x-component
    float y;    ///< y-component
    float z;    ///< z-component
} float3_t;

/** \struct float4_t
 *
 *  \brief Basic 4-vector data structure.
 */
typedef struct float4_t
{
    union
    {
        struct
        {
            float x;    ///< x-component
            float y;    ///< y-component
            float z;    ///< z-component
            float w;    ///< w-component
        };
        struct
        {
            float r;    ///< red-color-component
            float g;    ///< green-color-component
            float b;    ///< blue-color-component
            float a;    ///< alpha-transparency-component (pre-multiplied or not)
        };
        struct
        {
            float3_t v; ///< float3 vector
            float    t; ///< w-component
        };
    };
} float4_t;

/** \struct double2_t
 *
 *  \brief Double precision 2-vector data structure.
 */
typedef struct double2_t
{
    union {
        struct {
            double x;    ///< x-component
            double y;    ///< y-component
        };
        struct {
            double s;    ///< s-component
            double t;    ///< t-component
        };
    };
} double2_t;

/** \struct double3_t
 *
 *  \brief Double precision 3-vector data structure.
 */
typedef struct double3_t
{
    double x;    ///< x-component
    double y;    ///< y-component
    double z;    ///< z-component
} double3_t;

typedef struct double4_t
{
    union
    {
        struct
        {
            double x;    ///< x-component
            double y;    ///< y-component
            double z;    ///< z-component
            double w;    ///< w-component
        };
        struct
        {
            double r;    ///< red-color-component
            double g;    ///< green-color-component
            double b;    ///< blue-color-component
            double a;    ///< alpha-transparency-component (pre-multiplied or not)
        };
    };
} double4_t;

/** @struct float44_t
 *
 *  @brief Basic 4x4 matrix structure.
 */
typedef struct float44_t
{
    union
    {
        float m[4][4];                  ///< array[4][4] of componets
        struct
        {
            float m00, m01, m02, m03;
            float m10, m11, m12, m13;
            float m20, m21, m22, m23;
            float m30, m31, m32, m33;
        };
    };
} float44_t;

typedef struct double44_t
{
    union
    {
        double m[4][4];                  ///< array[4][4] of componets
        struct
        {
            double m00, m01, m02, m03;
            double m10, m11, m12, m13;
            double m20, m21, m22, m23;
            double m30, m31, m32, m33;
        };
    };
} double44_t;

typedef struct float33_t {
    union {
        float m[3][3];

        struct
        {
            float m00, m01, m02;
            float m10, m11, m12;
            float m20, m21, m22;
        };
    };
} float33_t;

typedef struct int2_t {
    int x;
    int y;
} int2_t;

static inline int2_t int2_from_float3_xy(float3_t a) { return (int2_t) {(int)a.x, (int)a.y}; }
static inline int2_t int2_from_float3_xz(float3_t a) { return (int2_t) {(int)a.x, (int)a.z}; }

typedef struct int3_t {
    int x;
    int y;
    int z;
} int3_t;

static inline int3_t int3_from_float3(float3_t a) { return (int3_t) {(int)a.x, (int)a.y, (int)a.z}; }

// declarations
static float3_t  const zero3_sc       = { 0.f, 0.f, 0.f };
static float4_t  const zero4_sc       = { 0.f, 0.f, 0.f, 0.f };
static float4_t  const zero4_point_sc = { 0.f, 0.f, 0.f, 1.f };
static float44_t const identity_sc    = { 1.f, 0.f, 0.f, 0.f,
                                          0.f, 1.f, 0.f, 0.f,
                                          0.f, 0.f, 1.f, 0.f,
                                          0.f, 0.f, 0.f, 1.f };

static float33_t const identity33_sc  = { 1.f, 0.f, 0.f,
                                          0.f, 1.f, 0.f,
                                          0.f, 0.f, 1.f };


// misc
#define EPSILON 0.000000954 /* near 1e-6 but as a division by 2^n */
#define EPSILON_BIG 0.000015259 /* near 1e-5 but as a division by 2^n */
#define PI      3.14159265358979323846
static inline long long snapf_eps(float a);
static inline bool is_zero( float x ) { return 0 == snapf_eps(x); }
static inline bool is_equal( float x, float y ) { return x == y || 0 == snapf_eps(x) - snapf_eps(y); }
static inline bool is_greater( float x, float y ) { return snapf_eps(x) > snapf_eps(y); }
static inline bool is_lesser( float x, float y ) { return snapf_eps(x) < snapf_eps(y); }
static inline bool is_greater_equal( float x, float y ) { return snapf_eps(x) >= snapf_eps(y); }
static inline bool is_lesser_equal( float x, float y ) { return snapf_eps(x) <= snapf_eps(y); }

static inline long long snapf_eps_coarse(float a);
static inline bool is_zero_coarse( float x ) { return 0 == snapf_eps_coarse(x); }
static inline bool is_equal_coarse( float x, float y ) { return x == y || 0 == snapf_eps_coarse(x) - snapf_eps_coarse(y); }
static inline bool is_greater_coarse( float x, float y ) { return snapf_eps_coarse(x) > snapf_eps_coarse(y); }
static inline bool is_lesser_coarse( float x, float y ) { return snapf_eps_coarse(x) < snapf_eps_coarse(y); }
static inline bool is_greater_equal_coarse( float x, float y ) { return snapf_eps_coarse(x) >= snapf_eps_coarse(y); }
static inline bool is_lesser_equal_coarse( float x, float y ) { return snapf_eps_coarse(x) <= snapf_eps_coarse(y); }

static inline long long snapf_eps_eps(float a, float eps);
static inline bool is_zero_eps( float x, float eps ) { return 0 == snapf_eps_eps(x, eps); }
static inline bool is_equal_eps( float x, float y, float eps ) { return 0 == snapf_eps_eps(x, eps) - snapf_eps_eps(y, eps); }
static inline bool is_greater_eps( float x, float y, float eps ) { return snapf_eps_eps(x, eps) > snapf_eps_eps(y, eps); }
static inline bool is_lesser_eps( float x, float y, float eps ) { return snapf_eps_eps(x, eps) < snapf_eps_eps(y, eps); }
static inline bool is_greater_equal_eps( float x, float y, float eps ) { return snapf_eps_eps(x, eps) >= snapf_eps_eps(y, eps); }
static inline bool is_lesser_equal_eps( float x, float y, float eps ) { return snapf_eps_eps(x, eps) <= snapf_eps_eps(y, eps); }

// signed angle difference
static inline float angle_difference(float a, float b) {
    return fmod(a - b + 3.0f*PI, 2.0f*PI) - PI;
}

// arithmetic mean of two unsigned numbers without risking overflow
#define UNSIGNED_AVERAGE(a,b) (((a) / 2) + ((b) / 2) + ((a) & (b) & 1))

// algorithms - interface
static float2_t m_float2( float3_t const );                               //!< float3 to float2 conversion (z omitted)
static float3_t m_float3( float4_t const );                               //!< float4 to float3 conversion
static float4_t m_float4( float3_t const, float const );                  //!< float3 to float4 conversion

static float3_t float3_float4(float4_t const p) { return (float3_t) { p.x / p.w, p.y / p.w, p.z / p.w }; }
static float33_t float33_float44(float44_t const m) { return (float33_t) { m.m00, m.m01, m.m02, m.m10, m.m11, m.m12, m.m20, m.m21, m.m22 }; }

static bool      equal_float2( float2_t const, float2_t const );
static bool      lesser_float2( float2_t const, float2_t const);
static bool      greater_float2( float2_t const, float2_t const);
static bool      lesser_equal_float2( float2_t const, float2_t const);
static bool      greater_equal_float2( float2_t const, float2_t const);
static float2_t  add_float2( float2_t const, float2_t const );
static float2_t  sub_float2( float2_t const, float2_t const );
static float2_t  scale_float2( float2_t const, float const );
static float     dot_float2 (float2_t const, float2_t const );
static float     length_float2(float2_t const);
static float2_t  lerp_float2(float2_t const a, float2_t const b, float const t);           //!< Linear interpolation from a (t = 0) to b (t = 1)
static unsigned long long morton_float2(float2_t a);

static float3_t  add_float3   ( float3_t const, float3_t const );
static float3_t  sub_float3   ( float3_t const, float3_t const );
static float3_t  mul_float3  ( float3_t const, float3_t const );
static float3_t  cdiv_float3  ( float3_t const, float3_t const );
static float3_t  scale_float3 ( float3_t const, float const );
static float3_t  normal_float3( float3_t const );
static float     dot_float3   ( float3_t const, float3_t const );               //!< dot product
static float3_t  cross_float3 ( float3_t const, float3_t const );               //!< cross product
static float     length_float3       ( float3_t const );
static inline float quadrance_float3   ( float3_t const ); //!< length squared, easier to compute than length
static float     manhattan_float3       ( float3_t const );
static float3_t  swizzle_float3      ( float3_t const, int const, int const, int const );
static float3_t  lerp_float3(float3_t const a, float3_t const b, float const t);           //!< Linear interpolation from a (t = 0) to b (t = 1)
static float3_t  midpoint_float3(float3_t const a, float3_t const b);
static float3_t  assign_values_float3(float x, float y, float z);                       //!< assign a variable value
static float3_t  mix_float3(float3_t a, float3_t b, float c);   //!< linear combination of both
static float3_t  mix3_float3(float3_t a, float3_t b, float3_t c, float3_t u);   //!< linear combination of three points using barycentric coordinates
static float3_t  abs_float3   ( float3_t const );

/** \brief Decompose a vector into parallel and perpendicular components.
 *
 * \note Pass NULL to output parameters if not needed.
 *
 * \param [in] reference    Reference vector. Must be normalized.
 * \param [in] v            Vector to decompose.
 * \param [out] parallel    Part of v parallel to the reference vector.
 * \param [out] parallel    Part of v perpendicular to the reference vector.
 */
static void  decompose_float3 (float3_t const reference, float3_t const v, float3_t *parallel, float3_t *perpendicular);

static double4_t double4_float3( float3_t  const a, double const w ) { return (double4_t) { a.x, a.y, a.z, w }; };
static double4_t double4_float4( float4_t  const a ) { return (double4_t) { a.x, a.y, a.z, a.w }; };
static float3_t  float3_double4( double4_t const a ) { return (float3_t) { a.x, a.y, a.z }; };
static float4_t  float4_double4( double4_t const a ) { return (float4_t) { a.x, a.y, a.z, a.w }; };

static double4_t  add_double4   ( double4_t const, double4_t const );
static double4_t  sub_double4   ( double4_t const, double4_t const );
static double4_t  mul_double4   ( double4_t const, double    const );
static double4_t  cmul_double4  ( double4_t const, double4_t const );
static double4_t  cdiv_double4  ( double4_t const, double4_t const );
static double     dot_double4   ( double4_t const, double4_t const );            //!< dot product of all components
static double     dot3_double4  ( double4_t const, double4_t const );            //!< dot product of x, y, z
static double4_t  cross3_double4( double4_t const, double4_t const );            //!< cross product of x,y,z components
static double4_t  normal_double4( double4_t const );
static double     length_double4( double4_t const );

/** \brief Oblique angle between two vectors.
 *
 * \param [in] a            Vector a.
 * \param [in] b            Vector b.
 *
 * \return Angle between the two vectors (from 0 to pi)
 */
static float oblique_angle_between_float3(float3_t a, float3_t b);

/** \brief Compute the circumcenter and radius of a triangle.
 *
 * \param [in] a            Vertex a.
 * \param [in] b            Vertex b.
 * \param [in] c            Vertex c.
 * \param [out] out_center  Circumcenter of the triangle.
 * \param [out] out_radius  Radius of the circumcircle.
 *
 * \return True if the vertices are not collinear.
 */
static bool triangle_circumcircle(float3_t a, float3_t b, float3_t c, float3_t *out_center, float *out_radius);

static float triangle_area(float3_t a, float3_t b, float3_t c);

/** \brief Compute the area of a voronoi cell of a triangle if possible.
 *
 * \note This is the area for the part of the triangle that is closest to point o,
 *       so the part between the vertex, the midpoints of the edges connecting to it,
 *       and the center of the triangle's circumcircle.
 * 
 *       For obtuse triangles (when the circumcircle's center is outside the triangle),
 *       instead of using the center it moves it onto the edge opposite to point o.
 */
static float mixed_voronoi_cell_area(float3_t o, float3_t a, float3_t b);

static float4_t plane_from_vector_and_point(float4_t vector, float3_t point){
    float4_t plane;
    plane.x = vector.x;
    plane.y = vector.y;
    plane.z = vector.z;

    plane.w = -1 *(vector.x * point.x + vector.y * point.y + vector.z * point.z);

    return plane;
}

static float3_t plane_point_float3(
    float2_t st_point,
    float3_t plane_origin,
    float3_t plane_s_axis,
    float3_t plane_t_axis
);

static bool      equal_float3            ( float3_t const, float3_t const );
static bool      equal_float3_ptr        ( float3_t const * const a, float3_t const * const b );
static bool      lesser_float3           ( float3_t const, float3_t const);
static bool      greater_float3          ( float3_t const, float3_t const);
static bool      lesser_equal_float3     ( float3_t const, float3_t const);
static bool      greater_equal_float3    ( float3_t const, float3_t const);
static bool      any_equal_float3        ( float3_t const, float3_t const );
static bool      any_lesser_float3       ( float3_t const, float3_t const);
static bool      any_greater_float3      ( float3_t const, float3_t const);
static bool      any_lesser_equal_float3 ( float3_t const, float3_t const);
static bool      any_greater_equal_float3( float3_t const, float3_t const);

static bool      equal_float3_eps            ( float3_t const, float3_t const, float eps);
static bool      equal_float3_ptr_eps        ( float3_t const * const a, float3_t const * const b, float eps);
static bool      lesser_float3_eps           ( float3_t const, float3_t const, float eps);
static bool      greater_float3_eps          ( float3_t const, float3_t const, float eps);
static bool      lesser_equal_float3_eps     ( float3_t const, float3_t const, float eps);
static bool      greater_equal_float3_eps    ( float3_t const, float3_t const, float eps);
static bool      any_equal_float3_eps        ( float3_t const, float3_t const, float eps);
static bool      any_lesser_float3_eps       ( float3_t const, float3_t const, float eps);
static bool      any_greater_float3_eps      ( float3_t const, float3_t const, float eps);
static bool      any_lesser_equal_float3_eps ( float3_t const, float3_t const, float eps);
static bool      any_greater_equal_float3_eps( float3_t const, float3_t const, float eps);

/** Project a vector A onto vector B.
 * 
 * \param [in] a            Vector A.
 * \param [in] b            Vector B.
 * \param [out] r           r = b = (a . b) / (b . b)
 *
 * \return True if the projection was possible (b is not zero)
 */
static inline bool project_float3(float3_t a, float3_t b, float3_t *r);

/** \brief Gram-Schmidt orthogonalization of a matrix, from first to last column (second index).
 */
static bool gram_schmidt_column_float33(float33_t m, float33_t *r);

/** \brief QR decomposition of a matrix.
 * 
 * \note A = Q*R, where Q is orthogonal and R is upper triangular.
 * 
 * \param [in] a            Matrix to decompose.
 * \param [out] q_out       Q matrix.
 * \param [out] r_out       R matrix.
 *
 * \return True if successful.
 */
static bool qr_float33(float33_t a, float33_t *q_out, float33_t *r_out);

/** \brief Compute eigenvalues of a matrix using QR decomposition.
 * 
 * \param [in] a            Matrix to compute eigenvalues of.
 * \param [in] epsilon      Epsilon for convergence.
 * \param [in] iter_limit   Maximum number of iterations.
 * \param [out] out         Array of eigenvalues, from largest to smallest.
 */
static bool eigenvalue_qr_float33(float33_t a, float epsilon, size_t iter_limit, float out_value[3]);

/** \brief Compute an eigenvector of a matrix corresponding to a given eigenvalue.
 * 
 * \note The input matrix must be normal, that is, m*m^T = m^T*m.
 *       This is true for symmetric matrices, and possibly others.
 *       This condition is not checked.
 *
 *       !!! normal matrix does not mean normalized matrix !!!
 *       The function has "symmetric" in its name to avoid confusion.
 * 
 * \note Result: v such that m*v = eigenvalue*v
 * 
 * \param [in] m                    Matrix to compute eigenvector of.
 * \param [in] eigenvalue           Eigenvalue to compute eigenvector for.
 * \param [out] out                 Output eigenvector(s). Allocate to at least 3 elements.
 * \param [out] out_multiplicity    How many eigenvectors were found for this eigenvalue.
 *
 * \return True if successful.
 */
static bool eigenvector_symmetric_float33(float33_t m, float eigenvalue, float3_t *out, int *out_multiplicity);

/** \brief Compute 3 eigenvectors of a matrix corresponding to the given eigenvalues.
 * 
 * \note The input matrix must be normal, that is, m*m^T = m^T*m.
 *       This is true for symmetric matrices, and possibly others.
 *       This condition is not checked.
 * 
 *       !!! normal matrix does not mean normalized matrix !!!
 *       The function has "symmetric" in its name to avoid confusion.
 */
static bool eigenvectors_symmetric_float33(float33_t m, float value[3], float3_t vector[3]);

/** \brief Compute eigenvalues and eigenvectors of a matrix.
 *
 * \note The input matrix must be normal, that is, m*m^T = m^T*m.
 *       This is true for symmetric matrices, and possibly others.
 *       This condition is not checked.
 * 
 *       !!! normal matrix does not mean normalized matrix !!!
 *       The function has "symmetric" in its name to avoid confusion.
 * 
 * \param [in] m            Matrix to compute eigenvalues and eigenvectors of.
 * \param [out] value       Array of eigenvalues, from smallest to largest (set to NULL if not needed)
 * \param [out] vector      Array of eigenvectors, for reach eigenvalue (set to NULL if not needed, but consider `eigenvalue_qr_float33` instead)
 *
 * \return True if successful.
 */
static bool eigen_symmetric_float33(float33_t m, float value[3], float3_t vector[3]);
 
static float4_t  transform_float44      ( float44_t const, float4_t  const );   //!< transform vector by the matrix (right side)
static float4_t  cotransform_float44    ( float44_t const, float4_t  const );   //!< transform vector by the matrix (left side)
static float3_t  transform_float33      ( float33_t const, float3_t  const );   //!< transform 3-vector by a 3x3 matrix (right side)
static float3_t  cotransform_float33    ( float33_t const, float3_t  const );   //!< transform 3-vector by a 3x3 matrix (left side)
static float33_t mul_float33            ( float33_t const, float33_t const );   //!< multiply matrix by matrix
static float44_t mul_float44            ( float44_t const, float44_t const );   //!< multiply matrix by matrix
static float44_t mul_scalar_float44     ( float44_t const, float     const );   //!< multiply matrix by scalar
static float44_t transpose_float44      ( float44_t const );                    //!< transpose matrix
static float33_t transpose_float33      ( float33_t const );                    //!< transpose a 3x3 matrix
static float     det_float44            ( float44_t const );                    //!< get determinant of the matrix
static float44_t invert_float44         ( float44_t const, float     const );   //!< invert matrix with its non-zero determinant
static float     det_float33            ( float33_t const );                    //!< get determinant of a 3x3 matrix
static float33_t invert_float33         ( float33_t const, float     const );   //!< invert a 3x3 matrix with its non-zero determinant
static float33_t rotation_float33       ( float angle_a );
static float33_t rotation_sincos_float33( float sin_angle, float cos_angle );
static float44_t rotation_float44       ( float3_t axis_a, float angle_a );
static float44_t rotation_sincos_float44( float3_t axis_a, float sin_angle_a, float cos_angle_a );
static float33_t translation_float33    ( float2_t vector );
static float33_t cotranslation_float33  ( float2_t vector );
static float44_t translation_float44    ( float3_t vector );
static float44_t cotranslation_float44  ( float3_t vector );
static float33_t scale_uniform_float33  ( float s );
static float33_t scale_float33          ( float2_t vector );
static float44_t scale_float44          ( float3_t vector );
static bool      equal_float44          ( float44_t const, float44_t const );
static float44_t decompose_float44      ( float44_t A );
static float44_t gram_schmidt_float44   ( float44_t A );
static float44_t add_float44            ( float44_t const a, float44_t const b);

static float frobenius_float33(float33_t m); /* usefull for retrieving scale from tranformation (divide by sqrt3()) */

/** \brief Compute a basis change matrix.
 *
 * \note
 *  * Vectors need not be orthonormal, if you want to have a weird basis.
 *  * Inverse matrix is computed as an intermediate step whether you want it or not.
 *  * If the function returns false, the output arguments are not touched.
 *
 * \param [in] a            X vector of the new basis.
 * \param [in] b            Y vector of the new basis.
 * \param [in] c            Z vector of the new basis.
 * \param [out] out         Output matrix (set to NULL if not needed)
 * \param [out] out_inv     Inverse matrix (set to NULL if not needed)
 *
 * \return True if the input vectors form a valid 3d basis.
 */
static bool basis_change_float44(float3_t a, float3_t b, float3_t c, float44_t *out, float44_t *out_inv);

/** \brief Compute a basis change matrix from three points.
 *
 * \note
 *  The computed matrix M is such that:
 *   det(M) = 1
 *   a * M = (a_x, a_y, 0)
 *   b * M = (b_x, 0,   0)
 *   c * M = (0,   0,   0)
 *  If the function returns false, m_out is not touched.
 *
 * \param [in] a            Point A.
 * \param [in] b            Point B.
 * \param [in] c            Point C.
 * \param [out] m_out       Basis change matrix.
 *
 * \return True on success.
 */
static bool basis_from_triangle(float3_t const a, float3_t const b, float3_t const c, float44_t *m_out);

/** \brief Compute a rotated and centered basis change matrix from three points.
 *
 * \note
 *  The computed matrix M is such that:
 *   det(M) = 1
 *   a * M = (a_x, a_y, 0)
 *   b * M = (b_x, b_y, 0)
 *   c * M = (c_x, c_y, 0)
 *   circumcenter(a,b,c) * M = (0,0,0)
 *  If the function returns false, m_out and m_inv_out are not touched.
 *
 * \param [in] a            Point A.
 * \param [in] b            Point B.
 * \param [in] c            Point C.
 * \param [in] rotation     Rotation around the circumcenter of the triangle.
 * \param [out] m_out       Basis change matrix.
 * \param [out] m_inv_out   Inverse basis change matrix.
 *
 * \return True on success.
 */
static bool centered_basis_from_triangle(float3_t const a, float3_t const b, float3_t const c, float rotation, float44_t *m_out, float44_t *m_inv_out);

static double44_t mul_double44       ( double44_t const, double44_t const );   //!< multiply matrix by matrix
static double44_t mul_scalar_double44( double44_t const, double     const );   //!< multiply matrix by scalar
static double44_t transpose_double44 ( double44_t const );                    //!< transpose matrix
static double     det_double44       ( double44_t const );                    //!< get determinant of the matrix
static double44_t invert_double44    ( double44_t const, double     const );   //!< invert matrix with its non-zero determinant
static bool       equal_double44     ( double44_t const, double44_t const );
static double44_t add_double44       ( double44_t const, double44_t const b);

/** \brief Compute a matrix that rotates unit vector a onto unit vector b.
 *
 * \note This kind of rotation is not unique. You can pre-multiply the result with any rotation
 *       around vector a, or post-multiply with any rotation around vector b,
 *       and you will still get a rotation that brings a onto b.
 *       This function just returns one of infinitely many such rotations.
 *
 * \param [in] a        Source unit (this is not checked) vector.
 * \param [in] b        Destination unit (this is not checked) vector.
 *
 * \return Matrix M such that Ma = b
 */
static float44_t rotation_of_vector_onto(float3_t const a, float3_t const b);

static float    line_point_distance( float3_t d /* line direction */, float3_t o /* line origin */, float3_t p /* point */);
static float    skew_lines_distance( float3_t da /* 1st line direction */, float3_t oa /* 1st line origin */, float3_t db/* 2nd line direction */, float3_t ob /* 2nd line origin */ );
static float    lines_distance     ( float3_t da /* 1st line direction */, float3_t oa /* 1st line origin */, float3_t db/* 2nd line direction */, float3_t ob /* 2nd line origin */ );
static float3_t lines_intersection ( float3_t da /* 1st line direction */, float3_t oa /* 1st line origin */, float3_t db/* 2nd line direction */, float3_t ob /* 2nd line origin */ );

/* Compute nearest points on skew lines.
    \param [in] da      Direction of line a.
    \param [in] oa      Origin of line a.
    \param [in] db      Direction of line b.
    \param [in] ob      Origin of line b.
    \param [out] ca     Point on line a nearest to line b.
    \param [out] cb     Point on line b nearest to line a.

    \return 0 on success, -1 when lines are parallel. On Success values are stored in ca and cb.
*/
static int skew_lines_nearest(float3_t da, float3_t oa,float3_t db, float3_t ob, float3_t *ca, float3_t *cb);

static void compose_axis_angle(
    float3_t axis_0_a, float first_angle_0_a,
    float3_t axis_1_a, float angle_1_a,
    float3_t *axis_out_a, float *angle_out_a
);

static float44_t orthographic   (float r, float l, float t, float b, float f, float n);

static float maxf  ( float, float );
static float minf  ( float, float );
static float clampf( float, float, float );
static float scalar_product_fun_for_matrices( float44_t matrix, int column, bool pointer, int previous );

static int   maxi( int, int );
static int   mini( int, int );
static int   clampi( int, int, int );

static bool is_power_of_2(unsigned int);
static unsigned int next_power_of_2(unsigned int);

static float2_t max_float2( float2_t const, float2_t const);
static float2_t min_float2( float2_t const, float2_t const);

static float3_t max_float3( float3_t const, float3_t const);
static float3_t min_float3( float3_t const, float3_t const);
static float3_t clamp_float3( float3_t const, float3_t const, float3_t const);
static float    max_of_float3( float3_t const a );
static float    min_of_float3( float3_t const a );

static float _signumf(float);

// debug
static void printf_matrix(float44_t matrix);
static void printf_matrix33(float33_t matrix);

// algorithms - implementation

static inline long long snapf_eps(float a) {
    return (long long)floor((double)a / (double)EPSILON);
}

static inline long long snapf_eps_coarse(float a) {
    return (long long)floor((double)a / (double)EPSILON_BIG);
}

static inline long long snapf_eps_eps(float a, float eps) {
    return (long long)floor((double)a / (double)eps);
}

static inline float2_t m_float2( float3_t v )
{
#ifdef __cplusplus
	float2_t ret;
	ret.x = v.x; ret.y = v.y;
	return ret;
#else
    return (float2_t){v.x, v.y};
#endif
}

static inline float3_t m_float3( float4_t v )
{
    return v.v;
}

static inline float4_t m_float4( float3_t a, float b )
{
    float4_t r;

    r.v = a;
    r.w = b;

    return r;
}

static bool inline equal_float2( float2_t const a, float2_t const b )
{
    return (0LL == snapf_eps(a.x) - snapf_eps(b.x)) &&
        (0LL == snapf_eps(a.y) - snapf_eps(b.y));
}

static bool lesser_float2(float2_t const a, float2_t const b) {
    return is_lesser(a.x, b.x) && is_lesser(a.y, b.y);
}

static bool greater_float2(float2_t const a, float2_t const b) {
    return is_greater(a.x, b.x) && is_greater(a.y, b.y);
}

static bool lesser_equal_float2(float2_t const a, float2_t const b) {
    return is_lesser_equal(a.x, b.x) && is_lesser_equal(a.y, b.y);
}

static bool greater_equal_float2(float2_t const a, float2_t const b) {
    return is_greater_equal(a.x, b.x) && is_greater_equal(a.y, b.y);
}

static inline float2_t  add_float2( float2_t const a, float2_t const b )
{
#ifdef __cplusplus
	float2_t ret;
	ret.x = a.x + b.x; ret.y = a.y + b.y;
	return ret;
#else
	return (float2_t) { a.x + b.x, a.y + b.y };
#endif
}

static inline float2_t  sub_float2( float2_t const a, float2_t const b )
{
#ifdef __cplusplus
	float2_t ret;
	ret.x = a.x - b.x; ret.y = a.y - b.y;
	return ret;
#else
	return (float2_t) { a.x - b.x, a.y - b.y };
#endif
}

static inline float2_t normal_float2(float2_t const v) {
	float l = sqrtf(v.x*v.x + v.y*v.y);
#ifdef __cplusplus
	float2_t ret;
	ret.x = v.x/l;
	ret.y = v.y/l;
	return ret;
#else
	return (float2_t) { .x = v.x/l, .y = v.y/l };
#endif
}

static inline float dot_float2(float2_t const a, float2_t const b) {
	return a.x*b.x + a.y*b.y;
}

static inline float2_t scale_float2(float2_t const a, float const s) {
#ifdef __cplusplus
	float2_t ret;
	ret.x = a.x*s;
	ret.y = a.y*s;
	return ret;
#else
	return (float2_t) { .x = a.x*s, .y = a.y*s };
#endif
}

static inline float length_float2(float2_t const v) {
	return sqrtf(v.x*v.x + v.y*v.y);
}

static inline float2_t  lerp_float2(float2_t const a, float2_t const b, float const t) {
	float2_t const v = sub_float2(b,a);
	return add_float2(a, scale_float2(v, t));
}

static inline float3_t  add_float3( float3_t const a, float3_t const b )
{
    float3_t r = { a.x + b.x, a.y + b.y, a.z + b.z };
    return r;
}
// a - b
static inline float3_t  sub_float3( float3_t const a, float3_t const b )
{
    float3_t r = { a.x - b.x, a.y - b.y, a.z - b.z };
    return r;
}

static float3_t  mul_float3( float3_t const a, float3_t const b )
{
    float3_t r = { a.x * b.x, a.y * b.y, a.z * b.z };
    return r;
}
static float3_t  cdiv_float3( float3_t const a, float3_t const b )
{
    float3_t r = { a.x / b.x, a.y / b.y, a.z / b.z };
    return r;
}

static inline float3_t scale_float3( float3_t const a, float const s )
{
    float3_t r = { a.x * s, a.y * s, a.z * s };
    return r;
}

static inline float3_t normal_float3( float3_t const a )
{
	float const dl = dot_float3(a, a);

	// left it EXACTLY as here, no EPSILON, dl may be realllly small
	if (0 == dl) {
#ifdef __cplusplus
		float3_t ret;
		ret.x = 0.f; ret.y = 0.f; ret.z = 0.f;
		return ret;
#else
		return (float3_t) { 0.f, 0.f, 0.f };
#endif
	}

    return scale_float3(a, 1.f/sqrtf(dl));
}

static inline float dot_float3( float3_t const a, float3_t const b )
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

static inline float3_t cross_float3( float3_t const a, float3_t const b )
{
    float3_t r;

    r.x = a.y * b.z - a.z * b.y;
    r.y = a.z * b.x - a.x * b.z;
    r.z = a.x * b.y - a.y * b.x;

    return r;
}

static float3_t assign_values_float3(float x, float y, float z){
    float3_t r;

    r.x = x;
    r.y = y;
    r.z = z;

    return r;
}

static float3_t  mix_float3(float3_t a, float3_t b, float c) {
    assert(c <= 1.0+EPSILON);
    assert(c >= 0.0-EPSILON);
    return add_float3(scale_float3(a, 1.0 - c), scale_float3(b, c));
}
static float3_t  abs_float3   ( float3_t const a ) {
    return (float3_t){ fabsf(a.x), fabsf(a.y), fabsf(a.z) };
}

static void decompose_float3 (float3_t const reference, float3_t const v, float3_t *parallel, float3_t *perpendicular) {
    float3_t const par = scale_float3(reference, dot_float3(reference, v));

    if (parallel) {
        *parallel = par;
    }

    if (perpendicular) {
        *perpendicular = sub_float3(v, par);
    }
}

static float3_t  mix3_float3(float3_t a, float3_t b, float3_t c, float3_t u) {
    return (float3_t) {
        a.x * u.x + b.x * u.y + c.x * u.z,
        a.y * u.x + b.y * u.y + c.y * u.z,
        a.z * u.x + b.z * u.y + c.z * u.z,
    };
}

static long long hilbert_float2(float2_t a) {
	union {
		unsigned i;
		float f;
	} u, v;

	u.f = a.x;
	v.f = a.y;

	int su = (int)(u.i) >> 31;
	int sv = (int)(v.i) >> 31;

	long long x = (((u.i & 0x7FFFFFFFULL) ^ su) - su) + 0x7FFFFFFFULL;
	long long y = (((v.i & 0x7FFFFFFFULL) ^ sv) - sv) + 0x7FFFFFFFULL;
	long long n = 1LL << (32LL);

	long long h = 0LL;
	long long p, q;
	for (long long b = n >> 1LL; b > 0; b >>= 1LL) {
		p = (x & b) > 0LL;
		q = (y & b) > 0LL;
		h += b * b * ((3LL * p) ^ q);
		if (q == 0LL) {
			if (p == 1LL) {
				x = n - 1LL - x;
				y = n - 1LL - y;
			}

			long long t = x;
			x = y;
			y = t;
		}
	}
	return h;
}

static unsigned long long morton_float2(float2_t a) {
    /*  Taking the fact that IEEE float may be sorted as 1-complement integer.
        For more info see: https://stackoverflow.com/questions/26856268/morton-index-from-2d-point-with-floats
    */
    union {
        unsigned int i;
        float f;
    } x, y;

    x.f = a.x;
    y.f = a.y;

    /* Sign mask. */
    unsigned long long sx = (unsigned long long)((int)(x.i) >> 31);
    unsigned long long sy = (unsigned long long)((int)(y.i) >> 31);

    /* Biasing. */
    unsigned long long xx = 0x80000000ULL - ((unsigned long long)x.i & sx & 0x7FFFFFFFULL) + ((unsigned long long)x.i & ~sx);
    unsigned long long yy = 0x80000000ULL - ((unsigned long long)y.i & sy & 0x7FFFFFFFULL) + ((unsigned long long)y.i & ~sy);

    /* Converting to morton code with bit-twiddling. Possibly some intrisic for
       permutation/shuffle may be thought for optimization. */
    xx = (xx | (xx << 16)) & 0x0000ffff0000ffffULL;
    yy = (yy | (yy << 16)) & 0x0000ffff0000ffffULL;

    xx = (xx | (xx <<  8)) & 0x00ff00ff00ff00ffULL;
    yy = (yy | (yy <<  8)) & 0x00ff00ff00ff00ffULL;

    xx = (xx | (xx <<  4)) & 0x0f0f0f0f0f0f0f0fULL;
    yy = (yy | (yy <<  4)) & 0x0f0f0f0f0f0f0f0fULL;

    xx = (xx | (xx <<  2)) & 0x3333333333333333ULL;
    yy = (yy | (yy <<  2)) & 0x3333333333333333ULL;

    xx = (xx | (xx <<  1)) & 0x5555555555555555ULL;
    yy = (yy | (yy <<  1)) & 0x5555555555555555ULL;

    return xx | (yy << 1);
}

static float quadrance_float3 ( float3_t const a ) {
    return a.x*a.x + a.y*a.y + a.z*a.z;
}

static inline float length_float3( float3_t const a )
{
    return sqrtf(quadrance_float3(a) );
}

static float manhattan_float3( float3_t const a ) {
    return fabsf(a.x) + fabsf(a.y) + fabsf(a.z);
}

static inline double4_t add_double4( double4_t const a, double4_t const b ) {
    return (double4_t) { a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w };
}

static inline double4_t sub_double4( double4_t const a, double4_t const b ) {
    return (double4_t) { a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w };
}

static inline double4_t mul_double4( double4_t const a, double const s ) {
    return (double4_t) { a.x * s, a.y * s, a.z * s, a.w * s };
}

static inline double4_t cmul_double4( double4_t const a, double4_t const b ) {
    return (double4_t) { a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w };
}

static inline double4_t cdiv_double4( double4_t const a, double4_t const b ){
    return (double4_t) { a.x / b.x, a.y / b.y, a.z / b.z, a.w / b.w };
}

static inline double dot_double4( double4_t const a, double4_t const b ) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

static inline double dot3_double4  ( double4_t const a, double4_t const b ) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

static inline double4_t cross3_double4( double4_t const a, double4_t const b ) {
    assert(0 == a.w);
    assert(0 == b.w);

    return (double4_t) {
        .x = a.y * b.z - a.z * b.y,
        .y = a.z * b.x - a.x * b.z,
        .z = a.x * b.y - a.y * b.x,
        .w = 0.0
    };
}

static inline double4_t normal_double4( double4_t const a ) {
    return mul_double4(a, 1.0 / length_double4(a));
}

static double length_double4( double4_t const a ) {
    return sqrt(dot_double4(a, a));
}

static inline float cos_between_floats3( float3_t const a, float3_t const b )
{
	return dot_float3(a,b) / ( length_float3(a) * length_float3(b) );
}

static float oblique_angle_between_float3(float3_t a, float3_t b) {
    return atan2f(length_float3(cross_float3(a, b)), dot_float3(a, b));
}

static bool triangle_circumcircle(float3_t a, float3_t b, float3_t c, float3_t *out_center, float *out_radius) {
    float3_t const ab = sub_float3(a, b);
    float3_t const bc = sub_float3(b, c);
    float3_t const ca = sub_float3(c, a);

    float const q_ab = quadrance_float3(ab);
    float const q_bc = quadrance_float3(bc);
    float const q_ca = quadrance_float3(ca);

    float const l_ab = sqrtf(q_ab);
    float const l_bc = sqrtf(q_bc);
    float const l_ca = sqrtf(q_ca);

    float const q_cross_ab_bc = quadrance_float3(cross_float3(ab, bc));

    if (is_zero_coarse(q_cross_ab_bc)) {
        return false;
    }

    if (out_radius) {
        *out_radius = 0.5f*l_ab*l_bc*l_ca / sqrtf(q_cross_ab_bc);
    }

    if (out_center) {
        float const scale = -0.5f/q_cross_ab_bc;

        float const u = scale*q_bc*dot_float3(ab, ca);
        float const v = scale*q_ca*dot_float3(ab, bc);
        float const w = scale*q_ab*dot_float3(ca, bc);

        float3_t const au = scale_float3(a, u);
        float3_t const bv = scale_float3(b, v);
        float3_t const cw = scale_float3(c, w);

        *out_center = add_float3(au, add_float3(bv, cw));
    }

    return true;
}

static float triangle_area(float3_t a, float3_t b, float3_t c) {
    return 0.5f * length_float3(cross_float3(sub_float3(b, a), sub_float3(c, a)));
}

static float mixed_voronoi_cell_area(float3_t o, float3_t a, float3_t b) {
    float const a_o = oblique_angle_between_float3(sub_float3(a, o), sub_float3(b, o));
    float const a_b = oblique_angle_between_float3(sub_float3(a, b), sub_float3(o, b));
    float const a_a = oblique_angle_between_float3(sub_float3(b, a), sub_float3(o, a));

    float3_t center;

    if ((a_o < M_PI_2) && (a_b < M_PI_2) && (a_a < M_PI_2)) {
        if (!triangle_circumcircle(o, a, b, &center, NULL)) {
            return NAN;
        }
    } else {
        if ((a_o > a_b) && (a_o > a_a)) {
            center = midpoint_float3(a, b);
        } else if ((a_a > a_b) && (a_a > a_o)) {
            center = midpoint_float3(o, b);
        } else {
            center = midpoint_float3(o, a);
        }
    }

    float3_t const m_oa = midpoint_float3(o, a);
    float3_t const m_ob = midpoint_float3(o, b);

    float const area_oa = triangle_area(o, m_oa, center);
    float const area_ob = triangle_area(o, m_ob, center);

    return area_oa + area_ob;
}

#if 0
static float mixed_voronoi_cell_area(float3_t o, float3_t a, float3_t b) {
    float3_t const ab = sub_float3(a, b);
    float3_t const ao = sub_float3(a, o);
    float3_t const bo = sub_float3(b, o);

    float const q_ab = quadrance_float3(ab);
    float const q_ao = quadrance_float3(ao);
    float const q_bo = quadrance_float3(bo);

    if ((q_ab + q_ao > q_bo) && (q_ab + q_bo > q_ao) && (q_ao + q_bo > q_ab)) {
        // triangle is acute
        float const angle_a = oblique_angle_between_float3(sub_float3(o, a), sub_float3(b, a));
        float const angle_b = oblique_angle_between_float3(sub_float3(o, b), sub_float3(a, b));

        float const sin_a = sinf(angle_a);
        float const cos_a = cosf(angle_a);
        float const sin_b = sinf(angle_b);
        float const cos_b = cosf(angle_b);

        return (q_ao*cos_a/sin_a + q_bo*cos_b/sin_b)/8.0f;
    } else {
        // triangle is obtuse
        float const angle = oblique_angle_between_float3(ao, bo);

        float const area = length_float3(cross_float3(ab, ao)) * 0.5f;

        // check if the angle at this vertex is obtuse
        if (angle > (M_PI_2 + EPSILON_BIG)) {
            return 0.5f*area;
        } else {
            return 0.25f*area;
        }
    }
}
#endif

static bool inline equal_float3( float3_t const a, float3_t const b )
{
    return (0LL == snapf_eps(a.x) - snapf_eps(b.x)) &&
        (0LL == snapf_eps(a.y) - snapf_eps(b.y)) &&
        (0LL == snapf_eps(a.z) - snapf_eps(b.z));
}

static bool inline equal_float3_ptr( float3_t const * const a, float3_t const * const b )
{
    return equal_float3(*a, *b);
}

static inline float3_t  lerp_float3(float3_t const a, float3_t const b, float const t) {
	float3_t const v = sub_float3(b,a);
	return add_float3(a, scale_float3(v, t));
}

static float3_t  midpoint_float3(float3_t const a, float3_t const b) {
    return (float3_t) {
        (a.x + b.x) * 0.5f,
        (a.y + b.y) * 0.5f,
        (a.z + b.z) * 0.5f
    };
}

static inline float4_t transform_float44( float44_t const a, float4_t const b )
{
    float4_t r;

    r.x = a.m00 * b.x + a.m01 * b.y + a.m02 * b.z + a.m03 * b.w;
    r.y = a.m10 * b.x + a.m11 * b.y + a.m12 * b.z + a.m13 * b.w;
    r.z = a.m20 * b.x + a.m21 * b.y + a.m22 * b.z + a.m23 * b.w;
    r.w = a.m30 * b.x + a.m31 * b.y + a.m32 * b.z + a.m33 * b.w;

    return r;
}

static float4_t  cotransform_float44( float44_t const a, float4_t  const b ) {
    return transform_float44(transpose_float44(a), b);
}

static inline float3_t transform_float33( float33_t const a, float3_t const b )
{
    float3_t r;

    r.x = a.m00 * b.x + a.m01 * b.y + a.m02 * b.z;
    r.y = a.m10 * b.x + a.m11 * b.y + a.m12 * b.z;
    r.z = a.m20 * b.x + a.m21 * b.y + a.m22 * b.z;

    return r;
}

static float3_t  cotransform_float33( float33_t const a, float3_t  const b ) {
    return transform_float33(transpose_float33(a), b);
}

static float33_t mul_float33( float33_t const a, float33_t const b) {
    float33_t r;

    r.m00 = a.m02*b.m20+a.m01*b.m10+a.m00*b.m00;
    r.m01 = a.m02*b.m21+a.m01*b.m11+a.m00*b.m01;
    r.m02 = a.m02*b.m22+a.m01*b.m12+a.m00*b.m02;

    r.m10 = a.m12*b.m20+a.m11*b.m10+a.m10*b.m00;
    r.m11 = a.m12*b.m21+a.m11*b.m11+a.m10*b.m01;
    r.m12 = a.m12*b.m22+a.m11*b.m12+a.m10*b.m02;

    r.m20 = a.m22*b.m20+a.m21*b.m10+a.m20*b.m00;
    r.m21 = a.m22*b.m21+a.m21*b.m11+a.m20*b.m01;
    r.m22 = a.m22*b.m22+a.m21*b.m12+a.m20*b.m02;

    return r;
}

static inline float44_t mul_float44( float44_t const a, float44_t const b )
{
    float44_t r;

    r.m00 = a.m03*b.m30+a.m02*b.m20+a.m01*b.m10+a.m00*b.m00;
    r.m01 = a.m03*b.m31+a.m02*b.m21+a.m01*b.m11+a.m00*b.m01;
    r.m02 = a.m03*b.m32+a.m02*b.m22+a.m01*b.m12+a.m00*b.m02;
    r.m03 = a.m03*b.m33+a.m02*b.m23+a.m01*b.m13+a.m00*b.m03;

    r.m10 = a.m13*b.m30+a.m12*b.m20+a.m11*b.m10+a.m10*b.m00;
    r.m11 = a.m13*b.m31+a.m12*b.m21+a.m11*b.m11+a.m10*b.m01;
    r.m12 = a.m13*b.m32+a.m12*b.m22+a.m11*b.m12+a.m10*b.m02;
    r.m13 = a.m13*b.m33+a.m12*b.m23+a.m11*b.m13+a.m10*b.m03;

    r.m20 = a.m23*b.m30+a.m22*b.m20+a.m21*b.m10+a.m20*b.m00;
    r.m21 = a.m23*b.m31+a.m22*b.m21+a.m21*b.m11+a.m20*b.m01;
    r.m22 = a.m23*b.m32+a.m22*b.m22+a.m21*b.m12+a.m20*b.m02;
    r.m23 = a.m23*b.m33+a.m22*b.m23+a.m21*b.m13+a.m20*b.m03;

    r.m30 = a.m33*b.m30+a.m32*b.m20+a.m31*b.m10+a.m30*b.m00;
    r.m31 = a.m33*b.m31+a.m32*b.m21+a.m31*b.m11+a.m30*b.m01;
    r.m32 = a.m33*b.m32+a.m32*b.m22+a.m31*b.m12+a.m30*b.m02;
    r.m33 = a.m33*b.m33+a.m32*b.m23+a.m31*b.m13+a.m30*b.m03;

    return r;
}

static inline float44_t mul_scalar_float44( float44_t const a, float const b )
{
    float44_t r;

    r.m00 = a.m00 * b;
    r.m01 = a.m01 * b;
    r.m02 = a.m02 * b;
    r.m03 = a.m03 * b;

    r.m10 = a.m10 * b;
    r.m11 = a.m11 * b;
    r.m12 = a.m12 * b;
    r.m13 = a.m13 * b;

    r.m20 = a.m20 * b;
    r.m21 = a.m21 * b;
    r.m22 = a.m22 * b;
    r.m23 = a.m23 * b;

    r.m30 = a.m30 * b;
    r.m31 = a.m31 * b;
    r.m32 = a.m32 * b;
    r.m33 = a.m33 * b;

    return r;
}

static inline float44_t transpose_float44( float44_t const a )
{
    float44_t r;

    r.m00 = a.m00;
    r.m01 = a.m10;
    r.m02 = a.m20;
    r.m03 = a.m30;

    r.m10 = a.m01;
    r.m11 = a.m11;
    r.m12 = a.m21;
    r.m13 = a.m31;

    r.m20 = a.m02;
    r.m21 = a.m12;
    r.m22 = a.m22;
    r.m23 = a.m32;

    r.m30 = a.m03;
    r.m31 = a.m13;
    r.m32 = a.m23;
    r.m33 = a.m33;

    return r;
}

static inline float33_t transpose_float33( float33_t const a )
{
    float33_t r;

    r.m00 = a.m00;
    r.m01 = a.m10;
    r.m02 = a.m20;

    r.m10 = a.m01;
    r.m11 = a.m11;
    r.m12 = a.m21;

    r.m20 = a.m02;
    r.m21 = a.m12;
    r.m22 = a.m22;

    return r;
}


static inline float det_float44( float44_t a )
{
    return  a.m00*(a.m11*(a.m22*a.m33-a.m23*a.m32)-a.m12*(a.m21*a.m33-a.m23*a.m31)+a.m13*(a.m21*a.m32-a.m22*a.m31))
           -a.m01*(a.m10*(a.m22*a.m33-a.m23*a.m32)-a.m12*(a.m20*a.m33-a.m23*a.m30)+a.m13*(a.m20*a.m32-a.m22*a.m30))
           +a.m02*(a.m10*(a.m21*a.m33-a.m23*a.m31)-a.m11*(a.m20*a.m33-a.m23*a.m30)+a.m13*(a.m20*a.m31-a.m21*a.m30))
           -a.m03*(a.m10*(a.m21*a.m32-a.m22*a.m31)-a.m11*(a.m20*a.m32-a.m22*a.m30)+a.m12*(a.m20*a.m31-a.m21*a.m30));
}

static inline float det_float33( float33_t a )
{
    return a.m00*(a.m11*a.m22 - a.m12*a.m21)
         - a.m01*(a.m10*a.m22 - a.m12*a.m20)
         - a.m02*(a.m10*a.m21 - a.m11*a.m20);
}

static inline float44_t invert_float44( float44_t const a, float const det )
{
    float44_t r;
    float idet;

    assert(isnormal(det));

    idet = 1.f / det;

    r.m00 = ( (a.m11*a.m22-a.m12*a.m21)*a.m33+(a.m13*a.m21-a.m11*a.m23)*a.m32+(a.m12*a.m23-a.m13*a.m22)*a.m31 ) * idet;
    r.m01 = ( (a.m02*a.m21-a.m01*a.m22)*a.m33+(a.m01*a.m23-a.m03*a.m21)*a.m32+(a.m03*a.m22-a.m02*a.m23)*a.m31 ) * idet;
    r.m02 = ( (a.m01*a.m12-a.m02*a.m11)*a.m33+(a.m03*a.m11-a.m01*a.m13)*a.m32+(a.m02*a.m13-a.m03*a.m12)*a.m31 ) * idet;
    r.m03 = ( (a.m02*a.m11-a.m01*a.m12)*a.m23+(a.m01*a.m13-a.m03*a.m11)*a.m22+(a.m03*a.m12-a.m02*a.m13)*a.m21 ) * idet;

    r.m10 = ( (a.m12*a.m20-a.m10*a.m22)*a.m33+(a.m10*a.m23-a.m13*a.m20)*a.m32+(a.m13*a.m22-a.m12*a.m23)*a.m30 ) * idet;
    r.m11 = ( (a.m00*a.m22-a.m02*a.m20)*a.m33+(a.m03*a.m20-a.m00*a.m23)*a.m32+(a.m02*a.m23-a.m03*a.m22)*a.m30 ) * idet;
    r.m12 = ( (a.m02*a.m10-a.m00*a.m12)*a.m33+(a.m00*a.m13-a.m03*a.m10)*a.m32+(a.m03*a.m12-a.m02*a.m13)*a.m30 ) * idet;
    r.m13 = ( (a.m00*a.m12-a.m02*a.m10)*a.m23+(a.m03*a.m10-a.m00*a.m13)*a.m22+(a.m02*a.m13-a.m03*a.m12)*a.m20 ) * idet;

    r.m20 = ( (a.m10*a.m21-a.m11*a.m20)*a.m33+(a.m13*a.m20-a.m10*a.m23)*a.m31+(a.m11*a.m23-a.m13*a.m21)*a.m30 ) * idet;
    r.m21 = ( (a.m01*a.m20-a.m00*a.m21)*a.m33+(a.m00*a.m23-a.m03*a.m20)*a.m31+(a.m03*a.m21-a.m01*a.m23)*a.m30 ) * idet;
    r.m22 = ( (a.m00*a.m11-a.m01*a.m10)*a.m33+(a.m03*a.m10-a.m00*a.m13)*a.m31+(a.m01*a.m13-a.m03*a.m11)*a.m30 ) * idet;
    r.m23 = ( (a.m01*a.m10-a.m00*a.m11)*a.m23+(a.m00*a.m13-a.m03*a.m10)*a.m21+(a.m03*a.m11-a.m01*a.m13)*a.m20 ) * idet;

    r.m30 = ( (a.m11*a.m20-a.m10*a.m21)*a.m32+(a.m10*a.m22-a.m12*a.m20)*a.m31+(a.m12*a.m21-a.m11*a.m22)*a.m30 ) * idet;
    r.m31 = ( (a.m00*a.m21-a.m01*a.m20)*a.m32+(a.m02*a.m20-a.m00*a.m22)*a.m31+(a.m01*a.m22-a.m02*a.m21)*a.m30 ) * idet;
    r.m32 = ( (a.m01*a.m10-a.m00*a.m11)*a.m32+(a.m00*a.m12-a.m02*a.m10)*a.m31+(a.m02*a.m11-a.m01*a.m12)*a.m30 ) * idet;
    r.m33 = ( (a.m00*a.m11-a.m01*a.m10)*a.m22+(a.m02*a.m10-a.m00*a.m12)*a.m21+(a.m01*a.m12-a.m02*a.m11)*a.m20 ) * idet;

    return r;
}

static inline float33_t invert_float33( float33_t const a, float const det ) {
    float33_t r;

    assert(isnormal(det));

    float const idet = 1.0f / det;

    r.m00 =  (a.m00*a.m11 - a.m12*a.m21) * idet;
    r.m01 = -(a.m01*a.m22 - a.m02*a.m21) * idet;
    r.m02 =  (a.m01*a.m12 - a.m02*a.m11) * idet;

    r.m10 = -(a.m10*a.m22 - a.m12*a.m20) * idet;
    r.m11 =  (a.m00*a.m22 - a.m02*a.m20) * idet;
    r.m12 = -(a.m00*a.m12 - a.m02*a.m10) * idet;

    r.m20 =  (a.m10*a.m21 - a.m11*a.m20) * idet;
    r.m21 = -(a.m00*a.m21 - a.m01*a.m20) * idet;
    r.m22 =  (a.m00*a.m11 - a.m01*a.m10) * idet;

    return r;
}

static inline float44_t rotation_common(float3_t axis_a, float const c, float const s, float const t)
{
#ifdef __cplusplus
	float44_t tmp;
	tmp.m00 = t*axis_a.x*axis_a.x + c; tmp.m01 = t*axis_a.x*axis_a.y - axis_a.z*s; tmp.m02 = t*axis_a.x*axis_a.z + axis_a.y*s; tmp.m03 = 0.f;
	tmp.m10 = t*axis_a.x*axis_a.y + axis_a.z*s; tmp.m11 = t*axis_a.y*axis_a.y + c; tmp.m12 = t*axis_a.y*axis_a.z - axis_a.x*s; tmp.m13 = 0.f;
	tmp.m20 = t*axis_a.x*axis_a.z - axis_a.y*s; tmp.m21 = t*axis_a.y*axis_a.z + axis_a.x*s; tmp.m22 = t*axis_a.z*axis_a.z + c; tmp.m23 = 0.f;
	tmp.m30 = 0.f; tmp.m31 = 0.f; tmp.m32 = 0.f; tmp.m33 = 1.f;
	return tmp;
#else
	return (float44_t)
	{
		t*axis_a.x*axis_a.x + c, t*axis_a.x*axis_a.y - axis_a.z*s, t*axis_a.x*axis_a.z + axis_a.y*s, 0.f,
		t*axis_a.x*axis_a.y + axis_a.z*s, t*axis_a.y*axis_a.y + c, t*axis_a.y*axis_a.z - axis_a.x*s, 0.f,
		t*axis_a.x*axis_a.z - axis_a.y*s, t*axis_a.y*axis_a.z + axis_a.x*s, t*axis_a.z*axis_a.z + c, 0.f,
		0.f, 0.f, 0.f, 1.f
	};
#endif
}

static inline float33_t rotation_float33       ( float angle_a ) {
    return rotation_sincos_float33(sinf(angle_a), cosf(angle_a));
}

static inline float33_t rotation_sincos_float33( float sin_angle, float cos_angle ) {
    return (float33_t) {
        cos_angle, -sin_angle, 0.0f,
        sin_angle, cos_angle, 0.0f,
        0.0f, 0.0f, 1.0f
    };
}

static inline float44_t rotation_float44( float3_t axis_a, float angle_a )
{
	float const c = cosf(angle_a);
	float const s = sinf(angle_a);
	float const t = 1 - c;

	return rotation_common(axis_a, c, s, t);
}

static inline float44_t rotation_sincos_float44( float3_t axis_a, float sin_angle_a, float cos_angle_a ) {
	float const c = cos_angle_a;
	float const s = sin_angle_a;
	float const t = 1 - c;

	return rotation_common(axis_a, c, s, t);
}

static inline float33_t cotranslation_float33    ( float2_t vector ) {
  return (float33_t) {
      1.0f, 0.0f, 0.0f,
      0.0f, 1.0f, 0.0f,
      vector.x, vector.y, 1.0f,
  };
}

static inline float33_t translation_float33  ( float2_t vector ) {
  return (float33_t) {
      1.0f, 0.0f, vector.x,
      0.0f, 1.0f, vector.y,
      0.0f, 0.0f, 1.0f
  };
}

static inline float44_t cotranslation_float44 ( float3_t vector ) {
  return (float44_t) {
      1.0f, 0.0f, 0.0f, 0.0f,
      0.0f, 1.0f, 0.0f, 0.0f,
      0.0f, 0.0f, 1.0f, 0.0f,
      vector.x, vector.y, vector.z, 1.0f
  };
}

static inline float44_t translation_float44 ( float3_t vector ) {
  return (float44_t) {
      1.0f, 0.0f, 0.0f, vector.x,
      0.0f, 1.0f, 0.0f, vector.y,
      0.0f, 0.0f, 1.0f, vector.z,
      0.0f, 0.0f, 0.0f, 1.0f
  };
}

static float33_t scale_uniform_float33  ( float s ) {
    return scale_float33((float2_t) { .x = s, .y = s });
}

static float33_t scale_float33( float2_t v ) {
    return (float33_t) {
        v.x, 0.0f, 0.0f,
        0.0f, v.y, 0.0f,
        0.0f, 0.0f, 1.0f
    };
}

static float44_t scale_float44( float3_t vector ) {
    return (float44_t) {
        vector.x, 0.0f,     0.0f,     0.0f,
        0.0f,     vector.y, 0.0f,     0.0f,
        0.0f,     0.0f,     vector.z, 0.0f,
        0.0f,     0.0f,     0.0f,     1.0f
    };
}

static bool equal_float44  ( float44_t const a, float44_t const b) {
	return is_equal(a.m00, b.m00) &&
	       is_equal(a.m01, b.m01) &&
	       is_equal(a.m02, b.m02) &&
	       is_equal(a.m03, b.m03) &&
	       is_equal(a.m10, b.m10) &&
	       is_equal(a.m11, b.m11) &&
	       is_equal(a.m12, b.m12) &&
	       is_equal(a.m13, b.m13) &&
	       is_equal(a.m20, b.m20) &&
	       is_equal(a.m21, b.m21) &&
	       is_equal(a.m22, b.m22) &&
	       is_equal(a.m23, b.m23) &&
	       is_equal(a.m30, b.m30) &&
	       is_equal(a.m31, b.m31) &&
	       is_equal(a.m32, b.m32) &&
	       is_equal(a.m33, b.m33);
}

static inline bool project_float3(float3_t a, float3_t b, float3_t *r) {
    float const s = dot_float3(b, b);
    if (s < EPSILON) {
        return false;
    }

    if (r) *r = scale_float3(b, (dot_float3(a, b) / s));
    return true;
}

static bool gram_schmidt_column_float33(float33_t m, float33_t *r) {
    // r can be NULL, in this case this function only checks if the Gram-Schmidt process is possible
    float3_t u[3], p;

    for (int i=0; i<3; ++i) {
        float3_t const a = (float3_t) { .x = m.m[0][i], .y = m.m[1][i], .z = m.m[2][i] };
        u[i] = a;

        for (int j=0; j<i; ++j) {
            if (!project_float3(a, u[j], &p)) {
                return false;
            }
            u[i] = sub_float3(u[i], p);
        }

        float const l = length_float3(u[i]);
        if (l < EPSILON) {
            return false;
        }

        float3_t const e = scale_float3(u[i], 1.0f / l);
        if (r) {
            r->m[0][i] = e.x;
            r->m[1][i] = e.y;
            r->m[2][i] = e.z;
        }
    }
    
    //printf_matrix33(m);

    return true;
}

static bool qr_float33(float33_t a, float33_t *q_out, float33_t *r_out) {
    float33_t q, r;

    if (!gram_schmidt_column_float33(a, &q)) {
        return false;
    }

    for (int i=0; i<3; ++i) {
        float3_t const ai = (float3_t) { .x = a.m[0][i], .y = a.m[1][i], .z = a.m[2][i] };
        for (int j=0; j<=i; ++j) {
            float3_t const ei = (float3_t) { .x = q.m[0][j], .y = q.m[1][j], .z = q.m[2][j] };
            r.m[j][i] = dot_float3(ai,ei);
        }
        for (int j=i+1; j<3; ++j) {
            r.m[j][i] = 0.0f;
        }
    }

    if (q_out) *q_out = q;
    if (r_out) *r_out = r;
    return true;
}

static inline float eigenvalue_qr_float33_helper_eval(float33_t const m) {
    float const a = m.m[1][0];
    float const b = m.m[2][0];
    float const c = m.m[2][1];

    return (a*a + b*b + c*c);
}

static bool eigenvalue_qr_float33(float33_t m, float epsilon, size_t iter_limit, float out[3]) {
    float33_t q, r;

    for (size_t i=0; i<iter_limit; ++i) {
        if (eigenvalue_qr_float33_helper_eval(m) < epsilon) {
            if (out) {
                float a = m.m[0][0];
                float b = m.m[1][1];
                float c = m.m[2][2];

                float temp;

                if (a < c) {
                    temp = a;
                    a = c;
                    c = temp;
                }

                if (a < b) {
                    temp = a;
                    a = b;
                    b = temp;
                }

                if (b < c) {
                    temp = b;
                    b = c;
                    c = temp;
                }

                out[0] = a;
                out[1] = b;
                out[2] = c;
            }

            return true;
        }

        if (!qr_float33(m, &q, &r)) {
            return false;
        }

        m = mul_float33(r, q);
    }

    return false;
}

static bool eigenvector_symmetric_float33(float33_t m, float eigenvalue, float3_t *out, int *out_multiplicity) {
    float33_t n = m;
    for (int i=0; i<3; ++i) {
        n.m[i][i] -= eigenvalue;
    }

    //printf_matrix33(m);
    //printf_matrix33(n);

    float3_t col0 = { n.m[0][0], n.m[1][0], n.m[2][0] };
    float3_t col1 = { n.m[0][1], n.m[1][1], n.m[2][1] };

    float3_t v = cross_float3(col0, col1);

    if (quadrance_float3(v) < EPSILON) {
        float3_t col2 = { n.m[0][2], n.m[1][2], n.m[2][2] };
        v = cross_float3(col0, col2);

        if (quadrance_float3(v) < EPSILON) {
            v = cross_float3(col1, col2);

            if (quadrance_float3(v) < EPSILON) {
                // eigenvalue is of multiplicity 2 or 3

                bool const col0_nonzero = n.m[0][0] != 0.0f || n.m[0][1] != 0.0f || n.m[0][2] != 0.0f;
                bool const col1_nonzero = n.m[0][0] != 0.0f || n.m[0][1] != 0.0f || n.m[0][2] != 0.0f;
                bool const col2_nonzero = n.m[0][0] != 0.0f || n.m[0][1] != 0.0f || n.m[0][2] != 0.0f;

                if (col0_nonzero) {
                    return false;
                } else if (col1_nonzero) {
                } else if (col2_nonzero) {
                    return false;
                } else {
                    if (out) {
                        // this means the M matrix is I * eigenvalue,
                        // so every vector is an eigenvector with that eigenvalue,
                        // so we choose any set that's orthonormal
                        out[0] = (float3_t) { 1.0f, 0.0f, 0.0f };
                        out[1] = (float3_t) { 0.0f, 1.0f, 0.0f };
                        out[2] = (float3_t) { 0.0f, 0.0f, 1.0f };
                    }
                    if (out_multiplicity) *out_multiplicity = 3;
                    return true;
                }
            }
        }
    }

    if (out) *out = normal_float3(v);
    if (out_multiplicity) *out_multiplicity = 1;
    return true;
}

static bool eigenvectors_symmetric_float33(float33_t m, float value[3], float3_t vector[3]) {
    float3_t temp[3];
    int ix = 0;

    while (ix < 3) {
        int multiplicity = 1;
        if (!eigenvector_symmetric_float33(m, value[ix], temp, &multiplicity)) {
            return false;
        }

        int const ix_next = ix + multiplicity;
        if (ix_next > 3) {
            return false;
        }

        memcpy(vector+ix, temp, sizeof(float3_t) * multiplicity);
        ix = ix_next;
	}

    return true;
}

static bool eigen_symmetric_float33(float33_t m, float value[3], float3_t vector[3]) {
    float val[3];
	if (!eigenvalue_qr_float33(m, 0.00001, 10000, val)) {
		return false;
	}

    if (value) {
        memcpy(value, val, sizeof(float) * 3);
    }

    float3_t vec[3];

    if (!eigenvectors_symmetric_float33(m, val, vec)) {
        return false;
    }

    if (vector) {
        memcpy(vector, vec, sizeof(float3_t) * 3);
    }

    return true;
}

/*
 * Commented out because "return Q, R" is not C syntax for returning both Q and R.
 * Should be rewritten to return a structure.
 * Also, it would be better to name it as decomposition_qr or something like this,
 * to indicate better what it is doing.
 *
static inline float44_t decomposition( float44_t A ) {
    float44_t Q;
    float44_t R;
    int i, j, k, n;
    n = 4;
    for (j = 0; j < n; j++) {
        for (i = 0; i < n; i++) {
            Q.m[i][j] = 0.0f;
        }
    }
    for (j = 0; j < n; j++) {
        for (i = 0; i < j; i++)
            for (k = 0; k < n; k++) R.m[i][j] = R.m[i][j] + A.m[k][j] * Q.m[k][i];
        for (k = 0; k < n; k++) {
            Q.m[k][j] = A.m[k][j];
            for (i = 0; i < j; i++) Q.m[k][j] = Q.m[k][j] - R.m[i][j] * Q.m[k][i];
        }
        for (k = 0; k < n; k++) R.m[j][j] = R.m[j][j] + Q.m[k][j] * Q.m[k][j];
        R.m[j][j] = sqrt(R.m[j][j]);
        for (k = 0; k < n; k++) Q.m[k][j] = Q.m[k][j]/R.m[j][j];
    }
    return Q, R;
}
*/

static inline float44_t gram_schmidt_float44( float44_t A ) {
    int       n = 4;
    bool      test = true;
	//    double    table[n]; unused?
    float     var = 0;
    float     varb = 0;
    int       temp = 1;
    int       i;
    int       j;
    float44_t B;
    for(i = 0; i < n; i++) {
        for(j = 0; j < n; j++) {
            B.m[i][j] = A.m[i][j];
        }
    }
    for(j = 1; j < n; j++) {
        do {
            var = scalar_product_fun_for_matrices(B,j,test, temp);
            test = false;
            varb = scalar_product_fun_for_matrices(B,j,test, temp);
            for(i = 0; i < n; i++) {
                B.m[i][j] = B.m[i][j] - (var/varb) * B.m[i][j-temp];
                test = true;
            }
            temp = temp-1;
        }
        while (temp > 0);
        temp = j;
        temp ++;
    }
    return B;
}

static float44_t add_float44( float44_t const a, float44_t const b) {
	float44_t result;
	for (int i=0; i<4; ++i) {
		for (int j=0; j<4; ++j) {
			result.m[i][j] = a.m[i][j] + b.m[i][j];
		}
	}
	return result;
}

static bool basis_change_float44( float3_t a, float3_t b, float3_t c, float44_t *out, float44_t *out_inv) {
    float44_t basis_change_inv;

    basis_change_inv.m[0][0] = a.x;
    basis_change_inv.m[0][1] = a.y;
    basis_change_inv.m[0][2] = a.z;
    basis_change_inv.m[0][3] = 0.0f;

    basis_change_inv.m[1][0] = b.x;
    basis_change_inv.m[1][1] = b.y;
    basis_change_inv.m[1][2] = b.z;
    basis_change_inv.m[1][3] = 0.0f;

    basis_change_inv.m[2][0] = c.x;
    basis_change_inv.m[2][1] = c.y;
    basis_change_inv.m[2][2] = c.z;
    basis_change_inv.m[2][3] = 0.0f;

    basis_change_inv.m[3][0] = 0.0f;
    basis_change_inv.m[3][1] = 0.0f;
    basis_change_inv.m[3][2] = 0.0f;
    basis_change_inv.m[3][3] = 1.0f;

    float const basis_change_inv_det = det_float44(basis_change_inv);

    if (!isnormal(basis_change_inv_det)) {
        return false;
    }

    if (out) {
        *out = invert_float44(basis_change_inv, basis_change_inv_det);
    }

    if (out_inv) {
        *out_inv = basis_change_inv;
    }

    return true;
}

static inline double44_t rows_double44(double4_t const a, double4_t const b, double4_t const c, double4_t const d) {
    return (double44_t) {
        .m00 = a.x, .m01 = a.y, .m02 = a.z, .m03 = a.w,
        .m10 = b.x, .m11 = b.y, .m12 = b.z, .m13 = b.w,
        .m20 = c.x, .m21 = c.y, .m22 = c.z, .m23 = c.w,
        .m30 = d.x, .m31 = d.y, .m32 = d.z, .m33 = d.w
    };
}

static inline double44_t cols_double44(double4_t const a, double4_t const b, double4_t const c, double4_t const d) {
    return (double44_t) {
        .m00 = a.x, .m10 = a.y, .m20 = a.z, .m30 = a.w,
        .m01 = b.x, .m11 = b.y, .m21 = b.z, .m31 = b.w,
        .m02 = c.x, .m12 = c.y, .m22 = c.z, .m32 = c.w,
        .m03 = d.x, .m13 = d.y, .m23 = d.z, .m33 = d.w
    };
}

static inline double4_t transform_double44( double44_t const a, double4_t const b ) {
    return (double4_t) {
        .x = a.m00 * b.x + a.m01 * b.y + a.m02 * b.z + a.m03 * b.w,
        .y = a.m10 * b.x + a.m11 * b.y + a.m12 * b.z + a.m13 * b.w,
        .z = a.m20 * b.x + a.m21 * b.y + a.m22 * b.z + a.m23 * b.w,
        .w = a.m30 * b.x + a.m31 * b.y + a.m32 * b.z + a.m33 * b.w
    };
}

static inline double4_t cotransform_double44( double44_t const a, double4_t const b ) {
    return (double4_t) {
        .x = a.m00 * b.x + a.m10 * b.y + a.m20 * b.z + a.m30 * b.w,
        .y = a.m01 * b.x + a.m11 * b.y + a.m21 * b.z + a.m31 * b.w,
        .z = a.m02 * b.x + a.m12 * b.y + a.m22 * b.z + a.m32 * b.w,
        .w = a.m03 * b.x + a.m13 * b.y + a.m23 * b.z + a.m33 * b.w
    };
}

static double44_t mul_double44( double44_t const a, double44_t const b) {
    return (double44_t) {
        .m00 = a.m03*b.m30+a.m02*b.m20+a.m01*b.m10+a.m00*b.m00,
        .m01 = a.m03*b.m31+a.m02*b.m21+a.m01*b.m11+a.m00*b.m01,
        .m02 = a.m03*b.m32+a.m02*b.m22+a.m01*b.m12+a.m00*b.m02,
        .m03 = a.m03*b.m33+a.m02*b.m23+a.m01*b.m13+a.m00*b.m03,

        .m10 = a.m13*b.m30+a.m12*b.m20+a.m11*b.m10+a.m10*b.m00,
        .m11 = a.m13*b.m31+a.m12*b.m21+a.m11*b.m11+a.m10*b.m01,
        .m12 = a.m13*b.m32+a.m12*b.m22+a.m11*b.m12+a.m10*b.m02,
        .m13 = a.m13*b.m33+a.m12*b.m23+a.m11*b.m13+a.m10*b.m03,

        .m20 = a.m23*b.m30+a.m22*b.m20+a.m21*b.m10+a.m20*b.m00,
        .m21 = a.m23*b.m31+a.m22*b.m21+a.m21*b.m11+a.m20*b.m01,
        .m22 = a.m23*b.m32+a.m22*b.m22+a.m21*b.m12+a.m20*b.m02,
        .m23 = a.m23*b.m33+a.m22*b.m23+a.m21*b.m13+a.m20*b.m03,

        .m30 = a.m33*b.m30+a.m32*b.m20+a.m31*b.m10+a.m30*b.m00,
        .m31 = a.m33*b.m31+a.m32*b.m21+a.m31*b.m11+a.m30*b.m01,
        .m32 = a.m33*b.m32+a.m32*b.m22+a.m31*b.m12+a.m30*b.m02,
        .m33 = a.m33*b.m33+a.m32*b.m23+a.m31*b.m13+a.m30*b.m03
    };
}
static double44_t mul_scalar_double44( double44_t const a, double const s ) {
    return (double44_t) {
        .m00 = a.m00 * s,
        .m01 = a.m10 * s,
        .m02 = a.m20 * s,
        .m03 = a.m30 * s,

        .m10 = a.m01 * s,
        .m11 = a.m11 * s,
        .m12 = a.m21 * s,
        .m13 = a.m31 * s,

        .m20 = a.m02 * s,
        .m21 = a.m12 * s,
        .m22 = a.m22 * s,
        .m23 = a.m32 * s,

        .m30 = a.m03 * s,
        .m31 = a.m13 * s,
        .m32 = a.m23 * s,
        .m33 = a.m33 * s
    };
}

static double44_t transpose_double44 ( double44_t const a ) {
    return (double44_t) {
        .m00 = a.m00,
        .m01 = a.m10,
        .m02 = a.m20,
        .m03 = a.m30,

        .m10 = a.m01,
        .m11 = a.m11,
        .m12 = a.m21,
        .m13 = a.m31,

        .m20 = a.m02,
        .m21 = a.m12,
        .m22 = a.m22,
        .m23 = a.m32,

        .m30 = a.m03,
        .m31 = a.m13,
        .m32 = a.m23,
        .m33 = a.m33
    };
}

static double     det_double44       ( double44_t const a ) {
    return  a.m00*(a.m11*(a.m22*a.m33-a.m23*a.m32)-a.m12*(a.m21*a.m33-a.m23*a.m31)+a.m13*(a.m21*a.m32-a.m22*a.m31))
           -a.m01*(a.m10*(a.m22*a.m33-a.m23*a.m32)-a.m12*(a.m20*a.m33-a.m23*a.m30)+a.m13*(a.m20*a.m32-a.m22*a.m30))
           +a.m02*(a.m10*(a.m21*a.m33-a.m23*a.m31)-a.m11*(a.m20*a.m33-a.m23*a.m30)+a.m13*(a.m20*a.m31-a.m21*a.m30))
           -a.m03*(a.m10*(a.m21*a.m32-a.m22*a.m31)-a.m11*(a.m20*a.m32-a.m22*a.m30)+a.m12*(a.m20*a.m31-a.m21*a.m30));
}
static double44_t invert_double44    ( double44_t const a, double const det) {
    assert(isnormal(det));

    double idet = 1.f / det;

    return (double44_t) {
        .m00 = ( (a.m11*a.m22-a.m12*a.m21)*a.m33+(a.m13*a.m21-a.m11*a.m23)*a.m32+(a.m12*a.m23-a.m13*a.m22)*a.m31 ) * idet,
        .m01 = ( (a.m02*a.m21-a.m01*a.m22)*a.m33+(a.m01*a.m23-a.m03*a.m21)*a.m32+(a.m03*a.m22-a.m02*a.m23)*a.m31 ) * idet,
        .m02 = ( (a.m01*a.m12-a.m02*a.m11)*a.m33+(a.m03*a.m11-a.m01*a.m13)*a.m32+(a.m02*a.m13-a.m03*a.m12)*a.m31 ) * idet,
        .m03 = ( (a.m02*a.m11-a.m01*a.m12)*a.m23+(a.m01*a.m13-a.m03*a.m11)*a.m22+(a.m03*a.m12-a.m02*a.m13)*a.m21 ) * idet,

        .m10 = ( (a.m12*a.m20-a.m10*a.m22)*a.m33+(a.m10*a.m23-a.m13*a.m20)*a.m32+(a.m13*a.m22-a.m12*a.m23)*a.m30 ) * idet,
        .m11 = ( (a.m00*a.m22-a.m02*a.m20)*a.m33+(a.m03*a.m20-a.m00*a.m23)*a.m32+(a.m02*a.m23-a.m03*a.m22)*a.m30 ) * idet,
        .m12 = ( (a.m02*a.m10-a.m00*a.m12)*a.m33+(a.m00*a.m13-a.m03*a.m10)*a.m32+(a.m03*a.m12-a.m02*a.m13)*a.m30 ) * idet,
        .m13 = ( (a.m00*a.m12-a.m02*a.m10)*a.m23+(a.m03*a.m10-a.m00*a.m13)*a.m22+(a.m02*a.m13-a.m03*a.m12)*a.m20 ) * idet,

        .m20 = ( (a.m10*a.m21-a.m11*a.m20)*a.m33+(a.m13*a.m20-a.m10*a.m23)*a.m31+(a.m11*a.m23-a.m13*a.m21)*a.m30 ) * idet,
        .m21 = ( (a.m01*a.m20-a.m00*a.m21)*a.m33+(a.m00*a.m23-a.m03*a.m20)*a.m31+(a.m03*a.m21-a.m01*a.m23)*a.m30 ) * idet,
        .m22 = ( (a.m00*a.m11-a.m01*a.m10)*a.m33+(a.m03*a.m10-a.m00*a.m13)*a.m31+(a.m01*a.m13-a.m03*a.m11)*a.m30 ) * idet,
        .m23 = ( (a.m01*a.m10-a.m00*a.m11)*a.m23+(a.m00*a.m13-a.m03*a.m10)*a.m21+(a.m03*a.m11-a.m01*a.m13)*a.m20 ) * idet,

        .m30 = ( (a.m11*a.m20-a.m10*a.m21)*a.m32+(a.m10*a.m22-a.m12*a.m20)*a.m31+(a.m12*a.m21-a.m11*a.m22)*a.m30 ) * idet,
        .m31 = ( (a.m00*a.m21-a.m01*a.m20)*a.m32+(a.m02*a.m20-a.m00*a.m22)*a.m31+(a.m01*a.m22-a.m02*a.m21)*a.m30 ) * idet,
        .m32 = ( (a.m01*a.m10-a.m00*a.m11)*a.m32+(a.m00*a.m12-a.m02*a.m10)*a.m31+(a.m02*a.m11-a.m01*a.m12)*a.m30 ) * idet,
        .m33 = ( (a.m00*a.m11-a.m01*a.m10)*a.m22+(a.m02*a.m10-a.m00*a.m12)*a.m21+(a.m01*a.m12-a.m02*a.m11)*a.m20 ) * idet
    };
}

static bool equal_double44( double44_t const a, double44_t const b ) {
    return is_equal(a.m00, b.m00) && /* TODO is this enough as it is float based equality? */
          is_equal(a.m01, b.m01) &&
          is_equal(a.m02, b.m02) &&
          is_equal(a.m03, b.m03) &&
          is_equal(a.m10, b.m10) &&
          is_equal(a.m11, b.m11) &&
          is_equal(a.m12, b.m12) &&
          is_equal(a.m13, b.m13) &&
          is_equal(a.m20, b.m20) &&
          is_equal(a.m21, b.m21) &&
          is_equal(a.m22, b.m22) &&
          is_equal(a.m23, b.m23) &&
          is_equal(a.m30, b.m30) &&
          is_equal(a.m31, b.m31) &&
          is_equal(a.m32, b.m32) &&
          is_equal(a.m33, b.m33);
}

static double44_t add_double44( double44_t const a, double44_t const b)  {
    return (double44_t) {
        .m00 = a.m00 + b.m00,
        .m01 = a.m01 + b.m01,
        .m02 = a.m02 + b.m02,
        .m03 = a.m03 + b.m03,

        .m10 = a.m10 + b.m10,
        .m11 = a.m11 + b.m11,
        .m12 = a.m12 + b.m12,
        .m13 = a.m13 + b.m13,

        .m20 = a.m20 + b.m20,
        .m21 = a.m21 + b.m21,
        .m22 = a.m22 + b.m22,
        .m23 = a.m23 + b.m23,

        .m30 = a.m30 + b.m30,
        .m31 = a.m31 + b.m31,
        .m32 = a.m32 + b.m32,
        .m33 = a.m33 + b.m33
    };
}

static float44_t rotation_of_vector_onto(float3_t const a, float3_t const b) {
	float c = dot_float3(a,b);

	if (is_equal(1.0f, c)) {
		/* If the vectors are already pointing in the same direction, do nothing. */
		return identity_sc;
	} else if (is_equal(-1.0f, c)) {
		/* If the vectors are pointing in opposite directions,
		 * return a rotation by 180 degrees about some axis perpendicular to a.
		 */

		float3_t const perp_a = (a.x || a.y) ? (float3_t) { .x = -a.y, .y = a.x, .z = 0.0f } : (float3_t) { .x = a.x, .y = 0.0f, .z = 0.0f};

		return rotation_float44(normal_float3(perp_a), M_PI);
	} else {
		/*
		 * The formula for rotation matrix that rotates unit vector A onto unit vector B
		 * (courtesy of linear algebra stack exchange) is:
		 *
		 * V = A x B                   // cross product
		 * s = ||V||                   // length_float3
		 * c = A * B                   // dot product
		 *
		 * (this means matrix, not determinant)
		 * [V]_x = |    0, -V_z,  V_y |
		 *         |  V_z,    0, -V_x |
		 *         | -V_y,  V_x,    0 |
		 * ^- skew-symmetric cross-product matrix of V
		 *
		 * R = I + [V]_x + ([V]_x)^2 * (1-c)/s^2
		 *
		 * where (1-c)/s^2 can be simplified to 1/(1+c)
		 * because c = cos(angle(a,b)) and s = sin(angle(a,b))
		 *
		 * Note: result is transposed to match conventions used in the rest of mathutils.
		 */

		float3_t const v = cross_float3(a,b);
		float44_t const u = (float44_t) {
			0.0f,   v.z,  -v.y, 0.0f,
			-v.z,   0.0f,  v.x, 0.0f,
			 v.y,  -v.x,  0.0f, 0.0f,
			0.0f,   0.0f, 0.0f, 0.0f
		};

		float44_t const u_squared = mul_scalar_float44(mul_float44(u,u), 1.0f/(1.0f + c));

		return add_float44(identity_sc, add_float44(u, u_squared));
	}
}


static float line_point_distance(float3_t d /* line direction */, float3_t o /* line origin */, float3_t p /* point */) {
    /* orthonormal coordinates assumption */
    return length_float3(cross_float3(sub_float3(p, o), d));
}

static float skew_lines_distance(float3_t da, float3_t oa, float3_t db, float3_t ob) {
    float3_t c = sub_float3  (ob, oa);
    float3_t axb = cross_float3(da, db);
    float    laxb = length_float3      (axb);

    /* Skew lines */
    if (laxb > EPSILON) {
        return fabsf(dot_float3(c, axb)) / laxb;
    }

    /* Parallel lines */
    return FLT_MAX;
}

static int skew_lines_nearest(float3_t da, float3_t oa,float3_t db, float3_t ob, float3_t *ca, float3_t *cb) {
    float3_t n = cross_float3(da, db);

    if (length_float3(n) > EPSILON_BIG) {
        float3_t na = cross_float3(da, n);
        float3_t nb = cross_float3(db, n);
        *ca = add_float3(oa, scale_float3(da, dot_float3(sub_float3(ob, oa), nb) / dot_float3(da, nb)));
        *cb = add_float3(ob, scale_float3(db, dot_float3(sub_float3(oa, ob), na) / dot_float3(db, na)));
        return 0;
    }

    return -1;
}

static float lines_distance(float3_t da, float3_t oa, float3_t db, float3_t ob) {
    float d = skew_lines_distance(da, oa, db, ob);

    /* Parallel lines */
    if (FLT_MAX == d) {
        return line_point_distance(da, oa, ob);
    }

    /* Skew lines */
    return d;
}

static float3_t lines_intersection (float3_t da, float3_t oa, float3_t db, float3_t ob) {
    float3_t const c     = sub_float3  (ob,  oa);
    float3_t const cxb   = cross_float3(c,   db);
    float3_t const axb   = cross_float3(da,  db);
    float    const laxb2 = dot_float3  (axb, axb);

    float const l = dot_float3(cxb, axb) / laxb2;

    return add_float3(scale_float3(da, l), oa);
}

static void compose_axis_angle(
    float3_t  axis_0_a,   float  angle_0_a,
    float3_t  axis_1_a,   float  angle_1_a,
    float3_t *axis_out_a, float *angle_out_a
) {
    float const c0 = cosf(angle_0_a * 0.5f);
    float const s0 = sinf(angle_0_a * 0.5f);
    float const c1 = cosf(angle_1_a * 0.5f);
    float const s1 = sinf(angle_1_a * 0.5f);

    float const c2  = c0 * c1 - s0 * s1 * dot_float3(axis_0_a, axis_1_a);

    float3_t const s2n = add_float3(
        add_float3(
            scale_float3(axis_0_a, s0 * c1),
            scale_float3(axis_1_a, c0 * s1)
        ),
        scale_float3(cross_float3(axis_0_a, axis_1_a), s0 * s1)
    );

    float const s2 = length_float3(s2n);

    *angle_out_a = atan2(s2, c2) * 2.0;

    if (is_zero(s2)) {
        *axis_out_a  = axis_0_a;
    } else {
        *axis_out_a = normal_float3(s2n);
    }

}

static float44_t orthographic(float r, float l, float t, float b, float f, float n) {
#ifdef __cplusplus
	float44_t tmp;
	tmp.m00 = 2.f / (r - l); tmp.m01 = 0.f; tmp.m02 = 0.f; tmp.m03 = (r + l) / (l - r);
	tmp.m10 = 0.f; tmp.m11 = 2.f / (t - b); tmp.m12 = 0.f; tmp.m13 = (t + b) / (b - t);
	tmp.m20 = 0.f; tmp.m21 = 0.f; tmp.m22 = 2.f / (n - f); tmp.m23 = (f + n) / (n - f);
	tmp.m30 = 0.f; tmp.m31 = 0.f; tmp.m32 = 0.f; tmp.m33 = 1.f;
	return transpose_float44(tmp);
#else
	return transpose_float44((float44_t)
	{
		2.f / (r-l), 0.f, 0.f, -(r+l) / (l-r),
		0.f, 2.f/(t-b), 0.f,   -(t+b) / (b-t),
		0.f, 0.f, -2.f/(n-f),   -(f+n) / (n-f),
		0.f, 0.f, 0.f, 1.f
	});
#endif
}

static inline float scalar_product_fun_for_matrices( float44_t matrix, int column, bool pointer, int previous ) {
    int   n = 4;
    int   i;
    float scalar_product;
    float sum = 0;
    if(pointer == 1) {
        for(i = 0; i < n; i++) {
            scalar_product = matrix.m[i][column]*matrix.m[i][column-previous];
            sum = sum + scalar_product;
        }
    }
    else {
        for(i = 0;i < n; i++) {
            scalar_product = matrix.m[i][column-previous] * matrix.m[i][column-previous];
            sum = sum + scalar_product;
        }
    }
    return sum;
}

static inline float maxf( float a, float b )
{
    return (a > b) ? a : b;
}
static inline float minf( float a, float b )
{
    return (a < b) ? a : b;
}
static inline int maxi( int a, int b )
{
    return (a > b) ? a : b;
}
static inline int mini( int a, int b )
{
    return (a < b) ? a : b;
}
static inline float clampf( float x, float a, float b )
{
    return minf(maxf(a, x), b);
}
static inline int clampi( int x, int a, int b )
{
    return mini(maxi(a, x), b);
}
static bool is_power_of_2(unsigned int n) {
    return (n & (n-1)) == 0;
}
static unsigned int next_power_of_2(unsigned int n) {
    unsigned int result = 1;
    while (result < n) {
        result <<= 1;
    }
    return result;
}
static inline float2_t max_float2( float2_t const a, float2_t const b )
{
    float2_t r = { maxf(a.x, b.x), maxf(a.y, b.y) };
    return r;
}
static inline float2_t min_float2( float2_t const a, float2_t const b )
{
    float2_t r = { minf(a.x, b.x), minf(a.y, b.y) };
    return r;
}
static inline float3_t max_float3( float3_t const a, float3_t const b )
{
    float3_t r = { maxf(a.x, b.x), maxf(a.y, b.y), maxf(a.z, b.z) };
    return r;
}
static inline float3_t min_float3( float3_t const a, float3_t const b )
{
    float3_t r = { minf(a.x, b.x), minf(a.y, b.y), minf(a.z, b.z) };
    return r;
}
static inline float max_of_float3( float3_t const a )
{
    return maxf(a.x, maxf(a.y, a.z));
}
static inline float min_of_float3( float3_t const a )
{
    return minf(a.x, minf(a.y, a.z));
}

static inline float3_t clamp_float3( float3_t const a, float3_t const b, float3_t const c )
{
    float3_t r = { clampf(a.x, b.x, c.x), clampf(a.y, b.y, c.y), clampf(a.z, b.z, c.z) };
    return r;
}

static inline float3_t swizzle_float3( float3_t const a, int const x_swizzle_a, int const y_swizzle_a, int const z_swizzle_a )
{
    float3_t r;
    float s[3] = { a.x, a.y, a.z };
    float t[3] = { 0.f, 0.f, 0.f };

    t[0] = s[x_swizzle_a];
    t[1] = s[y_swizzle_a];
    t[2] = s[z_swizzle_a];

    r.x = t[0];
    r.y = t[1];
    r.z = t[2];

    return r;
}

static inline float3_t plane_point_float3(
    float2_t st_point,
    float3_t plane_origin,
    float3_t plane_s_axis,
    float3_t plane_t_axis
) {
    return
        add_float3(
            add_float3(
                scale_float3(plane_s_axis, st_point.s),
                scale_float3(plane_t_axis, st_point.t)
            ),
            plane_origin
        );
}

static inline float _signumf(float a)
{
    if (a < 0.f) {
        return -1.f;
    }
    return 1.f;
}

static bool lesser_float3(float3_t const a, float3_t const b) {
    return is_lesser(a.x, b.x) && is_lesser(a.y, b.y) && is_lesser(a.z, b.z);
}

static bool greater_float3(float3_t const a, float3_t const b) {
    return is_greater(a.x, b.x) && is_greater(a.y, b.y) && is_greater(a.z, b.z);
}

static bool lesser_equal_float3(float3_t const a, float3_t const b) {
    return is_lesser_equal(a.x, b.x) && is_lesser_equal(a.y, b.y) && is_lesser_equal(a.z, b.z);
}

static bool greater_equal_float3(float3_t const a, float3_t const b) {
    return is_greater_equal(a.x, b.x) && is_greater_equal(a.y, b.y) && is_greater_equal(a.z, b.z);
}

static bool any_equal_float3(float3_t const a, float3_t const b) {
    return is_equal(a.x, b.x) || is_equal(a.y, b.y) || is_equal(a.z, b.z);
}

static bool any_lesser_float3(float3_t const a, float3_t const b) {
    return is_lesser(a.x, b.x) || is_lesser(a.y, b.y) || is_lesser(a.z, b.z);
}

static bool any_greater_float3(float3_t const a, float3_t const b) {
    return is_greater(a.x, b.x) || is_greater(a.y, b.y) || is_greater(a.z, b.z);
}

static bool any_lesser_equal_float3(float3_t const a, float3_t const b) {
    return is_lesser_equal(a.x, b.x) || is_lesser_equal(a.y, b.y) || is_lesser_equal(a.z, b.z);
}

static bool any_greater_equal_float3(float3_t const a, float3_t const b) {
    return is_greater_equal(a.x, b.x) || is_greater_equal(a.y, b.y) || is_greater_equal(a.z, b.z);
}

static bool      equal_float3_eps            ( float3_t const a, float3_t const b, float eps) {
    return is_equal_eps(a.x, b.x, eps) && is_equal_eps(a.y, b.y, eps) && is_equal_eps(a.z, b.z, eps);
}

static bool      equal_float3_ptr_eps        ( float3_t const * const a, float3_t const * const b, float eps) {
    return equal_float3_eps(*a, *b, eps);
}

static bool      lesser_float3_eps           ( float3_t const a, float3_t const b, float eps) {
    return is_lesser_eps(a.x, b.x, eps) && is_lesser_eps(a.y, b.y, eps) && is_lesser_eps(a.z, b.z, eps);
}

static bool      greater_float3_eps          ( float3_t const a, float3_t const b, float eps) {
    return is_greater_eps(a.x, b.x, eps) && is_greater_eps(a.y, b.y, eps) && is_greater_eps(a.z, b.z, eps);
}

static bool      lesser_equal_float3_eps     ( float3_t const a, float3_t const b, float eps) {
    return is_lesser_equal_eps(a.x, b.x, eps) && is_lesser_equal_eps(a.y, b.y, eps) && is_lesser_equal_eps(a.z, b.z, eps);
}

static bool      greater_equal_float3_eps    ( float3_t const a, float3_t const b, float eps) {
    return is_greater_equal_eps(a.x, b.x, eps) && is_greater_equal_eps(a.y, b.y, eps) && is_greater_equal_eps(a.z, b.z, eps);
}

static bool      any_equal_float3_eps        ( float3_t const a, float3_t const b, float eps) {
    return is_equal_eps(a.x, b.x, eps) || is_equal_eps(a.y, b.y, eps) || is_equal_eps(a.z, b.z, eps);
}

static bool      any_lesser_float3_eps       ( float3_t const a, float3_t const b, float eps) {
    return is_lesser_eps(a.x, b.x, eps) || is_lesser_eps(a.y, b.y, eps) || is_lesser_eps(a.z, b.z, eps);
}

static bool      any_greater_float3_eps      ( float3_t const a, float3_t const b, float eps) {
    return is_greater_eps(a.x, b.x, eps) || is_greater_eps(a.y, b.y, eps) || is_greater_eps(a.z, b.z, eps);
}

static bool      any_lesser_equal_float3_eps ( float3_t const a, float3_t const b, float eps) {
    return is_lesser_equal_eps(a.x, b.x, eps) || is_lesser_equal_eps(a.y, b.y, eps) || is_lesser_equal_eps(a.z, b.z, eps);
}

static bool      any_greater_equal_float3_eps( float3_t const a, float3_t const b, float eps) {
    return is_greater_equal_eps(a.x, b.x, eps) || is_greater_equal_eps(a.y, b.y, eps) || is_greater_equal_eps(a.z, b.z, eps);
}

static bool basis_from_triangle(float3_t const a, float3_t const b, float3_t const c, float44_t *m_out) {
    float3_t const u = normal_float3(sub_float3(b, c));
    float3_t const ba = sub_float3(b, a);
    float3_t const ba_proj_u = scale_float3(u, dot_float3(u, ba));
    float3_t const w = normal_float3(sub_float3(ba, ba_proj_u));
    float3_t const v = normal_float3(cross_float3(u, w));

    /*
    formated_log(
        "basis:\nx\ty\tz\n%.3f\t%.3f\t%.3f\n%.3f\t%.3f\t%.3f\n%.3f\t%.3f\t%.3f", u.x, u.y, u.z, w.x,
        w.y, w.z, v.x, v.y, v.z);
*/

    float44_t m = identity_sc;

    if (!basis_change_float44(u, w, v, &m, NULL)) {
        return false;
    }

    m = mul_float44(cotranslation_float44(scale_float3(c, -1.0f)), m);

#ifndef NDEBUG
    if (fabsf(det_float44(m) - 1.0f) > EPSILON_BIG) {
        return false;
    }
#endif // NDEBUG

    if (m_out) {
        *m_out = m;
    }

    return true;
}

static bool centered_basis_from_triangle(float3_t const a, float3_t const b, float3_t const c, float rotation, float44_t *m_out, float44_t *m_inv_out) {
    float44_t m = identity_sc;

    if (!basis_from_triangle(a, b, c, &m)) {
        return false;
    }

    float3_t const a_prim = m_float3(cotransform_float44(m, m_float4(a, 1.0f)));
    float3_t const b_prim = m_float3(cotransform_float44(m, m_float4(b, 1.0f)));
#ifndef NDEBUG
    float3_t const c_prim = m_float3(cotransform_float44(m, m_float4(c, 1.0f)));

/*
    formated_imp(
        "transformed:\nx\ty\tz\n%.3f\t%.3f\t%.3f\n%.3f\t%.3f\t%.3f\n%.3f\t%.3f\t%.3f", a_prim.x,
        a_prim.y, a_prim.z, b_prim.x, b_prim.y, b_prim.z, c_prim.x, c_prim.y, c_prim.z);
*/
#endif // NDEBUG

    assert(fabsf(a_prim.z) < EPSILON_BIG);
    assert((fabsf(b_prim.y) < EPSILON_BIG) && (fabsf(b_prim.z) < 1e-4));
#ifndef NDEBUG
    assert(
        (fabsf(c_prim.x) < 1e-4) && (fabsf(c_prim.y) < EPSILON_BIG) &&
        (fabsf(c_prim.z) < EPSILON_BIG));
#endif // NDEBUG

    float const center_d = -2.0f * b_prim.x * a_prim.y;
    float const a_prim_sq = (a_prim.x * a_prim.x + a_prim.y * a_prim.y) / center_d;
    float const b_prim_sq = (b_prim.x * b_prim.x) / center_d;

    float3_t const center = (float3_t){
        .x = -b_prim_sq * a_prim.y,
        .y = -a_prim_sq * b_prim.x + b_prim_sq * a_prim.x,
        .z = 0.0f,
    };

#ifndef NDEBUG
    assert(fabsf(det_float44(m) - 1.0f) < EPSILON_BIG);
#endif // NDEBUG

    /*
    float44_t const rot_m = mul_float44(
        cotranslation_float44(scale_float3(center, -1)),
            mul_float44(
                rotation_float44(
                    (float3_t){ .x = 0.0f, .y = 0.0f, .z = 1.0f }, rotation),
                cotranslation_float44(center)));
*/

    float44_t const rot_m = mul_float44(
                rotation_float44(
                    (float3_t){ .x = 0.0f, .y = 0.0f, .z = 1.0f }, rotation),
                cotranslation_float44(center));

    float44_t const m_inv = mul_float44(rot_m, invert_float44(m, 1.0f));

#ifndef NDEBUG
        assert(fabsf(det_float44(m_inv) - 1.0f) < EPSILON_BIG);
#endif // NDEBUG

    m = invert_float44(m_inv, 1.0f);

    if (m_out) {
        *m_out = m;
    }

    if (m_inv_out) {
        *m_inv_out = m_inv;
    }

    return true;
}

static float frobenius_float33( float33_t m) {
	double f = 0.0;
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			f += m.m[i][j] * m.m[i][j];

	return (float)sqrt(f);
}

#endif /* __common_mathutils_h */
