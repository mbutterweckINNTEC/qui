/**
 *  \brief  Quaternion algebra and conversion functions.
 *
 *  \file   mathutils.h
 *
 *  \author Marcin Mrotek
 *  \date   2017
 *
 *  \copyright Copyright (c) INNTEC.PL Sp. z o.o. All rights reserved.
 *             Prawa autorskie zarezerwowane - INNTEC.PL Sp. z o.o.
 *
 *  \note INNTEC.PL Sp. z o.o. Company Confidential Information.
 *        Informacje Poufne firmy INNTEC.PL Sp. z o.o.
 *
 */

#ifndef __COMMON_QUATERNION_H
#define __COMMON_QUATERNION_H

#include "mathutils.h"

/* Types */

/** \brief Quaternion structure.
 * \note Can be used as both (a+bi+cj+dk) and scalar-vector representations.
 * (I hope using unions for that isn't undefined behavior).
 */

typedef struct quaternion_t {
    union {
        struct {
            float scalar;
            float3_t vector;
        };
        struct {
            float a, b, c, d;
        };
    };
} quaternion_t;

/** \brief Axis-angle representation for rotations.
 */
typedef struct axis_angle_t {
    float3_t axis;            ///< Axis of rotation.
    float angle;              ///< Angle of rotation.
} axis_angle_t;

static inline bool equal_axis_angle(axis_angle_t a, axis_angle_t b);

/* Interface */

static inline bool equal_quaternion(quaternion_t const, quaternion_t const);
static inline float sign_of_variable(float x);

static inline quaternion_t add_quaternion(quaternion_t const, quaternion_t const);
static inline quaternion_t sub_quaternion(quaternion_t const, quaternion_t const);
static inline quaternion_t mul_quaternion(quaternion_t const, quaternion_t const);
static inline quaternion_t norm_quaternion(quaternion_t const);

static inline quaternion_t negative_angle_quaternion(quaternion_t const);

static inline float dot_quaternion(quaternion_t const, quaternion_t const);
static inline quaternion_t avg_quaternion(quaternion_t const, float, quaternion_t const, float);

static inline quaternion_t conjugate_quaternion(quaternion_t const);
static inline quaternion_t inverse_quaternion(quaternion_t const);
static inline float3_t     transform_quaternion(quaternion_t const, float3_t);  /* for unit quaternions */

static inline quaternion_t pure_quaternion(float3_t const v);
static inline quaternion_t axis_angle_to_quaternion(float3_t const axis, float angle);
static inline axis_angle_t quaternion_to_axis_angle(quaternion_t const);
static inline float44_t    quaternion_to_rotation_matrix(quaternion_t const);
static inline quaternion_t rotation_matrix_to_quaternion(float44_t const, bool *pointer);

static inline float inner_product_quaternion(quaternion_t const, quaternion_t const);

static inline quaternion_t rotate_vector_onto(float3_t const from, float3_t const to);

static inline quaternion_t quaternion_from_euler(float const angle_x, float const angle_y, float const angle_z);
static inline void quaternion_to_euler(quaternion_t q, float *angle_x, float *angle_y, float *angle_z);
/** /brief 				Rotating aabb with quaternion.
*
*	/note 				Function gets two float3_t arguments, which are the min and max of aabb, and one *					quaternion argument by which this function rotate aabb and change *					original float3_t values.
*
*	/param [in] aabb_min		A float3_t structure containing 3-dimensional coordinates of a point which is
* 					the closest corner of aabb.
*
*	/param [in] aabb_min		A float3_t structure containing 3-dimensional coordinates of a point which is
* 					the farthest corner of aabb.
*
*	/param [in] rotate_aabb		A quaternion_t structure containing values of rotation by which aabb is moved.
*
*
*	/return 			Function set values of aabb_min and aabb_max for minimal cube which contains *					rotated aabb.
*/
static inline void rotate_aabb(float3_t *aabb_min, float3_t *aabb_max, quaternion_t rotat_aabb);

/** \brief Swizzle quaternion axis.

    \note Note that input shall be in form: a != b,  a != c, b != c. Also a, b and c
          are in range [0..2] included. This is not checked. More (possibly) at swizzle_float3.

    \param [in] q       The quaternion.
    \param [in] a       What shall go into x component.
    \param [in] b       What shall go into y component.
    \param [in] c       What shall go into z component.

    \return Swizzled quaternion.
*/
static inline quaternion_t quaternion_swizzle(quaternion_t q, int a, int b, int c);

/* Implementations */

static inline bool equal_axis_angle(axis_angle_t a, axis_angle_t b) {
	return is_equal(a.axis.x, b.axis.x) && is_equal(a.axis.y, b.axis.y) && is_equal(a.axis.z, b.axis.z) && is_equal(a.angle, b.angle);
}

static inline bool equal_quaternion(quaternion_t const q, quaternion_t const p) {
    return is_equal(q.a, p.a) && is_equal(q.b, p.b) && is_equal(q.c, p.c) && is_equal(q.d, p.d);
}

static inline quaternion_t add_quaternion(quaternion_t const q, quaternion_t const p) {
    return (quaternion_t) {
        .a = q.a + p.a,
        .b = q.b + p.b,
        .c = q.c + p.c,
        .d = q.d + p.d,
    };
}

static inline quaternion_t sub_quaternion(quaternion_t const q, quaternion_t const p) {
    return (quaternion_t) {
        .a = q.a - p.a,
        .b = q.b - p.b,
        .c = q.c - p.c,
        .d = q.d - p.d,
    };
}

static inline quaternion_t mul_quaternion(quaternion_t const q, quaternion_t const p) {
    return (quaternion_t) {
        .a = (q.a*p.a - q.b*p.b - q.c*p.c - q.d*p.d),
        .b = (q.a*p.b + q.b*p.a + q.c*p.d - q.d*p.c),
        .c = (q.a*p.c - q.b*p.d + q.c*p.a + q.d*p.b),
        .d = (q.a*p.d + q.b*p.c - q.c*p.b + q.d*p.a),
    };
}

static inline float quadrance_quaternion(quaternion_t const q) {
    return q.a*q.a + q.b*q.b + q.c*q.c + q.d*q.d;
}

static inline quaternion_t norm_quaternion(quaternion_t const q) {
    float const dl = sqrtf(quadrance_quaternion(q));

    // left it EXACTLY as here, no EPSILON, dl may be realllly small
    if (isnormal(dl)) {
        float scale = 1.0f/dl;
        return (quaternion_t) {
            .a = q.a*scale,
            .b = q.b*scale,
            .c = q.c*scale,
            .d = q.d*scale,
        };
    } else {
        return (quaternion_t) {
            .a = 1.0f,
            .b = 0.0f,
            .c = 0.0f,
            .d = 0.0f,
        };
   }
}

static inline quaternion_t negative_angle_quaternion(quaternion_t const q) {
    return (quaternion_t){.a = -q.a, .b = q.b, .c = q.c, .d = q.d};
}

static inline float dot_quaternion(quaternion_t const q, quaternion_t const p) {
    return q.b * p.b + q.c * p.c + q.d * p.d;
}

static inline float inner_product_quaternion(quaternion_t const q, quaternion_t const p) {
    return q.a * p.a + q.b * p.b + q.c * p.c + q.d * p.d;
}

static inline quaternion_t avg_quaternion(quaternion_t const q1, float const w1, quaternion_t const q2, float const w2) {
    float dot = dot_quaternion(q1, q2);
    if (is_zero(dot) && is_equal(w1, w2)) {
        /* Not unique so do it other way, propably bad way */
        float _w1 = w1 / (w1+w2);
        float _w2 = w2 / (w1+w2);
        quaternion_t _avg = {
            .a = _w1 * q1.a + _w2 * q2.a,
            .b = _w1 * q1.b + _w2 * q2.b,
            .c = _w1 * q1.c + _w2 * q2.c,
            .d = _w1 * q1.d + _w2 * q2.d
        };
        return norm_quaternion(_avg);
    }

    float z = sqrtf(powf(w1-w2, 2.f) + 4.f*w1*w2*powf(dot, 2.f));
    float m = z*(w1+w2+z);
    float a = sqrtf((w1*(w1-w2+z))/m);
    float b = sqrtf((w2*(w2-w1+z))/m);
    float sgn = copysignf(1.f, dot);

    quaternion_t avg = {
        .a = a * q1.a + sgn * b * q2.a,
        .b = a * q1.b + sgn * b * q2.b,
        .c = a * q1.c + sgn * b * q2.c,
        .d = a * q1.d + sgn * b * q2.d
    };

    return norm_quaternion(avg);
}

static inline quaternion_t conjugate_quaternion(quaternion_t const q) {
    return (quaternion_t){q.a, -q.b, -q.c, -q.d};
}

static inline quaternion_t inverse_quaternion(quaternion_t const q) {
    float s = (q.a*q.a + q.b*q.b + q.c*q.c + q.d*q.d);
    quaternion_t _q = conjugate_quaternion(q);

    return (quaternion_t){_q.a / s, _q.b / s, _q.c / s, _q.d / s};
}

static inline float3_t transform_quaternion(quaternion_t const q, float3_t v) {
    return mul_quaternion(
        mul_quaternion(q, pure_quaternion(v)),
        conjugate_quaternion(q)
    ).vector;
}

static inline quaternion_t pure_quaternion(float3_t const v) {
    return (quaternion_t) {0.f, v};
}

static inline quaternion_t axis_angle_to_quaternion(float3_t const axis, float angle) {
    float const half_angle = 0.5f*angle;
    return (quaternion_t) {
        .scalar = cos(half_angle),
        .vector = scale_float3(axis, sin(half_angle)),
    };
}

static inline axis_angle_t quaternion_to_axis_angle(quaternion_t const q) {
    float const dl = length_float3(q.vector);
    if(0 == dl) {
        return (axis_angle_t) {
            .angle = 0.0f,
            .axis = (float3_t) { .x = 1.0f, .y = 0.0f, .z = 0.0f },
        };
    } else {
        return (axis_angle_t) {
            .angle = 2.0f*atan2(dl, q.scalar),
            .axis = scale_float3(q.vector, 1.0f/dl),
        };
    }
}

static inline float sign_of_variable(float x) {
    return (x >= 0.0f) ? +1.0f : -1.0f;
}

static inline float44_t quaternion_to_rotation_matrix(quaternion_t const q) {
    axis_angle_t const aa = quaternion_to_axis_angle(q);
    return rotation_float44(aa.axis, aa.angle);
}

static inline quaternion_t rotation_matrix_to_quaternion(float44_t const bb, bool *pointer){
    quaternion_t ww = {1,1,1,1};
    float44_t inv=transpose_float44(bb);
    if(equal_float44(inv,identity_sc) && is_equal(det_float44(bb),1)){
        float k0 = ( bb.m00 + bb.m11 + bb.m22 + 1.0f) / 4.0f;
        float k1 = ( bb.m00 - bb.m11 - bb.m22 + 1.0f) / 4.0f;
        float k2 = (-bb.m00 + bb.m11 - bb.m22 + 1.0f) / 4.0f;
        float k3 = (-bb.m00 - bb.m11 + bb.m22 + 1.0f) / 4.0f;

        if(k0 < 0.0f) {k0 = 0.0f;}
        if(k1 < 0.0f) {k1 = 0.0f;}
        if(k2 < 0.0f) {k2 = 0.0f;}
        if(k3 < 0.0f) {k3 = 0.0f;}

        k0 = sqrtf(k0);
        k1 = sqrtf(k1);
        k2 = sqrtf(k2);
        k3 = sqrtf(k3);

        if(k0 >= k1 && k0 >= k2 && k0 >= k3) {
            k0 *= +1.0f;
            k1 *= sign_of_variable(bb.m32 - bb.m23);
            k2 *= sign_of_variable(bb.m13 - bb.m31);
            k3 *= sign_of_variable(bb.m21 - bb.m12);
        } else if(k1 >= k0 && k1 >= k2 && k1 >= k3) {
            k0 *= sign_of_variable(bb.m32 - bb.m23);
            k1 *= +1.0f;
            k2 *= sign_of_variable(bb.m21 + bb.m12);
            k3 *= sign_of_variable(bb.m13 + bb.m31);
        } else if(k2 >= k0 && k2 >= k1 && k2 >= k3) {
            k0 *= sign_of_variable(bb.m13 - bb.m31);
            k1 *= sign_of_variable(bb.m21 + bb.m12);
            k2 *= +1.0f;
            k3 *= sign_of_variable(bb.m32 + bb.m23);
        } else if(k3 >= k0 && k3 >= k1 && k3 >= k2) {
            k0 *= sign_of_variable(bb.m21 - bb.m12);
            k1 *= sign_of_variable(bb.m31 + bb.m13);
            k2 *= sign_of_variable(bb.m32 + bb.m23);
            k3 *= +1.0f;
        } else {
            if (pointer) {
                *pointer = false;
            }
            return ww;
        }
        ww.a = k0;
        ww.b = k1;
        ww.c = k2;
        ww.d = k3;
        norm_quaternion(ww);
        if (pointer) {
            *pointer = true;
        }
    }
    return ww;
}

static inline quaternion_t quaternion_from_euler(float const angle_x, float const angle_y, float const angle_z) {
	float const cos_x = cos(0.5f*angle_x);
	float const sin_x = sin(0.5f*angle_x);
	float const cos_y = cos(0.5f*angle_y);
	float const sin_y = sin(0.5f*angle_y);
	float const cos_z = cos(0.5f*angle_z);
	float const sin_z = sin(0.5f*angle_z);

	return (quaternion_t) {
		.a = cos_x*cos_y*cos_z + sin_x*sin_y*sin_z,
		.b = sin_x*cos_y*cos_z + cos_x*sin_y*sin_z,
		.c = cos_x*sin_y*cos_z + sin_x*cos_y*sin_z,
		.d = cos_x*cos_y*sin_z + sin_x*sin_y*cos_z,
	};
}
static inline void quaternion_to_euler(quaternion_t q,float *angle_x, float *angle_y, float *angle_z) {
	if (angle_x) {
        *angle_x = atan2(2.f*(q.a*q.b + q.c*q.d), 1.f - 2.f *(q.b*q.b+q.c*q.c));
    }

    if (angle_y) {
        *angle_y = asin(2.f*(q.a*q.c-q.d*q.b));
    }

    if (angle_z) {
        *angle_z = atan2(2.f*(q.a*q.d + q.b*q.c), 1.f - 2.f *(q.c*q.c+q.d*q.d));
    }
}
static inline void rotate_aabb(float3_t *aabb_min, float3_t *aabb_max, quaternion_t rotat_aabb) {
		if (aabb_min && aabb_max) {
			float3_t vertices[8] = {
				{ aabb_min->x, aabb_min->y, aabb_min->z },
				{ aabb_min->x, aabb_min->y, aabb_max->z },

				{ aabb_min->x, aabb_max->y, aabb_min->z },
				{ aabb_min->x, aabb_max->y, aabb_max->z },

				{ aabb_max->x, aabb_min->y, aabb_min->z },
				{ aabb_max->x, aabb_min->y, aabb_max->z },

				{ aabb_max->x, aabb_max->y, aabb_min->z },
				{ aabb_max->x, aabb_max->y, aabb_max->z }
			};
			for (int i = 0; i < 8; i++) {
				vertices[i] = transform_quaternion(rotat_aabb, vertices[i]);
			}
            float3_t AABB_min = { FLT_MAX, FLT_MAX, FLT_MAX };
            float3_t AABB_max = { -FLT_MAX, -FLT_MAX, -FLT_MAX };
			for (int i = 1; i < 8; i++) {
				AABB_min = min_float3( AABB_min, vertices[i] );
				AABB_max = max_float3( AABB_max, vertices[i] );
			}
			*aabb_min = AABB_min;
			*aabb_max = AABB_max;
		}

}

static inline quaternion_t quaternion_swizzle(quaternion_t q, int a, int b, int c) {
    return (quaternion_t) {
        .scalar = q.scalar,
        .vector = swizzle_float3(q.vector, a, b, c)
    };
}

static inline quaternion_t rotate_vector_onto(float3_t const from, float3_t const to) {
    float3_t const a = normal_float3(from);
    float3_t const b = normal_float3(to);

    float3_t axis = cross_float3(a, b);
    float const sine = length_float3(axis);
    float const cosine = dot_float3(a, b);

    if (is_zero(sine)) {
        if (cosine > 0.0f) {
            return (quaternion_t) { 1.0f, 0.0f, 0.0f, 0.0f };
        } else {
            float3_t not_parallel = a;
            if (!is_zero(a.x)) {
                not_parallel.y += 1.0f;
            } else if (!is_zero(a.y)) {
                not_parallel.x += 1.0f;
            } else if (!is_zero(a.z)) {
                not_parallel.x += 1.0f;
            } else {
                // degenerate case, returning whatever.
                return (quaternion_t) { 1.0f, 0.0, 0.0f, 0.0f };
            }

            axis = cross_float3(not_parallel, a);

            return axis_angle_to_quaternion(normal_float3(axis), M_PI);
        }
    } else {
        return axis_angle_to_quaternion(scale_float3(axis, 1.0f/sine), atan2(sine, cosine));
    }
}

#endif // __COMMON_QUATERNION_H
