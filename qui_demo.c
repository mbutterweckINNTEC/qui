#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <unistd.h>

#include "mathutils.h"
#include "quaternion.h"

#include "qui_def.h"
#include "qui_in.h"
#include "qui_shdr.h"
#include "qui_util.h"

char *cube_vsh =
	"#version 440"				"\n"
	"in      vec4 v;"			"\n"
	"in      vec4 c;"			"\n"
	"uniform mat4 M;"			"\n"
	"out     vec4 k;"			"\n"
	"out     vec4 p;"			"\n"
	""					"\n"
	"void main() {"				"\n"
	"	p = M * v;"			"\n"
	"	gl_Position = p;"		"\n"
	"	k = c;"				"\n"
	"}"					"\n";

char *cube_fsh =
	"#version 440"					"\n"
	"in  vec4 k;"					"\n"
	"in  vec4 p;"					"\n"
	"out vec4 K;"					"\n"
	""						"\n"
	"void main() {"					"\n"
	"	vec3 s = dFdx(p.xyz);"			"\n"
	"	vec3 t = dFdy(p.xyz);"			"\n"
	"	vec3 n = normalize(cross(s, t));"	"\n"
	"	vec3 L = normalize(vec3(0.5,0.5,1));"	"\n"
	"	K = k * dot(n, L);"			"\n"
	"}"						"\n";


float cube_v[] = {
	 0.5f,  0.375f, -0.25f, 1.f,
	 0.5f, -0.375f, -0.25f, 1.f,
	 0.5f,  0.375f,  0.25f, 1.f,
	 0.5f, -0.375f,  0.25f, 1.f,
	-0.5f,  0.375f, -0.25f, 1.f,
	-0.5f, -0.375f, -0.25f, 1.f,
	-0.5f,  0.375f,  0.25f, 1.f,
	-0.5f, -0.375f,  0.25f, 1.f
};

int cube_i[] = {
	4, 2, 0,
	2, 7, 3,
	6, 5, 7,
	1, 7, 5,
	0, 3, 1,
	4, 1, 5,
	4, 6, 2,
	2, 6, 7,
	6, 4, 5,
	1, 3, 7,
	0, 2, 3,
	4, 0, 1
};

#define QUI_MAN_DRW_MX 12

enum {
	QUI_MAN_STTS_NIL,
	QUI_MAN_STTS_MOV_X,
	QUI_MAN_STTS_MOV_Y,
	QUI_MAN_STTS_MOV_Z,
	QUI_MAN_STTS_MOV_V,
	QUI_MAN_STTS_ROT_X,
	QUI_MAN_STTS_ROT_Y,
	QUI_MAN_STTS_ROT_Z,
	QUI_MAN_STTS_ROT_V,
	QUI_MAN_STTS_SCL		/* only uniform scalling */
};

struct qui_man {
	int bo, vao;

	int stts;
	int phi;
	float dlt;

	float3_t t;
	quaternion_t q;
	float s;
};

int qui_man_mk(struct qui_man *qm) {
	if (NULL == qm)
		return -1;

	memset(qm, 0, sizeof(struct qui_man));

	/* visual objects */
	glCreateBuffers(1, &qm->bo);

	if (!qm->bo)
		goto ouch;

	glNamedBufferStorage(qm->bo, QUI_MAN_SZ, NULL, GL_DYNAMIC_STORAGE_BIT | GL_MAP_WRITE_BIT);

	float3_t *v = glMapNamedBufferRange(qm->bo, 0, QUI_MAN_SZ, GL_MAP_WRITE_BIT);

	if (!v)
		goto ouch;

	/* x axis */
	*v++ = (float3_t) { QUI_MAN_L_XYZ, 0.f, 0.f };
	*v++ = QUI_MAN_X_CLR;
	*v++ = (float3_t) { 0.0f, 0.f, 0.f };
	*v++ = QUI_MAN_X_CLR;

	for (int i = 0; i < QUI_MAN_CRCL_X_N; ++i) {
		float phi = (float)i / (float)(QUI_MAN_CRCL_X_N - 1) * 2.0 * M_PI;
		*v++ = (float3_t) { 0.f, QUI_MAN_R_XYZ * cos(phi), QUI_MAN_R_XYZ * sin(phi) };
		*v++ = QUI_MAN_X_CLR;		
	}

	/* y axis */
	*v++ = (float3_t) { 0.0f, QUI_MAN_L_XYZ, 0.f };
	*v++ = QUI_MAN_Y_CLR;
	*v++ = (float3_t) { 0.0f, 0.f, 0.f };
	*v++ = QUI_MAN_Y_CLR;

	for (int i = 0; i < QUI_MAN_CRCL_Y_N; ++i) {
		float phi = (float)i / (float)(QUI_MAN_CRCL_Y_N - 1) * 2.0 * M_PI;
		*v++ = (float3_t) { QUI_MAN_R_XYZ * cos(phi), 0.f, QUI_MAN_R_XYZ * sin(phi) };
		*v++ = QUI_MAN_Y_CLR;		
	}

	/* z axis */
	*v++ = (float3_t) { 0.0f, 0.f, QUI_MAN_L_XYZ };
	*v++ = QUI_MAN_Z_CLR;
	*v++ = (float3_t) { 0.0f, 0.f, 0.f };
	*v++ = QUI_MAN_Z_CLR;

	for (int i = 0; i < QUI_MAN_CRCL_Z_N; ++i) {
		float phi = (float)i / (float)(QUI_MAN_CRCL_Z_N - 1) * 2.0 * M_PI;
		*v++ = (float3_t) { QUI_MAN_R_XYZ * cos(phi), QUI_MAN_R_XYZ * sin(phi), 0.f };
		*v++ = QUI_MAN_Z_CLR;		
	}

	/* view axis */
	*v++ = (float3_t) { 0.0f, 0.f, QUI_MAN_L_XYZ };
	*v++ = QUI_MAN_V_CLR;
	*v++ = (float3_t) { 0.0f, 0.f, 0.f };
	*v++ = QUI_MAN_V_CLR;

	for (int i = 0; i < QUI_MAN_CRCL_V_N; ++i) {
		float phi = (float)i / (float)(QUI_MAN_CRCL_V_N - 1) * 2.0 * M_PI;
		*v++ = (float3_t) { QUI_MAN_R_V * cos(phi), QUI_MAN_R_V * sin(phi), 0.f };
		*v++ = QUI_MAN_V_CLR;		
	}

	/* resize frame */
	*v++ = (float3_t) { QUI_MAN_S_XY, QUI_MAN_S_DXY, 0.f };
	*v++ = QUI_MAN_S_CLR;
	*v++ = (float3_t) { QUI_MAN_S_DXY, QUI_MAN_S_DXY, 0.f };
	*v++ = QUI_MAN_S_CLR;
	*v++ = (float3_t) { QUI_MAN_S_DXY, QUI_MAN_S_XY, 0.f };
	*v++ = QUI_MAN_S_CLR;

	*v++ = (float3_t) { -QUI_MAN_S_XY, QUI_MAN_S_DXY, 0.f };
	*v++ = QUI_MAN_S_CLR;
	*v++ = (float3_t) { -QUI_MAN_S_DXY, QUI_MAN_S_DXY, 0.f };
	*v++ = QUI_MAN_S_CLR;
	*v++ = (float3_t) { -QUI_MAN_S_DXY, QUI_MAN_S_XY, 0.f };
	*v++ = QUI_MAN_S_CLR;

	*v++ = (float3_t) { QUI_MAN_S_XY, -QUI_MAN_S_DXY, 0.f };
	*v++ = QUI_MAN_S_CLR;
	*v++ = (float3_t) { QUI_MAN_S_DXY, -QUI_MAN_S_DXY, 0.f };
	*v++ = QUI_MAN_S_CLR;
	*v++ = (float3_t) { QUI_MAN_S_DXY, -QUI_MAN_S_XY, 0.f };
	*v++ = QUI_MAN_S_CLR;

	*v++ = (float3_t) { -QUI_MAN_S_XY, -QUI_MAN_S_DXY, 0.f };
	*v++ = QUI_MAN_S_CLR;
	*v++ = (float3_t) { -QUI_MAN_S_DXY, -QUI_MAN_S_DXY, 0.f };
	*v++ = QUI_MAN_S_CLR;
	*v++ = (float3_t) { -QUI_MAN_S_DXY, -QUI_MAN_S_XY, 0.f };
	*v++ = QUI_MAN_S_CLR;


	glUnmapNamedBuffer(qm->bo);

	glGenVertexArrays(1, &qm->vao);

	if (!qm->vao)
		goto ouch;

	glBindVertexArray(qm->vao);
	glBindBuffer(GL_ARRAY_BUFFER, qm->bo);
	glVertexAttribPointer(0, 3, GL_FLOAT, 0, 24, 0);
	glEnableVertexAttribArray(0); 
	glVertexAttribPointer(1, 3, GL_FLOAT, 0, 24, (char*)12);
	glEnableVertexAttribArray(1); 
	glBindVertexArray(0);

	return 0;
ouch:
	if (qm->bo)
		glDeleteBuffers(1, &qm->bo);

	if (qm->vao)
		glDeleteVertexArrays(1, &qm->vao);

	memset(qm, 0, sizeof(struct qui_man));

	return -1;
}

enum {
	QUI_MAN_FLGS_AXIS_X = 0x1,
	QUI_MAN_FLGS_AXIS_Y = 0x2,
	QUI_MAN_FLGS_AXIS_Z = 0x4,
	QUI_MAN_FLGS_AXIS_V = 0x8,
	QUI_MAN_FLGS_AXIS_XYZV = 0xf,

	QUI_MAN_FLGS_CRCL_X = 0x10,
	QUI_MAN_FLGS_CRCL_Y = 0x20,
	QUI_MAN_FLGS_CRCL_Z = 0x40,
	QUI_MAN_FLGS_CRCL_V = 0x80,
	QUI_MAN_FLGS_CRCL_XYZV = 0xf0,
	QUI_MAN_FLGS_AXIS_CRCL = 0xff,

	QUI_MAN_FLGS_PIE_X = 0x100,
	QUI_MAN_FLGS_PIE_Y = 0x200,
	QUI_MAN_FLGS_PIE_Z = 0x400,
	QUI_MAN_FLGS_PIE_V = 0x800,
	QUI_MAN_FLGS_PIE_ALL = 0xf00,

	QUI_MAN_FLGS_FRM = 0xf000,	/* frame consist of 4 separate draw calls, that's why it spans 4 bits */
	QUI_MAN_FLGS_AXIS_CRCL_FRM = 0xf0ff,

	QUI_MAN_FLGS_BITS_N = 16
};

enum {
	QUI_MAN_OP_END,
	QUI_MAN_OP_PIE_STRT,
	QUI_MAN_OP_PIE_ANGL,
	QUI_MAN_OP_LN_WDTH,
	QUI_MAN_OP_SCL
};

#define QUI_MAN_DRW_PASS_N 4

int qui_man_drw(struct qui_man *qm, struct qui_shdr *qs, float44_t P, float44_t V, float44_t W, int op[], int flgs) {
	static int const sv_def[] = {
		/* x axis */ 0,
		/* y axis */ QUI_MAN_PLN_X_N,
		/* z axis */ QUI_MAN_PLN_X_N + QUI_MAN_PLN_Y_N,
		/* v axis */ QUI_MAN_PLN_X_N + QUI_MAN_PLN_Y_N + QUI_MAN_PLN_Z_N,
		/* x crcl */ QUI_MAN_AXIS_X_N,
		/* y crcl */ QUI_MAN_PLN_X_N + QUI_MAN_AXIS_Y_N,
		/* z crcl */ QUI_MAN_PLN_X_N + QUI_MAN_PLN_Y_N + QUI_MAN_AXIS_Z_N,
		/* v crcl */ QUI_MAN_PLN_X_N + QUI_MAN_PLN_Y_N + QUI_MAN_PLN_Z_N + QUI_MAN_AXIS_V_N,
		/* x disc */ QUI_MAN_AXIS_X_N - 1,
		/* y disc */ QUI_MAN_PLN_X_N + QUI_MAN_AXIS_Y_N - 1,
		/* z disc */ QUI_MAN_PLN_X_N + QUI_MAN_PLN_Y_N + QUI_MAN_AXIS_Z_N - 1,
		/* v disc */ QUI_MAN_PLN_X_N + QUI_MAN_PLN_Y_N + QUI_MAN_PLN_Z_N + QUI_MAN_AXIS_V_N - 1,
		/* resize */ QUI_MAN_PLN_X_N + QUI_MAN_PLN_Y_N + QUI_MAN_PLN_Z_N + QUI_MAN_PLN_V_N,
		/* resize */ QUI_MAN_PLN_X_N + QUI_MAN_PLN_Y_N + QUI_MAN_PLN_Z_N + QUI_MAN_PLN_V_N + 3,
		/* resize */ QUI_MAN_PLN_X_N + QUI_MAN_PLN_Y_N + QUI_MAN_PLN_Z_N + QUI_MAN_PLN_V_N + 6,
		/* resize */ QUI_MAN_PLN_X_N + QUI_MAN_PLN_Y_N + QUI_MAN_PLN_Z_N + QUI_MAN_PLN_V_N + 9
	};
	static int const nv_def[] = {
		/* x axis */ QUI_MAN_AXIS_X_N,
		/* y axis */ QUI_MAN_AXIS_Y_N,
		/* z axis */ QUI_MAN_AXIS_Z_N,
		/* v axis */ QUI_MAN_AXIS_V_N,
		/* x crcl */ QUI_MAN_CRCL_X_N,
		/* y crcl */ QUI_MAN_CRCL_Y_N,
		/* z crcl */ QUI_MAN_CRCL_Z_N,
		/* v crcl */ QUI_MAN_CRCL_V_N,
		/* x disc */ 2,
		/* y disc */ 2,
		/* z disc */ 2,
		/* v disc */ 2,
		/* resize */ 3,
		/* resize */ 3,
		/* resize */ 3,
		/* resize */ 3,
	};
	static int const pass_def[QUI_MAN_DRW_PASS_N] = { 0x0077, 0xf080, 0x0700, 0x0800 };

	int sv[QUI_MAN_DRW_PASS_N][QUI_MAN_DRW_MX];
	int nv[QUI_MAN_DRW_PASS_N][QUI_MAN_DRW_MX];
	int iv[QUI_MAN_DRW_PASS_N] = { 0 };
	int pass = 0, dsc_strt = 0, dsc_angl = 360, lnwdth = 0;

	float44_t M = mul_float44(V, P);

	if (!qm || !qs)
		return -1;

	if (op) {
		for (int i = 0; op[i]; ++i) {
			switch (op[i++]) {
			case QUI_MAN_OP_PIE_STRT: dsc_strt = op[i]; break;
			case QUI_MAN_OP_PIE_ANGL: dsc_angl = op[i]; break;
			case QUI_MAN_OP_LN_WDTH: lnwdth = op[i]; break;
			};
		}
	}

	W = mul_float44(W, P);

	for (int i = 0; i < QUI_MAN_FLGS_BITS_N; ++i) {
		for (int j = 0; j < QUI_MAN_DRW_PASS_N; ++j) {
			if (flgs & 1 << i && pass_def[j] & 1 << i) {
				sv[j][iv[j]] = sv_def[i];
				nv[j][iv[j]] = nv_def[i];
				if (QUI_MAN_FLGS_PIE_ALL & 1 << i) {
					nv[j][iv[j]] += dsc_angl;
				}
				iv[j]++;
				pass |= 1 << j;
			}
		}
	}

	if (0 == pass)
		return 0;

	glEnable(GL_BLEND);
	glBlendEquation(GL_FUNC_ADD);
	glBlendFunc(GL_CONSTANT_COLOR, GL_ONE_MINUS_CONSTANT_COLOR );
	glBlendColor(0.75, 0.75, 0.75, 0.75);

	glUseProgram(qs->po);
	glBindVertexArray(qm->vao);

	if (pass & 3 && lnwdth)
		glLineWidth(lnwdth);

	if (pass & 1) {
		glUniformMatrix4fv(qs->M, 1, 0, &M.m[0][0]);
		glMultiDrawArrays(GL_LINE_STRIP, sv[0], nv[0], iv[0]);
	}

	if (pass & 2) {
		glUniformMatrix4fv(qs->M, 1, 0, &W.m[0][0]);
		glMultiDrawArrays(GL_LINE_STRIP, sv[1], nv[1], iv[1]);
	}

	if (pass & 3 && lnwdth)
		glLineWidth(1);

	if (pass & 4) {
		if (dsc_strt) {
			quaternion_t q;
			if (flgs & QUI_MAN_FLGS_PIE_X)
				q = axis_angle_to_quaternion((float3_t) {1.f, 0.f, 0.f }, -dsc_strt * M_PI / 180.0);
			if (flgs & QUI_MAN_FLGS_PIE_Y)
				q = axis_angle_to_quaternion((float3_t) { 0.f, 1.f, 0.f }, dsc_strt * M_PI / 180.0);
			if (flgs & QUI_MAN_FLGS_PIE_Z)
				q = axis_angle_to_quaternion((float3_t) { 0.f, 0.f, 1.f }, -dsc_strt * M_PI / 180.0);

			M = mul_float44(mul_float44(quaternion_to_rotation_matrix(q), V), P);
		}
		glUniformMatrix4fv(qs->M, 1, 0, &M.m[0][0]);
		glMultiDrawArrays(GL_TRIANGLE_FAN, sv[2], nv[2], iv[2]);
	}

	if (pass & 8) {
		if (dsc_strt) {
			quaternion_t q = axis_angle_to_quaternion((float3_t) { 0.f, 0.f, 1.f }, -dsc_strt * M_PI / 180.0);
			M = mul_float44(quaternion_to_rotation_matrix(q), W);
		} else {
			M = W;
		}
		glUniformMatrix4fv(qs->M, 1, 0, &M.m[0][0]);
		glMultiDrawArrays(GL_TRIANGLE_FAN, sv[3], nv[3], iv[3]);
	}

	glDisable(GL_BLEND);

	return 0;
}

#define QUI_MAN_EPS 0.01

enum {
	QUI_MAN_NIL,
	QUI_MAN_ACT,
	QUI_MAN_SET
};

int qui_man(struct qui_man *qm, struct qui_shdr *qs, struct qui_in *qi, float44_t P, float44_t V, float3_t *mt, quaternion_t *mq, float *ms) {
	int rstts = qm->stts ? QUI_MAN_ACT :QUI_MAN_NIL;

	float3_t mt_ = *mt;

	if (qm->stts == QUI_MAN_STTS_MOV_X || qm->stts == QUI_MAN_STTS_MOV_Y || qm->stts == QUI_MAN_STTS_MOV_Z) {
		mt_ = qm->t;
	}
	float44_t T =  {
		*ms, 0.f, 0.f, 0.f,
		0.f, *ms, 0.f, 0.f,
		0.f, 0.f, *ms, 0.f,
		mt_.x, mt_.y, mt_.z, 1.f
	};

	V = mul_float44(T, V);

	float44_t PV = mul_float44(V, P);
	float detPV = det_float44(PV);
	float44_t iPV = invert_float44(PV, detPV);
	float3_t p = float3_float4(cotransform_float44(iPV, (float4_t){ qi->p.x, qi->p.y, 0.f, 1.f }));
	float3_t d = normal_float3(m_float3(cotransform_float44(iPV, (float4_t){ 0.f, 0.f, 1.f, 0.f })));
	float fV = frobenius_float33(float33_float44(V)) / sqrt(3);
	float44_t W = {
		fV, 0.f, 0.f, 0.f,
		0.f, fV, 0.f, 0.f,
		0.f, 0.f, fV, 0.f,
		V.m30, V.m31, V.m32, 1.f
	};
	float44_t PW = (mul_float44(W, P));
	float detPW = det_float44(PW);
	float44_t iPW = invert_float44(PW, detPW);
	float2_t pv = m_float2(float3_float4(cotransform_float44(iPW, (float4_t){ qi->p.x, qi->p.y, 0.f, 1.f })));

	int op[16];
	float ds = 1.f;
	float dlt = qm->dlt;
	int hvrd = 0;
	int dphi = 0;

	/* input handling */
	if (QUI_MAN_STTS_NIL == qm->stts) {
		float l, nl;
		int stts = 0, phi = 0, nphi = 0;

		l = qui_ray_xcrcl_(p, d, &phi);
		stts = QUI_MAN_STTS_ROT_X;

		nl = qui_ray_ycrcl_(p, d, &nphi);
		if (nl < l) {
			l = nl;
			phi = nphi;
			stts = QUI_MAN_STTS_ROT_Y;
		}

		nl = qui_ray_zcrcl_(p, d, &nphi);
		if (nl < l) {
			l = nl;
			phi = nphi;
			stts = QUI_MAN_STTS_ROT_Z;
		}

		nl = qui_ray_vcrcl_(pv, &nphi);
		if (nl < l) {
			l = nl;
			phi = nphi;
			stts = QUI_MAN_STTS_ROT_V;
		}

		nl = qui_ray_crnr_(pv);
		if (nl < l) {
			l = nl;
			ds = length_float2(pv) / sqrtf(2.f) / QUI_MAN_S_DXY;
			stts = QUI_MAN_STTS_SCL;
		}

		nl = qui_ray_seg_dst(p, d, (float3_t){ 0.f, 0.f, 0.f }, (float3_t){ QUI_MAN_L_XYZ, 0.f, 0.f });
		if (nl < l) {
			l = nl;
			qm->dlt = qui_ray_ln_near(p, d, (float3_t){ 0.f, 0.f, 0.f }, (float3_t){ QUI_MAN_L_XYZ, 0.f, 0.f }).x;
			stts = QUI_MAN_STTS_MOV_X;
		}

		nl = qui_ray_seg_dst(p, d, (float3_t){ 0.f, 0.f, 0.f }, (float3_t){ 0.f, QUI_MAN_L_XYZ, 0.f });
		if (nl < l) {
			l = nl;
			qm->dlt = qui_ray_ln_near(p, d, (float3_t){ 0.f, 0.f, 0.f }, (float3_t){ 0.f, QUI_MAN_L_XYZ, 0.f }).y;
			stts = QUI_MAN_STTS_MOV_Y;
		}

		nl = qui_ray_seg_dst(p, d, (float3_t){ 0.f, 0.f, 0.f }, (float3_t){ 0.f, 0.f, QUI_MAN_L_XYZ });
		if (nl < l) {
			l = nl;
			qm->dlt = qui_ray_ln_near(p, d, (float3_t){ 0.f, 0.f, 0.f }, (float3_t){ 0.f, 0.f, QUI_MAN_L_XYZ }).z;
			stts = QUI_MAN_STTS_MOV_Z;
		}

		if (l * fV < QUI_MAN_EPS) {
			if (qi->rls & QUI_IN_LMB) {
				qm->phi = phi;
				qm->stts = stts;
				qm->q = *mq;
				qm->t = *mt;
				qm->s = *ms;
				rstts = QUI_MAN_ACT;
			} else {
				hvrd = stts;
				ds = 1.f;
			}
		} else {
			ds = 1.f;
		}
	} else {
		int phi = 0;
		switch(qm->stts) {
		case QUI_MAN_STTS_ROT_X:
			qui_ray_xcrcl_(p, d, &phi);
			dphi = phi - qm->phi;
			if (dphi < 0) dphi += 360;
			break;
		case QUI_MAN_STTS_ROT_Y:
			qui_ray_ycrcl_(p, d, &phi);
			dphi = phi - qm->phi;
			if (dphi < 0) dphi += 360;
			break;
		case QUI_MAN_STTS_ROT_Z:
			qui_ray_zcrcl_(p, d, &phi);
			dphi = phi - qm->phi;
			if (dphi < 0) dphi += 360;
			break;
		case QUI_MAN_STTS_ROT_V:
			qui_ray_vcrcl_(pv, &phi);
			dphi = phi - qm->phi;
			if (dphi < 0) dphi += 360;
			break;
		case QUI_MAN_STTS_SCL:
			ds = length_float2(pv) / sqrtf(2.f) / QUI_MAN_S_DXY;
			break;
		case QUI_MAN_STTS_MOV_X:
			dlt = qui_ray_ln_near(p, d, (float3_t){ 0.f, 0.f, 0.f }, (float3_t){ QUI_MAN_L_XYZ, 0.f, 0.f }).x;
			break;
		case QUI_MAN_STTS_MOV_Y:
			dlt = qui_ray_ln_near(p, d, (float3_t){ 0.f, 0.f, 0.f }, (float3_t){ 0.f, QUI_MAN_L_XYZ, 0.f }).y;
			break;
		case QUI_MAN_STTS_MOV_Z:
			dlt = qui_ray_ln_near(p, d, (float3_t){ 0.f, 0.f, 0.f }, (float3_t){ 0.f, 0.f, QUI_MAN_L_XYZ }).z;
			break;
		};

		if (qi->rls & QUI_IN_LMB) {
			rstts = QUI_MAN_SET;
		}

		if (qi->rls & QUI_IN_ESC) {
			switch(qm->stts) {
			case QUI_MAN_STTS_ROT_X:
			case QUI_MAN_STTS_ROT_Y:
			case QUI_MAN_STTS_ROT_Z:
			case QUI_MAN_STTS_ROT_V:
				*mq = qm->q;
				break;
			case QUI_MAN_STTS_SCL:
				*ms = qm->s;
				break;
			case QUI_MAN_STTS_MOV_X:
			case QUI_MAN_STTS_MOV_Y:
			case QUI_MAN_STTS_MOV_Z:
				*mt = qm->t;
				break;
			};

			qm->stts = QUI_MAN_STTS_NIL;
			rstts = QUI_MAN_NIL;
		}
	}

	W = (float44_t) {
		fV * ds, 0.f, 0.f, 0.f,
		0.f, fV * ds, 0.f, 0.f,
		0.f, 0.f, fV * ds, 0.f,
		V.m30, V.m31, V.m32, 1.f
	};

	switch(qm->stts) {
	case QUI_MAN_STTS_MOV_X:
		qui_man_drw(qm, qs, P, V, W, (int[]){QUI_MAN_OP_LN_WDTH, 2, QUI_MAN_OP_END}, QUI_MAN_FLGS_AXIS_X);
		break;
	case QUI_MAN_STTS_MOV_Y:
		qui_man_drw(qm, qs, P, V, W, (int[]){QUI_MAN_OP_LN_WDTH, 2, QUI_MAN_OP_END}, QUI_MAN_FLGS_AXIS_Y);
		break;
	case QUI_MAN_STTS_MOV_Z:
		qui_man_drw(qm, qs, P, V, W, (int[]){QUI_MAN_OP_LN_WDTH, 2, QUI_MAN_OP_END}, QUI_MAN_FLGS_AXIS_Z);
		break;
	case QUI_MAN_STTS_NIL:
		qui_man_drw(qm, qs, P, V, W, (int[]){QUI_MAN_OP_LN_WDTH, 2, QUI_MAN_OP_END}, QUI_MAN_FLGS_AXIS_CRCL_FRM);
		break;
	case QUI_MAN_STTS_ROT_X:
		qui_man_drw(qm, qs, P, V, W, (int[]){QUI_MAN_OP_PIE_STRT, qm->phi, QUI_MAN_OP_PIE_ANGL, dphi, QUI_MAN_OP_LN_WDTH, 2, QUI_MAN_OP_END}, QUI_MAN_FLGS_CRCL_X | QUI_MAN_FLGS_PIE_X);
		break;
	case QUI_MAN_STTS_ROT_Y:
		qui_man_drw(qm, qs, P, V, W, (int[]){QUI_MAN_OP_PIE_STRT, qm->phi, QUI_MAN_OP_PIE_ANGL, dphi, QUI_MAN_OP_LN_WDTH, 2, QUI_MAN_OP_END}, QUI_MAN_FLGS_CRCL_Y | QUI_MAN_FLGS_PIE_Y);
		break;
	case QUI_MAN_STTS_ROT_Z:
		qui_man_drw(qm, qs, P, V, W, (int[]){QUI_MAN_OP_PIE_STRT, qm->phi, QUI_MAN_OP_PIE_ANGL, dphi, QUI_MAN_OP_LN_WDTH, 2, QUI_MAN_OP_END}, QUI_MAN_FLGS_CRCL_Z | QUI_MAN_FLGS_PIE_Z);
		break;
	case QUI_MAN_STTS_ROT_V:
		qui_man_drw(qm, qs, P, V, W, (int[]){QUI_MAN_OP_PIE_STRT, qm->phi, QUI_MAN_OP_PIE_ANGL, dphi, QUI_MAN_OP_LN_WDTH, 2, QUI_MAN_OP_END}, QUI_MAN_FLGS_CRCL_V | QUI_MAN_FLGS_PIE_V);
		break;
		case QUI_MAN_STTS_SCL:
		qui_man_drw(qm, qs, P, V, W, (int[]){QUI_MAN_OP_LN_WDTH, 2, QUI_MAN_OP_END}, QUI_MAN_FLGS_FRM);
		break;
	};

	switch(hvrd) {
	case QUI_MAN_STTS_NIL:
		break;
	case QUI_MAN_STTS_MOV_X:
		qui_man_drw(qm, qs, P, V, W, (int[]){QUI_MAN_OP_LN_WDTH, 4, QUI_MAN_OP_END}, QUI_MAN_FLGS_AXIS_X);
		break;
	case QUI_MAN_STTS_MOV_Y:
		qui_man_drw(qm, qs, P, V, W, (int[]){QUI_MAN_OP_LN_WDTH, 4, QUI_MAN_OP_END}, QUI_MAN_FLGS_AXIS_Y);
		break;
	case QUI_MAN_STTS_MOV_Z:
		qui_man_drw(qm, qs, P, V, W, (int[]){QUI_MAN_OP_LN_WDTH, 4, QUI_MAN_OP_END}, QUI_MAN_FLGS_AXIS_Z);
		break;
	case QUI_MAN_STTS_ROT_X:
		qui_man_drw(qm, qs, P, V, W, (int[]){QUI_MAN_OP_LN_WDTH, 4, QUI_MAN_OP_END}, QUI_MAN_FLGS_CRCL_X);
		break;
	case QUI_MAN_STTS_ROT_Y:
		qui_man_drw(qm, qs, P, V, W, (int[]){QUI_MAN_OP_LN_WDTH, 4, QUI_MAN_OP_END}, QUI_MAN_FLGS_CRCL_Y);
		break;
	case QUI_MAN_STTS_ROT_Z:
		qui_man_drw(qm, qs, P, V, W, (int[]){QUI_MAN_OP_LN_WDTH, 4, QUI_MAN_OP_END}, QUI_MAN_FLGS_CRCL_Z);
		break;
	case QUI_MAN_STTS_ROT_V:
		qui_man_drw(qm, qs, P, V, W, (int[]){QUI_MAN_OP_LN_WDTH, 4, QUI_MAN_OP_END}, QUI_MAN_FLGS_CRCL_V);
		break;
		case QUI_MAN_STTS_SCL:
		qui_man_drw(qm, qs, P, V, W, (int[]){QUI_MAN_OP_LN_WDTH, 4, QUI_MAN_OP_END}, QUI_MAN_FLGS_FRM);
		break;
	};

	switch(qm->stts) {
	case QUI_MAN_STTS_NIL:
		break;
	case QUI_MAN_STTS_ROT_X:
		*mq = mul_quaternion(qm->q, axis_angle_to_quaternion((float3_t) {1.f, 0.f, 0.f }, -dphi * M_PI / 180.0));
		break;
	case QUI_MAN_STTS_ROT_Y:
		*mq = mul_quaternion(qm->q, axis_angle_to_quaternion((float3_t) { 0.f, 1.f, 0.f }, dphi * M_PI / 180.0));
		break;
	case QUI_MAN_STTS_ROT_Z:
		*mq = mul_quaternion(qm->q, axis_angle_to_quaternion((float3_t) { 0.f, 0.f, 1.f }, -dphi * M_PI / 180.0));
		break;
	case QUI_MAN_STTS_ROT_V:
		*mq = mul_quaternion(qm->q, axis_angle_to_quaternion(d, dphi * M_PI / 180.0));
		break;
	case QUI_MAN_STTS_MOV_X:
		*mt = add_float3(qm->t, (float3_t){ dlt - qm->dlt, 0.f, 0.f });
		break;
	case QUI_MAN_STTS_MOV_Y:
		*mt = add_float3(qm->t, (float3_t){ 0.f, dlt - qm->dlt, 0.f });
		break;
	case QUI_MAN_STTS_MOV_Z:
		*mt = add_float3(qm->t, (float3_t){  0.f, 0.f, dlt - qm->dlt });
		break;
	};

	*ms *= ds;

	if (rstts == QUI_MAN_SET)
			qm->stts = QUI_MAN_STTS_NIL;

	return rstts;
}

float qui_scrll;
float qui_scrll_sgn = 1.f;

void qui_scrll_cb(GLFWwindow *wndw, double x, double y) {
	qui_scrll = y;
}


int main(int argc, char *argv[]) {
	GLFWwindow *wndw;
	int w = 1024, h = 1024, qui_po, cube_po, cube_bo[2], cube_vao, M_loc_cube, qui_bo, qui_vao, M_loc_qui;
	float bg[4] = { 170.f/255.f, 170.f/255.f, 170.f/255.f, 0};
	float44_t M = identity_sc;
	float44_t Q = identity_sc;
	float44_t V = (float44_t){
				1.f, 0.f, 0.f, 0.f,
				0.f, 1.f, 0.f, 0.f,
				0.f, 0.f, 1.f, 0.f,
				0.5, 0.3333, 0.0, 1.f
			};
	float44_t VM = identity_sc;
	float44_t P = orthographic(1.f, -1.f, 1.f, -1.f, -1.f, 1.f);
	float44_t PV = identity_sc;
	float44_t PVM = identity_sc;
	float44_t S = identity_sc;
	float44_t T = identity_sc;

	quaternion_t q = { 1.0 };
	float angl = 0.0, ar;
	double cx, cy, cl;
	struct qui_shdr qs = {};
	struct qui_in qi = {};
	struct qui_man qm = {};

	glfwInit();

	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 5);
	glfwWindowHint(GLFW_DOUBLEBUFFER, GLFW_TRUE);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	wndw = glfwCreateWindow(w, h, "Quaternion UI prototype", NULL, NULL);
//	glfwSetScrollCallback(wndw, sqrrl);
	glfwMakeContextCurrent(wndw);
	glewExperimental = 1;
	glewInit();

	glfwSetScrollCallback(wndw, qui_scrll_cb);

	cube_po = shdr_mk(cube_vsh, cube_fsh);

	qui_shdr_mk(&qs);
	qui_man_mk(&qm);

	glCreateBuffers(2, cube_bo);
	glNamedBufferStorage(cube_bo[0], sizeof(cube_v), cube_v, 0);
	glNamedBufferStorage(cube_bo[1], sizeof(cube_i), cube_i, 0);
	glGenVertexArrays(1, &cube_vao);
	glBindVertexArray(cube_vao);
	glBindBuffer(GL_ARRAY_BUFFER, cube_bo[0]);
	glVertexAttribPointer(0, 4, GL_FLOAT, 0, 0, 0);
	glEnableVertexAttribArray(0); 
	glVertexAttribPointer(0, 4, GL_FLOAT, 0, 0, 0);
	glVertexAttrib4f(1, 1,1,1,1);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cube_bo[1]);
	glBindVertexArray(0);

	glClearColor(bg[0], bg[1], bg[2], bg[3]);
	M_loc_cube = glGetUniformLocation(cube_po, "M");

	glDepthFunc(GL_LESS);

	quaternion_t mq = { 1, 0, 0, 0 };
	float3_t mt = {0, 0, 0};
	float ms = 1.f;

	while(!glfwWindowShouldClose(wndw)) {
		glfwGetWindowSize(wndw, &w, &h);
	
		ar = (float)w / (float)h;

		glViewport(0, 0, w, h);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		if (glfwGetKey(wndw, GLFW_KEY_UP) == GLFW_PRESS) {
			if (glfwGetKey(wndw, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS) {
				V = mul_float44(V, (float44_t){ 1, 0, 0, 0,    0, 1, 0, 0,   0, 0, 1.0, 0,   0, 0.125, 0, 1});
			} else {
				V = mul_float44(V, quaternion_to_rotation_matrix((quaternion_t){cos(0.03125),-sin(0.03125), 0, 0 }));
			}
		}

		if (glfwGetKey(wndw, GLFW_KEY_DOWN) == GLFW_PRESS) {
			if (glfwGetKey(wndw, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS) {
				V = mul_float44(V, (float44_t){ 1, 0, 0, 0,    0, 1, 0, 0,   0, 0, 1.0, 0,   0, -0.125, 0, 1});
			} else {
				V = mul_float44(V, quaternion_to_rotation_matrix((quaternion_t){cos(0.03125), sin(0.03125), 0, 0 }));
			}
		}

		if (glfwGetKey(wndw, GLFW_KEY_LEFT) == GLFW_PRESS) {
			if (glfwGetKey(wndw, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS) {
				V = mul_float44(V, (float44_t){ 1, 0, 0, 0,    0, 1, 0, 0,   0, 0, 1.0, 0, -0.125, 0, 0, 1});
			} else {
				V = mul_float44(V, quaternion_to_rotation_matrix((quaternion_t){cos(0.03125), 0,-sin(0.03125), 0 }));
			}
		}

		if (glfwGetKey(wndw, GLFW_KEY_RIGHT) == GLFW_PRESS) {
			if (glfwGetKey(wndw, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS) {
				V = mul_float44(V, (float44_t){ 1, 0, 0, 0,    0, 1, 0, 0,   0, 0, 1.0, 0,  0.125, 0, 0, 1});
			} else {
				V = mul_float44(V, quaternion_to_rotation_matrix((quaternion_t){cos(0.03125), 0, sin(0.03125), 0 }));
			}
		}

		if (glfwGetKey(wndw, GLFW_KEY_KP_ADD) == GLFW_PRESS) {
			V = mul_float44(V, (float44_t){ 2, 0, 0, 0,    0, 2, 0, 0,   0, 0, 2.0, 0,   0, 0, 0, 1});
		}

		if (glfwGetKey(wndw, GLFW_KEY_KP_SUBTRACT) == GLFW_PRESS) {
			V = mul_float44(V, (float44_t){ 0.5, 0, 0, 0,    0, 0.5, 0, 0,   0, 0, 0.5, 0,   0, 0, 0, 1});
		}

		if (glfwGetKey(wndw, GLFW_KEY_X) == GLFW_PRESS) {
			q.b = sqrt(1.0 - q.a*q.a);
			q.c = 0.0;
			q.d = 0.0;
		}

		if (glfwGetKey(wndw, GLFW_KEY_Y) == GLFW_PRESS) {
			q.c = sqrt(1.0 - q.a*q.a);
			q.b = 0.0;
			q.d = 0.0;
		}

		if (glfwGetKey(wndw, GLFW_KEY_Z) == GLFW_PRESS) {
			q.d = sqrt(1.0 - q.a*q.a);
			q.c = 0.0;
			q.b = 0.0;
		}

		if (glfwGetMouseButton(wndw, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS) {
			qui_in_prss(&qi, QUI_IN_LMB);
		} else {
			qui_in_rls(&qi, QUI_IN_LMB);
		}
		
		if (glfwGetKey(wndw, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
			qui_in_prss(&qi, QUI_IN_ESC);
		} else {
			qui_in_rls(&qi, QUI_IN_ESC);
		}

		{
			double x = 0.0, y = 0.0;
			glfwGetCursorPos(wndw, &x, &y);
			qui_in_mv(&qi, (float2_t){ (x / (double)w * 2.0 - 1.0), 1.0 - 2.0 * y / (double)h });
		}

		if (qui_scrll) {
			qui_in_scrll(&qi, qui_scrll);
			qui_scrll = 0.0;
		}


		P.m[0][0] = 1.f / ar;
	//	P.m[3][2] = 0.5f;	/* adding perspective */

		PV = mul_float44(V, P);

		glEnable(GL_DEPTH_TEST);
		glUseProgram(cube_po);
		glUniformMatrix4fv(M_loc_cube, 1, 0, &PVM.m[0][0]);
		glBindVertexArray(cube_vao);
		glDrawElements(GL_TRIANGLES, sizeof(cube_i) / sizeof(int), GL_UNSIGNED_INT, 0);

		if (qui_man(&qm, &qs, &qi, P, V, &mt, &mq, &ms)) {
			printf("set\n");
		}

		S = (float44_t){
			ms, 0.f, 0.f, 0.f,
			0.f, ms, 0.f, 0.f,
			0.f, 0.f, ms, 0.f,
			0.f, 0.f, 0.f, 1.f
		};
		T = (float44_t){
			1.f, 0.f, 0.f, 0.f,
			0.f, 1.f, 0.f, 0.f,
			0.f, 0.f, 1.f, 0.f,
			mt.x, mt.y, mt.z, 1.f
		};
		Q = quaternion_to_rotation_matrix(mq);
		M = mul_float44(mul_float44(S, Q), T);
		VM = mul_float44(M, V);
		PVM = mul_float44(VM, P);

		glfwSwapBuffers(wndw);
		glfwWaitEventsTimeout(1.0 / 30.0);

		qui_in_nxt(&qi);

		angl += 0.005;
	}

	glfwTerminate();
	return 0;
}
