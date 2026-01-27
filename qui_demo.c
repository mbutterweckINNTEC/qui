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

static int shdr_mk(char  const *vshsrc, char  const *fshsrc) {
        char log[4096];
        int po = 0, vso = 0, fso = 0, stts = 0;

        vso = glCreateShader(GL_VERTEX_SHADER);
        if (!vso)
                goto end;

        glShaderSource(vso, 1, &vshsrc, NULL);
        glCompileShader(vso);

        fso = glCreateShader(GL_FRAGMENT_SHADER);
        if (!fso)
                goto end;

        glShaderSource(fso, 1, (char const * const*)&fshsrc, NULL);
        glCompileShader(fso);

        glGetShaderiv(vso, GL_COMPILE_STATUS, &stts);
        if (!stts) {
                glGetShaderInfoLog(vso, 4096, NULL, log);
                printf("vsh: %d\n", stts);
                puts(log);
                goto end;
        }

        glGetShaderiv(fso, GL_COMPILE_STATUS, &stts);
        if (!stts) {
                glGetShaderInfoLog(fso, 4096, NULL, log);
                printf("fsh: %d\n", stts);
                puts(log);
                goto end;
        }

        po = glCreateProgram();
        if (!po)
                goto end;

        glAttachShader(po, vso);
        glAttachShader(po, fso);
        glLinkProgram(po);

        glGetProgramiv(po, GL_LINK_STATUS, &stts);
        if (!stts) {
                glGetProgramInfoLog(po, 4096, NULL, log);
                printf("po: %d\n", stts);
                puts(log);
                goto end;
        }

end:
        if (!stts || 0 == po) {
                glDeleteProgram(po);
                po = 0;
        }

        if (vso)
                glDeleteShader(vso);

        if (fso)
                glDeleteShader(fso);

        return po;
}

char *qui_vsh =
	"#version 440"				"\n"
	"in      vec4 v;"			"\n"
	"in      vec4 c;"			"\n"
	"uniform mat4 M;"			"\n"
	"out     vec4 k;"			"\n"
	""					"\n"
	"void main() {"				"\n"
	"	gl_Position = M * v;"		"\n"
	"	k = c;"				"\n"
	"}"					"\n";

char *qui_fsh =
	"#version 440"				"\n"
	"in  vec4 k;"				"\n"
	"out vec4 K;"				"\n"
	""					"\n"
	"void main() {"				"\n"
	"	K = k;"				"\n"
	"}"					"\n";

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

static float3_t rnbw(float t) {
	int s = t;
	int r = t - s;
	switch(s) {
	case 0: return (float3_t) { 1.f, r, 0.f };
	case 1: return (float3_t) { 1.f - r, 1.f, 0.f };
	case 2: return (float3_t) { 0.f, 1.f, r };
	case 3: return (float3_t) { 0.f, 1.f - r, 1.f };
	case 4: return (float3_t) { r, 1.f, 0.f };
	case 5: return (float3_t) { 1.f, 0.f, 1.f - r };
	};

	return (float3_t) { 0.f, 0.f, 0.f };
}

struct qui_shdr {
	int po;

	/* uniform locations */
	int M;
};

int qui_shdr_mk(struct qui_shdr *qs) {
	if (!qs)
		return -1;

	qs->po = shdr_mk(qui_vsh, qui_fsh);

	if (!qs->po)
		return -1;

	qs->M = glGetUniformLocation(qs->po, "M");

	if (qs->M == -1)
		return -1;

	return 0;
}

enum {
	QUI_IN_LMB = 0x1,
	QUI_IN_X = 0x2,
	QUI_IN_Y = 0x4,
	QUI_IN_Z = 0x8,
	QUI_IN_0 = 0x10,
	QUI_IN_1 = 0x20,
	QUI_IN_2 = 0x40,
	QUI_IN_3 = 0x80,
	QUI_IN_4 = 0x100,
	QUI_IN_5 = 0x200,
	QUI_IN_6 = 0x400,
	QUI_IN_7 = 0x800,
	QUI_IN_8 = 0x1000,
	QUI_IN_9 = 0x2000,
	QUI_IN_DOT = 0x4000,
	QUI_IN_MINUS = 0x8000,
	QUI_IN_ESC = 0x10000,

	QUI_IN_ALL = ~0
};

struct qui_in {
	int prss;
	int rls;

	float2_t p;
	float2_t d;

	float s;
};

static int qui_in_prss(struct qui_in *qi, int bttn);
static int qui_in_rls(struct qui_in *qin, int bttn);
static int qui_in_mv(struct qui_in *qi, float2_t p);
static int qui_in_scrll(struct qui_in *qi, float scrll);
static int qui_in_nxt(struct qui_in *qi);

static int qui_in_prss(struct qui_in *qi, int bttn) {
	if (!qi)
		return -1;

	qi->prss |= bttn;
	qi->rls &=~ bttn;

	return 0;
}

static int qui_in_rls(struct qui_in *qi, int bttn) {
	if (!qi)
		return -1;

	qi->rls |= qi->prss & bttn;
	qi->prss &=~ bttn;

	return 0;
}

static int qui_in_mv(struct qui_in *qi, float2_t p) {
	if (!qi)
		return -1;

	qi->d = add_float2(qi->d, sub_float2(p, qi->p));
	qi->p = p;

	return 0;
}

static int qui_in_scrll(struct qui_in *qi, float scrll) {
	if (!qi)
		return -1;

	qi->s += scrll;

	return 0;
}

static int qui_in_nxt(struct qui_in *qi) {
	if (!qi)
		return -1;

	qi->rls = 0;
	qi->d = (float2_t){ 0.f, 0.f };
	qi->s = 0.f;

	return 0;
}

#define QUI_MAN_AXIS_X_N 2
#define QUI_MAN_CRCL_X_N 360
#define QUI_MAN_AXIS_Y_N 2
#define QUI_MAN_CRCL_Y_N 360
#define QUI_MAN_AXIS_Z_N 2
#define QUI_MAN_CRCL_Z_N 360
#define QUI_MAN_AXIS_V_N 2
#define QUI_MAN_CRCL_V_N 360

#define QUI_MAN_PLN_X_N (QUI_MAN_AXIS_X_N + QUI_MAN_CRCL_X_N)
#define QUI_MAN_PLN_Y_N (QUI_MAN_AXIS_Y_N + QUI_MAN_CRCL_Y_N)
#define QUI_MAN_PLN_Z_N (QUI_MAN_AXIS_Z_N + QUI_MAN_CRCL_Z_N)
#define QUI_MAN_PLN_V_N (QUI_MAN_AXIS_Z_N + QUI_MAN_CRCL_Z_N)

#define QUI_MAN_RSZ_N 12

#define QUI_MAN_ALL_N (QUI_MAN_PLN_X_N + QUI_MAN_PLN_Y_N + QUI_MAN_PLN_Z_N + QUI_MAN_PLN_V_N + QUI_MAN_RSZ_N)

#define QUI_MAN_ATTR_N 2
#define QUI_MAN_SZ (QUI_MAN_ALL_N * QUI_MAN_ATTR_N * sizeof(float3_t))

#define QUI_MAN_X_CLR (float3_t) { 1.0f, 0.5f, 0.5f };
#define QUI_MAN_Y_CLR (float3_t) { 0.5f, 1.0f, 0.5f };
#define QUI_MAN_Z_CLR (float3_t) { 0.5f, 0.5f, 1.0f };
#define QUI_MAN_V_CLR (float3_t) { 1.0f, 0.75f, 0.5f };	/* view axis */
#define QUI_MAN_S_CLR (float3_t) { 0.5f, 0.75f, 1.0f };	/* resize frame */

#define QUI_MAN_R_XYZ 0.875f
#define QUI_MAN_R_V 1.f
#define QUI_MAN_L_XYZ 0.75f

#define QUI_MAN_S_XY 0.875f
#define QUI_MAN_S_DXY 0.9375f

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
	QUI_MAN_STTS_SCL_U		/* only uniform scalling */
};

struct qui_man {
	int bo, vao;

	int stts;
	int phi, dphi;
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
	QUI_MAN_OP_LN_WDTH
};

#define QUI_MAN_DRW_PASS_N 4

int qui_man_drw(struct qui_man *qm, struct qui_shdr *qs, float44_t P, float44_t V, int op[], int flgs) {
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

	//float dV = cbrt(det_float33(float33_float44(V)));
//	float dV = sqrtf(V.m00 * V.m11 - V.m01 * V.m10);
	float sV = frobenius_float33(float33_float44(V)) / sqrt(3);

	float44_t DP = {
		sV, 0.f, 0.f, 0.f,
		0.f, sV, 0.f, 0.f,
		0.f, 0.f, 1.f, 0.f,
		0.f, 0.f, 0.f, 1.f
	};

	DP = mul_float44(DP, P);

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

	//glDisable(GL_DEPTH_TEST);
	glUseProgram(qs->po);
	glBindVertexArray(qm->vao);

	if (pass & 3 && lnwdth)
		glLineWidth(lnwdth);

	if (pass & 1) {
		glUniformMatrix4fv(qs->M, 1, 0, &M.m[0][0]);
		glMultiDrawArrays(GL_LINE_STRIP, sv[0], nv[0], iv[0]);
	}

	if (pass & 2) {
		glUniformMatrix4fv(qs->M, 1, 0, &DP.m[0][0]);
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
			M = mul_float44(quaternion_to_rotation_matrix(q), P);
		} else {
			M = P;
		}
		glUniformMatrix4fv(qs->M, 1, 0, &M.m[0][0]);
		glMultiDrawArrays(GL_TRIANGLE_FAN, sv[3], nv[3], iv[3]);
	}

//	glEnable(GL_DEPTH_TEST);

	return 0;
}

static inline float ray_pnt_dst(float3_t ro, float3_t rd, float3_t p) {
	float t = dot_float3(rd, sub_float3(p, ro)) / dot_float3(rd, rd);
	return length_float3(sub_float3(p, add_float3(ro, scale_float3(rd, t))));
}

#define QUI_MAN_EPS 0.01

/* todo:@michal: make analytic function? */

static float qui_ray_xcrcl_(float3_t p, float3_t d, int *out) {
	float3_t v;
	float l = FLT_MAX, vr, phi;
	for (int i = 0; i < QUI_MAN_CRCL_X_N; ++i) {
		phi = (float)i / (float)(QUI_MAN_CRCL_X_N - 1) * 2.0 * M_PI;
		v = (float3_t) { 0.f, QUI_MAN_R_XYZ * cos(phi), QUI_MAN_R_XYZ * sin(phi) };
		vr = ray_pnt_dst(p, d, v);

		if (vr < l) {
			l = vr;
			*out = i;
		}
	}
	return l;
}

static float qui_ray_ycrcl_(float3_t p, float3_t d, int *out) {
	float3_t v;
	float l = FLT_MAX, vr, phi;
	for (int i = 0; i < QUI_MAN_CRCL_Y_N; ++i) {
		phi = (float)i / (float)(QUI_MAN_CRCL_Y_N - 1) * 2.0 * M_PI;
		v = (float3_t) { QUI_MAN_R_XYZ * cos(phi), 0.f, QUI_MAN_R_XYZ * sin(phi) };
		vr = ray_pnt_dst(p, d, v);

		if (vr < l) {
			l = vr;
			*out = i;
		}
	}
	return l;
}

static float qui_ray_zcrcl_(float3_t p, float3_t d, int *out) {
	float3_t v;
	float l = FLT_MAX, vr, phi;
	for (int i = 0; i < QUI_MAN_CRCL_Z_N; ++i) {
		phi = (float)i / (float)(QUI_MAN_CRCL_Z_N - 1) * 2.0 * M_PI;
		v = (float3_t) { QUI_MAN_R_XYZ * cos(phi), QUI_MAN_R_XYZ * sin(phi), 0.f };
		vr = ray_pnt_dst(p, d, v);

		if (vr < l) {
			l = vr;
			*out = i;
		}
	}
	return l;
}

static float qui_ray_vcrcl_(float2_t p, int *out) {
	float2_t v;
	float l = FLT_MAX, vr, phi;
	for (int i = 0; i < QUI_MAN_CRCL_V_N; ++i) {
		phi = (float)i / (float)(QUI_MAN_CRCL_V_N - 1) * 2.0 * M_PI;
		v = (float2_t) { QUI_MAN_R_V * cos(phi), QUI_MAN_R_V * sin(phi) };
		vr = length_float2(sub_float2(v, p));	

		if (vr < l) {
			l = vr;
			*out = i;
		}
	}

	printf ("vl = %f, phi = %d\n", l, *out);
	return l;
}

enum {
	QUI_MAN_NIL,
	QUI_MAN_MOV,
	QUI_MAN_ROT,
	QUI_MAN_SCL
};

int qui_man(struct qui_man *qm, struct qui_shdr *qs, struct qui_in *qi, float44_t P, float44_t V, quaternion_t *out) {
	float44_t PV = (mul_float44(V, P));
	float detPV = det_float44(PV);
	float44_t iPV = invert_float44(PV, detPV);
	float3_t p = float3_float4(cotransform_float44(iPV, (float4_t){ qi->p.x, qi->p.y, 0.f, 1.f }));
	float3_t d = normal_float3(m_float3(cotransform_float44(iPV, (float4_t){ 0.f, 0.f, 1.f, 0.f })));
	float2_t pv = (float2_t){ qi->p.x * P.m11 / P.m00, qi->p.y };
	int op[16], stts = 0;

	/* input handling */
	if (QUI_MAN_STTS_NIL == qm->stts) {
		if (qi->rls & QUI_IN_LMB) {
			float l, nl, scl = frobenius_float33(float33_float44(V))/ sqrt(3.0);
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

			if (l * scl < QUI_MAN_EPS) {
				qm->phi = phi;
				qm->stts = stts;
			}
		}
	} else {
		int phi = 0;
		switch(qm->stts) {
		case QUI_MAN_STTS_ROT_X:
			qui_ray_xcrcl_(p, d, &phi);
			break;
		case QUI_MAN_STTS_ROT_Y:
			qui_ray_ycrcl_(p, d, &phi);
			break;
		case QUI_MAN_STTS_ROT_Z:
			qui_ray_zcrcl_(p, d, &phi);
			break;
		case QUI_MAN_STTS_ROT_V:
			qui_ray_vcrcl_(pv, &phi);
			break;
		};

		printf("angles %d %d", qm->phi, phi);

		qm->dphi = phi - qm->phi;

		if (qm->dphi < 0)
			qm->dphi += 360;
		printf(" %d\n", qm->dphi);

		if (qi->rls & QUI_IN_LMB) {
			stts = qm->stts;
			qm->stts = QUI_MAN_STTS_NIL;
		}

		if (qi->rls & QUI_IN_ESC) {
			stts = QUI_MAN_STTS_NIL;
			qm->stts = QUI_MAN_STTS_NIL;
		}
	}

	switch(qm->stts) {
	case QUI_MAN_STTS_NIL:
		qui_man_drw(qm, qs, P, V, (int[]){QUI_MAN_OP_LN_WDTH, 2, QUI_MAN_OP_END}, QUI_MAN_FLGS_AXIS_CRCL_FRM);
		break;
	case QUI_MAN_STTS_ROT_X:
		qui_man_drw(qm, qs, P, V, (int[]){QUI_MAN_OP_PIE_STRT, qm->phi, QUI_MAN_OP_PIE_ANGL, qm->dphi, QUI_MAN_OP_LN_WDTH, 2, QUI_MAN_OP_END}, QUI_MAN_FLGS_CRCL_X | QUI_MAN_FLGS_PIE_X);
		break;
	case QUI_MAN_STTS_ROT_Y:
		qui_man_drw(qm, qs, P, V, (int[]){QUI_MAN_OP_PIE_STRT, qm->phi, QUI_MAN_OP_PIE_ANGL, qm->dphi, QUI_MAN_OP_LN_WDTH, 2, QUI_MAN_OP_END}, QUI_MAN_FLGS_CRCL_Y | QUI_MAN_FLGS_PIE_Y);
		break;
	case QUI_MAN_STTS_ROT_Z:
		qui_man_drw(qm, qs, P, V, (int[]){QUI_MAN_OP_PIE_STRT, qm->phi, QUI_MAN_OP_PIE_ANGL, qm->dphi, QUI_MAN_OP_LN_WDTH, 2, QUI_MAN_OP_END}, QUI_MAN_FLGS_CRCL_Z | QUI_MAN_FLGS_PIE_Z);
		break;
	case QUI_MAN_STTS_ROT_V:
		qui_man_drw(qm, qs, P, V, (int[]){QUI_MAN_OP_PIE_STRT, qm->phi, QUI_MAN_OP_PIE_ANGL, qm->dphi, QUI_MAN_OP_LN_WDTH, 2, QUI_MAN_OP_END}, QUI_MAN_FLGS_CRCL_V | QUI_MAN_FLGS_PIE_V);
		break;
	};

	switch(qm->stts | stts) {
	case QUI_MAN_STTS_NIL:
		*out = (quaternion_t) { 1.f, 0.f, 0.f, 0.f };
		break;
	case QUI_MAN_STTS_ROT_X:
		*out = axis_angle_to_quaternion((float3_t) {1.f, 0.f, 0.f }, -qm->dphi * M_PI / 180.0);
		break;
	case QUI_MAN_STTS_ROT_Y:
		*out = axis_angle_to_quaternion((float3_t) { 0.f, 1.f, 0.f }, qm->dphi * M_PI / 180.0);
		break;
	case QUI_MAN_STTS_ROT_Z:
		*out = axis_angle_to_quaternion((float3_t) { 0.f, 0.f, 1.f }, -qm->dphi * M_PI / 180.0);
		break;
	case QUI_MAN_STTS_ROT_V:
		*out = axis_angle_to_quaternion(d, qm->dphi * M_PI / 180.0);
		break;
	};

	return stts;
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
	float44_t V = identity_sc;
	float44_t VM = identity_sc;
	float44_t P = orthographic(1.f, -1.f, 1.f, -1.f, -1.f, 1.f);
	float44_t PV = identity_sc;
	float44_t PVM = identity_sc;

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

	while(!glfwWindowShouldClose(wndw)) {
		glfwGetWindowSize(wndw, &w, &h);
	
		ar = (float)w / (float)h;

		glViewport(0, 0, w, h);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		if (glfwGetKey(wndw, GLFW_KEY_UP) == GLFW_PRESS) {
			V = mul_float44(V, quaternion_to_rotation_matrix((quaternion_t){cos(0.03125),-sin(0.03125), 0, 0 }));
		}

		if (glfwGetKey(wndw, GLFW_KEY_DOWN) == GLFW_PRESS) {
			V = mul_float44(V, quaternion_to_rotation_matrix((quaternion_t){cos(0.03125), sin(0.03125), 0, 0 }));
		}

		if (glfwGetKey(wndw, GLFW_KEY_LEFT) == GLFW_PRESS) {
			V = mul_float44(V, quaternion_to_rotation_matrix((quaternion_t){cos(0.03125), 0,-sin(0.03125), 0 }));
		}

		if (glfwGetKey(wndw, GLFW_KEY_RIGHT) == GLFW_PRESS) {
			V = mul_float44(V, quaternion_to_rotation_matrix((quaternion_t){cos(0.03125), 0, sin(0.03125), 0 }));
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

		quaternion_t obr;

		if (qui_man(&qm, &qs, &qi, (P),(V), &obr)) {
			printf("set\n");
			M = mul_float44(M, quaternion_to_rotation_matrix(obr));
			obr = (quaternion_t){ 1.f, 0.f, 0.f, 0.f };
		}

		Q = quaternion_to_rotation_matrix(obr);
		VM = mul_float44(M, mul_float44(Q, V));
		PVM = mul_float44(VM, P);

		glEnable(GL_DEPTH_TEST);
		glUseProgram(cube_po);
		glUniformMatrix4fv(M_loc_cube, 1, 0, &PVM.m[0][0]);
		glBindVertexArray(cube_vao);
		glDrawElements(GL_TRIANGLES, sizeof(cube_i) / sizeof(int), GL_UNSIGNED_INT, 0);


		glfwSwapBuffers(wndw);
		glfwWaitEventsTimeout(1.0 / 30.0);

		qui_in_nxt(&qi);

		angl += 0.005;
	}

	glfwTerminate();
	return 0;
}
