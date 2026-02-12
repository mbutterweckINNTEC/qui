struct qui_man {
	int stts;
	int phi;
	float t0, dt, ds;
	int dphi;

	float3_t t;
	quaternion_t q;
	float s;

	int mod;
};

enum {
	QUI_MAN_NIL,
	QUI_MAN_ACT,
	QUI_MAN_SET
};


int qui_man(struct qui_man *qm, float3_t *mt, quaternion_t *mq, float *ms);

int qui_man_mk();
int qui_man_rm();

#ifdef QUI_IMPL

int qui_man_vbo, qui_man_vao;

enum {
	QUI_MAN_MOD_MOUSE,
	QUI_MAN_MOD_KEY
};

/* IMPLEMENTATION */

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

int qui_man_mk() {
	/* visual objects */
	glCreateBuffers(1, &qui_man_vbo);

	if (!qui_man_vbo)
		goto ouch;

	glNamedBufferStorage(qui_man_vbo, QUI_MAN_SZ, NULL, GL_DYNAMIC_STORAGE_BIT | GL_MAP_WRITE_BIT);

	float3_t *v = glMapNamedBufferRange(qui_man_vbo, 0, QUI_MAN_SZ, GL_MAP_WRITE_BIT);

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


	glUnmapNamedBuffer(qui_man_vbo);

	glGenVertexArrays(1, &qui_man_vao);

	if (!qui_man_vao)
		goto ouch;

	glBindVertexArray(qui_man_vao);
	glBindBuffer(GL_ARRAY_BUFFER, qui_man_vbo);
	glVertexAttribPointer(0, 3, GL_FLOAT, 0, 24, 0);
	glEnableVertexAttribArray(0); 
	glVertexAttribPointer(1, 3, GL_FLOAT, 0, 24, (char*)12);
	glEnableVertexAttribArray(1); 
	glBindVertexArray(0);

	return 0;
ouch:
	qui_man_rm();

	return -1;
}

int qui_man_rm() {
	if (qui_man_vbo)
		glDeleteBuffers(1, &qui_man_vbo);

	if (qui_man_vao)
		glDeleteVertexArrays(1, &qui_man_vao);

	return 0;
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

int qui_man_drw(float44_t P, float44_t V, float44_t W, int op[], int flgs) {
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
	glBlendFunc(GL_CONSTANT_COLOR, GL_ONE_MINUS_CONSTANT_COLOR);
	glBlendColor(0.75, 0.75, 0.75, 0.75);

	glUseProgram(qui_shdr_po);
	glBindVertexArray(qui_man_vao);

	if (pass & 3 && lnwdth)
		glLineWidth(lnwdth);

	if (pass & 1) {
		glUniformMatrix4fv(qui_shdr_M, 1, 0, &M.m[0][0]);
		glMultiDrawArrays(GL_LINE_STRIP, sv[0], nv[0], iv[0]);
	}

	if (pass & 2) {
		glUniformMatrix4fv(qui_shdr_M, 1, 0, &W.m[0][0]);
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
		glUniformMatrix4fv(qui_shdr_M, 1, 0, &M.m[0][0]);
		glMultiDrawArrays(GL_TRIANGLE_FAN, sv[2], nv[2], iv[2]);
	}

	if (pass & 8) {
		if (dsc_strt) {
			quaternion_t q = axis_angle_to_quaternion((float3_t) { 0.f, 0.f, 1.f }, -dsc_strt * M_PI / 180.0);
			M = mul_float44(quaternion_to_rotation_matrix(q), W);
		} else {
			M = W;
		}
		glUniformMatrix4fv(qui_shdr_M, 1, 0, &M.m[0][0]);
		glMultiDrawArrays(GL_TRIANGLE_FAN, sv[3], nv[3], iv[3]);
	}

	glDisable(GL_BLEND);

	return 0;
}

#define QUI_MAN_EPS 0.01

int qui_man(struct qui_man *qm, float3_t *mt, quaternion_t *mq, float *ms) {
	float44_t P = qui_mtrx_top(QUI_MTRX_P);
	float44_t V = qui_mtrx_top(QUI_MTRX_V);

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
	float3_t p = float3_float4(cotransform_float44(iPV, (float4_t){ qui_in.p.x, qui_in.p.y, 0.f, 1.f }));
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
	float2_t pv = m_float2(float3_float4(cotransform_float44(iPW, (float4_t){ qui_in.p.x, qui_in.p.y, 0.f, 1.f })));
	float44_t U = {
		1.f, 0.f, 0.f, 0.f,
		0.f, 1.f, 0.f, 0.f,
		0.f, 0.f, 1.f, 0.f,
		V.m30, V.m31, V.m32, 1.f
	};


	int op[16];
	float ds = 1.f;
	int hvrd = 0;


//	qm->dt = qm->t0;

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
			qm->t0 = qui_ray_ln_near(p, d, (float3_t){ 0.f, 0.f, 0.f }, (float3_t){ QUI_MAN_L_XYZ, 0.f, 0.f }).x;
			stts = QUI_MAN_STTS_MOV_X;
		}

		nl = qui_ray_seg_dst(p, d, (float3_t){ 0.f, 0.f, 0.f }, (float3_t){ 0.f, QUI_MAN_L_XYZ, 0.f });
		if (nl < l) {
			l = nl;
			qm->t0 = qui_ray_ln_near(p, d, (float3_t){ 0.f, 0.f, 0.f }, (float3_t){ 0.f, QUI_MAN_L_XYZ, 0.f }).y;
			stts = QUI_MAN_STTS_MOV_Y;
		}

		nl = qui_ray_seg_dst(p, d, (float3_t){ 0.f, 0.f, 0.f }, (float3_t){ 0.f, 0.f, QUI_MAN_L_XYZ });
		if (nl < l) {
			l = nl;
			qm->t0 = qui_ray_ln_near(p, d, (float3_t){ 0.f, 0.f, 0.f }, (float3_t){ 0.f, 0.f, QUI_MAN_L_XYZ }).z;
			stts = QUI_MAN_STTS_MOV_Z;
		}

		if (l * fV < QUI_MAN_EPS) {
			if (qui_in.rls & QUI_IN_LMB) {
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

		if (qm->mod == QUI_MAN_MOD_MOUSE) {
			switch(qm->stts) {
			case QUI_MAN_STTS_ROT_X:
				qui_ray_xcrcl_(p, d, &phi);
				qm->dphi = phi - qm->phi;
				if (qm->dphi < 0) qm->dphi += 360;
				break;
			case QUI_MAN_STTS_ROT_Y:
				qui_ray_ycrcl_(p, d, &phi);
				qm->dphi = phi - qm->phi;
				if (qm->dphi < 0) qm->dphi += 360;
				break;
			case QUI_MAN_STTS_ROT_Z:
				qui_ray_zcrcl_(p, d, &phi);
				qm->dphi = phi - qm->phi;
				if (qm->dphi < 0) qm->dphi += 360;
				break;
			case QUI_MAN_STTS_ROT_V:
				qui_ray_vcrcl_(pv, &phi);
				qm->dphi = phi - qm->phi;
				if (qm->dphi < 0) qm->dphi += 360;
				break;
			case QUI_MAN_STTS_SCL:
				ds = length_float2(pv) / sqrtf(2.f) / QUI_MAN_S_DXY;
				break;
			case QUI_MAN_STTS_MOV_X:
				qm->dt = qui_ray_ln_near(p, d, (float3_t){ 0.f, 0.f, 0.f }, (float3_t){ QUI_MAN_L_XYZ, 0.f, 0.f }).x;
				break;
			case QUI_MAN_STTS_MOV_Y:
				qm->dt = qui_ray_ln_near(p, d, (float3_t){ 0.f, 0.f, 0.f }, (float3_t){ 0.f, QUI_MAN_L_XYZ, 0.f }).y;
				break;
			case QUI_MAN_STTS_MOV_Z:
				qm->dt = qui_ray_ln_near(p, d, (float3_t){ 0.f, 0.f, 0.f }, (float3_t){ 0.f, 0.f, QUI_MAN_L_XYZ }).z;
				break;
			};

			if (qui_in.rls & QUI_IN_LMB) {
				rstts = QUI_MAN_SET;
			}
		}

		if (qui_in.rls & QUI_IN_ESC) {
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

			qm->mod = QUI_MAN_MOD_MOUSE;
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
		qui_man_drw(P, V, W, (int[]){QUI_MAN_OP_LN_WDTH, 2, QUI_MAN_OP_END}, QUI_MAN_FLGS_AXIS_X);
		break;
	case QUI_MAN_STTS_MOV_Y:
		qui_man_drw(P, V, W, (int[]){QUI_MAN_OP_LN_WDTH, 2, QUI_MAN_OP_END}, QUI_MAN_FLGS_AXIS_Y);
		break;
	case QUI_MAN_STTS_MOV_Z:
		qui_man_drw(P, V, W, (int[]){QUI_MAN_OP_LN_WDTH, 2, QUI_MAN_OP_END}, QUI_MAN_FLGS_AXIS_Z);
		break;
	case QUI_MAN_STTS_NIL:
		qui_man_drw(P, V, W, (int[]){QUI_MAN_OP_LN_WDTH, 2, QUI_MAN_OP_END}, QUI_MAN_FLGS_AXIS_CRCL_FRM);
		break;
	case QUI_MAN_STTS_ROT_X:
		qui_man_drw(P, V, W, (int[]){QUI_MAN_OP_PIE_STRT, qm->phi, QUI_MAN_OP_PIE_ANGL, qm->dphi, QUI_MAN_OP_LN_WDTH, 2, QUI_MAN_OP_END}, QUI_MAN_FLGS_CRCL_X | QUI_MAN_FLGS_PIE_X);
		break;
	case QUI_MAN_STTS_ROT_Y:
		qui_man_drw(P, V, W, (int[]){QUI_MAN_OP_PIE_STRT, qm->phi, QUI_MAN_OP_PIE_ANGL, qm->dphi, QUI_MAN_OP_LN_WDTH, 2, QUI_MAN_OP_END}, QUI_MAN_FLGS_CRCL_Y | QUI_MAN_FLGS_PIE_Y);
		break;
	case QUI_MAN_STTS_ROT_Z:
		qui_man_drw(P, V, W, (int[]){QUI_MAN_OP_PIE_STRT, qm->phi, QUI_MAN_OP_PIE_ANGL, qm->dphi, QUI_MAN_OP_LN_WDTH, 2, QUI_MAN_OP_END}, QUI_MAN_FLGS_CRCL_Z | QUI_MAN_FLGS_PIE_Z);
		break;
	case QUI_MAN_STTS_ROT_V:
		qui_man_drw(P, V, W, (int[]){QUI_MAN_OP_PIE_STRT, qm->phi, QUI_MAN_OP_PIE_ANGL, qm->dphi, QUI_MAN_OP_LN_WDTH, 2, QUI_MAN_OP_END}, QUI_MAN_FLGS_CRCL_V | QUI_MAN_FLGS_PIE_V);
		break;
		case QUI_MAN_STTS_SCL:
		qui_man_drw(P, V, W, (int[]){QUI_MAN_OP_LN_WDTH, 2, QUI_MAN_OP_END}, QUI_MAN_FLGS_FRM);
		break;
	};

	switch(hvrd) {
	case QUI_MAN_STTS_NIL:
		break;
	case QUI_MAN_STTS_MOV_X:
		qui_man_drw(P, V, W, (int[]){QUI_MAN_OP_LN_WDTH, 4, QUI_MAN_OP_END}, QUI_MAN_FLGS_AXIS_X);
		break;
	case QUI_MAN_STTS_MOV_Y:
		qui_man_drw(P, V, W, (int[]){QUI_MAN_OP_LN_WDTH, 4, QUI_MAN_OP_END}, QUI_MAN_FLGS_AXIS_Y);
		break;
	case QUI_MAN_STTS_MOV_Z:
		qui_man_drw(P, V, W, (int[]){QUI_MAN_OP_LN_WDTH, 4, QUI_MAN_OP_END}, QUI_MAN_FLGS_AXIS_Z);
		break;
	case QUI_MAN_STTS_ROT_X:
		qui_man_drw(P, V, W, (int[]){QUI_MAN_OP_LN_WDTH, 4, QUI_MAN_OP_END}, QUI_MAN_FLGS_CRCL_X);
		break;
	case QUI_MAN_STTS_ROT_Y:
		qui_man_drw(P, V, W, (int[]){QUI_MAN_OP_LN_WDTH, 4, QUI_MAN_OP_END}, QUI_MAN_FLGS_CRCL_Y);
		break;
	case QUI_MAN_STTS_ROT_Z:
		qui_man_drw(P, V, W, (int[]){QUI_MAN_OP_LN_WDTH, 4, QUI_MAN_OP_END}, QUI_MAN_FLGS_CRCL_Z);
		break;
	case QUI_MAN_STTS_ROT_V:
		qui_man_drw(P, V, W, (int[]){QUI_MAN_OP_LN_WDTH, 4, QUI_MAN_OP_END}, QUI_MAN_FLGS_CRCL_V);
		break;
		case QUI_MAN_STTS_SCL:
		qui_man_drw(P, V, W, (int[]){QUI_MAN_OP_LN_WDTH, 4, QUI_MAN_OP_END}, QUI_MAN_FLGS_FRM);
		break;
	};

	float3_t bg;
	char *nm = "";
	char *unt = "";
	switch(qm->stts) {
	case QUI_MAN_STTS_NIL:
		break;
	case QUI_MAN_STTS_ROT_X:
		nm = "∢x";
		unt = "°";
		bg = QUI_MAN_X_CLR;
		break;
	case QUI_MAN_STTS_MOV_X:
		nm = "δx";
		unt = "mm";
		bg = QUI_MAN_X_CLR;
		break;
	case QUI_MAN_STTS_ROT_Y:
		nm = "∢y";
		unt = "°";
		bg = QUI_MAN_Y_CLR;
		break;
	case QUI_MAN_STTS_MOV_Y:
		nm = "δy";
		unt = "mm";
		bg = QUI_MAN_Y_CLR;
		break;
	case QUI_MAN_STTS_ROT_Z:
		nm = "∢z";
		unt = "°";
		bg = QUI_MAN_Z_CLR;
		break;
	case QUI_MAN_STTS_MOV_Z:
		nm = "δz";
		unt = "mm";
		bg = QUI_MAN_Z_CLR;
		break;
	case QUI_MAN_STTS_ROT_V:
		nm = "∢v";
		unt = "°";
		bg = QUI_MAN_V_CLR;
		break;
	case QUI_MAN_STTS_SCL:
		nm = "↔";
		bg = QUI_MAN_S_CLR;
		break;
	};

	int val_flgs = qm->mod == QUI_MAN_MOD_MOUSE ? QUI_VAL_FLGS_RST : 0;

	qui_mtrx_psh(QUI_MTRX_V, U);

	float44_t M = {
		1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0,
		-0.75*fV, -0.75*fV, 0, 1
	};

	switch(qm->stts) {
	case QUI_MAN_STTS_NIL:
		break;
	case QUI_MAN_STTS_ROT_X:
	case QUI_MAN_STTS_ROT_Y:
	case QUI_MAN_STTS_ROT_Z:
	case QUI_MAN_STTS_ROT_V:
		if (qui_val_i(M, bg, nm, unt, &qm->dphi, val_flgs)) {
			while (qm->dphi < 0) qm->dphi += 360;
			while (360 <= qm->dphi) qm->dphi -= 360;
			qm->mod = QUI_MAN_MOD_KEY;
		}
		break;
	case QUI_MAN_STTS_MOV_X:
	case QUI_MAN_STTS_MOV_Y:
	case QUI_MAN_STTS_MOV_Z:
		if (qui_val_f(M, bg, nm, unt, &qm->dt, val_flgs)) {
			qm->mod = QUI_MAN_MOD_KEY;
		}
		break;
	case QUI_MAN_STTS_SCL:
		if (qm->mod == QUI_MAN_MOD_MOUSE) {
			qm->ds = *ms * ds / qm->s;
		}
		switch (qui_val_f(M, bg, nm, unt, &qm->ds, val_flgs)) {
		case QUI_VAL_RET_NIL: break;
		case QUI_VAL_RET_SET:
			if (qm->ds < 0.001)
				qm->ds = 1.f;
		case QUI_VAL_RET_ED:
			qm->mod = QUI_MAN_MOD_KEY;

			if (qm->ds < 0.f)
				qm->ds *= -1.f;

			ds = qm->s * fmax(0.001, qm->ds) / *ms;
			break;
		};
		break;
	};

	qui_mtrx_pop(QUI_MTRX_V);

	switch(qm->stts) {
	case QUI_MAN_STTS_NIL:
		break;
	case QUI_MAN_STTS_ROT_X:
		*mq = mul_quaternion(qm->q, axis_angle_to_quaternion((float3_t) {1.f, 0.f, 0.f }, -qm->dphi * M_PI / 180.0));
		break;
	case QUI_MAN_STTS_ROT_Y:
		*mq = mul_quaternion(qm->q, axis_angle_to_quaternion((float3_t) { 0.f, 1.f, 0.f }, qm->dphi * M_PI / 180.0));
		break;
	case QUI_MAN_STTS_ROT_Z:
		*mq = mul_quaternion(qm->q, axis_angle_to_quaternion((float3_t) { 0.f, 0.f, 1.f }, -qm->dphi * M_PI / 180.0));
		break;
	case QUI_MAN_STTS_ROT_V:
		*mq = mul_quaternion(qm->q, axis_angle_to_quaternion(d, qm->dphi * M_PI / 180.0));
		break;
	case QUI_MAN_STTS_MOV_X:
		*mt = add_float3(qm->t, (float3_t){ qm->dt - qm->t0, 0.f, 0.f });
		break;
	case QUI_MAN_STTS_MOV_Y:
		*mt = add_float3(qm->t, (float3_t){ 0.f, qm->dt - qm->t0, 0.f });
		break;
	case QUI_MAN_STTS_MOV_Z:
		*mt = add_float3(qm->t, (float3_t){  0.f, 0.f, qm->dt - qm->t0 });
		break;
	};

	if (qui_in.rls & QUI_IN_RET) {
		rstts = QUI_MAN_SET;
	}

	*ms *= ds;

	if (rstts == QUI_MAN_SET) {
			qm->stts = QUI_MAN_STTS_NIL;
			qm->mod = QUI_MAN_MOD_MOUSE;
	}

	return rstts;
}

#endif /* QUI_IMPL */