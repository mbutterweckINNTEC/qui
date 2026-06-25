#ifndef QUI_VAL_H
#define QUI_VAL_H

enum {
	QUI_VAL_FLGS_RST = 0x1,		/* reset value on edit */
	QUI_VAL_FLGS_CNST = 0x2		/* immutable */
};

enum {
	QUI_VAL_RET_NIL,
	QUI_VAL_RET_ED,
	QUI_VAL_RET_SET
};

int qui_val_i(float44_t M, float3_t bg_, char *nm, char *unt, int *val, int flgs);
int qui_val_f(float44_t M, float3_t bg_, char *nm, char *unt, float *val, int flgs);
int qui_val_g(float44_t M, float3_t bg_, char *nm, char *unt, double *val, int flgs);

int qui_val_hvr(float44_t M);

/* Implementation */
#ifdef QUI_IMPL

static float const qui_val_du = 0.05273437f;
static float const qui_val_scl = 0.03125;
static float3_t const qui_val_mv = { 0.5 * qui_val_scl, 0.02 / 0.0625 * qui_val_scl };

static float3_t qui_val_nm_ngon[] = {
	{ -1.f, 0.f, 0.f },
	{-0.5f, 0.8660254038f, 0.f },
	{ 0.5f, 0.8660254038f, 0.f },
	{ 1.0, 0.0, 0.f},
	{ 0.5f,-0.8660254038f, 0.f },
	{-0.5f,-0.8660254038f, 0.f },
};

static float3_t qui_val_unt_ngon[] = {
	{ -1.f + 6.375f + 0.35f, 0.f           - 2 * 0.8660254038f - 2 * 0.08660254038f, 0.f },
	{-0.5f + 6.375f + 0.35f, 0.8660254038f - 2 * 0.8660254038f - 2 * 0.08660254038f, 0.f },
	{ 0.5f + 6.375f + 0.35f, 0.8660254038f - 2 * 0.8660254038f - 2 * 0.08660254038f, 0.f },
	{ 1.0f + 6.375f + 0.35f, 0.0f          - 2 * 0.8660254038f - 2 * 0.08660254038f, 0.f },
	{ 0.5f + 6.375f + 0.35f,-0.8660254038f - 2 * 0.8660254038f - 2 * 0.08660254038f, 0.f },
	{-0.5f + 6.375f + 0.35f,-0.8660254038f - 2 * 0.8660254038f - 2 * 0.08660254038f, 0.f },
};

static float3_t qui_val_val_ngon[] = {
	{ -1.f + 1.6875f, 0.f           - 0.8660254038f - 0.08660254038f, 0.f },
	{-0.5f + 1.6875f, 0.8660254038f - 0.8660254038f - 0.08660254038f, 0.f },
	{ 3.875f + 1.6875f, 0.8660254038f - 0.8660254038f - 0.08660254038f, 0.f },
	{ 4.375f + 1.6875f, 0.0f          - 0.8660254038f - 0.08660254038f, 0.f },
	{ 3.875f + 1.6875f,-0.8660254038f - 0.8660254038f - 0.08660254038f, 0.f },
	{-0.5f + 1.6875f,-0.8660254038f - 0.8660254038f - 0.08660254038f, 0.f },
};

int qui_val_hvr(float44_t M_) {
	float44_t P = qui_mtrx_top(QUI_MTRX_P);
	float44_t V = qui_mtrx_top(QUI_MTRX_V);
	float44_t S = {
		qui_val_scl, 0, 0, 0,
		0, qui_val_scl, 0, 0,
		0, 0, 1, 0,
		- 3.f * qui_val_du, 0, 0, 1
        };
	float44_t M = mul_float44(S, M_);
	float44_t PVM = mul_float44(M, mul_float44(V, P));

	float4_t x[3] = {
		{ 1.5 + 0.1f, - 0.8660254038f - 0.08660254038f, 0.f, 1.f },
		{ 3.0 + 1.5 + 0.1f, - 0.8660254038f - 0.08660254038f, 0.f, 1.f },
		{ 2.0 + 1.5 + 0.1f, - 0.8660254038f - 0.08660254038f, 0.f, 1.f }
	};

	float r = frobenius_float33(float33_float44(PVM));

	for (int i = 0; i < 3; ++i) {
		float3_t p = float3_float4(cotransform_float44(PVM, x[i]));
		float dx = p.x - qui_in.p.x;
		float dy = p.y - qui_in.p.y;

		float d = sqrt(dx * dx + dy * dy);

		if (d < r) {
			qui_in.st |= QUI_IN_ST_HVRD;
			return 1;
		}
	}

	return 0;
}

static int qui_val_drw(float44_t M_, float3_t clr, char *nm, char *unt, char *val) {
	float4_t bg = m_float4(clr, 1.f);
	float4_t fg = m_float4(mix_float3(clr, (float3_t){1,1,1}, 0.75), 1.f);
	float4_t fgv = m_float4(mix_float3(clr, (float3_t){0,0,0}, 0.75), 1.f);

	float44_t S = {
		qui_val_scl, 0, 0, 0,
		0, qui_val_scl, 0, 0,
		0, 0, 1, 0,
		- 3.f * qui_val_du, 0, 0, 1
	};

	float44_t Z = {
		qui_val_scl, 0, 0, 0,
		0, qui_val_scl, 0, 0,
		0, 0, 1, 0,
		-qui_val_mv.x - 3.f * qui_val_du, -qui_val_mv.y, 0, 1
	};

	float44_t M = mul_float44(S, M_);

	float44_t T = mul_float44(Z, M_);
		
	float44_t T2 = mul_float44(
		(float44_t) {
			1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 1, 0,
			2.f, -1, 0, 1
		},
		T
	);

	float44_t T3 = mul_float44(
		(float44_t) {
			1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 1, 0,
			6.75, -2, 0, 1
		},
		T
	);

	if (nm) {
		qui_ngon(6, qui_val_nm_ngon, M, bg);
		qui_txt(nm, T, fg);
	}
	qui_ngon(6, qui_val_val_ngon, M, bg);
	qui_txt(val, T2, fgv);
	if (unt) {
		qui_ngon(6, qui_val_unt_ngon, M, bg);
		qui_txt(unt, T3, fg);
	}

	return 0;
}

int qui_val_i(float44_t M, float3_t bg_, char *nm, char *unt, int *val, int flgs) {
	char sval[64];

	sprintf(sval, "%4d", *val);
	qui_val_drw(M, bg_, nm, unt, sval);

	if (flgs & QUI_VAL_FLGS_CNST)
		return QUI_VAL_RET_NIL;

	if (qui_in.st == QUI_IN_ST_CNSMD)
		return QUI_VAL_RET_NIL;

	qui_tip_sgnl |= QUI_TIP_SGNL_FCS & qui_tip_msk;

	if (qui_in.rls & QUI_IN_RET) {
		qui_in.st |= QUI_IN_ST_CNSMD;
		return QUI_VAL_RET_SET;
	}

	if (qui_in.rls & QUI_IN_BCK) {
		sprintf(sval, "%d", *val);
		int sl = strlen(sval);

		if (sl) {
			sval[--sl] = '\0';

			if (0 == sl) {
				*val = 0;
				return QUI_VAL_RET_ED;
			}
				
			if (1 == sscanf(sval, "%d", val)) {
				return QUI_VAL_RET_ED;
			}
		}
		qui_in.st |= QUI_IN_ST_CNSMD;
	}

	if (qui_in.rls & QUI_IN_NUM) {
		memset(sval, 0, 64);
		if (*val && (flgs & QUI_VAL_FLGS_RST) == 0) {
			sprintf(sval, "%d", *val);
		}
		int sl = strlen(sval);

		switch (qui_in.rls & QUI_IN_NUM) {
		case QUI_IN_0: sval[sl] = '0'; break;
		case QUI_IN_1: sval[sl] = '1'; break;
		case QUI_IN_2: sval[sl] = '2'; break;
		case QUI_IN_3: sval[sl] = '3'; break;
		case QUI_IN_4: sval[sl] = '4'; break;
		case QUI_IN_5: sval[sl] = '5'; break;
		case QUI_IN_6: sval[sl] = '6'; break;
		case QUI_IN_7: sval[sl] = '7'; break;
		case QUI_IN_8: sval[sl] = '8'; break;
		case QUI_IN_9: sval[sl] = '9'; break;
		};

		if (1 == sscanf(sval, "%d", val)) {
			qui_in.st |= QUI_IN_ST_CNSMD;
			return QUI_VAL_RET_ED;
		}
	}

	if (qui_in.rls & QUI_IN_MINUS) {
		*val *= -1;
		qui_in.st |= QUI_IN_ST_CNSMD;
		return QUI_VAL_RET_ED;
	}
	return QUI_VAL_RET_NIL;
}

char *qui_val_hot_nm;
char *qui_val_hot_unt;
int qui_val_hot_prc;

int qui_val_is_hot(char *nm, char *unt) {
	int r = 0;
	if (nm && !qui_val_hot_nm || !nm && qui_val_hot_nm)
		r = -1;
	if (unt && !qui_val_hot_unt || !unt && qui_val_hot_unt)
		r = -1;
	if (nm && qui_val_hot_nm)
		r |= strcmp(nm, qui_val_hot_nm);
	if (unt && qui_val_hot_unt)
		r |= strcmp(unt, qui_val_hot_unt);
	return !r;
}

void qui_val_f2a(char *dst, float f, char *nm, char *unt) {
	int p = 0, l;

	if (qui_val_is_hot(nm, unt))
		p = qui_val_hot_prc;

	l = sprintf(dst, "%.3f", f);

	if (10.f <= fabs(f)) {
		dst[l-1] = '0';
		p = p < 3 ? p : 2;
	}
	if (100.f <= fabs(f)) {
		dst[l-2] = '0';
		p = p < 2 ? p : 1;
	}
	if (1000.f <= fabs(f)) {
		dst[l-3] = '0';
		p = 0;
	}

	if (dst[l-1] == '0')
		dst[l-1] = '\0';
	else
		return;

	if (dst[l-2] == '0' && p < 3)
		dst[l-2] = '\0';
	else
		return;

	if (dst[l-3] == '0' && p < 2)
		dst[l-3] = '\0';
	else
		return;

	if (p == 0)
		dst[l-4] = '\0';
}

int qui_val_f(float44_t M, float3_t bg_, char *nm, char *unt, float *val, int flgs) {
	char sval[64];

	qui_val_f2a(sval, *val, nm, unt);
	qui_val_drw(M, bg_, nm, unt, sval);

	if (flgs & QUI_VAL_FLGS_CNST)
		return QUI_VAL_RET_NIL;

	if (qui_in.st == QUI_IN_ST_CNSMD)
		return QUI_VAL_RET_NIL;

	qui_tip_sgnl |= QUI_TIP_SGNL_FCS & qui_tip_msk;

	if (qui_in.rls & QUI_IN_RET) {
		qui_val_hot_nm = 0;
		qui_val_hot_unt = 0;

		qui_in.st |= QUI_IN_ST_CNSMD;

		return QUI_VAL_RET_SET;
	}

	if (qui_in.rls & QUI_IN_BCK) {
		int sl = strlen(sval);

		if (sl) {
			sval[--sl] = '\0';

			if (0 == sl) {
				*val = 0;
				return 1;
			}
			*val = atof(sval);

			int prc = 0;
			char *a = strchr(sval, '.');
			if (a) {
				prc = 1;
				if (*++a) {
					++prc;
					if (a && *++a) {
						++prc;
					}
				}
			}
			qui_val_hot_nm = nm;
			qui_val_hot_unt = unt;
			qui_val_hot_prc = prc;


			qui_in.st |= QUI_IN_ST_CNSMD;
			return QUI_VAL_RET_ED;
		}
	}

	if (qui_in.rls & (QUI_IN_NUM | QUI_IN_DOT)) {
		if (flgs & QUI_VAL_FLGS_RST || *val ==  0 && !qui_val_is_hot(nm, unt)) {
			memset(sval, 0, 64);
		}
		int sl = strlen(sval);

		switch (qui_in.rls & (QUI_IN_NUM | QUI_IN_DOT)) {
		case QUI_IN_0: sval[sl++] = '0'; break;
		case QUI_IN_1: sval[sl++] = '1'; break;
		case QUI_IN_2: sval[sl++] = '2'; break;
		case QUI_IN_3: sval[sl++] = '3'; break;
		case QUI_IN_4: sval[sl++] = '4'; break;
		case QUI_IN_5: sval[sl++] = '5'; break;
		case QUI_IN_6: sval[sl++] = '6'; break;
		case QUI_IN_7: sval[sl++] = '7'; break;
		case QUI_IN_8: sval[sl++] = '8'; break;
		case QUI_IN_9: sval[sl++] = '9'; break;
		case QUI_IN_DOT:
			if (0 == sl) {
				sval[sl++] = '0';
			}
			sval[sl++] = '.';
			break;
		};

		sval[sl] = '\0';
		*val = atof(sval);

		int prc = 0;
		char *a = strchr(sval, '.');
		if (a) {
			prc = 1;
			if (*++a) {
				++prc;
				if (a && *++a) {
					++prc;
				}
			}
		}
		qui_val_hot_nm = nm;
		qui_val_hot_unt = unt;
		qui_val_hot_prc = prc;

		qui_in.st |= QUI_IN_ST_CNSMD;
		return QUI_VAL_RET_ED;
	}

	if (qui_in.rls & QUI_IN_MINUS) {
		*val *= -1;
		qui_val_hot_nm = nm;
		qui_val_hot_unt = unt;
		qui_in.st |= QUI_IN_ST_CNSMD;
		return QUI_VAL_RET_ED;
	}
	return QUI_VAL_RET_NIL;
}

int qui_val_g(float44_t M, float3_t bg, char *nm, char *unt, double *val, int flgs) {
	float v = *val;
	int r = qui_val_f(M, bg, nm, unt, &v, flgs);
	*val = v;
	return r;
}
#endif /* QUI_IMPL */
#endif /* QUI_VAL_H */
