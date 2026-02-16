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

/* Implementation */

static float const qui_val_scl = 0.03125;
static float2_t const qui_val_mv = {0.035 / 0.0625 * qui_val_scl, 0.02 / 0.0625 * qui_val_scl };

static float2_t qui_val_nm_ngon[] = {
	{ -1.f, 0.f },
	{-0.5f, 0.8660254038f },
	{ 0.5f, 0.8660254038f },
	{ 1.0, 0.0},
	{ 0.5f,-0.8660254038f },
	{-0.5f,-0.8660254038f },
};

static float2_t qui_val_unt_ngon[] = {
	{ -1.f + 6.f + 2 * 0.1f, 0.f           - 2 * 0.8660254038f - 2 * 0.08660254038f },
	{-0.5f + 6.f +2 *  0.1f, 0.8660254038f - 2 * 0.8660254038f - 2 * 0.08660254038f },
	{ 0.5f + 6.f + 2 * 0.1f, 0.8660254038f - 2 * 0.8660254038f - 2 * 0.08660254038f },
	{ 1.0f + 6.f + 2 * 0.1f, 0.0f          - 2 * 0.8660254038f - 2 * 0.08660254038f },
	{ 0.5f + 6.f + 2 * 0.1f,-0.8660254038f - 2 * 0.8660254038f - 2 * 0.08660254038f },
	{-0.5f + 6.f + 2 * 0.1f,-0.8660254038f - 2 * 0.8660254038f - 2 * 0.08660254038f },
};

static float2_t qui_val_val_ngon[] = {
	{ -1.f + 1.5 + 0.1f, 0.f           - 0.8660254038f - 0.08660254038f },
	{-0.5f + 1.5 + 0.1f, 0.8660254038f - 0.8660254038f - 0.08660254038f },
	{ 3.5f + 1.5 + 0.1f, 0.8660254038f - 0.8660254038f - 0.08660254038f },
	{ 4.0f + 1.5 + 0.1f, 0.0f          - 0.8660254038f - 0.08660254038f },
	{ 3.5f + 1.5 + 0.1f,-0.8660254038f - 0.8660254038f - 0.08660254038f },
	{-0.5f + 1.5 + 0.1f,-0.8660254038f - 0.8660254038f - 0.08660254038f },
};

static int qui_val_drw(float44_t M_, float3_t clr, char *nm, char *unt, char *val) {
	float4_t bg = m_float4(clr, 1.f);
	float4_t fg = m_float4(mix_float3(clr, (float3_t){1,1,1}, 0.75), 1.f);
	float4_t fgv = m_float4(mix_float3(clr, (float3_t){0,0,0}, 0.75), 1.f);

	float44_t S = {
		qui_val_scl, 0, 0, 0,
		0, qui_val_scl, 0, 0,
		0, 0, 1, 0,
		-qui_val_scl*5, 0, 0, 1
	};

	float44_t Z = {
		qui_val_scl, 0, 0, 0,
		0, qui_val_scl, 0, 0,
		0, 0, 1, 0,
		-qui_val_mv.x-qui_val_scl*5, -qui_val_mv.y, 0, 1
	};

	float44_t M = mul_float44(S, M_);

	float44_t T = mul_float44(Z, M_);
		
	float44_t T2 = mul_float44(
		(float44_t) {
			1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 1, 0,
			2.0+0.1, -1, 0, 1
		},
		T
	);

	float44_t T3 = mul_float44(
		(float44_t) {
			1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 1, 0,
			6.2, -2, 0, 1
		},
		T
	);

	qui_ngon(6, qui_val_nm_ngon, M, bg);
	qui_txt(nm, T, fg);
	qui_ngon(6, qui_val_val_ngon, M, bg);
	qui_txt(val, T2, fgv);
	qui_ngon(6, qui_val_unt_ngon, M, bg);
	qui_txt(unt, T3, fg);

	return 0;
}

int qui_val_i(float44_t M, float3_t bg_, char *nm, char *unt, int *val, int flgs) {
	char sval[64];

	sprintf(sval, "%4d", *val);
	qui_val_drw(M, bg_, nm, unt, sval);

	if (flgs & QUI_VAL_FLGS_CNST)
		return QUI_VAL_RET_NIL;

	qui_tip_sgnl |= QUI_TIP_SGNL_FCS & qui_tip_msk;

	if (qui_in.rls & QUI_IN_RET)
		return QUI_VAL_RET_SET;

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

		if (1 == sscanf(sval, "%d", val))
			return QUI_VAL_RET_ED;
	}

	if (qui_in.rls & QUI_IN_MINUS) {
		*val *= -1;
		return QUI_VAL_RET_ED;
	}
	return QUI_VAL_RET_NIL;
}

void qui_val_f2a(char *dst, float f) {
	int p, l;
	union {
		float f;
		unsigned b;
	} u;

	u.f = f;
	p = u.b & 3u;
	u.b &=~ 3u;

	l = sprintf(dst, "%.3f", u.f);

	if (10.f <= fabs(u.f)) {
		dst[l-1] = '0';
		p = p < 3 ? p : 2;
	}
	if (100.f <= fabs(u.f)) {
		dst[l-2] = '0';
		p = p < 2 ? p : 1;
	}
	if (1000.f <= fabs(u.f)) {
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

	qui_val_f2a(sval, *val);
	qui_val_drw(M, bg_, nm, unt, sval);

	if (flgs & QUI_VAL_FLGS_CNST)
		return QUI_VAL_RET_NIL;

	qui_tip_sgnl |= QUI_TIP_SGNL_FCS & qui_tip_msk;

	if (qui_in.rls & QUI_IN_RET) {
		union {
			float f;
			unsigned b;
		} u;
		u.f = *val;
		u.b &=~ 3;
		*val = u.f;

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

			union {
				float f;
				unsigned b;
			} u;

			sscanf(sval, "%f", &u.f);

			u.b &=~ 3;
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
			u.b |= prc;

			*val = u.f;
			return QUI_VAL_RET_ED;
		}
	}

	if (qui_in.rls & (QUI_IN_NUM | QUI_IN_DOT)) {
		if (flgs & QUI_VAL_FLGS_RST || *val ==  0) {
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

		union {
			float f;
			unsigned b;
		} u;

		sscanf(sval, "%f", &u.f);

		u.b &=~ 3;
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
		u.b |= prc;

		*val = u.f;

		return QUI_VAL_RET_ED;
	}

	if (qui_in.rls & QUI_IN_MINUS) {
		*val *= -1;
		return QUI_VAL_RET_ED;
	}
	return QUI_VAL_RET_NIL;
}