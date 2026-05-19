#ifndef QUI_MSR_H
#define QUI_MSR_H

enum {
	QUI_MSR_MNR,	/* minor tic length: double */
	QUI_MSR_MYR,	/* mayor tic length: double */
	QUI_MSR_CLR,	/* color: float4_t */
	QUI_MSR_UNT,	/* unit name: char */
	QUI_MSR_SNP,	/* snap distance: double */
	QUI_MSR_ARRWHD,	/* arrow head length: double*/
	QUI_MSR_FNTH,	/* arrow head length: double*/
};

int qui_msr_set(int param, ...);

int qui_msr(float3_t r, float3_t o);

#ifdef QUI_IMPL

#include <stdarg.h>

enum {
	QUI_MSR_ST_IDL,
	QUI_MSR_ST_RDY,
	QUI_MSR_ST_LN,
	QUI_MSR_ST_ANG,
};

#define QUI_MSR_N 8

int qui_msr_st;

float3_t qui_msr_p[QUI_MSR_N][3];
int qui_msr_n;
int qui_msr_pn[QUI_MSR_N];

float qui_msr_mnr = 1.f, qui_msr_myr = 10.f, qui_msr_snp = 0.0625f, qui_msr_arrwhd = 2.5, qui_msr_fnth = 3.5;
float4_t qui_msr_clr = { 1.f, 1.f, 1.f, 1.f };
char *qui_msr_unt = "mm";

void qui_msr_flp(float3_t *a, float3_t *b) {
	float44_t P = qui_mtrx_top(QUI_MTRX_P);
	float44_t V = qui_mtrx_top(QUI_MTRX_V);
	float44_t PV = mul_float44(V, P);
	float3_t a_ = float3_float4(cotransform_float44(PV, m_float4(*a, 1.f)));
	float3_t b_ = float3_float4(cotransform_float44(PV, m_float4(*b, 1.f)));
	float3_t t;
	float dx = (b_.x - a_.x);
	float dy = (b_.y - a_.y);

	if (fabs(dy) < fabs(dx)) {
		if (dx < 0) {
			t = *a;
			*a = *b;
			*b = t;
		}
	} else {
		if (a_.x < 0.f || b_.x < 0.f) {
			if (0 > dy) {
				t = *a;
				*a = *b;
				*b = t;
			}
		} else {
			if (0 < dy) {
				t = *a;
				*a = *b;
				*b = t;
			}

		}
	}
}

int qui_msr_drw_ln(float3_t p_[2], float3_t n) {
	float3_t p[2];
	memcpy(p, p_, sizeof(p));
	qui_msr_flp(p, p+1);

	qui_lns(2, p, 1, identity_sc, qui_msr_clr);

	float3_t l = sub_float3(p[1], p[0]);
	float3_t s = normal_float3(l);
	float3_t t = normal_float3(cross_float3(s, n));
	s = normal_float3(cross_float3(n, t));
	s = scale_float3(s, qui_msr_arrwhd);
	t = scale_float3(t, qui_msr_arrwhd / 3.f);
	float3_t aha[3] = { p[0], add_float3(p[0], sub_float3(s, t)), add_float3(p[0], add_float3(s, t)) };

	qui_ngon(3, aha, identity_sc, qui_msr_clr);

	float3_t ahb[3] = { p[1], sub_float3(p[1], sub_float3(s, t)), sub_float3(p[1], add_float3(s, t)) };
	qui_ngon(3, ahb, identity_sc, qui_msr_clr);

	float3_t c = scale_float3(add_float3(p[0], p[1]), 0.5f);

	char a[32];
	sprintf(a, "%.3f %s", length_float3(l), qui_msr_unt);
	s = scale_float3(normal_float3(s), qui_msr_fnth);
	t = scale_float3(normal_float3(t), -qui_msr_fnth);
	n = scale_float3(n, qui_msr_fnth);
	c = add_float3(c, scale_float3(normal_float3(s), -0.5f * qui_msr_fnth * qui_txt_len(a)));
	c = add_float3(c, scale_float3(normal_float3(t), qui_msr_fnth));

	float44_t T = {
		s.x, s.y, s.z, 0.f,
		t.x, t.y, t.z, 0.f,
		n.x, n.y, n.z, 0.f,
		c.x, c.y, c.z, 1.f
	};

	qui_txt(a, T, qui_msr_clr);

	return 0;
}

int qui_msr_drw_rlr(float3_t pa, float3_t pb, float3_t n) {
	float3_t d = sub_float3(pb, pa);
	float l = length_float3(d);
	float3_t t = normal_float3(cross_float3(n, d));
	float mnr = qui_msr_mnr;
	float myr = qui_msr_myr;
	int mnrn = (int)(l / mnr), myrn = (int)(l / myr) + 1, pn = (mnrn + myrn + 1) << 1, off = 0;
	float3_t *p = alloca(pn * sizeof(float3_t));
	float3_t *p_ = p;
	float3_t s;
	int tic = 5;

	if (NULL == p)
		return -1;

	*p++ = (float3_t){ 0.f, 0.f };
	*p++ = (float3_t){ l, 0.f };

	for (int i = 1; i <= mnrn; ++i) {
		*p++ = (float3_t) { i * mnr };
		*p++ = (float3_t) { i * mnr, qui_msr_arrwhd };
	}

	for (int i = 0; i <= myrn; ++i) {
		*p++ = (float3_t) { i * myr };
		*p++ = (float3_t) { i * myr, qui_msr_arrwhd * 2.f };
	}

	float3_t d_ = normal_float3(d);
	float3_t rn = normal_float3(cross_float3(t, d_));
	float44_t B = {
		d_.x, d_.y, d_.z, 0.f,
		t.x, t.y, t.z, 0.f,
		rn.x, rn.y, rn.z, 0.f,
		pa.x, pa.y, pa.z, 1.f
	};

	qui_lns(pn, p_, 1, B, qui_msr_clr);

	for (int i = 0; i <= myrn; i += tic) {
		char a[32];
		if (myrn < i + tic)
			sprintf(a, "%.1f %s", myr * i, qui_msr_unt);
		else
			sprintf(a, "%.1f", myr * i);
		float44_t T = {
			qui_msr_fnth, 0.f, 0.f, 0.f,
			0.f, qui_msr_fnth, 0.f, 0.f,
			0.f, 0.f, qui_msr_fnth, 0.f,
			myr * i-0.5*qui_msr_fnth*qui_txt_len(a), 2.0 * qui_msr_arrwhd, 0.f, 1.f
		};
		qui_txt(a, mul_float44(T, B), qui_msr_clr);
	}

	return 0;
}

int qui_msr_drw_ang(float3_t p[3]) {
	qui_plln(3, p, 1, identity_sc, qui_msr_clr);

	if (!equal_float3(p[1], p[2])) {
		float3_t u_ = normal_float3(sub_float3(p[0], p[1]));
		float3_t v_ = normal_float3(sub_float3(p[2], p[1]));
		float cang = dot_float3(u_, v_);
		float R = (2.f - sqrt(0.5f + 0.5f * cang));
		float3_t u = scale_float3(u_, 2.f * qui_msr_fnth);
		float3_t v = scale_float3(v_, 2.f * qui_msr_fnth);
		float3_t r = scale_float3(normal_float3(add_float3(u, v)), 2.f * R * qui_msr_fnth);
		float3_t c[3] = {
			add_float3(p[1], u),
			add_float3(p[1], r),
			add_float3(p[1], v),
		};

		qui_bzr(3, c, 1, identity_sc, qui_msr_clr);
		/* todo: value */
	}

	return 0;
}

int qui_msr_drw_prtrctr(float3_t p[]) {
	qui_plln(3, p, 1, identity_sc, qui_msr_clr);

	if (!equal_float3(p[1], p[2])) {
		float3_t u_ = normal_float3(sub_float3(p[0], p[1]));
		float3_t v_ = normal_float3(sub_float3(p[2], p[1]));
		float cang = dot_float3(u_, v_);
		float R = (2.f - sqrt(0.5f + 0.5f * cang));
		float3_t u = scale_float3(u_, 2.f * qui_msr_fnth);
		float3_t v = scale_float3(v_, 2.f * qui_msr_fnth);
		float3_t r = scale_float3(normal_float3(add_float3(u, v)), 2.f * R * qui_msr_fnth);
		float3_t c[3] = {
			add_float3(p[1], u),
			add_float3(p[1], r),
			add_float3(p[1], v),
		};

		qui_bzr(3, c, 1, identity_sc, qui_msr_clr);
		/* todo: scale */
	}

	return 0;
}

int qui_msr_set(int param, ...) {
	va_list ap;
	int err = 0;

	va_start(ap, param);

	switch(param) {
	case QUI_MSR_MNR:	qui_msr_mnr = va_arg(ap, double);	break;
	case QUI_MSR_MYR:	qui_msr_myr = va_arg(ap, double);	break;
	case QUI_MSR_CLR:	qui_msr_clr = va_arg(ap, float4_t);	break;
	case QUI_MSR_UNT:	qui_msr_unt = va_arg(ap, char*);	break;
	case QUI_MSR_SNP:	qui_msr_snp = va_arg(ap, double);	break;
	case QUI_MSR_ARRWHD:	qui_msr_arrwhd = va_arg(ap, double);	break;
	case QUI_MSR_FNTH:	qui_msr_fnth = va_arg(ap, double);	break;
	default:
		err = -1;
	};

	va_end(ap);

	return err;
}

int qui_msr_drw_x(float3_t o) {
	float44_t P = qui_mtrx_top(QUI_MTRX_P);
	float44_t V = qui_mtrx_top(QUI_MTRX_V);
	float44_t PV = mul_float44(V, P);
	float3_t x = float3_float4(cotransform_float44(PV, m_float4(o, 1.f)));
	float3_t p[4] = {
		{ x.x - qui_msr_fnth, x.y, x.z }, { x.x + qui_msr_fnth, x.y, x.z },
		{ x.x, x.y - qui_msr_fnth, x.y }, { x.x, x.y + qui_msr_fnth, x.z },
	};

	qui_mtrx_psh(QUI_MTRX_P, identity_sc);
	qui_mtrx_psh(QUI_MTRX_V, identity_sc);
	qui_lns(4, p, 1, identity_sc, qui_msr_clr);
	qui_mtrx_pop(QUI_MTRX_V);
	qui_mtrx_pop(QUI_MTRX_P);
	return 0;
}

float3_t qui_msr_qr(float3_t o) {
	float3_t p = { qui_in.p.x, qui_in.p.y, 0.f };
	int qrt, qri;

	if (0 > (qri = qui_qr_find(QUI_QR_FND_BST, o, qui_msr_snp, &qrt)))
		goto end;

	p = qui_qr[qrt][qri + qui_qr_s[qrt] - 1];
end:
	return p;
}

int qui_msr(float3_t r, float3_t o) {
	switch(qui_msr_st) {
	case QUI_MSR_ST_IDL:
		if (qui_in.rls & QUI_IN_M) {
			qui_msr_st = QUI_MSR_ST_RDY;
			break;
		}
		return 0;
	case QUI_MSR_ST_RDY:
		qui_msr_drw_x(o);
		if (qui_msr_n < QUI_MSR_N) {
			if (qui_in.rls & QUI_IN_LMB) {
				qui_msr_pn[qui_msr_n] = 1;
				qui_msr_p[qui_msr_n][0] = qui_msr_qr(o);
				qui_msr_st = QUI_MSR_ST_LN;
			}
			if (qui_in.rls & QUI_IN_RMB) {
				qui_msr_pn[qui_msr_n] = 1;
				qui_msr_p[qui_msr_n][0] = qui_msr_qr(o);
				qui_msr_st = QUI_MSR_ST_ANG;
			}
		}
		if (qui_in.rls & QUI_IN_BCK && qui_msr_n)
			--qui_msr_n;
		if (qui_in.rls & QUI_IN_ESC)
			qui_msr_st = QUI_MSR_ST_IDL;
		break;
	case QUI_MSR_ST_LN:
		if (qui_in.rls & QUI_IN_LMB) {
			qui_msr_p[qui_msr_n][qui_msr_pn[qui_msr_n]++] = qui_msr_qr(o);
			qui_msr_n++;
			qui_msr_st = QUI_MSR_ST_RDY;
		}
		if (qui_in.rls & QUI_IN_BCK && qui_msr_pn[qui_msr_n])
			--qui_msr_pn[qui_msr_n];
		if (qui_in.rls & QUI_IN_ESC)
			qui_msr_st = QUI_MSR_ST_RDY;
		break;
	case QUI_MSR_ST_ANG:
		if (qui_in.rls & QUI_IN_RMB) {
			qui_msr_p[qui_msr_n][qui_msr_pn[qui_msr_n]++] = qui_msr_qr(o);

			if (qui_msr_pn[qui_msr_n] == 3) {
				++qui_msr_n;
				qui_msr_st = QUI_MSR_ST_RDY;
			}
		}
		if (qui_in.rls & QUI_IN_BCK && qui_msr_pn[qui_msr_n])
			--qui_msr_pn[qui_msr_n];
		if (qui_in.rls & QUI_IN_ESC)
			qui_msr_st = QUI_MSR_ST_RDY;
		break;
	};

	float3_t n = normal_float3(r);//qui_msr_vwn();

	for (int i = 0; i < qui_msr_n; ++i) {
		switch (qui_msr_pn[i]) {
		case 2:
			qui_msr_drw_ln(qui_msr_p[i], n);
			break; 
		case 3:
			qui_msr_drw_ang(qui_msr_p[i]);
			break;
		};
	}

	if (qui_msr_st == QUI_MSR_ST_LN && qui_msr_pn[qui_msr_n]) {
		float3_t pcrsr = qui_msr_qr(o);
		qui_msr_drw_rlr(qui_msr_p[qui_msr_n][0], pcrsr, n);
	}

	if (qui_msr_st == QUI_MSR_ST_ANG && qui_msr_pn[qui_msr_n] == 1) {
		float3_t p[3] = { qui_msr_p[qui_msr_n][0], qui_msr_qr(o) };
		p[2] = p[1];
		qui_msr_drw_prtrctr(p);
	}

	if (qui_msr_st == QUI_MSR_ST_ANG && qui_msr_pn[qui_msr_n] == 2) {
		float3_t p[3] = { qui_msr_p[qui_msr_n][0], qui_msr_p[qui_msr_n][1], qui_msr_qr(o) };
		qui_msr_drw_prtrctr(p);
	}

	return 1;
}

#endif /* QUI_IMPL */
#endif /* QUI_MSR_H */
