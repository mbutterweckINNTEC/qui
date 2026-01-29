
static float3_t qui_rnbw(float t) {
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

static inline float qui_ray_pnt_dst(float3_t ro, float3_t rd, float3_t p) {
	float t = dot_float3(rd, sub_float3(p, ro)) / dot_float3(rd, rd);
	return length_float3(sub_float3(p, add_float3(ro, scale_float3(rd, t))));
}

/* todo:@michal: make analytic function? */

static float qui_ray_xcrcl_(float3_t p, float3_t d, int *out) {
	float3_t v;
	float l = FLT_MAX, vr, phi;
	for (int i = 0; i < QUI_MAN_CRCL_X_N; ++i) {
		phi = (float)i / (float)(QUI_MAN_CRCL_X_N - 1) * 2.0 * M_PI;
		v = (float3_t) { 0.f, QUI_MAN_R_XYZ * cos(phi), QUI_MAN_R_XYZ * sin(phi) };
		vr = qui_ray_pnt_dst(p, d, v);

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
		vr = qui_ray_pnt_dst(p, d, v);

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
		vr = qui_ray_pnt_dst(p, d, v);

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

	return l;
}

static float qui_ray_crnr_(float2_t p) {
	static float2_t const P[12] = {
		{  QUI_MAN_S_DXY,  QUI_MAN_S_DXY },
		{  QUI_MAN_S_XY,   QUI_MAN_S_DXY },
		{  QUI_MAN_S_DXY,  QUI_MAN_S_XY  },
		{ -QUI_MAN_S_DXY,  QUI_MAN_S_DXY },
		{ -QUI_MAN_S_XY,   QUI_MAN_S_DXY },
		{ -QUI_MAN_S_DXY,  QUI_MAN_S_XY  },
		{  QUI_MAN_S_DXY, -QUI_MAN_S_DXY },
		{  QUI_MAN_S_XY,  -QUI_MAN_S_DXY },
		{  QUI_MAN_S_DXY, -QUI_MAN_S_XY  },
		{ -QUI_MAN_S_DXY, -QUI_MAN_S_DXY },
		{ -QUI_MAN_S_XY,  -QUI_MAN_S_DXY },
		{ -QUI_MAN_S_DXY, -QUI_MAN_S_XY  }
	};

	float l = FLT_MAX, l_;

	for (int i = 0; i < 12; ++i) {
		l_ = length_float2(sub_float2(P[i], p));

		if (l_ < l) {
			l = l_;
		}
	}

	return l;
}

static float3_t qui_ray_ln_near(float3_t ro, float3_t rd, float3_t lo, float3_t ld) {
	float3_t n = cross_float3(rd, cross_float3(ld, rd));

	return add_float3(lo, scale_float3(ld, dot_float3(sub_float3(ro, lo), n) / dot_float3(ld, n)));
}

static float qui_ray_seg_dst(float3_t ro, float3_t rd, float3_t a, float3_t b) {
	rd = normal_float3(rd);
	float3_t abd = sub_float3(b, a);
	float lab = length_float3(abd);

	float3_t c = qui_ray_ln_near(ro, rd, a, abd);

	if (length_float3(sub_float3(c, a)) < lab && length_float3(sub_float3(c, b)) < lab) {
		float3_t n = normal_float3(cross_float3(abd, rd));
		return fabs(dot_float3(n, sub_float3(ro, a)));
	}

	return fmin(
		length_float3(sub_float3(a, add_float3(ro, scale_float3(rd, dot_float3(rd, sub_float3(a, ro)))))),
		length_float3(sub_float3(b, add_float3(ro, scale_float3(rd, dot_float3(rd, sub_float3(b, ro))))))
	);
}
