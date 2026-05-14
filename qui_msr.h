#ifndef QUI_MSR_H
#define QUI_MSR_H

int qui_msr2(float3_t pa, float3_t pb, float3_t n, float mnr, float myr, float4_t c);

#ifdef QUI_IMPL

int qui_msr2(float3_t pa, float3_t pb, float3_t n, float mnr, float myr, float4_t c) {
	float3_t d = sub_float3(pb, pa);
	float l = length_float3(d);
	float3_t t = normal_float3(cross_float3(n, d));
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
		*p++ = (float3_t) { i * mnr, 0.0625 };
	}

	for (int i = 0; i <= myrn; ++i) {
		*p++ = (float3_t) { i * myr };
		*p++ = (float3_t) { i * myr, 0.125 };
	}

	float3_t d_ = normal_float3(d);
	float3_t rn = normal_float3(cross_float3(t, d_));
	float44_t B = {
		d_.x, d_.y, d_.z, 0.f,
		t.x, t.y, t.z, 0.f,
		rn.x, rn.y, rn.z, 0.f,
		0.f, 0.f, 0.f, 1.f
	};

	qui_lns(pn, p_, B, c);

	for (int i = 0; i <= myrn; i += tic) {
		char a[32];
		sprintf(a, "%.1f", myr * i);
		float44_t T = {
			0.0625f, 0.f, 0.f, 0.f,
			0.f, 0.0625f, 0.f, 0.f,
			0.f, 0.f, 0.0625f, 0.f,
			myr * i-0.5*0.0625*qui_txt_len(a), 0.125f, 0.f, 1.f
		};
		qui_txt(a, mul_float44(T, B), c);
	}

	return 0;
}

#endif /* QUI_IMPL */
#endif /* QUI_MSR_H */
