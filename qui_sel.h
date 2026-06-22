#ifndef QUI_SEL_H
#define QUI_SEL_H

enum {
	QUI_SEL_RET_NIL,
	QUI_SEL_RET_CHNGD,
};

int qui_sel(char *it[], int n, int *i, float44_t M, float3_t bg);

/* Implementation */
#ifdef QUI_IMPL

static float const qui_sel_scl = 0.03125;

#define QUI_SEL_W 4

static float3_t qui_sel_ngon[] = {
	{ -1.f, 0.f, 0.f },
	{-0.5f, 0.8660254038f, 0.f },
	{ QUI_SEL_W, 0.8660254038f, 0.f },
	{ 0.5f + QUI_SEL_W, 0.f, 0.f },
	{ QUI_SEL_W,-0.8660254038f, 0.f },
	{-0.5f,-0.8660254038f, 0.f }
};

static float3_t qui_sel_larrw_ngon[] = {
	{ -1.f, 0.f, 0.f },
	{-0.5f, 0.8660254038f, 0.f },
	{ 0.f, 0.8660254038f, 0.f },
	{-0.5f, 0.f, 0.f },
	{ 0.f,-0.8660254038f, 0.f },
	{-0.5f,-0.8660254038f, 0.f }
};

static float3_t qui_sel_rarrw_ngon[] = {
	{-0.5f + QUI_SEL_W, 0.8660254038f, 0.f },
	{ QUI_SEL_W, 0.8660254038f, 0.f },
	{ 0.5f + QUI_SEL_W, 0.f, 0.f },
	{ QUI_SEL_W,-0.8660254038f, 0.f },
	{-0.5f + QUI_SEL_W,-0.8660254038f, 0.f },
	{ QUI_SEL_W, 0.f, 0.f }
};

void qui_sel_drw(char *it, float44_t M_, float4_t bg, float4_t lbc, float4_t rbc, float4_t tc) {
	char txt[32] = "";
	int len;

	strncpy(txt, it, 31);
	len = strlen(txt);


	float44_t S = {
		qui_sel_scl, 0, 0, 0,
		0, qui_sel_scl, 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1
	};
	float44_t M = mul_float44(S, M_);
	float44_t T = mul_float44(
		(float44_t) {
			1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 1, 0,
			0, -0.5f, 0, 1
		}, M
	);

	while (len && QUI_SEL_W + 1 < qui_txt_len(txt))
		txt[--len] = '\0';

	qui_ngon(6, qui_sel_ngon, M, bg);
	qui_txt(txt, T, tc);
	qui_ngon(6, qui_sel_larrw_ngon, M, lbc);
	qui_ngon(6, qui_sel_rarrw_ngon, M, rbc);
}

int qui_sel_hvr(float44_t M_) {
	float44_t P = qui_mtrx_top(QUI_MTRX_P);
	float44_t V = qui_mtrx_top(QUI_MTRX_V);
	float44_t S = {
		qui_sel_scl, 0, 0, 0,
		0, qui_sel_scl, 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1
	};
	float44_t M = mul_float44(S, M_);
	float44_t PVM = mul_float44(M, mul_float44(V, P));
	float detPVM = det_float44(PVM);

	if (detPVM) {
		float44_t iPVM = invert_float44(PVM, detPVM);

		float3_t p = float3_float4(cotransform_float44(iPVM, (float4_t){ qui_in.p.x, qui_in.p.y, 0.f, 1.f }));
		float3_t r = normal_float3(m_float3(cotransform_float44(iPVM, (float4_t){ 0.f, 0.f, 1.f, 0.f })));
		float d = qui_ray_pnt_dst(p, r, (float3_t){ 0.f, 0.f, 0.f });

		if (d < 0.8660254038f) {
			qui_tip_sgnl |= QUI_TIP_SGNL_FCS;
			qui_in.st = QUI_IN_ST_HVRD;

			return 1;
		}

		r = normal_float3(m_float3(cotransform_float44(iPVM, (float4_t){ 0.f, 0.f, 1.f, 0.f })));
		d = qui_ray_pnt_dst(p, r, (float3_t){ QUI_SEL_W - 0.5f, 0.f, 0.f });

		if (d < 0.8660254038f) {
			qui_tip_sgnl |= QUI_TIP_SGNL_FCS;
			qui_in.st = QUI_IN_ST_HVRD;

			return 2;
		}
	}

	return 0;
}

int qui_sel(char *it[], int n, int *i, float44_t M, float3_t clr) {
	float4_t bg = m_float4(clr, 1.f);
	float4_t bch = m_float4(mix_float3(clr, (float3_t){1,1,1}, 0.75), 1.f);
	float4_t bcm = m_float4(mix_float3(clr, (float3_t){1,1,1}, 0.5), 1.f);
	float4_t bcl = m_float4(mix_float3(clr, (float3_t){1,1,1}, 0.25), 1.f);
	float4_t tc = m_float4(mix_float3(clr, (float3_t){0,0,0}, 0.75), 1.f);
	float4_t lbc = bcl, rbc = bcl, hbc = bcm;

	int hvr = qui_sel_hvr(M);

	if (qui_in.prss & QUI_IN_LMB)
		hbc = bch;
	if (hvr == 1)
		lbc = hbc;
	if (hvr == 2)
		rbc = hbc;

	qui_sel_drw(it[*i % n], M, bg, lbc, rbc, tc);

	if (qui_in.rls & QUI_IN_LMB) {
		if (hvr == 1) {
			*i = *i ? *i - 1 : n - 1;
			return QUI_SEL_RET_CHNGD;
		}

		if (hvr == 2) {
			*i = (*i + 1) % n;
			return QUI_SEL_RET_CHNGD;
		}
	}

	return QUI_SEL_RET_NIL;
}

#endif /* QUI_IMPL */
#endif /* QUI_SEL_H */
