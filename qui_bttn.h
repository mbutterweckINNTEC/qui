
enum {
	QUI_BTTN_NIL,	/* always zero */
	QUI_BTTN_L,		/* left click */
	QUI_BTTN_R		/* right click */
};

int qui_bttn(char *nm, float44_t M, float3_t clr);

#ifdef QUI_IMPL

static float const qui_bttn_scl = 0.03125;
static float2_t const qui_bttn_mv = {0.035 / 0.0625 * qui_bttn_scl, 0.02 / 0.0625 * qui_bttn_scl };

static float2_t qui_bttn_ngon[] = {
	{ -1.f, 0.f },
	{-0.5f, 0.8660254038f },
	{ 0.5f, 0.8660254038f },
	{ 1.0, 0.0},
	{ 0.5f,-0.8660254038f },
	{-0.5f,-0.8660254038f },
};

int qui_bttn(char *nm, float44_t M, float3_t clr) {
	int ret = QUI_BTTN_NIL;
	float4_t bg = m_float4(clr, 1.f);
	float4_t fg = m_float4(mix_float3(clr, (float3_t){1,1,1}, 0.75), 1.f);

	float44_t S = {
		qui_bttn_scl, 0, 0, 0,
		0, qui_bttn_scl, 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1
	};
	float44_t N = mul_float44(S, M);

	/* input */
	float44_t P = qui_mtrx_top(QUI_MTRX_P);
	float44_t V = qui_mtrx_top(QUI_MTRX_V);
	float44_t PVN = mul_float44(N, mul_float44(V, P));
	float detPVN = det_float44(PVN);

	if (detPVN) {
		float44_t iPVN = invert_float44(PVN, detPVN);

		float3_t p = float3_float4(cotransform_float44(iPVN, (float4_t){ qui_in.p.x, qui_in.p.y, 0.f, 1.f }));
		float3_t r = normal_float3(m_float3(cotransform_float44(iPVN, (float4_t){ 0.f, 0.f, 1.f, 0.f })));

		float d = qui_ray_pnt_dst(p, r, (float3_t){ 0.f, 0.f, 0.f });

		if (d < 0.8660254038f) {
			bg = m_float4(mix_float3(clr, (float3_t){ 0.5f, 0.5f, 0.5f }, 0.5), 1.f);

			if (qui_in.prss & QUI_IN_LMB) {
				bg = m_float4(mix_float3(clr, (float3_t){ 0.25f, 0.25f, 0.25f }, 0.5), 1.f);
			}

			if (qui_in.prss & QUI_IN_RMB) {
				bg = fg;
				fg = m_float4(mix_float3(clr, (float3_t){ 0.25f, 0.25f, 0.25f }, 0.5), 1.f);
			}

			if (qui_in.rls & QUI_IN_LMB) {
				bg = m_float4(mix_float3(clr, (float3_t){ 0.25f, 0.25f, 0.25f }, 0.5), 1.f);
				fg = bg;
				ret = QUI_BTTN_L;
			}

			if (qui_in.rls & QUI_IN_RMB) {
				bg = fg;
				ret = QUI_BTTN_R;
			}
		}
	}

	/* draw */

	float44_t Z = {
		qui_bttn_scl, 0, 0, 0,
		0, qui_bttn_scl, 0, 0,
		0, 0, 1, 0,
		-qui_bttn_mv.x, -qui_bttn_mv.y, 0, 1
	};

	float44_t T = mul_float44(Z, M);

	qui_ngon(6, qui_bttn_ngon, N, bg);
	qui_txt(nm, T, fg);

	return ret;
}

#endif /* QUI_IMPL */