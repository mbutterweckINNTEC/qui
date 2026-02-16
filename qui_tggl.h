
int qui_tggl(char *nm, int *flgs, int flg, float44_t M, float3_t clr);

#ifdef QUI_IMPL

static float const qui_tggl_scl = 0.03125;
static float2_t const qui_tggl_mv = {0.035 / 0.0625 * qui_tggl_scl, 0.02 / 0.0625 * qui_tggl_scl };

static float2_t qui_tggl_ngon[] = {
	{ -2.f, 0.f },
	{-1.5f, 0.8660254038f },
	{ 1.5f, 0.8660254038f },
	{ 2.0, 0.0},
	{ 1.5f,-0.8660254038f },
	{-1.5f,-0.8660254038f }
};

static float2_t qui_tggl_ins_ngon[] = {
	{-0.25f, 0.f },
	{ 0.125f, 0.6495190528f },
	{ 1.375f, 0.6495190528f },
	{ 1.75f, 0.0},
	{ 1.375f,-0.6495190528f },
	{ 0.125f,-0.6495190528f },
};

static float2_t qui_tggl_knob_off_ngon[] = {
	{-0.25f, 0.f },
	{ 0.125f, 0.6495190528f },
	{ 0.875f, 0.6495190528f },
	{ 1.25f, 0.0},
	{ 0.875f,-0.6495190528f },
	{ 0.125f,-0.6495190528f },
};

static float2_t qui_tggl_knob_on_ngon[] = {
	{ 0.25f, 0.f },
	{ 0.625f, 0.6495190528f },
	{ 1.375f, 0.6495190528f },
	{ 1.75f, 0.0},
	{ 1.375f,-0.6495190528f },
	{ 0.625f,-0.6495190528f },
};



int qui_tggl(char *nm, int *flgs, int flg, float44_t M, float3_t clr) {
	int ret = 0;
	float4_t bg = m_float4(clr, 1.f);
	float4_t ig = m_float4(mix_float3(clr, (float3_t){ 0.75f, 0.75f, 0.75f }, 0.5), 1.f);
	float4_t fg = m_float4(mix_float3(clr, (float3_t){1,1,1}, 0.75), 1.f);

	float4_t kc_on = m_float4(mix_float3(clr, (float3_t){1,1,1}, 0.75), 1.f);
	float4_t kc_off = m_float4(mix_float3(clr, (float3_t){0.875f,0.875f,0.875f}, 0.75), 1.f);

	float44_t S = {
		qui_tggl_scl, 0, 0, 0,
		0, qui_tggl_scl, 0, 0,
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

		float d = fminf(
			qui_ray_pnt_dst(p, r, (float3_t){ 0.f, 0.f, 0.f }),
			fminf(
				qui_ray_pnt_dst(p, r, (float3_t){ -1.f, 0.f, 0.f }),
				qui_ray_pnt_dst(p, r, (float3_t){ 1.f, 0.f, 0.f })
			)
		);

		if (d < 0.8660254038f) {
			qui_tip_sgnl |= QUI_TIP_SGNL_FCS & qui_tip_msk;

			bg = m_float4(mix_float3(clr, (float3_t){ 0.5f, 0.5f, 0.5f }, 0.5), 1.f);

			if (qui_in.prss & QUI_IN_LMB) {
				bg = m_float4(mix_float3(clr, (float3_t){ 0.25f, 0.25f, 0.25f }, 0.5), 1.f);
			}

			if (qui_in.rls & QUI_IN_LMB) {
				bg = m_float4(mix_float3(clr, (float3_t){ 0.25f, 0.25f, 0.25f }, 0.5), 1.f);
				fg = bg;
				*flgs ^= flg;
				ret = 1;
			}
		}
	}

	/* draw */

	float44_t Z = {
		qui_tggl_scl, 0, 0, 0,
		0, qui_tggl_scl, 0, 0,
		0, 0, 1, 0,
	 	- 2.5 * qui_tggl_mv.x, -qui_tggl_mv.y, 0, 1
	};

	float44_t T = mul_float44(Z, M);

	qui_ngon(6, qui_tggl_ngon, N, bg);
	qui_txt(nm, T, fg);
	qui_ngon(6, qui_tggl_ins_ngon, N, ig);

	if (*flgs & flg)
		qui_ngon(6, qui_tggl_knob_on_ngon, N, kc_on);
	else
		qui_ngon(6, qui_tggl_knob_off_ngon, N, kc_off);

	return ret;
}

#endif /* QUI_IMPL */