static float2_t qui_val_nm_ngon[] = {
	{ -1.f, 0.f },
	{-0.5f, 0.8660254038f },
	{ 0.5f, 0.8660254038f },
	{ 1.0, 0.0},
	{ 0.5f,-0.8660254038f },
	{-0.5f,-0.8660254038f },
};

static float2_t qui_val_unt_ngon[] = {
	{ -1.f + 5.f + 2 * 0.1f, 0.f           - 2 * 0.8660254038f - 2 * 0.08660254038f },
	{-0.5f + 5.f +2 *  0.1f, 0.8660254038f - 2 * 0.8660254038f - 2 * 0.08660254038f },
	{ 0.5f + 5.f + 2 * 0.1f, 0.8660254038f - 2 * 0.8660254038f - 2 * 0.08660254038f },
	{ 1.0f + 5.f + 2 * 0.1f, 0.0f          - 2 * 0.8660254038f - 2 * 0.08660254038f },
	{ 0.5f + 5.f + 2 * 0.1f,-0.8660254038f - 2 * 0.8660254038f - 2 * 0.08660254038f },
	{-0.5f + 5.f + 2 * 0.1f,-0.8660254038f - 2 * 0.8660254038f - 2 * 0.08660254038f },
};

static float2_t qui_val_val_ngon[] = {
	{ -1.f + 1.5 + 0.1f, 0.f           - 0.8660254038f - 0.08660254038f },
	{-0.5f + 1.5 + 0.1f, 0.8660254038f - 0.8660254038f - 0.08660254038f },
	{ 2.5f + 1.5 + 0.1f, 0.8660254038f - 0.8660254038f - 0.08660254038f },
	{ 3.0f + 1.5 + 0.1f, 0.0f          - 0.8660254038f - 0.08660254038f },
	{ 2.5f + 1.5 + 0.1f,-0.8660254038f - 0.8660254038f - 0.08660254038f },
	{-0.5f + 1.5 + 0.1f,-0.8660254038f - 0.8660254038f - 0.08660254038f },
};

static int qui_val_i(struct qui_ctx *qc, float44_t PV, float2_t p, float3_t bg_, char *nm, char *unt, int *i) {
	float4_t bg = m_float4(bg_, 1.f);
	float4_t fg = m_float4(mix_float3(bg_, (float3_t){1,1,1}, 0.75), 1.f);
	float4_t fgv = m_float4(mix_float3(bg_, (float3_t){0,0,0}, 0.75), 1.f);

	float uisc = 0.03125;
	float2_t uit = {0.035 / 0.0625 * uisc, 0.02 / 0.0625 * uisc };
	float44_t PWM = mul_float44(
		(float44_t) {
			uisc, 0, 0, 0,
			0, uisc, 0, 0,
			0, 0, 1, 0,
			p.x-uisc*5, p.y, 0, 1
		},
		PV
	);
	float44_t PWT = mul_float44(
		(float44_t) {
			uisc, 0, 0, 0,
			0, uisc, 0, 0,
			0, 0, 1, 0,
			p.x-uit.x-uisc*5,p.y-uit.y, 0, 1
		},
		PV
	);
		
	float44_t PWT2 = mul_float44(
		(float44_t) {
			1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 1, 0,
			2.0+0.1, -1, 0, 1
		},
		PWT
	);

	float44_t PWT3 = mul_float44(
		(float44_t) {
			1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 1, 0,
			5.5+0.1, -2, 0, 1
		},
		PWT
	);

	char sphi[64];
	sprintf(sphi, "%4d", *i);
	qui_ngon(qc, 6, qui_val_nm_ngon, PWM, bg);
	qui_txt(qc, nm, PWT, fg);
	qui_ngon(qc, 6, qui_val_val_ngon, PWM, bg);
	qui_txt(qc, sphi, PWT2, fgv);
	qui_ngon(qc, 6, qui_val_unt_ngon, PWM, bg);
	qui_txt(qc, unt, PWT3, fg);

	return 0;
}