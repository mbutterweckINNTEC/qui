#define QUI_MAN_AXIS_X_N 2
#define QUI_MAN_CRCL_X_N 360
#define QUI_MAN_AXIS_Y_N 2
#define QUI_MAN_CRCL_Y_N 360
#define QUI_MAN_AXIS_Z_N 2
#define QUI_MAN_CRCL_Z_N 360
#define QUI_MAN_AXIS_V_N 2
#define QUI_MAN_CRCL_V_N 360

#define QUI_MAN_PLN_X_N (QUI_MAN_AXIS_X_N + QUI_MAN_CRCL_X_N)
#define QUI_MAN_PLN_Y_N (QUI_MAN_AXIS_Y_N + QUI_MAN_CRCL_Y_N)
#define QUI_MAN_PLN_Z_N (QUI_MAN_AXIS_Z_N + QUI_MAN_CRCL_Z_N)
#define QUI_MAN_PLN_V_N (QUI_MAN_AXIS_Z_N + QUI_MAN_CRCL_Z_N)

#define QUI_MAN_RSZ_N 12

#define QUI_MAN_ALL_N (QUI_MAN_PLN_X_N + QUI_MAN_PLN_Y_N + QUI_MAN_PLN_Z_N + QUI_MAN_PLN_V_N + QUI_MAN_RSZ_N)

#define QUI_MAN_ATTR_N 2
#define QUI_MAN_SZ (QUI_MAN_ALL_N * QUI_MAN_ATTR_N * sizeof(float3_t))

#define QUI_MAN_X_CLR (float3_t) { 1.0f, 0.5f, 0.5f }
#define QUI_MAN_Y_CLR (float3_t) { 0.5f, 1.0f, 0.5f }
#define QUI_MAN_Z_CLR (float3_t) { 0.5f, 0.5f, 1.0f }
#define QUI_MAN_V_CLR (float3_t) { 1.0f, 0.75f, 0.5f }	/* view axis */
#define QUI_MAN_S_CLR (float3_t) { 0.5f, 0.75f, 1.0f }	/* resize frame */

#define QUI_MAN_R_XYZ 0.875f
#define QUI_MAN_R_V 1.f
#define QUI_MAN_L_XYZ 0.75f

#define QUI_MAN_S_XY 0.875f
#define QUI_MAN_S_DXY 0.9375f

