#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <unistd.h>

#include "mathutils.h"
#include "quaternion.h"

#include "qui_def.h"
#include "qui_in.h"
#include "qui_shdr.h"
#include "qui_ctx.h"
#include "qui_util.h"
#include "qui_man.h"
#include "qui_fnt.h"
#include "qui_txt.h"

char *cube_vsh =
	"#version 440"				"\n"
	"in      vec4 v;"			"\n"
	"in      vec4 c;"			"\n"
	"uniform mat4 M;"			"\n"
	"out     vec4 k;"			"\n"
	"out     vec4 p;"			"\n"
	""					"\n"
	"void main() {"				"\n"
	"	p = M * v;"			"\n"
	"	gl_Position = p;"		"\n"
	"	k = c;"				"\n"
	"}"					"\n";

char *cube_fsh =
	"#version 440"					"\n"
	"in  vec4 k;"					"\n"
	"in  vec4 p;"					"\n"
	"out vec4 K;"					"\n"
	""						"\n"
	"void main() {"					"\n"
	"	vec3 s = dFdx(p.xyz);"			"\n"
	"	vec3 t = dFdy(p.xyz);"			"\n"
	"	vec3 n = normalize(cross(s, t));"	"\n"
	"	vec3 L = normalize(vec3(0.5,0.5,1));"	"\n"
	"	K = k * dot(n, L);"			"\n"
	"}"						"\n";


float cube_v[] = {
	 0.5f,  0.375f, -0.25f, 1.f,
	 0.5f, -0.375f, -0.25f, 1.f,
	 0.5f,  0.375f,  0.25f, 1.f,
	 0.5f, -0.375f,  0.25f, 1.f,
	-0.5f,  0.375f, -0.25f, 1.f,
	-0.5f, -0.375f, -0.25f, 1.f,
	-0.5f,  0.375f,  0.25f, 1.f,
	-0.5f, -0.375f,  0.25f, 1.f
};

int cube_i[] = {
	4, 2, 0,
	2, 7, 3,
	6, 5, 7,
	1, 7, 5,
	0, 3, 1,
	4, 1, 5,
	4, 6, 2,
	2, 6, 7,
	6, 4, 5,
	1, 3, 7,
	0, 2, 3,
	4, 0, 1
};

float qui_scrll;
float qui_scrll_sgn = 1.f;

void qui_scrll_cb(GLFWwindow *wndw, double x, double y) {
	qui_scrll = y;
}

int main(int argc, char *argv[]) {
	GLFWwindow *wndw;
	int w = 1024, h = 1024, qui_po, cube_po, cube_bo[2], cube_vao, M_loc_cube, qui_bo, qui_vao, M_loc_qui;
	float bg[4] = { 170.f/255.f, 170.f/255.f, 170.f/255.f, 0};
	float44_t M = identity_sc;
	float44_t Q = identity_sc;
	float44_t V = (float44_t){
				1.f, 0.f, 0.f, 0.f,
				0.f, 1.f, 0.f, 0.f,
				0.f, 0.f, 1.f, 0.f,
				0.5, 0.3333, 0.0, 1.f
			};
	float44_t VM = identity_sc;
	float44_t P = orthographic(1.f, -1.f, 1.f, -1.f, -1.f, 1.f);
	float44_t PV = identity_sc;
	float44_t PVM = identity_sc;
	float44_t S = identity_sc;
	float44_t T = identity_sc;

	quaternion_t q = { 1.0 };
	float angl = 0.0, ar;
	double cx, cy, cl;
	struct qui_ctx qc = {};
	struct qui_in qi = {};
	struct qui_man qm = {};
	static struct qui_fnt fnt = {};

	glfwInit();

	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 5);
	glfwWindowHint(GLFW_DOUBLEBUFFER, GLFW_TRUE);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	wndw = glfwCreateWindow(w, h, "Quaternion UI prototype", NULL, NULL);
//	glfwSetScrollCallback(wndw, sqrrl);
	glfwMakeContextCurrent(wndw);
	glewExperimental = 1;
	glewInit();

	glfwSetScrollCallback(wndw, qui_scrll_cb);

	cube_po = shdr_mk(cube_vsh, cube_fsh);

	qui_ctx_mk(&qc);

	glCreateBuffers(2, cube_bo);
	glNamedBufferStorage(cube_bo[0], sizeof(cube_v), cube_v, 0);
	glNamedBufferStorage(cube_bo[1], sizeof(cube_i), cube_i, 0);
	glGenVertexArrays(1, &cube_vao);
	glBindVertexArray(cube_vao);
	glBindBuffer(GL_ARRAY_BUFFER, cube_bo[0]);
	glVertexAttribPointer(0, 4, GL_FLOAT, 0, 0, 0);
	glEnableVertexAttribArray(0); 
	glVertexAttribPointer(0, 4, GL_FLOAT, 0, 0, 0);
	glVertexAttrib4f(1, 1,1,1,1);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cube_bo[1]);
	glBindVertexArray(0);

	glClearColor(bg[0], bg[1], bg[2], bg[3]);
	M_loc_cube = glGetUniformLocation(cube_po, "M");

	glDepthFunc(GL_LESS);

	quaternion_t mq = { 1, 0, 0, 0 };
	float3_t mt = {0, 0, 0};
	float ms = 1.f;

	qui_fnt_ld(&fnt, "./iosevka.obj");

	glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
	qui_fnt(&qc, &fnt);

	while(!glfwWindowShouldClose(wndw)) {
		glfwGetWindowSize(wndw, &w, &h);
	
		ar = (float)w / (float)h;

		glViewport(0, 0, w, h);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		if (glfwGetKey(wndw, GLFW_KEY_UP) == GLFW_PRESS) {
			if (glfwGetKey(wndw, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS) {
				V = mul_float44(V, (float44_t){ 1, 0, 0, 0,    0, 1, 0, 0,   0, 0, 1.0, 0,   0, 0.125, 0, 1});
			} else {
				V = mul_float44(V, quaternion_to_rotation_matrix((quaternion_t){cos(0.03125),-sin(0.03125), 0, 0 }));
			}
		}

		if (glfwGetKey(wndw, GLFW_KEY_DOWN) == GLFW_PRESS) {
			if (glfwGetKey(wndw, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS) {
				V = mul_float44(V, (float44_t){ 1, 0, 0, 0,    0, 1, 0, 0,   0, 0, 1.0, 0,   0, -0.125, 0, 1});
			} else {
				V = mul_float44(V, quaternion_to_rotation_matrix((quaternion_t){cos(0.03125), sin(0.03125), 0, 0 }));
			}
		}

		if (glfwGetKey(wndw, GLFW_KEY_LEFT) == GLFW_PRESS) {
			if (glfwGetKey(wndw, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS) {
				V = mul_float44(V, (float44_t){ 1, 0, 0, 0,    0, 1, 0, 0,   0, 0, 1.0, 0, -0.125, 0, 0, 1});
			} else {
				V = mul_float44(V, quaternion_to_rotation_matrix((quaternion_t){cos(0.03125), 0,-sin(0.03125), 0 }));
			}
		}

		if (glfwGetKey(wndw, GLFW_KEY_RIGHT) == GLFW_PRESS) {
			if (glfwGetKey(wndw, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS) {
				V = mul_float44(V, (float44_t){ 1, 0, 0, 0,    0, 1, 0, 0,   0, 0, 1.0, 0,  0.125, 0, 0, 1});
			} else {
				V = mul_float44(V, quaternion_to_rotation_matrix((quaternion_t){cos(0.03125), 0, sin(0.03125), 0 }));
			}
		}

		if (glfwGetKey(wndw, GLFW_KEY_KP_ADD) == GLFW_PRESS) {
			V = mul_float44(V, (float44_t){ 2, 0, 0, 0,    0, 2, 0, 0,   0, 0, 2.0, 0,   0, 0, 0, 1});
		}

		if (glfwGetKey(wndw, GLFW_KEY_KP_SUBTRACT) == GLFW_PRESS) {
			V = mul_float44(V, (float44_t){ 0.5, 0, 0, 0,    0, 0.5, 0, 0,   0, 0, 0.5, 0,   0, 0, 0, 1});
		}

		if (glfwGetKey(wndw, GLFW_KEY_X) == GLFW_PRESS) {
			q.b = sqrt(1.0 - q.a*q.a);
			q.c = 0.0;
			q.d = 0.0;
		}

		if (glfwGetKey(wndw, GLFW_KEY_Y) == GLFW_PRESS) {
			q.c = sqrt(1.0 - q.a*q.a);
			q.b = 0.0;
			q.d = 0.0;
		}

		if (glfwGetKey(wndw, GLFW_KEY_Z) == GLFW_PRESS) {
			q.d = sqrt(1.0 - q.a*q.a);
			q.c = 0.0;
			q.b = 0.0;
		}

		if (glfwGetMouseButton(wndw, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS) {
			qui_in_prss(&qi, QUI_IN_LMB);
		} else {
			qui_in_rls(&qi, QUI_IN_LMB);
		}
		
		if (glfwGetKey(wndw, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
			qui_in_prss(&qi, QUI_IN_ESC);
		} else {
			qui_in_rls(&qi, QUI_IN_ESC);
		}

		{
			double x = 0.0, y = 0.0;
			glfwGetCursorPos(wndw, &x, &y);
			qui_in_mv(&qi, (float2_t){ (x / (double)w * 2.0 - 1.0), 1.0 - 2.0 * y / (double)h });
		}

		if (qui_scrll) {
			qui_in_scrll(&qi, qui_scrll);
			qui_scrll = 0.0;
		}


		P.m[0][0] = 1.f / ar;
	//	P.m[3][2] = 0.5f;	/* adding perspective */

		PV = mul_float44(V, P);

		glEnable(GL_DEPTH_TEST);
		glUseProgram(cube_po);
		glUniformMatrix4fv(M_loc_cube, 1, 0, &PVM.m[0][0]);
		glBindVertexArray(cube_vao);
		glDrawElements(GL_TRIANGLES, sizeof(cube_i) / sizeof(int), GL_UNSIGNED_INT, 0);

		if (qui_man(&qc, &qm, &qi, P, V, &mt, &mq, &ms)) {
		}

		S = (float44_t){
			ms, 0.f, 0.f, 0.f,
			0.f, ms, 0.f, 0.f,
			0.f, 0.f, ms, 0.f,
			0.f, 0.f, 0.f, 1.f
		};
		T = (float44_t){
			1.f, 0.f, 0.f, 0.f,
			0.f, 1.f, 0.f, 0.f,
			0.f, 0.f, 1.f, 0.f,
			mt.x, mt.y, mt.z, 1.f
		};
		Q = quaternion_to_rotation_matrix(mq);
		M = mul_float44(mul_float44(S, Q), T);
		VM = mul_float44(M, V);
		PVM = mul_float44(VM, P);

		qui_txt(&qc, "haha!", PV);

		glfwSwapBuffers(wndw);
		glfwWaitEventsTimeout(1.0 / 30.0);

		qui_in_nxt(&qi);

		angl += 0.005;
	}

	glfwTerminate();
	return 0;
}
