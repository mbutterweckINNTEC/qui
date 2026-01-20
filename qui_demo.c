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

static int shdr_mk(char  const *vshsrc, char  const *fshsrc) {
        char log[4096];
        int po = 0, vso = 0, fso = 0, stts = 0;

        vso = glCreateShader(GL_VERTEX_SHADER);
        if (!vso)
                goto end;

        glShaderSource(vso, 1, &vshsrc, NULL);
        glCompileShader(vso);

        fso = glCreateShader(GL_FRAGMENT_SHADER);
        if (!fso)
                goto end;

        glShaderSource(fso, 1, (char const * const*)&fshsrc, NULL);
        glCompileShader(fso);

        glGetShaderiv(vso, GL_COMPILE_STATUS, &stts);
        if (!stts) {
                glGetShaderInfoLog(vso, 4096, NULL, log);
                printf("vsh: %d\n", stts);
                puts(log);
                goto end;
        }

        glGetShaderiv(fso, GL_COMPILE_STATUS, &stts);
        if (!stts) {
                glGetShaderInfoLog(fso, 4096, NULL, log);
                printf("fsh: %d\n", stts);
                puts(log);
                goto end;
        }

        po = glCreateProgram();
        if (!po)
                goto end;

        glAttachShader(po, vso);
        glAttachShader(po, fso);
        glLinkProgram(po);

        glGetProgramiv(po, GL_LINK_STATUS, &stts);
        if (!stts) {
                glGetProgramInfoLog(po, 4096, NULL, log);
                printf("po: %d\n", stts);
                puts(log);
                goto end;
        }

end:
        if (!stts || 0 == po) {
                glDeleteProgram(po);
                po = 0;
        }

        if (vso)
                glDeleteShader(vso);

        if (fso)
                glDeleteShader(fso);

        return po;
}

char *qui_vsh =
	"#version 440"				"\n"
	"in      vec4 v;"			"\n"
	"in      vec4 c;"			"\n"
	"uniform mat4 M;"			"\n"
	"out     vec4 k;"			"\n"
	""					"\n"
	"void main() {"				"\n"
	"	gl_Position = M * v;"		"\n"
	"	k = c;"				"\n"
	"}"					"\n";

char *qui_fsh =
	"#version 440"				"\n"
	"in  vec4 k;"				"\n"
	"out vec4 K;"				"\n"
	""					"\n"
	"void main() {"				"\n"
	"	K = k;"				"\n"
	"}"					"\n";

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

#define QUI_CRCL_N 64 
#define QUI_CRSS_N 8
#define QUI_VCTR_N 2 //4 
float qui_vrtx[(QUI_CRCL_N + QUI_CRSS_N + QUI_VCTR_N) * 8];

int main(int argc, char *argv[]) {
	GLFWwindow *wndw;
	int w = 1024, h = 1024, qui_po, cube_po, cube_bo[2], cube_vao, M_loc_cube, qui_bo, qui_vao, M_loc_qui;
	float bg[4] = { 170.f/255.f, 170.f/255.f, 170.f/255.f, 0};
	float44_t M = identity_sc;
	float44_t V = identity_sc;
	float44_t VM = identity_sc;
	float44_t P = identity_sc;
	float44_t PV = identity_sc;
	float44_t PVM = identity_sc;
	quaternion_t q = { 1.0 };
	float angl = 0.0, ar;
	double cx, cy, cl;

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

	qui_po = shdr_mk(qui_vsh, qui_fsh);
	cube_po = shdr_mk(cube_vsh, cube_fsh);

	int j = 0;
	for (int i = 0; i < QUI_CRCL_N; ++i) {
		angl = (float)i / (float)(QUI_CRCL_N - 1) * 2.0 * M_PI;
		qui_vrtx[j++] = cos(angl);
		qui_vrtx[j++] = sin(angl);
		qui_vrtx[j++] = 0.0;
		qui_vrtx[j++] = 1.0;

		float t_ = 2 * (float)i / (float)(QUI_CRCL_N - 1);
		float t	= 6.f * (t_ - (int)t_);
		int s = t;
		int f = t - s;
		switch(s) {
		case 0:
			qui_vrtx[j++] = 1.f;
			qui_vrtx[j++] = f;
			qui_vrtx[j++] = 0.f;
			break;
		case 1:
			qui_vrtx[j++] = 1.f - f;
			qui_vrtx[j++] = 1.f;
			qui_vrtx[j++] = 0.f;
			break;
		case 2:
			qui_vrtx[j++] = 0.f;
			qui_vrtx[j++] = 1.f;
			qui_vrtx[j++] = f;
			break;
		case 3:
			qui_vrtx[j++] = 0.f;
			qui_vrtx[j++] = 1.f - f;
			qui_vrtx[j++] = 1.f;
			break;
		case 4:
			qui_vrtx[j++] = f;
			qui_vrtx[j++] = 1.f;
			qui_vrtx[j++] = 0.f;
			break;
		case 5:
			qui_vrtx[j++] = 1.f;
			qui_vrtx[j++] = 0.f;
			qui_vrtx[j++] = 1.f - f;
			break;
		};

		qui_vrtx[j++] = 1.0;
	}

	qui_vrtx[j++] = 0.f;
	qui_vrtx[j++] = 0.f;
	qui_vrtx[j++] = 0.f;
	qui_vrtx[j++] = 1.f;
	qui_vrtx[j++] = 0.f;
	qui_vrtx[j++] = 0.f;
	qui_vrtx[j++] = 0.f;
	qui_vrtx[j++] = 1.f;

	qui_vrtx[j++] = 1.f;
	qui_vrtx[j++] = 0.f;
	qui_vrtx[j++] = 0.f;
	qui_vrtx[j++] = 1.f;
	qui_vrtx[j++] = 1.f;
	qui_vrtx[j++] = 1.f;
	qui_vrtx[j++] = 1.f;
	qui_vrtx[j++] = 1.f;

	qui_vrtx[j++] = 0.f;
	qui_vrtx[j++] = 0.f;
	qui_vrtx[j++] = 0.f;
	qui_vrtx[j++] = 1.f;
	qui_vrtx[j++] = 0.f;
	qui_vrtx[j++] = 0.f;
	qui_vrtx[j++] = 0.f;
	qui_vrtx[j++] = 1.f;

	qui_vrtx[j++] =-1.f;
	qui_vrtx[j++] = 0.f;
	qui_vrtx[j++] = 0.f;
	qui_vrtx[j++] = 1.f;
	qui_vrtx[j++] = 1.f;
	qui_vrtx[j++] = 1.f;
	qui_vrtx[j++] = 1.f;
	qui_vrtx[j++] = 1.f;

	qui_vrtx[j++] = 0.f;
	qui_vrtx[j++] = 0.f;
	qui_vrtx[j++] = 0.f;
	qui_vrtx[j++] = 1.f;
	qui_vrtx[j++] = 0.f;
	qui_vrtx[j++] = 0.f;
	qui_vrtx[j++] = 0.f;
	qui_vrtx[j++] = 1.f;

	qui_vrtx[j++] = 0.f;
	qui_vrtx[j++] = 1.f;
	qui_vrtx[j++] = 0.f;
	qui_vrtx[j++] = 1.f;
	qui_vrtx[j++] = 1.f;
	qui_vrtx[j++] = 1.f;
	qui_vrtx[j++] = 1.f;
	qui_vrtx[j++] = 1.f;

	qui_vrtx[j++] = 0.f;
	qui_vrtx[j++] = 0.f;
	qui_vrtx[j++] = 0.f;
	qui_vrtx[j++] = 1.f;
	qui_vrtx[j++] = 0.f;
	qui_vrtx[j++] = 0.f;
	qui_vrtx[j++] = 0.f;
	qui_vrtx[j++] = 1.f;

	qui_vrtx[j++] = 0.f;
	qui_vrtx[j++] =-1.f;
	qui_vrtx[j++] = 0.f;
	qui_vrtx[j++] = 1.f;
	qui_vrtx[j++] = 1.f;
	qui_vrtx[j++] = 1.f;
	qui_vrtx[j++] = 1.f;
	qui_vrtx[j++] = 1.f;

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
	
	glCreateBuffers(1, &qui_bo);
	glNamedBufferStorage(qui_bo, sizeof(qui_vrtx), qui_vrtx, GL_DYNAMIC_STORAGE_BIT | GL_MAP_WRITE_BIT);
	glGenVertexArrays(1, &qui_vao);
	glBindVertexArray(qui_vao);
	glBindBuffer(GL_ARRAY_BUFFER, qui_bo);
	glVertexAttribPointer(0, 4, GL_FLOAT, 0, 32, 0);
	glEnableVertexAttribArray(0); 
	glVertexAttribPointer(1, 4, GL_FLOAT, 0, 32, (char*)16);
	glEnableVertexAttribArray(1); 
	glBindVertexArray(0);

	glClearColor(bg[0], bg[1], bg[2], bg[3]);
	M_loc_cube = glGetUniformLocation(cube_po, "M");
	M_loc_qui = glGetUniformLocation(qui_po, "M");


	glDepthFunc(GL_LESS);

	while(!glfwWindowShouldClose(wndw)) {
		glfwGetWindowSize(wndw, &w, &h);
	
		ar = (float)w / (float)h;

		glViewport(0, 0, w, h);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		if (glfwGetKey(wndw, GLFW_KEY_UP) == GLFW_PRESS) {
			V = mul_float44(V, quaternion_to_rotation_matrix((quaternion_t){cos(0.03125),-sin(0.03125), 0, 0 }));
		}

		if (glfwGetKey(wndw, GLFW_KEY_DOWN) == GLFW_PRESS) {
			V = mul_float44(V, quaternion_to_rotation_matrix((quaternion_t){cos(0.03125), sin(0.03125), 0, 0 }));
		}

		if (glfwGetKey(wndw, GLFW_KEY_LEFT) == GLFW_PRESS) {
			V = mul_float44(V, quaternion_to_rotation_matrix((quaternion_t){cos(0.03125), 0,-sin(0.03125), 0 }));
		}

		if (glfwGetKey(wndw, GLFW_KEY_RIGHT) == GLFW_PRESS) {
			V = mul_float44(V, quaternion_to_rotation_matrix((quaternion_t){cos(0.03125), 0, sin(0.03125), 0 }));
		}

		if (glfwGetMouseButton(wndw, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS) {
			glfwGetCursorPos(wndw, &cx, &cy);
			float sx = 2.0 * cx / w - 1.0;
		     	float sy = 1.0 - 2.0 * cy / h;
			float ma = PV.m[0][0];
			float mb = PV.m[1][0];
			float mc = PV.m[0][1];
			float md = PV.m[1][1];
			float mdet = ma*md-mb*mc;
			float ima = md / mdet;
			float imb = -mb / mdet;
			float imc = -mc / mdet;
			float imd = ma / mdet;
			cx = ima * sx + imb * sy;
			cy = imc * sx + imd * sy;
//			cx = (2.0 * cx / w - 1.0) * ar;
///			cy = 1.0 - 2.0 * cy / h;
			cl = sqrt(cx*cx+cy*cy);
			if (cl <= 1.0) {
				printf("%f, %f, %f\n", cx, cy, cl);
			//	angl = M_PI * cl * 0.5;
			//	q.a = cos(angl);
			//	q.b = cx / cl * sin(angl);
			//	q.c = cy / cl * sin(angl);

				q.b = cx; 
				q.c = cy;
				q.a = sqrt(1.0 - cl*cl);
				q.d = 0.0;

				float *v = glMapNamedBufferRange(qui_bo, (QUI_CRCL_N + QUI_CRSS_N) * 8 * sizeof(float), QUI_VCTR_N * 8 * sizeof(float), GL_MAP_WRITE_BIT);

				if (v) {
					*v++ = 0.f;
					*v++ = 0.f;
					*v++ = 0.f;
					*v++ = 1.f;
					*v++ = cl; 
					*v++ = cl;
					*v++ = cl;
					*v++ = 1.f;

					*v++ = cx;
					*v++ = cy;
					*v++ = 0.f;
					*v++ = 1.f;
					*v++ = cl;
					*v++ = cl;
					*v++ = cl;
					*v++ = 1.f;
/*
					*v++ = 0.f;
					*v++ = 0.f;
					*v++ = 0.f;
					*v++ = 1.f;
					*v++ = cl; 
					*v++ = cl;
					*v++ = cl;
					*v++ = 1.f;

					*v++ = q.a;
					*v++ = q.b;
					*v++ = q.c;
					*v++ = 1.f;
					*v++ = cl;
					*v++ = cl;
					*v++ = cl;
					*v++ = 1.f;
*/

					glUnmapNamedBuffer(qui_bo);
				}
			}

		}

		M = quaternion_to_rotation_matrix(q);
		VM = mul_float44(M, V);

		P.m[0][0] = 1.f / ar;

		PV = mul_float44(V, P);
		PVM = mul_float44(VM, P);

		glEnable(GL_DEPTH_TEST);
		glUseProgram(cube_po);
		glUniformMatrix4fv(M_loc_cube, 1, 0, &PVM.m[0][0]);
		glBindVertexArray(cube_vao);
		glDrawElements(GL_TRIANGLES, sizeof(cube_i)/ sizeof(int), GL_UNSIGNED_INT, 0);
		
//		glDisable(GL_DEPTH_TEST);
		glUseProgram(qui_po);
		glUniformMatrix4fv(M_loc_qui, 1, 0, &PV.m[0][0]);
		glBindVertexArray(qui_vao);
		glDrawArrays(GL_LINE_STRIP, 0, QUI_CRCL_N);
		glDrawArrays(GL_LINES, QUI_CRCL_N, QUI_CRSS_N);
		glUniformMatrix4fv(M_loc_qui, 1, 0, &PV.m[0][0]);
		glDrawArrays(GL_LINES, QUI_CRCL_N + QUI_CRSS_N, QUI_VCTR_N);

		glfwSwapBuffers(wndw);
		glfwWaitEventsTimeout(1.0 / 30.0);

		angl += 0.005;
	}

	glfwTerminate();
	return 0;
}
