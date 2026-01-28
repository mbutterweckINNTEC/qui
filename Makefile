CC = gcc
CFLGS = -ggdb -O1
LFLGS = -lglfw -lGLEW -lEGL -lGL -lGLU -lOpenGL -lm

qui_demo: qui_demo.c
	$(CC) $(CFLGS) $^ -o $@ $(LFLGS)

