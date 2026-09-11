/* Minimal GL/gl.h stub.
 *
 * CUDA's cudaGL.h includes <GL/gl.h> for the graphics-interop entry points,
 * and Roihu's GPU nodes carry the GLX headers but not Mesa's gl.h, so the
 * include fails before anything of ours is compiled. TurboGAP uses no GL
 * interop at all -- the header arrives only because hop's runtime shim pulls
 * in cuda_gl_interop.h -- so the four typedefs cudaGL.h needs are enough.
 *
 * Delete this the day Roihu ships mesa headers, or the day hop stops
 * including cuda_gl_interop.h.
 */
#ifndef TURBOGAP_GL_STUB_H
#define TURBOGAP_GL_STUB_H

typedef unsigned int GLenum;
typedef unsigned int GLuint;
typedef int GLint;
typedef int GLsizei;
typedef unsigned char GLboolean;
typedef signed char GLbyte;
typedef short GLshort;
typedef unsigned char GLubyte;
typedef unsigned short GLushort;
typedef unsigned long GLulong;
typedef float GLfloat;
typedef float GLclampf;
typedef double GLdouble;
typedef double GLclampd;
typedef void GLvoid;

#endif
