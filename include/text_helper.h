// Std. Includes
#include <iostream>
#include <map>
#include <string>
// GLEW
#define GLEW_STATIC
#include <GL/glew.h>
#include <GL/gl.h>
// #include <GL/glx.h>
#include <GL/glu.h>
#include <GL/glext.h>
// GLM
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
// FreeType
#include <ft2build.h>
#include FT_FREETYPE_H
// GL includes


void init_text_helper(GLuint va, GLuint vb);
void RenderText(std::string text, GLfloat scale, glm::vec3 color); 
