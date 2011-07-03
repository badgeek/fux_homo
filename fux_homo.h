 /* ------------------------------------------------------------------
  * GEM - Graphics Environment for Multimedia
  *
  *  Copyright (c) 2004 tigital@mac.com
  *  For information on usage and redistribution, and for a DISCLAIMER
  *  OF ALL WARRANTIES, see the file, "GEM.LICENSE.TERMS"
  *
  * ------------------------------------------------------------------
  */

#ifndef INCLUDE_GEM_GLMULTMATRIXF_H_
#define INCLUDE_GEM_GLMULTMATRIXF_H_

#include "GemGLBase.h"
#include <math.h>

/*
 CLASS
	fux_homo
 KEYWORDS
	openGL	0
 DESCRIPTION
	wrapper for the openGL-function
	"glMultMatrixf( GLfloat matrix)"
 */

struct ofPoint
{
	float x;
	float y;
};

class GEM_EXTERN fux_homo : public GemGLBase
{
	CPPEXTERN_HEADER(fux_homo, GemGLBase)

	public:
	  // Constructor
	  fux_homo (t_float);	// CON
	  void findHomography(ofPoint src[4], ofPoint dst[4], float homography[16]);
	  void gaussian_elimination(float *input, int n);
	
	protected:
	  // Destructor
	  virtual ~fux_homo ();
      // check extensions
      virtual bool isRunnable(void);

	  // Do the rendering
	  virtual void	render (GemState *state);

	  // variables
	  GLfloat	m_matrix[16];		// VAR
	  virtual void	srcMess(int argc, t_atom*argv);	// FUN
	  virtual void	destMess(int argc, t_atom*argv);// FUN


	private:

	  //we need some inlets
	  t_inlet *src_inlet;
	  t_inlet *dest_inlet;
	
	  //static member functions
	  static void	 srcMessCallback (void*, t_symbol*, int, t_atom*);
	  static void	 destMessCallback (void*, t_symbol*, int, t_atom*);
	
	  //quad var
	  float quad_src[8];
	  float quad_dest[8];
	
};
#endif // for header file
