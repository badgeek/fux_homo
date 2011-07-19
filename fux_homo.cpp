////////////////////////////////////////////////////////
//
// GEM - Graphics Environment for Multimedia
//
// Implementation file
//
// Copyright (c) 2004 tigital@mac.com
//  For information on usage and redistribution, and for a DISCLAIMER
//  *  OF ALL WARRANTIES, see the file, "GEM.LICENSE.TERMS"
//
////////////////////////////////////////////////////////

#include "fux_homo.h"

CPPEXTERN_NEW_WITH_ONE_ARG ( fux_homo , t_floatarg, A_DEFFLOAT )

/////////////////////////////////////////////////////////
//
// fux_homo
//
/////////////////////////////////////////////////////////
// Constructor
//
fux_homo :: fux_homo	(t_floatarg arg0=0)
{
	src_inlet = inlet_new(this->x_obj, &this->x_obj->ob_pd, &s_float, gensym("src"));
	//dest_inlet = inlet_new(this->x_obj, &this->x_obj->ob_pd, &s_float, gensym("dest"));	
}
/////////////////////////////////////////////////////////
// Destructor
//
fux_homo :: ~fux_homo () {
	inlet_free(src_inlet);
	//inlet_free(dest_inlet);
}

//////////////////
// extension check
bool fux_homo :: isRunnable(void) {
  if(GLEW_VERSION_1_1)return true;
  error("your system does not support OpenGL-1.1");
  return false;
}


/////////////////////////////////////////////////////////
// Render
//
void fux_homo :: render(GemState *state) {
	glMultMatrixf (m_matrix);
}

/////////////////////////////////////////////////////////
// Variables
//
void fux_homo :: srcMess (int argc, t_atom*argv) {	// FUN
	if(argc!=16){
		error("need 16 src + dest xy coordinate");
		return;
		}
	int i;
	for (i=0;i<16;i++) {
	  quad_src[i]=(GLfloat)atom_getfloat(argv+i);
	}
	
	//findHomography();
	
	ofPoint src[4];
	ofPoint dest[4];
	
	src[0].x = -2;//quad_src[0];
	src[0].y = -2;//quad_src[1];
	
	src[1].x = 2;//quad_src[2];
	src[1].y = -2;//quad_src[3];
	
	src[2].x = 2;//quad_src[4];
	src[2].y = 2;//quad_src[5];
	
	src[3].x = -2;//quad_src[6];
	src[3].y = 2;//quad_src[7];
	
	//dest
	
	dest[0].x = quad_src[0];
	dest[0].y = quad_src[1];
	
	dest[1].x = quad_src[2];
	dest[1].y = quad_src[3];
	
	dest[2].x = quad_src[4];
	dest[2].y = quad_src[5];
	
	dest[3].x = quad_src[6];
	dest[3].y = quad_src[7];
	
	
	float matrix_homo[16];
	
	
	findHomography(src, dest, matrix_homo);
	
	for (i=0;i<16;i++) {
	  m_matrix[i] = (GLfloat) matrix_homo[i];
	//  post("matrix[%i] = %f", i, matrix_homo[i]);
	}

	setModified();
}
/////////////////////////////////////////////////////////
// Variables
//
void fux_homo :: destMess (int argc, t_atom*argv) {	// FUN
	if(argc!=8){
		error("need 8 xy coordinate");
		return;
		}
	int i;
	for (i=0;i<8;i++) {
	  quad_dest[i]=(GLfloat)atom_getfloat(argv+i);
	}
	setModified();
}
/////////////////////////////////////////////////////////
// homography functions
//

void fux_homo :: gaussian_elimination(float *input, int n)
{
    // ported to c from pseudocode in
    // http://en.wikipedia.org/wiki/Gaussian_elimination

    float * A = input;
    int i = 0;
    int j = 0;
    int m = n-1;
    while (i < m && j < n)
    {
        // Find pivot in column j, starting in row i:
        int maxi = i;
        for(int k = i+1; k<m; k++)
        {
            if(fabs(A[k*n+j]) > fabs(A[maxi*n+j]))
            {
                maxi = k;
            }
        }
        if (A[maxi*n+j] != 0)
        {
            //swap rows i and maxi, but do not change the value of i
            if(i!=maxi)
                for(int k=0; k<n; k++)
                {
                    float aux = A[i*n+k];
                    A[i*n+k]=A[maxi*n+k];
                    A[maxi*n+k]=aux;
                }
            //Now A[i,j] will contain the old value of A[maxi,j].
            //divide each entry in row i by A[i,j]
            float A_ij=A[i*n+j];
            for(int k=0; k<n; k++)
            {
                A[i*n+k]/=A_ij;
            }
            //Now A[i,j] will have the value 1.
            for(int u = i+1; u< m; u++)
            {
                //subtract A[u,j] * row i from row u
                float A_uj = A[u*n+j];
                for(int k=0; k<n; k++)
                {
                    A[u*n+k]-=A_uj*A[i*n+k];
                }
                //Now A[u,j] will be 0, since A[u,j] - A[i,j] * A[u,j] = A[u,j] - 1 * A[u,j] = 0.
            }

            i++;
        }
        j++;
    }

    //back substitution
    for(int i=m-2; i>=0; i--)
    {
        for(int j=i+1; j<n-1; j++)
        {
            A[i*n+m]-=A[i*n+j]*A[j*n+m];
            //A[i*n+j]=0;
        }
    }
}

void fux_homo :: findHomography(ofPoint src[4], ofPoint dst[4], float homography[16])
{

    // create the equation system to be solved
    //
    // from: Multiple View Geometry in Computer Vision 2ed
    //       Hartley R. and Zisserman A.
    //
    // x' = xH
    // where H is the homography: a 3 by 3 matrix
    // that transformed to inhomogeneous coordinates for each point
    // gives the following equations for each point:
    //
    // x' * (h31*x + h32*y + h33) = h11*x + h12*y + h13
    // y' * (h31*x + h32*y + h33) = h21*x + h22*y + h23
    //
    // as the homography is scale independent we can let h33 be 1 (indeed any of the terms)
    // so for 4 points we have 8 equations for 8 terms to solve: h11 - h32
    // after ordering the terms it gives the following matrix
    // that can be solved with gaussian elimination:

    float P[8][9]=
    {
        {-src[0].x, -src[0].y, -1,   0,   0,  0, src[0].x*dst[0].x, src[0].y*dst[0].x, -dst[0].x }, // h11
        {  0,   0,  0, -src[0].x, -src[0].y, -1, src[0].x*dst[0].y, src[0].y*dst[0].y, -dst[0].y }, // h12

        {-src[1].x, -src[1].y, -1,   0,   0,  0, src[1].x*dst[1].x, src[1].y*dst[1].x, -dst[1].x }, // h13
        {  0,   0,  0, -src[1].x, -src[1].y, -1, src[1].x*dst[1].y, src[1].y*dst[1].y, -dst[1].y }, // h21

        {-src[2].x, -src[2].y, -1,   0,   0,  0, src[2].x*dst[2].x, src[2].y*dst[2].x, -dst[2].x }, // h22
        {  0,   0,  0, -src[2].x, -src[2].y, -1, src[2].x*dst[2].y, src[2].y*dst[2].y, -dst[2].y }, // h23

        {-src[3].x, -src[3].y, -1,   0,   0,  0, src[3].x*dst[3].x, src[3].y*dst[3].x, -dst[3].x }, // h31
        {  0,   0,  0, -src[3].x, -src[3].y, -1, src[3].x*dst[3].y, src[3].y*dst[3].y, -dst[3].y }, // h32
    };

    gaussian_elimination(&P[0][0],9);

    // gaussian elimination gives the results of the equation system
    // in the last column of the original matrix.
    // opengl needs the transposed 4x4 matrix:
    float aux_H[]= { P[0][8],P[3][8],0,P[6][8], // h11  h21 0 h31
                     P[1][8],P[4][8],0,P[7][8], // h12  h22 0 h32
                     0      ,      0,0,0,       // 0    0   0 0
                     P[2][8],P[5][8],0,1
                   };      // h13  h23 0 h33

    for(int i=0; i<16; i++) homography[i] = aux_H[i];
}


/////////////////////////////////////////////////////////
// static member functions
//

void fux_homo :: obj_setupCallback(t_class *classPtr) {
	 class_addmethod(classPtr, (t_method)&fux_homo::srcMessCallback, gensym("src"), A_GIMME, A_NULL);
	 //class_addmethod(classPtr, (t_method)&fux_homo::destMessCallback, gensym("dest"), A_GIMME, A_NULL);
}

void fux_homo :: srcMessCallback (void* data, t_symbol*,int argc, t_atom*argv){
	GetMyClass(data)->srcMess( argc, argv);
}

void fux_homo :: destMessCallback (void* data, t_symbol*,int argc, t_atom*argv){
	GetMyClass(data)->destMess( argc, argv);
}
