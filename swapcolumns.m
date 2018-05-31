/* swapcolumns.c  Swaps columns of X in-place.
 *                Relies on user to make sure X is not shared with another variable
 *                since there is no checking for this in the code.
 * Syntax:  swapcolumns(X,i,j)
 *          X = a matrix (any standard class ... no objects)
 *          i,j = columns to swap
 * Programmer:  James Tursa
 * Date:        June 24, 2016
 */
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    unsigned char b;
    unsigned char *data, *target, *source;
    size_t M, N, i, j, k, c, bytes;

      if( nrhs == 0 ) {
          mexPrintf("swapcolumns -> Swaps columns of X in-place.\n");
          mexPrintf("Relies on user to make sure X is not shared with another variable\n");
          mexPrintf("since there is no checking for this in the code.\n");
          mexPrintf("Syntax:  swapcolumns(X,i,j)\n");
          mexPrintf("         X = a matrix (any standard class ... no objects)\n");
          mexPrintf("         i,j = columns to swap\n");
          return;
      }
      if( nrhs != 3 ) {
          mexErrMsgTxt("Need exactly 3 inputs");
      }
      if( nlhs > 0 ) {
          mexErrMsgTxt("Too many outputs\n");
      }
      if( !(mxIsNumeric(prhs[0]) || mxIsChar(prhs[0]) || mxIsLogical(prhs[0]) || 
            mxIsCell(prhs[0])) ) {
          mexErrMsgTxt("1st argument needs to be standard class");
      }
      M = mxGetM(prhs[0]);
      N = mxGetN(prhs[0]);
      if( M == 0 || N == 0 ) {
          return;
      }
      i = (size_t) mxGetScalar(prhs[1]);
      j = (size_t) mxGetScalar(prhs[2]);
      if( i > N || j > N ) {
          mexErrMsgTxt("Column index(es) are too large");
      }
      if( i == 0 || j == 0 ) {
          mexErrMsgTxt("Column index(es) cannot be 0");
      }
      data = mxGetData(prhs[0]);
      bytes = M * mxGetElementSize(prhs[0]);
      c = 2;
      while( c-- ) {
          target = data + (i-1) * bytes;
          source = data + (j-1) * bytes;
          for( k=0; k<bytes; k++ ) {
              b = *target;
              *target++ = *source;
              *source++ = b;
          }
          if( mxIsNumeric(prhs[0]) && mxIsComplex(prhs[0]) ) {
              data = mxGetImagData(prhs[0]);
          } else {
              break;
          }
      }
  }