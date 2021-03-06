////
//  DLib: A simple image processing library.
//
//  David Crandall, 2003-2005
//  crandall@cs.cornell.edu
//
//  Please do not redistribute this code.
//
//
//
//
#ifndef __DRAWTEXT_H_
#define __DRAWTEXT_H_

#include <math.h>
#include "SImage.h"
#include "DrawText.h"


/*  GIMP header image file format (RGB): /Users/dave/maya/test_imgs/letters2.h  */

static int width = 564;
static int height = 13;

/*  Call this macro repeatedly.  After each use, the pixel data can be extracted  */

#define HEADER_PIXEL(data,pixel) {\
  pixel[0] = (((data[0] - 33) << 2) | ((data[1] - 33) >> 4)); \
  pixel[1] = ((((data[1] - 33) & 0xF) << 4) | ((data[2] - 33) >> 2)); \
  pixel[2] = ((((data[2] - 33) & 0x3) << 6) | ((data[3] - 33))); \
  data += 4; \
}
static const char *header_data =
	"````````````````````````````````````````````````````````````````"
	"````````````````!!!!````````````````````````````````````````````"
	"````````````````````````````````````````````````````````!!!!````"
	"````!!!!````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````!!!!````````````````````````````````````````````"
	"````````````````````````!!!!````````````````!!!!````!!!!````````"
	"````````!!!!````!!!!````````!!!!!!!!!!!!````````````!!!!````````"
	"````!!!!````!!!!!!!!````````````````````````!!!!````````````````"
	"````!!!!````````````````!!!!````````````````````!!!!````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````!!!!````````!!!!"
	"!!!!!!!!````````````````!!!!````````````````!!!!!!!!!!!!````````"
	"````!!!!!!!!!!!!````````````````````!!!!````````!!!!!!!!!!!!!!!!"
	"!!!!````````````!!!!!!!!````````!!!!!!!!!!!!!!!!!!!!````````!!!!"
	"!!!!!!!!````````````!!!!!!!!!!!!````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````!!!!!!!!!!!!````````````````"
	"!!!!!!!!````````````````!!!!````````````!!!!!!!!!!!!!!!!````````"
	"````!!!!!!!!!!!!````````!!!!!!!!!!!!````````````!!!!!!!!!!!!!!!!"
	"!!!!````!!!!!!!!!!!!!!!!!!!!````````!!!!!!!!!!!!````````!!!!````"
	"````````!!!!````````!!!!!!!!!!!!````````````!!!!!!!!!!!!!!!!````"
	"!!!!````````````!!!!````!!!!````````````````````!!!!````````````"
	"!!!!````!!!!````````````!!!!````````!!!!!!!!!!!!````````!!!!!!!!"
	"!!!!!!!!````````````!!!!!!!!!!!!````````!!!!!!!!!!!!!!!!````````"
	"````!!!!!!!!!!!!````````!!!!!!!!!!!!!!!!!!!!````!!!!````````````"
	"!!!!````!!!!````````````!!!!````!!!!````````````!!!!````!!!!````"
	"````````!!!!````!!!!````````````!!!!````!!!!!!!!!!!!!!!!!!!!````"
	"````````!!!!!!!!!!!!````````!!!!````````````````````!!!!!!!!!!!!"
	"````````````````!!!!````````````````````````````````````````````"
	"````````````````````````````````````````!!!!````````````````````"
	"````````````````````````````````````````!!!!````````````````````"
	"````````````````!!!!!!!!````````````````````````````````!!!!````"
	"````````````````````````!!!!````````````````````````!!!!````````"
	"!!!!````````````````````````````!!!!````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````!!!!````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````!!!!!!!!````````````!!!!````````````!!!!!!!!````````"
	"````````````````````````````````````````!!!!````````````````!!!!"
	"````!!!!````````````````!!!!````!!!!````!!!!````!!!!````!!!!````"
	"!!!!````!!!!````````!!!!!!!!````````!!!!````````````````````!!!!"
	"````````````````!!!!````````````````````````!!!!````````!!!!````"
	"!!!!````!!!!````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"!!!!````!!!!````````````!!!!````````!!!!!!!!````````````!!!!````"
	"````````!!!!````!!!!````````````!!!!````````````!!!!!!!!````````"
	"!!!!````````````````````````!!!!````````````````!!!!````````````"
	"!!!!````!!!!````````````!!!!````!!!!````````````!!!!````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````!!!!````````````"
	"!!!!````````!!!!````````!!!!````````````!!!!````````````!!!!````"
	"````````!!!!````!!!!````````````!!!!````!!!!````````!!!!````````"
	"!!!!````````````````````!!!!````````````````````!!!!````````````"
	"!!!!````!!!!````````````!!!!````````````!!!!````````````````````"
	"````````!!!!````!!!!````````!!!!````````!!!!````````````````````"
	"!!!!!!!!````!!!!!!!!````!!!!!!!!````````!!!!````!!!!````````````"
	"!!!!````!!!!````````````!!!!````!!!!````````````!!!!````!!!!````"
	"````````!!!!````!!!!````````````!!!!````````````!!!!````````````"
	"!!!!````````````!!!!````!!!!````````````!!!!````!!!!````````````"
	"!!!!````!!!!````````````!!!!````!!!!````````````!!!!````````````"
	"````````!!!!````````````!!!!````````````````!!!!````````````````"
	"````````````!!!!````````````````!!!!````````````````````````````"
	"````````````````!!!!````````````````````````````````````!!!!````"
	"````````````````````````````````````````````````````````!!!!````"
	"````````````````````````````!!!!````````````````````````````````"
	"````````!!!!````````````````````````````````````````````````````"
	"````````````````!!!!````````````````````````````!!!!````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````!!!!````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````!!!!````````````````````!!!!````````````"
	"````````!!!!````````````````````````````````````````````!!!!````"
	"````````````!!!!````!!!!````````!!!!!!!!!!!!!!!!!!!!!!!!!!!!````"
	"!!!!````````````!!!!````!!!!````!!!!````!!!!````````````````````"
	"````````````!!!!````````````````!!!!````````````````````````!!!!"
	"````````````!!!!!!!!!!!!````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````!!!!````````!!!!````````!!!!!!!!````!!!!````!!!!````"
	"````````````````````````!!!!````````````````````!!!!````````!!!!"
	"````!!!!````````!!!!````````````````````!!!!````````````````````"
	"````````````````!!!!````!!!!````````````!!!!````!!!!````````````"
	"!!!!````````````!!!!!!!!````````````````!!!!!!!!````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````!!!!````!!!!````````!!!!!!!!````````````!!!!````"
	"````````!!!!````````````!!!!````!!!!````````````````````!!!!````"
	"````````!!!!````!!!!````````````````````!!!!````````````````````"
	"!!!!````````````````````!!!!````````````!!!!````````````!!!!````"
	"````````````````````````!!!!````!!!!````!!!!````````````!!!!````"
	"````````````````!!!!!!!!````!!!!!!!!````!!!!!!!!````````!!!!````"
	"!!!!````````````!!!!````!!!!````````````!!!!````!!!!````````````"
	"!!!!````!!!!````````````!!!!````!!!!````````````````````````````"
	"!!!!````````````!!!!````````````!!!!````!!!!````````````!!!!````"
	"!!!!````!!!!````!!!!````````!!!!````!!!!````````````!!!!````!!!!"
	"````````````````````!!!!````````````````!!!!````````````````````"
	"!!!!````````````````````````!!!!````````````!!!!````!!!!````````"
	"````````````````````````````````````!!!!````````````!!!!!!!!!!!!"
	"!!!!````!!!!````!!!!!!!!````````````!!!!!!!!!!!!````````````!!!!"
	"!!!!!!!!!!!!````````!!!!!!!!!!!!````````!!!!!!!!!!!!!!!!````````"
	"````!!!!!!!!!!!!!!!!````!!!!````!!!!!!!!````````````!!!!!!!!````"
	"````````````!!!!!!!!!!!!````````!!!!````````!!!!````````````````"
	"!!!!````````````!!!!!!!!````!!!!````````!!!!````!!!!!!!!````````"
	"````!!!!!!!!!!!!````````!!!!````!!!!!!!!````````````!!!!!!!!!!!!"
	"!!!!````!!!!!!!!````!!!!!!!!````````!!!!!!!!!!!!````````!!!!!!!!"
	"!!!!!!!!````````!!!!````````````!!!!````!!!!````````````!!!!````"
	"!!!!````````````!!!!````!!!!````````````!!!!````!!!!````````````"
	"!!!!````!!!!!!!!!!!!!!!!!!!!````````````!!!!````````````````````"
	"!!!!````````````````````!!!!````````````````````````````````````"
	"````````!!!!````````````````````````````````````````````!!!!````"
	"!!!!````````!!!!!!!!````````````````!!!!````!!!!````````````!!!!"
	"````````````````````````````````````````````!!!!````````````````"
	"````````````````!!!!````````!!!!!!!!!!!!````````````````!!!!````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````!!!!````````!!!!````!!!!````!!!!````"
	"````````!!!!````````````````````````!!!!````````````````!!!!!!!!"
	"````````!!!!````````!!!!````````!!!!!!!!!!!!!!!!````````!!!!!!!!"
	"!!!!!!!!````````````````````!!!!````````````!!!!!!!!!!!!````````"
	"!!!!````````````!!!!````````````!!!!!!!!````````````````!!!!!!!!"
	"````````````````````````!!!!!!!!````````````````````````!!!!!!!!"
	"````````````````````````````!!!!````````!!!!````!!!!````!!!!````"
	"````!!!!````!!!!````````!!!!!!!!!!!!!!!!````````!!!!````````````"
	"````````!!!!````````````!!!!````!!!!!!!!!!!!!!!!````````!!!!!!!!"
	"!!!!!!!!````````!!!!````````````````````!!!!!!!!!!!!!!!!!!!!````"
	"````````!!!!````````````````````````````!!!!````!!!!!!!!````````"
	"````````!!!!````````````````````!!!!````!!!!````!!!!````!!!!````"
	"!!!!````!!!!````!!!!````````````!!!!````!!!!````````````!!!!````"
	"!!!!````````````!!!!````!!!!````````````!!!!````````!!!!!!!!!!!!"
	"````````````````!!!!````````````!!!!````````````!!!!````!!!!````"
	"````````!!!!````!!!!````!!!!````!!!!````````````!!!!````````````"
	"````!!!!````!!!!````````````````!!!!````````````````````!!!!````"
	"````````````````!!!!````````````````````````!!!!````````````!!!!"
	"````!!!!````````````````````````````````````````````````````````"
	"!!!!````````````!!!!````!!!!!!!!````````!!!!````!!!!````````````"
	"!!!!````!!!!````````````!!!!````!!!!````````````!!!!````````!!!!"
	"````````````````!!!!````````````!!!!````!!!!!!!!````````!!!!````"
	"````````!!!!````````````````````````!!!!````````!!!!````!!!!````"
	"````````````````!!!!````````````!!!!````!!!!````!!!!````!!!!!!!!"
	"````````!!!!````!!!!````````````!!!!````!!!!!!!!````````!!!!````"
	"!!!!````````````!!!!````````!!!!!!!!````````````!!!!````````````"
	"!!!!````````!!!!````````````````!!!!````````````!!!!````!!!!````"
	"````````!!!!````!!!!````````````!!!!````````!!!!````!!!!````````"
	"!!!!````````````!!!!````````````````!!!!````````````````!!!!````"
	"````````````````!!!!````````````````````!!!!````````````````!!!!"
	"!!!!````````!!!!````````!!!!````````````````````````````````````"
	"````````!!!!````!!!!````````````!!!!!!!!````````````````!!!!````"
	"!!!!````!!!!````!!!!````````````````````````````````````````!!!!"
	"````````````````````````````````!!!!````!!!!````!!!!````!!!!````"
	"````````!!!!````````````````````````````````````````````````````"
	"````````````````````````````````````````!!!!````````````!!!!!!!!"
	"````````!!!!````````````!!!!````````````````````!!!!````````````"
	"````````````````!!!!````!!!!!!!!!!!!!!!!!!!!````!!!!````````````"
	"!!!!````!!!!````````````!!!!````````````````!!!!````````!!!!````"
	"````````!!!!````````!!!!!!!!!!!!!!!!````````````````````````````"
	"````````````````````````````````!!!!!!!!````````!!!!!!!!!!!!!!!!"
	"!!!!````````````!!!!!!!!````````````````!!!!````````````!!!!````"
	"!!!!````!!!!````````!!!!````!!!!````````!!!!````````````!!!!````"
	"!!!!````````````````````!!!!````````````!!!!````!!!!````````````"
	"````````!!!!````````````````````!!!!````````!!!!!!!!````!!!!````"
	"````````!!!!````````````!!!!````````````````````````````!!!!````"
	"!!!!!!!!````````````````!!!!````````````````````!!!!````````````"
	"!!!!````!!!!````!!!!````!!!!````!!!!````````````!!!!````!!!!!!!!"
	"!!!!!!!!````````!!!!````````````!!!!````!!!!!!!!!!!!!!!!````````"
	"````````````````!!!!````````````!!!!````````````!!!!````````````"
	"!!!!````````!!!!````!!!!````````!!!!````!!!!````!!!!````````````"
	"!!!!````````````````````!!!!````````````````````!!!!````````````"
	"````````!!!!````````````````````````!!!!````````````````````!!!!"
	"````````!!!!````````````!!!!````````````````````````````````````"
	"````````````````!!!!````````````!!!!````!!!!````````````!!!!````"
	"!!!!````````````````````!!!!````````````!!!!````!!!!!!!!!!!!!!!!"
	"!!!!````````!!!!````````````````!!!!````````````!!!!````!!!!````"
	"````````!!!!````````````!!!!````````````````````````!!!!````````"
	"!!!!!!!!````````````````````````!!!!````````````!!!!````!!!!````"
	"!!!!````!!!!````````````!!!!````!!!!````````````!!!!````!!!!````"
	"````````!!!!````!!!!````````````!!!!````````!!!!````````````````"
	"````!!!!!!!!````````````````!!!!````````````````!!!!````````````"
	"!!!!````````!!!!````!!!!````````!!!!````!!!!````!!!!````````````"
	"!!!!````````````````!!!!````!!!!````````````````!!!!````````````"
	"````````!!!!````````````````````!!!!````````````````````!!!!````"
	"````````!!!!````````!!!!!!!!````````````!!!!````````````````````"
	"````````````````!!!!!!!!!!!!!!!!!!!!!!!!````````!!!!````!!!!````"
	"````!!!!````!!!!````!!!!!!!!````````!!!!````!!!!````````````````"
	"````````````!!!!````````````````````````````````!!!!````````````"
	"!!!!````````````!!!!!!!!!!!!!!!!!!!!````````````````````````````"
	"!!!!!!!!!!!!!!!!!!!!````````````````````````````````!!!!````````"
	"````````!!!!````````````!!!!````````````!!!!````````````````!!!!"
	"````````````````````````````````!!!!````````````````!!!!````````"
	"````````````````!!!!````!!!!````````````!!!!````````````!!!!````"
	"````````!!!!````````````!!!!````````````````````!!!!````````````"
	"````````````````````````````````````````!!!!!!!!````````````````"
	"````````````````````````````````````````!!!!!!!!````````!!!!````"
	"````````!!!!````````!!!!!!!!````!!!!!!!!!!!!!!!!!!!!````!!!!````"
	"````````!!!!````!!!!````````````````````!!!!````````````!!!!````"
	"!!!!````````````````````!!!!````````````````````!!!!````````````"
	"!!!!````!!!!````````````!!!!````````````!!!!````````````!!!!````"
	"````````!!!!````!!!!````!!!!````````````!!!!````````````````````"
	"!!!!````````````!!!!````!!!!````````!!!!!!!!````!!!!````````````"
	"!!!!````!!!!````````````````````!!!!````````````!!!!````!!!!````"
	"!!!!````````````````````````````!!!!````````````!!!!````````````"
	"!!!!````````````!!!!````````!!!!````!!!!````````````!!!!````!!!!"
	"````````````!!!!````!!!!````````````````!!!!````````````````!!!!"
	"````````````````````````!!!!````````````````````````!!!!````````"
	"````````````!!!!````````````````````````````````````````````````"
	"````````````````````````````````!!!!````````````!!!!````!!!!````"
	"````````!!!!````!!!!````````````````````!!!!````````````!!!!````"
	"!!!!````````````````````````!!!!````````````````!!!!````````````"
	"!!!!````!!!!````````````!!!!````````````!!!!````````````````````"
	"````!!!!````````!!!!````!!!!````````````````````!!!!````````````"
	"!!!!````!!!!````!!!!````!!!!````````````!!!!````!!!!````````````"
	"!!!!````!!!!````````````!!!!````!!!!````````````!!!!````````!!!!"
	"````````````````````````````!!!!````````````!!!!````````````````"
	"!!!!````````````!!!!````````!!!!````!!!!````````!!!!````!!!!````"
	"!!!!````````````!!!!````````````````!!!!````!!!!````````````!!!!"
	"````````````````!!!!!!!!````````````````````````!!!!````````````"
	"````````````!!!!!!!!````````````````````````````````````````````"
	"````````````````````````````````````!!!!````!!!!````````!!!!````"
	"!!!!````!!!!````!!!!````````!!!!````!!!!!!!!````````````!!!!````"
	"````````````````````````````!!!!````````````````````````````````"
	"!!!!````````````````````````````````````!!!!````````````````````"
	"!!!!!!!!````````````````````````````````````````!!!!!!!!````````"
	"````!!!!````````````````!!!!````````````!!!!````````````!!!!````"
	"````````!!!!````````````````````!!!!````````````!!!!````````````"
	"````!!!!````````````````````````!!!!````!!!!````````````!!!!````"
	"````````!!!!````````````!!!!````````````!!!!````````````````!!!!"
	"````````````````!!!!!!!!````````````````!!!!!!!!````````````````"
	"!!!!!!!!````````!!!!!!!!!!!!!!!!!!!!````````````!!!!!!!!````````"
	"````````````````````````````!!!!````````````````!!!!````````````"
	"!!!!````!!!!````````````!!!!````!!!!````````````!!!!````!!!!````"
	"````!!!!````````!!!!````````````````````!!!!````````````````````"
	"!!!!````````````!!!!````!!!!````````````!!!!````````````!!!!````"
	"````````!!!!````````````!!!!````!!!!````````!!!!````````!!!!````"
	"````````````````!!!!````````````!!!!````!!!!````````!!!!!!!!````"
	"!!!!````````````!!!!````!!!!````````````````````!!!!````````````"
	"!!!!````!!!!````````!!!!````````!!!!````````````!!!!````````````"
	"!!!!````````````!!!!````````````!!!!````````````!!!!````````````"
	"````!!!!````!!!!````````!!!!````````````!!!!````````````!!!!````"
	"````````!!!!````````````````````````````!!!!````````````````````"
	"````````!!!!````````````````!!!!````````````````````````````````"
	"````````````````````````````````````````````````!!!!````````!!!!"
	"!!!!````!!!!````````````!!!!````!!!!````````````!!!!````!!!!````"
	"````!!!!!!!!````!!!!````````````!!!!````````!!!!````````````````"
	"!!!!````````!!!!!!!!````!!!!````````````!!!!````````````!!!!````"
	"````````````````````!!!!````````!!!!````````!!!!````````````````"
	"!!!!````````````!!!!````!!!!````!!!!````!!!!````````````!!!!````"
	"!!!!````````````!!!!````!!!!````````````!!!!````!!!!````````!!!!"
	"!!!!````````!!!!````````````````!!!!````````````!!!!````````!!!!"
	"````````!!!!````!!!!````````!!!!!!!!````````````!!!!````````````"
	"````!!!!````!!!!````````````!!!!````!!!!````````````````!!!!````"
	"````````!!!!````````````````````````````!!!!````````````````````"
	"!!!!````````````````````!!!!````````````````````````````````````"
	"````````!!!!````````````````````````````````````````!!!!````!!!!"
	"````````````!!!!!!!!!!!!````````!!!!````````````!!!!````````!!!!"
	"!!!!!!!!````!!!!````````````````````````````````!!!!````````````"
	"````````````!!!!````````````````````````````````````````!!!!````"
	"````````````````!!!!!!!!````````````````````````````````````````"
	"!!!!!!!!````````!!!!````````````````````````!!!!!!!!!!!!````````"
	"!!!!!!!!!!!!!!!!!!!!````!!!!!!!!!!!!!!!!!!!!````````!!!!!!!!!!!!"
	"````````````````````!!!!````````!!!!!!!!!!!!!!!!````````````!!!!"
	"!!!!!!!!````````````````!!!!````````````````!!!!!!!!!!!!````````"
	"````!!!!!!!!````````````````````!!!!!!!!````````````````!!!!!!!!"
	"````````````````````````!!!!!!!!````````````````````````!!!!!!!!"
	"````````````````````````!!!!````````````````````!!!!!!!!!!!!````"
	"!!!!````````````!!!!````!!!!!!!!!!!!!!!!````````````!!!!!!!!!!!!"
	"````````!!!!!!!!!!!!````````````!!!!!!!!!!!!!!!!!!!!````!!!!````"
	"````````````````````!!!!!!!!!!!!````````!!!!````````````!!!!````"
	"````!!!!!!!!!!!!````````````!!!!!!!!!!!!````````!!!!````````````"
	"!!!!````!!!!!!!!!!!!!!!!!!!!````!!!!````````````!!!!````!!!!````"
	"````````!!!!````````!!!!!!!!!!!!````````!!!!````````````````````"
	"````!!!!!!!!!!!!````````!!!!````````````!!!!````````!!!!!!!!!!!!"
	"````````````````!!!!````````````````!!!!!!!!!!!!````````````````"
	"!!!!````````````````!!!!````!!!!````````!!!!````````````!!!!````"
	"````````!!!!````````````!!!!!!!!!!!!!!!!!!!!````````````!!!!````"
	"````````````````````````!!!!````````````````!!!!````````````````"
	"````````````````````````````````````````````````````````````````"
	"````!!!!!!!!````!!!!````!!!!!!!!!!!!!!!!````````````!!!!!!!!!!!!"
	"````````````!!!!!!!!````!!!!````````!!!!!!!!!!!!````````````!!!!"
	"````````````````````!!!!!!!!````!!!!````!!!!````````````!!!!````"
	"````````!!!!!!!!````````````````````!!!!````````!!!!````````````"
	"!!!!````````````!!!!!!!!````````!!!!````!!!!````!!!!````!!!!````"
	"````````!!!!````````!!!!!!!!!!!!````````!!!!!!!!!!!!!!!!````````"
	"````!!!!!!!!````!!!!````!!!!!!!!!!!!````````````````!!!!!!!!!!!!"
	"````````````````!!!!!!!!````````````!!!!!!!!````!!!!````````````"
	"!!!!````````````````!!!!````!!!!````````!!!!````````````!!!!````"
	"````````!!!!````````````!!!!!!!!!!!!!!!!!!!!````````````!!!!````"
	"````````````````!!!!````````````````````!!!!````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````!!!!````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"!!!!````````````````````````!!!!````````````````````````````````"
	"````````````````````````````````````!!!!````````````````````````"
	"````````````````````````````````!!!!````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````!!!!````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````!!!!````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````!!!!````````````````````````````````````````````````!!!!"
	"````````````````````````````````!!!!!!!!!!!!!!!!!!!!````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````!!!!````````````"
	"````````````````````````````````````````````````````!!!!````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````!!!!````"
	"````````````````````````````````!!!!````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````!!!!````````````````````````````````````````"
	"````````!!!!````````````````````!!!!````````````````````!!!!````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````!!!!````````````````!!!!````````````````````"
	"````````````````````````````````````````````````!!!!````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````!!!!````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````!!!!!!!!````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````!!!!````````````````````````````````````"
	"````````````!!!!````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````!!!!!!!!!!!!"
	"````````````````````````````````````````````````````````!!!!!!!!"
	"!!!!````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````!!!!````````````````````````````````````!!!!````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````!!!!````````````````````````````"
	"````````````````````````!!!!````````````````````````````````````"
	"````````!!!!````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````!!!!````````!!!!````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````!!!!!!!!!!!!````````````"
	"````````````````````!!!!!!!!!!!!````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````!!!!!!!!````````````"
	"````````````````!!!!!!!!````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````";



static const int numbers[10][120] = {
    {0, 0, 0, 1, 1, 0, 0, 0,
     0, 0, 1, 1, 1, 1, 0, 0,
     0, 1, 0, 0, 0, 0, 1, 0,
     0, 1, 0, 0, 0, 0, 1, 0,
     0, 1, 0, 0, 0, 0, 1, 0,
     0, 1, 0, 0, 0, 0, 1, 0,
     0, 1, 0, 0, 0, 0, 1, 0,
     0, 1, 0, 0, 0, 0, 1, 0,
     0, 1, 0, 0, 0, 0, 1, 0,
     0, 1, 0, 0, 0, 0, 1, 0,
     0, 1, 0, 0, 0, 0, 1, 0,
     0, 0, 1, 1, 1, 1, 0, 0,
     0, 0, 0, 1, 1, 0, 0, 0},

    {0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 1, 0, 0, 0,
     0, 0, 0, 0, 1, 0, 0, 0,
     0, 0, 0, 0, 1, 0, 0, 0,
     0, 0, 0, 0, 1, 0, 0, 0,
     0, 0, 0, 0, 1, 0, 0, 0,
     0, 0, 0, 0, 1, 0, 0, 0,
     0, 0, 0, 0, 1, 0, 0, 0,
     0, 0, 0, 0, 1, 0, 0, 0,
     0, 0, 0, 0, 1, 0, 0, 0,
     0, 0, 0, 0, 1, 0, 0, 0,
     0, 0, 0, 0, 1, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0},

    {0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 1, 1, 1, 0, 0, 0,
     0, 1, 0, 0, 0, 1, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0,
     0, 0, 0, 0, 1, 0, 0, 0,
     0, 0, 0, 1, 0, 0, 0, 0,
     0, 0, 1, 0, 0, 0, 0, 0,
     0, 1, 0, 0, 0, 0, 0, 0,
     0, 1, 0, 0, 0, 0, 0, 0,
     0, 1, 1, 1, 1, 1, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0},

    {0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 1, 1, 1, 0, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0,
     0, 0, 0, 0, 1, 0, 0, 0,
     0, 0, 0, 1, 0, 0, 0, 0,
     0, 0, 0, 0, 1, 0, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0,
     0, 1, 0, 0, 0, 1, 0, 0,
     0, 0, 1, 1, 1, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0},

    {0, 0, 0, 0, 0, 0, 0, 0,
     0, 1, 0, 0, 0, 1, 0, 0,
     0, 1, 0, 0, 0, 1, 0, 0,
     0, 1, 0, 0, 0, 1, 0, 0,
     0, 1, 0, 0, 0, 1, 0, 0,
     0, 1, 0, 0, 0, 1, 0, 0,
     0, 1, 1, 1, 1, 1, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0},

    {0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 1, 1, 1, 1, 1, 0,
     0, 0, 1, 0, 0, 0, 0, 0,
     0, 0, 1, 0, 0, 0, 0, 0,
     0, 0, 1, 0, 0, 0, 0, 0,
     0, 0, 1, 1, 1, 1, 1, 0,
     0, 0, 0, 0, 0, 0, 1, 0,
     0, 0, 0, 0, 0, 0, 1, 0,
     0, 0, 0, 0, 0, 0, 1, 0,
     0, 0, 0, 0, 0, 0, 1, 0,
     0, 0, 0, 0, 0, 0, 1, 0,
     0, 0, 1, 1, 1, 1, 1, 0,
     0, 0, 0, 0, 0, 0, 0, 0},

    {0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 1, 1, 1, 1, 0, 0,
     0, 0, 1, 0, 0, 0, 0, 0,
     0, 0, 1, 0, 0, 0, 0, 0,
     0, 0, 1, 0, 0, 0, 0, 0,
     0, 0, 1, 0, 0, 0, 0, 0,
     0, 0, 1, 1, 1, 1, 1, 0,
     0, 0, 1, 0, 0, 0, 1, 0,
     0, 0, 1, 0, 0, 0, 1, 0,
     0, 0, 1, 0, 0, 0, 1, 0,
     0, 0, 1, 0, 0, 0, 1, 0,
     0, 0, 1, 1, 1, 1, 1, 0,
     0, 0, 0, 0, 0, 0, 0, 0},

    {0, 0, 0, 0, 0, 0, 0, 0,
     0, 1, 1, 1, 1, 1, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0},

    {0, 0, 0, 0, 0, 0, 0, 0,
     0, 1, 1, 1, 1, 1, 1, 0,
     0, 1, 0, 0, 0, 0, 1, 0,
     0, 1, 0, 0, 0, 0, 1, 0,
     0, 1, 0, 0, 0, 0, 1, 0,
     0, 1, 0, 0, 0, 0, 1, 0,
     0, 1, 1, 1, 1, 1, 1, 0,
     0, 1, 0, 0, 0, 0, 1, 0,
     0, 1, 0, 0, 0, 0, 1, 0,
     0, 1, 0, 0, 0, 0, 1, 0,
     0, 1, 0, 0, 0, 0, 1, 0,
     0, 1, 1, 1, 1, 1, 1, 0,
     0, 0, 0, 0, 0, 0, 0, 0},

    {0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 1, 1, 1, 1, 1, 0,
     0, 1, 0, 0, 0, 0, 1, 0,
     0, 1, 0, 0, 0, 0, 1, 0,
     0, 1, 0, 0, 0, 0, 1, 0,
     0, 1, 0, 0, 0, 0, 1, 0,
     0, 1, 1, 1, 1, 1, 1, 0,
     0, 0, 0, 0, 0, 0, 1, 0,
     0, 0, 0, 0, 0, 0, 1, 0,
     0, 0, 0, 0, 0, 0, 1, 0,
     0, 0, 0, 0, 0, 0, 1, 0,
     0, 0, 0, 0, 0, 0, 1, 0,
     0, 0, 0, 0, 0, 0, 0, 0}};



static const int chr_rows=13;
static const int chr_cols=6;

static const int num_rows=13;
static const int num_cols=8;

static int first=1;

void draw_text(SDoublePlane &img, const char *str, int row, int col, int value, int scale)
{
  static SDoublePlane characters(height, width);

  if(first)
    {
      const char *data = header_data;
      for(int i=0; i<height; i++)
	for(int j=0; j<width; j++)
	  {
	    char pixel[4];

	    HEADER_PIXEL(data, pixel);

	    characters[i][j] = pixel[0];
	  }

      first=0;
    }

  for(size_t i=0; i<strlen(str); i++)
    {
      int char_no = str[i] - 33;

      if(char_no < 0) 
	{
	  col += chr_cols * scale;
	  continue;
	}

      if(row >= 0 && col >= 0 && row + chr_rows*scale < img.rows() && col + chr_cols*scale < img.cols())
	for(int i=0; i<chr_rows; i++)
	  for(int j=0; j<chr_cols; j++)
	    for(int ii=0; ii<scale; ii++)
	      for(int jj=0; jj<scale; jj++)
		{
		  img[row+i*scale+ii][col+j*scale+jj]=characters[i][char_no*chr_cols+j];
		}

      col = (col + chr_cols*scale);
    }
}

#endif
