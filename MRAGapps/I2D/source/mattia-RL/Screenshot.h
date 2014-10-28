/*
 *  screenshot.c
 *
 */
#pragma once

#include "RL_Environment.h"

#ifdef _RL_VIZ

#ifdef _WIN32
#include "windows.h"
#include "GL/gl.h"
#endif
  
#ifdef __APPLE__
#include "GLUT/glut.h"
#endif 

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
//#include <string>

//#ifndef linux
#define BYTE_SWAP(x)  x  
//#else
//#define BYTE_SWAP(x) x = ((x) >> 8) + ((x) << 8)
//#endif
#pragma pack(1)
typedef struct
{
    GLbyte	identsize;              // Size of ID field that follows header (0)
    GLbyte	colorMapType;           // 0 = None, 1 = paletted
    GLbyte	imageType;              // 0 = none, 1 = indexed, 2 = rgb, 3 = grey, +8=rle
    unsigned short	colorMapStart;          // First colour map entry
    unsigned short	colorMapLength;         // Number of colors
    unsigned char 	colorMapBits;   // bits per palette entry
    unsigned short	xstart;                 // image x origin
    unsigned short	ystart;                 // image y origin
    unsigned short	width;                  // width in pixels
    unsigned short	height;                 // height in pixels
    GLbyte	bits;                   // bits per pixel (8 16, 24, 32)
    GLbyte	descriptor;             // image descriptor
} TGAHEADER;
#pragma pack(8)

inline GLint gltWriteTGA(const char *szFileName, int nSizeX = 0, int nSizeY = 0)
{
    FILE *pFile;                // File pointer
    TGAHEADER tgaHeader;		// TGA file header
    unsigned long lImageSize;   // Size in bytes of image
    GLbyte	*pBits = NULL;      // Pointer to bits
    GLint iViewport[4];         // Viewport in pixels
	GLint nImageSize[2];
    //GLint lastBuffer;          // Storage for the current read buffer setting
    
	bool bUseViewport = (nSizeX == 0 && nSizeY == 0);
	
	if (bUseViewport)
	{
		glGetIntegerv(GL_VIEWPORT, iViewport);
		nImageSize[0] =  iViewport[2];
		nImageSize[1] =  iViewport[3];
	}
	else
	{
		nImageSize[0] = nSizeX;
		nImageSize[1] = nSizeY;
	}
	// Get the viewport dimensions
	
//	printf("Screenshot: viewport: %d %d %d %d\n", iViewport[0], iViewport[1], iViewport[2], iViewport[3]);
	
    // How big is the image going to be (targas are tightly packed)
	lImageSize = nImageSize[0]* nImageSize[1]*4;	
	
    // Allocate block. If this doesn't work, go home
    pBits = (GLbyte *)malloc(lImageSize);
    if(pBits == NULL)
        return 0;
	
    // Read bits from color buffer
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
	glPixelStorei(GL_PACK_ROW_LENGTH, 0);
	glPixelStorei(GL_PACK_SKIP_ROWS, 0);
	glPixelStorei(GL_PACK_SKIP_PIXELS, 0);
    
    // Get the current read buffer setting and save it. Switch to
    // the front buffer and do the read operation. Finally, restore
    // the read buffer state
   // glGetIntegerv(GL_READ_BUFFER, &lastBuffer);
    //glReadBuffer(GL_BACK);

	glReadPixels(0, 0, nImageSize[0], nImageSize[1], GL_BGRA_EXT, GL_UNSIGNED_BYTE, pBits);
		
   // glReadBuffer(lastBuffer);
    
    // Initialize the Targa header
    tgaHeader.identsize = 0;
    tgaHeader.colorMapType = 0;
    tgaHeader.imageType = 2;
    tgaHeader.colorMapStart = 0;
    tgaHeader.colorMapLength = 0;
    tgaHeader.colorMapBits = 0;
    tgaHeader.xstart = 0;
    tgaHeader.ystart = 0;
    tgaHeader.width = nImageSize[0];
    tgaHeader.height = nImageSize[1];
    tgaHeader.bits = 32;
    tgaHeader.descriptor = 8;
    
    // Do byte swap for big vs little endian
#ifdef __APPLE__
    BYTE_SWAP(tgaHeader.colorMapStart);
    BYTE_SWAP(tgaHeader.colorMapLength);
    BYTE_SWAP(tgaHeader.xstart);
    BYTE_SWAP(tgaHeader.ystart);
    BYTE_SWAP(tgaHeader.width);
    BYTE_SWAP(tgaHeader.height);
#endif
    
    // Attempt to open the file
    pFile = fopen(szFileName, "wb");
    if(pFile == NULL)
	{
        free(pBits);    // Free buffer and return error
        return 0;
	}
	
    // Write the header
    fwrite(&tgaHeader, sizeof(TGAHEADER), 1, pFile);
    
    // Write the image data
    fwrite(pBits, lImageSize, 1, pFile);
	
    // Free temporary buffer and close the file
    free(pBits);    
    fclose(pFile);
    
    // Success!
    return 1;
}

#endif
