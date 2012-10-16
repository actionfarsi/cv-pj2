///////////////////////////////////////////////////////////////////////////
//
// NAME
//  BlendImages.cpp -- blend together a set of overlapping images
//
// DESCRIPTION
//  This routine takes a collection of images aligned more or less horizontally
//  and stitches together a mosaic.
//
//  The images can be blended together any way you like, but I would recommend
//  using a soft halfway blend of the kind Steve presented in the first lecture.
//
//  Once you have blended the images together, you should crop the resulting
//  mosaic at the halfway points of the first and last image.  You should also
//  take out any accumulated vertical drift using an affine warp.
//  Lucas-Kanade Taylor series expansion of the registration error.
//
// SEE ALSO
//  BlendImages.h       longer description of parameters
//
// Copyright ?Richard Szeliski, 2001.  See Copyright.h for more details
// (modified for CSE455 Winter 2003 and CS4670 Fall 2012)
//
///////////////////////////////////////////////////////////////////////////

#include "ImageLib/ImageLib.h"
#include "BlendImages.h"
#include <float.h>
#include <math.h>

#define MAX(x,y) (((x) < (y)) ? (y) : (x))
#define MIN(x,y) (((x) < (y)) ? (x) : (y))


// Return the closest integer to x, rounding up
static int iround(double x) {
	if (x < 0.0) {
		return (int) (x - 0.5);
	} else {
		return (int) (x + 0.5);
	}
}

/******************* TO DO *********************
* AccumulateBlend:
*	INPUT:
*		img: a new image to be added to acc
*		acc: portion of the accumulated image where img is to be added
*      M: transformation matrix for computing a bounding box
*		blendWidth: width of the blending function (horizontal hat function;
*	    try other blending functions for extra credit)
*	OUTPUT:
*		add a weighted copy of img to the subimage specified in acc
*		the first 3 band of acc records the weighted sum of pixel colors
*		the fourth band of acc records the sum of weight
*/
static void AccumulateBlend(CByteImage& img, CFloatImage& acc, CTransform3x3 M, float blendWidth)
{
	// BEGIN TODO
	// Fill in this routine
	//repeat of blend function to find max_x,max_y etc
	CShape sh        = img.Shape();
	int width        = sh.width;
	int height       = sh.height;
	int nBands       = sh.nBands;
	CShape sha        = acc.Shape();

	int widtha=sha.width;
	int heighta=sha.height;
	//compute image bounding box
	float min_x = FLT_MAX, min_y = FLT_MAX;
	float max_x = 0, max_y = 0;


	CByteImage& img_it = img;
	CTransform3x3 &T= M;
	CVector3 corners[4];

	//corner pixels
	//bottom left corner
	corners[0][0] = 0.0;
	corners[0][1] = 0.0;
	corners[0][2] = 1.0;

	//bottom right corner
	corners[1][0] = width - 1;
	corners[1][1] = 0.0;
	corners[1][2] = 1.0;

	//top left corner
	corners[2][0] = 0.0;
	corners[2][1] = height - 1;
	corners[2][2] = 1.0;

	//top right corner
	corners[3][0] = width - 1;
	corners[3][1] = height - 1;
	corners[3][2] = 1.0;
	for(int k=0; k<4; k++){
		//apply transform to corner
		corners[k]=T*corners[k];

		//divide by 
		corners[k][0] /= corners[0][2];
		corners[k][1] /= corners[0][2];

	}




	// *** BEGIN TODO #1 ***
	// add some code here to update min_x, ..., max_y
	for (int j = 0; j < 4; j++)
	{
		//find the transformed corner position that has the min and max vales
		CVector3 c = corners[j];
		//corner[0] is x
		min_x = float(MIN(min_x, c[0]));
		max_x = float(MAX(max_x, c[0]));
		//corner[1] is y
		min_y = float(MIN(min_y, c[1]));
		max_y = float(MAX(max_y, c[1]));
	}



	CTransform3x3 M_inv=M.Inverse(); //get inverse transform
	//now we go through every pixel in img, map it with M then put it in the appropriate place in acc
	for(int i=0; i<widtha; i++){
		for(int j=0; j<heighta; j++){
			CVector3 pdes,psrc;
			//vector for acc
			pdes[0]=i; //x
			pdes[1]=j;  //y
			pdes[2]=1.0; //weight at 1 for img

			//apply inverse transform, becomes vector for img
			psrc=M_inv*pdes;

			//divide by new weight
			float x=(float) (psrc[0]/psrc[2]);		
			float y=(float) (psrc[1]/psrc[2]);	
			//edge detect
			if (y < 0.0 || y >= height-1 || x<0.0 || x >= width-1){
				continue;
			}

			double weight=1.0;


			//distance from bounding box
			//get the min
			//dist from right bound
			int dist = i - min_x;
			//dist from left bound
			dist = MIN(dist, max_x - i);
			//weight is equivalent to the distance over blendwidth
			if (dist<=(blendWidth)){
				weight = dist/(blendWidth);
			}



			//skip if black pixel
			if(img.Pixel(floor(x),floor(y),0)==0 ||img.Pixel(floor(x),floor(y),1)==0 ||img.Pixel(floor(x),floor(y),2)==0){
				continue;
			}



			//accumulate the linear interpolated pixels
			acc.Pixel(i,j,0)+=(float)(weight*img.PixelLerp(x,y,0));

			acc.Pixel(i,j,1)+=(float)(weight*img.PixelLerp(x,y,1));
			acc.Pixel(i,j,2)+=(float)(weight*img.PixelLerp(x,y,2));

			acc.Pixel(i,j,3)+=(float)weight;



			//weight for now set to 1, no actual blending yet

		}
	}
	// END TODO
}


/******************* TO DO *********************
* NormalizeBlend:
*	INPUT:
*		acc: input image whose alpha channel (4th channel) contains
*		     normalizing weight values
*		img: where output image will be stored
*	OUTPUT:
*		normalize r,g,b values (first 3 channels) of acc and store it into img
*/
static void NormalizeBlend(CFloatImage& acc, CByteImage& img)
{
	// BEGIN TODO
	// fill in this routine..
	CShape sh        = acc.Shape();
	int width        = sh.width;
	int height       = sh.height;
	int nBands       = sh.nBands;
	for(int i=0; i<width; i++){
		for(int j=0; j<height; j++){

			double weight=acc.Pixel(i,j,3);
			if(weight>0){
				//divide by weights
				img.Pixel(i,j,0)=acc.Pixel(i,j,0)/weight;
				img.Pixel(i,j,1)=acc.Pixel(i,j,1)/weight;
				img.Pixel(i,j,2)=acc.Pixel(i,j,2)/weight;
				//set alpha to opaque
				img.Pixel(i,j,3)=255;

			}else{
				img.Pixel(i,j,0)=0;
				img.Pixel(i,j,1)=0;
				img.Pixel(i,j,2)=0;
			}
		}
	}

	// END TODO
}



/******************* TO DO *********************
* BlendImages:
*	INPUT:
*		ipv: list of input images and their relative positions in the mosaic
*		blendWidth: width of the blending function
*	OUTPUT:
*		create & return final mosaic by blending all images
*		and correcting for any vertical drift
*/
CByteImage BlendImages(CImagePositionV& ipv, float blendWidth)
{
	// Assume all the images are of the same shape (for now)
	CByteImage& img0 = ipv[0].img;
	CShape sh        = img0.Shape();
	int width        = sh.width;
	int height       = sh.height;
	int nBands       = sh.nBands;

	int n = ipv.size();
	if (n == 0) return CByteImage(0,0,1);

	bool is360 = false;

	if (ipv[0].imgName == ipv[n-1].imgName)
		is360 = true;

	// Compute the bounding box for the mosaic
	float min_x = FLT_MAX, min_y = FLT_MAX;
	float max_x = 0, max_y = 0;
	int i;
	float firstcorner,lastcorner;
	//go through each image, find corners
	for (i = 0; i < n; i++)
	{
		CByteImage& img_it = ipv[i].img;
		CTransform3x3 &T= ipv[i].position;
		CVector3 corners[4];

		//corner pixels
		//bottom left corner
		corners[0][0] = 0.0;
		corners[0][1] = 0.0;
		corners[0][2] = 1.0;

		//bottom right corner
		corners[1][0] = width - 1;
		corners[1][1] = 0.0;
		corners[1][2] = 1.0;

		//top left corner
		corners[2][0] = 0.0;
		corners[2][1] = height - 1;
		corners[2][2] = 1.0;

		//top right corner
		corners[3][0] = width - 1;
		corners[3][1] = height - 1;
		corners[3][2] = 1.0;
		for(int k=0; k<4; k++){
			//apply transform to corner
			corners[k]=T*corners[k];

			//divide by weight
			corners[k][0] /= corners[0][2];
			corners[k][1] /= corners[0][2];

		}

		//code for drift, keep track of first and last image top left corner
		if(i==0){
		firstcorner=corners[0][1];
		}
		if(i==n-1){
		lastcorner=corners[0][1];
		}

		// *** BEGIN TODO #1 ***
		// add some code here to update min_x, ..., max_y

		for (int j = 0; j < 4; j++)
		{
			//find the transformed corner position that has the min and max vales
			CVector3 c = corners[j];
			//corner[0] is x
			min_x = float(MIN(min_x, c[0]));
			max_x = float(MAX(max_x, c[0]));
			//corner[1] is y
			min_y = float(MIN(min_y, c[1]));
			max_y = float(MAX(max_y, c[1]));
		}

	}
	//

	//how to apply offset? max_x-min_x equal real width of image
	printf("w:%d, h:%d,mw:%lf,mh:%lf,maxw:%lf,maxh:%lf",width,height,min_x,min_y,max_x,max_y);
	// Create a floating point accumulation image
	CShape mShape((int)(ceil(max_x) - floor(min_x)),
		(int)(ceil(max_y) - floor(min_y)), nBands + 1);
	CFloatImage accumulator(mShape);
	accumulator.ClearPixels();

	double x_init, x_final;
	double y_init, y_final;

	// Add in all of the images
	for (i = 0; i < n; i++) {
		// Compute the sub-image involved
		CTransform3x3 &M = ipv[i].position;
		CTransform3x3 M_t = CTransform3x3::Translation(-min_x, -min_y) * M; //isnt this the offset?
		CByteImage& img = ipv[i].img;

		// Perform the accumulation
		AccumulateBlend(img, accumulator, M_t, blendWidth);

		if (i == 0) {
			CVector3 p;
			p[0] = 0.5 * width;
			p[1] = 0.0;
			p[2] = 1.0;

			p = M_t * p;
			x_init = p[0];
			y_init = p[1];
		} else if (i == n - 1) {
			CVector3 p;
			p[0] = 0.5 * width;
			p[1] = 0.0;
			p[2] = 1.0;

			p = M_t * p;
			x_final = p[0];
			y_final = p[1];
		}
	}

	// Normalize the results
	mShape = CShape((int)(ceil(max_x) - floor(min_x)),
		(int)(ceil(max_y) - floor(min_y)), nBands);

	CByteImage compImage(mShape);
	NormalizeBlend(accumulator, compImage);
	bool debug_comp = false;
	if (debug_comp)
		WriteFile(compImage, "tmp_comp.tga");

	// Allocate the final image shape
	int outputWidth;

	bool crop = false;  // set to true to crop
	if (crop) {
		outputWidth = mShape.width - width;
	} else {
		outputWidth = mShape.width;
	}

	CShape cShape(outputWidth, mShape.height, nBands);

	CByteImage croppedImage(cShape);

	// Compute the affine translations
	CTransform3x3 A;

	// BEGIN TODO
	// fill in appropriate entries in A to trim the left edge and
	// to take out the vertical drift if this is a 360 panorama
	// (i.e. is360 is true)
	if(is360 !=true){
		A[0][0] = 1;
		A[0][1] = 0;
		A[0][2] = 0;

		A[1][0] = 0; 
		A[1][1] = 1;
		A[1][2] = 0; 

		A[2][0] = 0;
		A[2][1] = 0;
		A[2][2] = 1;
	}else{
		//take out drift if 360 panorama
		A[0][0] = 1;
		A[0][1] = 0;
		A[0][2] = width/2; // x translation have width of image, this is the crop
		ipv[0].
		A[1][0] = firstcorner-lastcorner; 
		A[1][1] = 1;
		A[1][2] = 0; // y translation

		A[2][0] = 0;
		A[2][1] = 0;
		A[2][2] = 1;
	}
	// END TODO 

	// Warp and crop the composite
	WarpGlobal(compImage, croppedImage, A, eWarpInterpLinear);

	return croppedImage;
}
