//
//  kimage.cpp
//  ktools
//
//  Created by Joohwi Lee on 1/22/14.
//
//

#include "kimage.h"
#include <itkSignedDanielssonDistanceMapImageFilter.h>
#include <itkInvertIntensityImageFilter.h>


void SetArrayTuple(vtkDataArray* a, int i, VectorType v) {
    for (int j = 0; j < a->GetNumberOfComponents(); j++) {
        a->SetComponent(i, j, v[j]);
    }
}

void SetArrayTuple(vtkDataArray* a, int i, double v) {
    a->SetTuple1(i, v);
}



/// @brief Compute the distance map for inside and outside both. The input image must be binary. Otherwise, it will not produce a correct result.
VectorImageType::Pointer ComputeDistanceMap(MaskImageType::Pointer binaryMap) {
    cout << "Computing distance map ..." << flush;

    typedef itk::SignedDanielssonDistanceMapImageFilter<MaskImageType, ImageType> SignedDistanceMapFilterType;

    // construct signed distance filter
    SignedDistanceMapFilterType::Pointer distmapFilter = SignedDistanceMapFilterType::New();
    distmapFilter->SetInput(binaryMap);
    distmapFilter->InsideIsPositiveOff();
    distmapFilter->UseImageSpacingOn();
    distmapFilter->Update();
    SignedDistanceMapFilterType::OutputImagePointer distmap = distmapFilter->GetDistanceMap();


    // invert the mask image
    typedef itk::InvertIntensityImageFilter<MaskImageType> InvertFilterType;
    InvertFilterType::Pointer invertFilter = InvertFilterType::New();
    invertFilter->SetInput(binaryMap);
    invertFilter->SetMaximum(1);
    invertFilter->Update();

    // construct signed distance filter
    SignedDistanceMapFilterType::Pointer distmapFilterInv = SignedDistanceMapFilterType::New();
    distmapFilterInv->SetInput(invertFilter->GetOutput());
    distmapFilterInv->InsideIsPositiveOff();
    distmapFilterInv->UseImageSpacingOn();
    distmapFilterInv->Update();
    SignedDistanceMapFilterType::OutputImagePointer distmapInv = distmapFilterInv->GetDistanceMap();

//    pi::ImageIO<MaskImageType> mio;
//    mio.WriteImage("a.nii.gz", invertFilter->GetOutput());
//    pi::ImageIO<ImageType> rio;
//    rio.WriteImage("b.nii.gz", distmapInv);

    // merge two distance maps
    /** Pointer Type for the vector distance image */
    SignedDistanceMapFilterType::VectorImageType::Pointer distanceOffsetImage = distmapFilter->GetVectorDistanceMap();
    SignedDistanceMapFilterType::VectorImageType::Pointer distanceOffsetImageInv = distmapFilterInv->GetVectorDistanceMap();

    pi::ImageIO<VectorImageType> io;
    VectorImageType::Pointer distanceVectorImage = io.NewImageS<MaskImageType>(binaryMap);

    const int nPixels = binaryMap->GetPixelContainer()->Size();

    SignedDistanceMapFilterType::VectorImageType::PixelType* offsetPointer = distanceOffsetImage->GetBufferPointer();
    SignedDistanceMapFilterType::VectorImageType::PixelType* offsetPointerInv = distanceOffsetImageInv->GetBufferPointer();
    MaskImageType::PixelType* maskImagePointer = binaryMap->GetBufferPointer();

    VectorType* vectorPointer = distanceVectorImage->GetBufferPointer();

    for (int i = 0; i < nPixels; i++) {
        if (maskImagePointer[i] > 0) {
            for (int k = 0; k < 3; k++) {
                vectorPointer[i][k] = offsetPointerInv[i][k];
            }
        } else {
            for (int k = 0; k < 3; k++) {
                vectorPointer[i][k] = offsetPointer[i][k];
            }
        }
    }
    return distanceVectorImage;
}
