//
//  kimage.h
//  ktools
//
//  Created by Joohwi Lee on 1/22/14.
//
//

#ifndef __ktools__kimage__
#define __ktools__kimage__

#include <iostream>
#include <string>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkImageData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataArray.h>

#include "piImageIO.h"

typedef itk::Image<float,3> ImageType;
typedef itk::Image< itk::Vector<double,3>, 3 > VectorImageType;
typedef VectorImageType::PixelType VectorType;
typedef itk::Image<short,3> MaskImageType;


VectorImageType::Pointer ComputeDistanceMap(MaskImageType::Pointer maskImage);

void SetArrayTuple(vtkDataArray* a, int i, VectorType v);
void SetArrayTuple(vtkDataArray* a, int i, double v);

template <class X>
void ConvertImageT(std::string& imageFile, vtkImageData* imgData, const char* attrName, int numberOfComponents, MaskImageType::Pointer maskImage) {
    pi::ImageIO<X> itkIO;
    typename X::Pointer srcImg = itkIO.ReadImage(imageFile.c_str());
    typename X::SizeType srcSize = srcImg->GetRequestedRegion().GetSize();
    typename X::PointType srcOrigin = srcImg->GetOrigin();
    typename X::SpacingType srcSpacing = srcImg->GetSpacing();
    
    imgData->SetOrigin(srcOrigin[0], srcOrigin[1], srcOrigin[2]);
    imgData->SetSpacing(srcSpacing[0], srcSpacing[1], srcSpacing[2]);
    imgData->SetDimensions(srcSize[0], srcSize[1], srcSize[2]);
    
    const int nPoints = srcImg->GetRequestedRegion().GetNumberOfPixels();
    
    vtkDoubleArray* attr = vtkDoubleArray::New();
    attr->SetNumberOfComponents(numberOfComponents);
    attr->SetName(attrName);
    attr->SetNumberOfTuples(nPoints);
    
    vtkPointData* pdata = imgData->GetPointData();
    
    switch (numberOfComponents) {
        case 1:
            pdata->SetScalars(attr);
            break;
        case 3:
            pdata->SetVectors(attr);
            break;
        case 9:
            pdata->SetTensors(attr);
            break;
        default:
            pdata->AddArray(attr);
            break;
    }
    
    int cnt = 0;
#pragma omp parallel for
    for (unsigned int z = 0; z < srcSize[2]; z++) {
        for (unsigned int y = 0; y < srcSize[1]; y++) {
            for (unsigned int x = 0; x < srcSize[0]; x++) {
                typename X::IndexType idx;
                idx[0] = x;
                idx[1] = y;
                idx[2] = z;
                typename X::PixelType v = srcImg->GetPixel(idx);
                /// If the mask image is not provided, copy the entire pixels
                /// Otherwise, copy only the masked region will be copied.
                if (maskImage.IsNull() || maskImage->GetPixel(idx) > 0) {
                    SetArrayTuple(attr, cnt, v);
                } else {
                    if (numberOfComponents == 1) {
                        SetArrayTuple(attr, cnt, 0);
                    } else if (numberOfComponents == 3){
                        VectorType zero;
                        zero.Fill(0);
                        SetArrayTuple(attr, cnt, zero);
                    }
                }
                cnt ++;
            }
        }
    }
}



/// @brief Convert a vector-valued itk image to vtkUnstructuredGrid
template <class X, class Y>
void ConvertVectorImageT(std::string& imageFile, vtkUnstructuredGrid* imgData, typename Y::Pointer maskImage, const char* attrName, int numberOfComponents) {
    pi::ImageIO<X> itkIO;
    typename X::Pointer srcImg = itkIO.ReadImage(imageFile.c_str());


    /// - Create an instance for the output grid
    vtkUnstructuredGrid* imageData = imgData;
    vtkPointData* pdata = imageData->GetPointData();

    /// - Create a point set to store valid points
    vtkPoints* pointSet = vtkPoints::New();

    /// - Create an array to store the pixel data
    vtkDoubleArray* attr = vtkDoubleArray::New();
    attr->SetNumberOfComponents(numberOfComponents);
    attr->SetName(attrName);


    std::vector<double> tuple;
    tuple.resize(numberOfComponents);

    /// - Loop over the entire pixel of the mask
    itk::ImageRegionIteratorWithIndex<Y> iter(maskImage, maskImage->GetBufferedRegion());
    for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter) {
        if (iter.Get() > 0) {

            typename X::PointType point;
            srcImg->TransformIndexToPhysicalPoint(iter.GetIndex(), point);

            /// - Add a point
            pointSet->InsertNextPoint(point[0], point[1], point[2]);

            typename X::PixelType pixel = srcImg->GetPixel(iter.GetIndex());
            if (numberOfComponents > 1) {
                for (int j = 0; j < numberOfComponents; j++) {
                    tuple[j] = -pixel[j];
                }
                cout << pixel << endl;
            }

            /// - Add a pixel value (a scalar or a vector)
            attr->InsertNextTupleValue(&tuple[0]);
        }
    }

    switch (numberOfComponents) {
        case 1:
            pdata->SetScalars(attr);
            break;
        case 3:
            pdata->SetVectors(attr);
            break;
        case 9:
            pdata->SetTensors(attr);
            break;
        default:
            pdata->AddArray(attr);
            break;
    }

    imageData->SetPoints(pointSet);
}


/// @brief Convert a scalar-valued itk image to vtkUnstructuredGrid
template <class X, class Y>
void ConvertImageT(std::string& imageFile, vtkUnstructuredGrid* imgData, typename Y::Pointer maskImage, const char* attrName, int numberOfComponents) {
    pi::ImageIO<X> itkIO;
    typename X::Pointer srcImg = itkIO.ReadImage(imageFile.c_str());


    /// - Create an instance for the output grid
    vtkUnstructuredGrid* imageData = vtkUnstructuredGrid::New();
    vtkPointData* pdata = imageData->GetPointData();

    /// - Create a point set to store valid points
    vtkPoints* pointSet = vtkPoints::New();

    /// - Create an array to store the pixel data
    vtkDoubleArray* attr = vtkDoubleArray::New();
    attr->SetNumberOfComponents(numberOfComponents);
    attr->SetName(attrName);


    std::vector<double> tuple;
    tuple.resize(numberOfComponents);

    /// - Loop over the entire pixel of the mask
    itk::ImageRegionIteratorWithIndex<Y> iter(maskImage, maskImage->GetBufferedRegion());
    for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter) {
        if (iter.Get() > 0) {
            typename X::PointType point;
            srcImg->TransformIndexToPhysicalPoint(iter.GetIndex(), point);

            /// - Add a point
            pointSet->InsertNextPoint(point[0], point[1], point[2]);

            typename X::PixelType pixel = srcImg->GetPixel(iter.GetIndex());
            
            /// - Add a pixel value (a scalar or a vector)
            attr->InsertNextValue(pixel);
        }
    }

    switch (numberOfComponents) {
        case 1:
            pdata->SetScalars(attr);
            break;
        case 3:
            pdata->SetVectors(attr);
            break;
        case 9:
            pdata->SetTensors(attr);
            break;
        default:
            pdata->AddArray(attr);
            break;
    }
    
    imageData->SetPoints(pointSet);
}

#endif /* defined(__ktools__kimage__) */
