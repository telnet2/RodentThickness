/**

Copyright (c) 2011 Ingo Berg 
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

*/


/*=========================================================================

 Program:   Insight Segmentation & Registration Toolkit
 Module:    $RCSfile: itkPixelMathImageFilter.h,v $
 Language:  C++
 Date:      $Date: 2008-10-07 17:31:02 $
 Version:   $Revision: 1.25 $

 Copyright (c) Insight Software Consortium. All rights reserved.
 See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.

 =========================================================================*/
#ifndef __itkPixelMathImageFilter_h
#define __itkPixelMathImageFilter_h

#include "itkInPlaceImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"

#define NUM_MAX_VARS 5

#include "string"

/** \class PixelMathImageFilter
 * \brief Implements pixel-wise generic operation on one image.
 *
 * This class is parameterized over the type of the input image and
 * the type of the output image.  It is also parameterized by the
 * operation to be applied, using a Functor style.
 *
 * PixelMathImageFilter allows the output dimension of the filter
 * to be larger than the input dimension. Thus subclasses of the
 * PixelMathImageFilter (like the CastImageFilter) can be used
 * to promote a 2D image to a 3D image, etc.
 *
 * \sa BinaryFunctorImageFilter TernaryFunctorImageFilter
 *
 * \ingroup   IntensityImageFilters     Multithreaded
 */
using namespace itk;

template <class TInputImage, class TOutputImage >
class ITK_EXPORT PixelMathImageFilter : public InPlaceImageFilter<TInputImage,TOutputImage>
{
public:
    /** Standard class typedefs. */
    typedef PixelMathImageFilter                       Self;
    typedef InPlaceImageFilter<TInputImage,TOutputImage>  Superclass;
    typedef SmartPointer<Self>                            Pointer;
    typedef SmartPointer<const Self>                      ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(PixelMathImageFilter, InPlaceImageFilter);

    /** Some typedefs. */
    typedef TInputImage                              InputImageType;
    typedef typename    InputImageType::ConstPointer InputImagePointer;
    typedef typename    InputImageType::RegionType   InputImageRegionType;
    typedef typename    InputImageType::PixelType    InputImagePixelType;

    typedef TOutputImage                             OutputImageType;
    typedef typename     OutputImageType::Pointer    OutputImagePointer;
    typedef typename     OutputImageType::RegionType OutputImageRegionType;
    typedef typename     OutputImageType::PixelType  OutputImagePixelType;

    std::string& GetEquation() { return m_Equation; };
    const std::string& GetEquation() const { return m_Equation; }

	void SetEquation(const std::string& equation) {
        if (m_Equation != equation) {
            m_Equation = equation;
            this->Modified();
		}
	}

protected:
    PixelMathImageFilter();
    virtual ~PixelMathImageFilter() {};

    /** PixelMathImageFilter can produce an image which is a different
     * resolution than its input image.  As such, PixelMathImageFilter
     * needs to provide an implementation for
     * GenerateOutputInformation() in order to inform the pipeline
     * execution model.  The original documentation of this method is
     * below.
     *
     * \sa ProcessObject::GenerateOutputInformaton()  */
    virtual void GenerateOutputInformation();

    /** PixelMathImageFilter can be implemented as a multithreaded filter.
     * Therefore, this implementation provides a ThreadedGenerateData() routine
     * which is called for each processing thread. The output image data is
     * allocated automatically by the superclass prior to calling
     * ThreadedGenerateData().  ThreadedGenerateData can only write to the
     * portion of the output image specified by the parameter
     * "outputRegionForThread"
     *
     * \sa ImageToImageFilter::ThreadedGenerateData(),
     *     ImageToImageFilter::GenerateData()  */
    void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                              ThreadIdType threadId );

    void ComputeIntensityStatistics();

private:
    PixelMathImageFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    std::string m_Equation;
    std::vector<InputImagePixelType> m_MinIntensities;
    std::vector<InputImagePixelType> m_MaxIntensities;
};

#ifndef ITK_MANUAL_INSTANTIATION
#include "PixelMathImageFilter.hxx"
#endif

#endif
