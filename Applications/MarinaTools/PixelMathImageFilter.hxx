/*=========================================================================

 Program:   Insight Segmentation & Registration Toolkit
 Module:    $RCSfile: itkPixelMathImageFilter.txx,v $
 Language:  C++
 Date:      $Date: 2009-10-28 03:37:14 $
 Version:   $Revision: 1.34 $

 Copyright (c) Insight Software Consortium. All rights reserved.
 See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.

 =========================================================================*/
#ifndef __itkPixelMathImageFilter_txx
#define __itkPixelMathImageFilter_txx

#include "PixelMathImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include <itkStatisticsImageFilter.h>
#include "itkProgressReporter.h"
#include "iostream"

#include "muParser.h"


/**
 * Constructor
 */
template <class TInputImage, class TOutputImage >
PixelMathImageFilter<TInputImage,TOutputImage>
::PixelMathImageFilter()
{
    this->SetNumberOfRequiredInputs( 1 );
    this->InPlaceOff();
}

/**
 * PixelMathImageFilter can produce an image which is a different resolution
 * than its input image.  As such, PixelMathImageFilter needs to provide an
 * implementation for GenerateOutputInformation() in order to inform
 * the pipeline execution model.  The original documentation of this
 * method is below.
 *
 * \sa ProcessObject::GenerateOutputInformaton()
 */
template <class TInputImage, class TOutputImage>
void
PixelMathImageFilter<TInputImage,TOutputImage>
::GenerateOutputInformation()
{
    // do not call the superclass' implementation of this method since
    // this filter allows the input the output to be of different dimensions

    // get pointers to the input and output
    typename Superclass::OutputImagePointer      outputPtr = this->GetOutput();
    typename Superclass::InputImageConstPointer  inputPtr  = this->GetInput();

    if ( !outputPtr || !inputPtr)
    {
        return;
    }

    // Set the output image largest possible region.  Use a RegionCopier
    // so that the input and output images can be different dimensions.
    OutputImageRegionType outputLargestPossibleRegion;
    this->CallCopyInputRegionToOutputRegion(outputLargestPossibleRegion,
                                            inputPtr->GetLargestPossibleRegion());
    outputPtr->SetLargestPossibleRegion( outputLargestPossibleRegion );

    // Set the output spacing and origin
    const ImageBase<Superclass::InputImageDimension> *phyData;

    phyData
    = dynamic_cast<const ImageBase<Superclass::InputImageDimension>*>(this->GetInput());

    if (phyData)
    {
        // Copy what we can from the image from spacing and origin of the input
        // This logic needs to be augmented with logic that select which
        // dimensions to copy
        unsigned int i, j;
        const typename InputImageType::SpacingType&
        inputSpacing = inputPtr->GetSpacing();
        const typename InputImageType::PointType&
        inputOrigin = inputPtr->GetOrigin();
        const typename InputImageType::DirectionType&
        inputDirection = inputPtr->GetDirection();

        typename OutputImageType::SpacingType outputSpacing;
        typename OutputImageType::PointType outputOrigin;
        typename OutputImageType::DirectionType outputDirection;

        // copy the input to the output and fill the rest of the
        // output with zeros.
        for (i=0; i < Superclass::InputImageDimension; ++i)
        {
            outputSpacing[i] = inputSpacing[i];
            outputOrigin[i] = inputOrigin[i];
            for (j=0; j < Superclass::OutputImageDimension; j++)
            {
                if (j < Superclass::InputImageDimension)
                {
                    outputDirection[j][i] = inputDirection[j][i];
                }
                else
                {
                    outputDirection[j][i] = 0.0;
                }
            }
        }
        for (; i < Superclass::OutputImageDimension; ++i)
        {
            outputSpacing[i] = 1.0;
            outputOrigin[i] = 0.0;
            for (j=0; j < Superclass::OutputImageDimension; j++)
            {
                if (j == i)
                {
                    outputDirection[j][i] = 1.0;
                }
                else
                {
                    outputDirection[j][i] = 0.0;
                }
            }
        }

        // set the spacing and origin
        outputPtr->SetSpacing( outputSpacing );
        outputPtr->SetOrigin( outputOrigin );
        outputPtr->SetDirection( outputDirection );
        outputPtr->SetNumberOfComponentsPerPixel( // propagate vector length info
                                                 inputPtr->GetNumberOfComponentsPerPixel());
    }
    else
    {
        // pointer could not be cast back down
        itkExceptionMacro(<< "itk::PixelMathImageFilter::GenerateOutputInformation "
                          << "cannot cast input to "
                          << typeid(ImageBase<Superclass::InputImageDimension>*).name() );
    }



    ComputeIntensityStatistics();
}


/**
 * ThreadedGenerateData Performs the pixel-wise addition
 */
template <class TInputImage, class TOutputImage >
void
PixelMathImageFilter<TInputImage,TOutputImage>
::ThreadedGenerateData( const OutputImageRegionType &outputRegionForThread,
                       ThreadIdType threadId)
{
    int nInputs = static_cast<int>(this->GetNumberOfInputs());
    InputImagePointer inputPtr[NUM_MAX_VARS];
    ImageRegionConstIterator<TInputImage>* inputIt[NUM_MAX_VARS];

    try {
        OutputImagePointer outputPtr = this->GetOutput(0);

        // Define the portion of the input to walk for this thread, using
        // the CallCopyOutputRegionToInputRegion method allows for the input
        // and output images to be different dimensions
        InputImageRegionType inputRegionForThread;
        this->CallCopyOutputRegionToInputRegion(inputRegionForThread, outputRegionForThread);

        // Define the iterators
        for (int i = 0; i < nInputs; i++) {
            inputPtr[i] = this->GetInput(i);
            inputIt[i] = new ImageRegionConstIterator<TInputImage>(inputPtr[i], inputRegionForThread);
            inputIt[i]->GoToBegin();
        }

        ImageRegionIterator<TOutputImage> outputIt(outputPtr, outputRegionForThread);
        ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());
        outputIt.GoToBegin();

		mu::Parser mathParser;
		double m_Vars[NUM_MAX_VARS + 2 * NUM_MAX_VARS];
        for (int i = 0; i < NUM_MAX_VARS; i++) {
            m_Vars[NUM_MAX_VARS + i * 2] = m_MinIntensities[i];
            m_Vars[NUM_MAX_VARS + i * 2 + 1] = m_MaxIntensities[i];
        }
		const char* VarNames[NUM_MAX_VARS + 2 * NUM_MAX_VARS] = { "A", "B", "C", "D", "E", "AMIN", "AMAX", "BMIN", "BMAX", "CMIN", "CMAX", "DMIN", "DMAX", "EMIN", "EMAX"  };
		for (int i = 0; i < NUM_MAX_VARS + 2 * NUM_MAX_VARS; i++) {
			mathParser.DefineVar(VarNames[i], &m_Vars[i]);
		}
		mathParser.SetExpr(m_Equation);
		while( !(inputIt[0]->IsAtEnd()) ) {
			for (int i = 0; i < nInputs; i++) {
				m_Vars[i] = static_cast<double>(inputIt[i]->Get());
				++(*inputIt[i]);
			}
			typename TOutputImage::PixelType o = 
            static_cast<typename TOutputImage::PixelType>(mathParser.Eval());
			outputIt.Set(o);
			++outputIt;
			progress.CompletedPixel();  // potential exception thrown here
		}
	} catch (mu::Parser::exception_type &e) {
		std::cout << e.GetMsg() << std::endl;
	} catch (itk::ExceptionObject& e) {
        std::cout << e.what() << endl;
    } catch (...) {
        std::cout << "Unknown error" << endl;
    }
	for (int i = 0;i < nInputs; i++) {
		delete inputIt[i];
	}
}


/**
 * Compute Intensity Statistics such as min and max
 */
template <class TInputImage, class TOutputImage >
void
PixelMathImageFilter<TInputImage,TOutputImage>::
ComputeIntensityStatistics() {
    const int nInputs = static_cast<int>(this->GetNumberOfInputs());

    m_MinIntensities.resize(nInputs);
    m_MaxIntensities.resize(nInputs);

    typedef itk::StatisticsImageFilter<TInputImage> StatisticsFilterType;

    for (int i = 0; i < nInputs; i++) {
        typename StatisticsFilterType::Pointer filter = StatisticsFilterType::New();
        filter->SetInput(this->GetInput(i));
        filter->Update();
        m_MinIntensities[i] = filter->GetMinimum();
        m_MaxIntensities[i] = filter->GetMaximum();
    }
}
#endif
