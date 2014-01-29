#include "piOptions.h"
#include "piImageIO.h"
#include "itkImage.h"
#include "PixelMathImageFilter.h"

using namespace std;
using namespace pi;

#ifndef PIXEL_TYPE
#define PIXEL_TYPE float
#endif

typedef itk::Image<PIXEL_TYPE,3> ImageType;
typedef itk::Image<PIXEL_TYPE,2> Image2D;
ImageIO<ImageType> io;
std::vector<ImageType::Pointer> inputImages;

int main(int argc, char* argv[]) {
    Options argParser;
    argParser.addOption("-e", "The equation to compute each output pixel.", "-e (A+B)", SO_REQ_SEP);
    argParser.addOption("-o", "output filename (the same data type with the last input)", "-o output.nrrd", SO_REQ_SEP);
    argParser.addOption("-h", "print this message", SO_NONE);

    StringVector args = argParser.ParseOptions(argc, argv, NULL);
    string eq = argParser.GetString("-e");
    string outputFilename = argParser.GetString("-o");

    if (argParser.GetBool("-h") || eq == "" || args.size() == 0 || outputFilename == "") {
        cout << "## kcalc usage \n"
            "\tkcalc [-e equation] [-o output-file] input1:A input2:B ...\n\n"
            "The kcalc performs a pixel-wise arithmetic. The pixel value of each input image is given as variables, A,B,C,D, and E. Several functions implemented in [MuParser](http://muparser.beltoforion.de/) includes +,-,*,/, and ? as well as trigonometric functions.\n\n"
            "Also, there are the min, max values of each input image for the use of scaling and other purposes, which are given as AMIN, AMAX, BMIN, BMAX, and etc.\n\n"
            "Note that the output data type is the same with the last input file. The order of images may produce different results, if images with different types are used.\n\n"
            "Some examples are:\n"
            "* **Addition**: kcalc -e \"(A+B)\" input1.nrrd input2.nrrd -o output.nrrd\n"
            "* **Averaging**: kcalc -e \"(A+B)/2\" input1.nrrd input2.nrrd -o output.nrrd\n"
            "* **Thresholding**: kcalc -e \"(A>10?1:0)\" input.nrrd -o output.nrrd\n"
            "* **Scaling**: -e (A-AMIN)/AMAX*255\n"
            "* **Masking**: -e (A==8?B:0)\n"
            "* ...\n\n"
            "### Options\n";
        argParser.PrintUsage();
        cout << endl;
        return 0;
    }

    if (argParser.GetBool("-2")) {
        cout << "Working on 2D images" << endl;
        ImageIO<Image2D> io2;
        ImageInfo lastImageInfo;
        PixelMathImageFilter<Image2D, Image2D>::Pointer pixelFilter = PixelMathImageFilter<Image2D, Image2D>::New();
        pixelFilter->SetEquation(eq);
        for (int i = 0; i < args.size(); i++) {
            pixelFilter->PushBackInput(io2.ReadCastedImage(args[i], lastImageInfo));
        }
        try {
            pixelFilter->Update();
        } catch (itk::ExceptionObject& e) {
            cout << e.what() << endl;
        }
        io2.WriteCastedImage(outputFilename, pixelFilter->GetOutput(), lastImageInfo.componenttype);
    } else {
        ImageInfo lastImageInfo;
        PixelMathImageFilter<ImageType, ImageType>::Pointer pixelFilter = PixelMathImageFilter<ImageType, ImageType>::New();
        pixelFilter->SetEquation(eq);
        for (int i = 0; i < args.size(); i++) {
            pixelFilter->PushBackInput(io.ReadCastedImage(args[i], lastImageInfo));
        }
        try {
            pixelFilter->Update();
        } catch (itk::ExceptionObject& e) {
            cout << e.what() << endl;
        }

        io.WriteCastedImage(outputFilename, pixelFilter->GetOutput(), lastImageInfo.componenttype);
    }
    return 0;
}
