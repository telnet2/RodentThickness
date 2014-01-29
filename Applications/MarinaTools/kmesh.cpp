
//
//  kmesh.cpp
//  ktools
//
//  Created by Joohwi Lee on 12/5/13.
//
//

#include "kmesh.h"

#include <set>
#include <iostream>

#include "piOptions.h"
#include "vtkio.h"

#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkIdList.h>
#include <vtkFloatArray.h>
#include <vtkMath.h>
#include <vtkAppendPolyData.h>
#include <vtkXMLImageDataReader.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkStreamTracer.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkLine.h>
#include <vtkPolyLine.h>
#include <vtkCleanPolyData.h>
#include <vtkMath.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkCellLocator.h>
#include <vtkKdTreePointLocator.h>
#include <vtkModifiedBSPTree.h>
#include <vtkCurvatures.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkMetaImageWriter.h>
#include <vtkImageStencil.h>


#include <itkImage.h>
#include <itkVectorNearestNeighborInterpolateImageFunction.h>
#include <itkEllipseSpatialObject.h>
#include <itkSpatialObjectToImageFilter.h>
#include <itkNearestNeighborInterpolateImageFunction.h>

#include "piImageIO.h"
#include "kimage.h"
#include "kstreamtracer.h"



#include <vnl/vnl_vector.h>


using namespace std;
using namespace pi;

/// @brief Convert a cartesian coordinate point to a spherical coordinate point
/// @param v a input point in the cartesian coordinate
/// @param phi output parameter for phi
/// @param theta output parameter for theta
static void cart2sph(const float *v, float *phi, float *theta) {
    // phi: azimuth, theta: elevation
    float d = v[0] * v[0] + v[1] * v[1];
    *phi = (d == 0) ? 0: atan2(v[1], v[0]);
    *theta = (v[2] == 0) ? 0: atan2(v[2], sqrt(d));
}


/// @brief Compute the factorial from x to y
static double factorial(double x, double y) {
    double f = 1;
    for (; x >= y; x--) {
        f *= x;
    }
    return f;
}

/// @brief Compute the basis function for spherical harmonics
static void basis(int degree, float *p, float *Y) {
    // real spherical harmonics basis functions
    // polar coordinate
    float phi, theta;
    cart2sph(p, &phi, &theta);
    theta = M_PI_2 - theta;  // convert to interval [0, PI]
    float *Pm = new float[degree + 1];

    // square root of 2
    const float sqr2 = sqrt(2.0f);

    const int yCount = (degree - 1) * (degree - 1);

    for (int l = 0; l <= degree; l++)
    {
        // legendre part
//        Series::legendre(l, cos(theta), Pm);
        float lconstant = sqrt((2 * l + 1) / (4 * M_PI));

        int center = (l + 1) * (l + 1) - l - 1;

        Y[center] = lconstant * Pm[0];

        for (int m = 1; m <= l; m++)
        {
            const double f = factorial(l+m, l-m+1);
            float precoeff = lconstant * (float)sqrt(1 / f);

            if (m % 2 == 1) precoeff = -precoeff;
            Y[center + m] = sqr2 * precoeff * Pm[m] * cos(m * phi);
            Y[center - m] = sqr2 * precoeff * Pm[m] * sin(m * phi);
        }
    }

    delete [] Pm;
}


// append polydatas into one
void runAppendData(Options& opts, StringVector& args) {
    vtkAppendPolyData* appender = vtkAppendPolyData::New();
    vtkIO io;

    // read all files
    for (int i = 0; i < args.size(); i++) {
        vtkPolyData* data = io.readFile(args[i]);
        appender->AddInput(data);
    }

    appender->Update();
    io.writeFile(opts.GetString("-appendData"), appender->GetOutput());

    return;
}

// add scalar value to a mesh
void runImportScalars(Options& opts, StringVector& args) {
    cout << "importing scalars from " << args[0] << " to " << args[2] << endl;
    vtkIO io;

    // read polydata
    vtkPolyData* poly = io.readFile(args[0]);
    if (poly == NULL) {
        cout << "can't read " << args[0] << endl;
        return;
    }
    vtkFloatArray* scalar = vtkFloatArray::New();
    scalar->SetNumberOfValues(poly->GetNumberOfPoints());
    scalar->SetName(opts.GetString("-scalarName").c_str());

    ifstream file(args[1].c_str());
    for (int i = 0; i < poly->GetNumberOfPoints() && !file.eof(); i++) {
        string line;
        file >> line;
        cout << line << endl;
        if (line == "nan") {
            scalar->SetValue(i, NAN);
        } else {
            scalar->SetValue(i, atof(line.c_str()));
        }

    }

    poly->GetPointData()->AddArray(scalar);

    io.writeFile(args[2], poly);
}


/// @brief export scalar values to a text file
void runExportScalars(Options& opts, StringVector& args) {
    vtkIO io;

    /// - Read polydata
    vtkPolyData* poly = io.readFile(args[0]);

    /// - Check if file is loaded
    if (poly == NULL) {
        cout << "can't read file: " << args[0];
        return;
    }

    bool isPointData = true;
    /// - Find the scalar attribute given as '-scalarName'
    vtkDataArray* scalar = poly->GetPointData()->GetScalars(opts.GetString("-scalarName").c_str());
    /// - Check if the scalar exists
    if (scalar == NULL) {
        /// - Try with cell arrays
        scalar = poly->GetCellData()->GetScalars(opts.GetString("-scalarName").c_str());
        if (scalar == NULL) {
            cout << "can't find the scalar attribute: " << opts.GetString("-scalarName") << endl;
            return;
        }
        isPointData = false;
    }

    ofstream file(args[1].c_str());

    int nScalars = scalar->GetNumberOfTuples();
    for (int i = 0; i < nScalars; i++) {
        file << scalar->GetTuple1(i) << endl;
    }
    file.close();
}


/// @brief Copy a scalar list to another object
void runCopyScalars(Options& opts, StringVector& args) {
    vtkIO vio;
    string inputModelFile = args[0];
    string inputModelFile2 = args[1];

    vtkPolyData* inputModel = vio.readFile(inputModelFile);
    vtkPolyData* inputModel2 = vio.readFile(inputModelFile2);
    vtkDataArray* scalars = inputModel->GetPointData()->GetScalars(opts.GetString("-scalarName").c_str());
    inputModel2->GetPointData()->AddArray(scalars);
    vio.writeFile(args[2], inputModel2);
}



/// @brief Perform smoothing on manifold by iterative averaging as used in FreeSurfer
void runScalarSmoothing(Options& opts, StringVector& args) {
    // number of iterations and sigma affects the smoothing results
    vtkIO io;
    vtkMath* math = vtkMath::New();

    // sigma for smoothing
    double sigma2 = opts.GetStringAsReal("-sigma", 1);
    sigma2 *= sigma2;

    // number of iterations
    int numIters = opts.GetStringAsInt("-iter", 1);


    // read polydata
    vtkPolyData* poly = io.readFile(args[0]);

    // access data
    string scalarName = opts.GetString("-scalarName");
    vtkDataArray* scalars = poly->GetPointData()->GetScalars(scalarName.c_str());
    if (scalars == NULL) {
        cout << "can't find scalars: " << scalarName << endl;
        return;
    }

    // copy the scalars to iteratively apply smoothing
    vtkFloatArray* data = vtkFloatArray::New();
    data->DeepCopy(scalars);

    // prepare new data array
    vtkFloatArray* newData = vtkFloatArray::New();

    string outputScalarName = opts.GetString("-outputScalarName", "smoothed_" + scalarName);
    newData->SetName(outputScalarName.c_str());
    newData->SetNumberOfTuples(data->GetNumberOfTuples());
    poly->GetPointData()->AddArray(newData);


    // check if the scalar array exists
    if (data == NULL) {
        cout << "can't access scalar array: " << scalarName << endl;
        return;
    }

    // iterate over all points
    vtkIdList* cellIds = vtkIdList::New();
    vtkIdList* ptIds = vtkIdList::New();
    std::set<int> ptSet;

    // build cells
    poly->BuildCells();
    poly->BuildLinks();

    for (int n = 0; n < numIters; n++) {
        cout << "Iter: " << n << endl;
        for (int i = 0; i < poly->GetNumberOfPoints(); i++) {
            double center[3];
            poly->GetPoint(i, center);

            // collect neighbor cells
            ptSet.clear();

            cellIds->Reset();
            poly->GetPointCells(i, cellIds);

            // iterate over neighbor cells
            for (int j = 0; j < cellIds->GetNumberOfIds(); j++) {
                int cellId = cellIds->GetId(j);
                ptIds->Reset();

                // collect cell points
                poly->GetCellPoints(cellId, ptIds);

                // iterate over all cell points
                for (int k = 0; k < ptIds->GetNumberOfIds(); k++) {
                    int ptId = ptIds->GetId(k);
                    ptSet.insert(ptId);
                }
            }

            // iterate over all neighbor points
            std::set<int>::iterator iter = ptSet.begin();

            // compute weight
            vnl_vector<float> weights;
            weights.set_size(ptSet.size());

            for (int j = 0; iter != ptSet.end(); iter++, j++) {
                int ptId = *iter;
                double neighbor[3];
                poly->GetPoint(ptId, neighbor);

                double dist2 = math->Distance2BetweenPoints(center, neighbor);

                // apply the heat kernel with the sigma
                weights[j] = exp(-dist2/sigma2);
            }

            // add one for the center
            double weightSum = weights.sum() + 1;
            weights /= weightSum;


            // iterate over neighbors and compute weighted sum
            double smoothedValue = data->GetTuple1(i) / weightSum;
            iter = ptSet.begin();

            // compute the weighted averge
            for (uint j = 0; j < ptSet.size(); j++, iter++) {
                int ptId = *iter;
                int value = data->GetTuple1(ptId);
                smoothedValue += (value * weights[j]);
            }
            newData->SetTuple1(i, smoothedValue);
        }


        // prepare next iteration by copying newdata to data
        data->DeepCopy(newData);
    }

    // write to file
    io.writeFile(args[1], poly);
}



/// @brief Convert an ITK image to a VTKImageData
void runConvertITK2VTI(Options& opts, StringVector& args) {
    if (args.size() < 2) {
        cout << "requires input-image-file and output-vti-file" << endl;
        return;
    }

    int attrDim = opts.GetStringAsInt("-attrDim", 1);
    string scalarName = opts.GetString("-scalarName", "Intensity");
    string maskImageFile = opts.GetString("-maskImage");

    MaskImageType::Pointer maskImage;
    if (maskImageFile != "") {
        ImageIO<MaskImageType> io;
        maskImage = io.ReadCastedImage(maskImageFile);
    }

    string input = args[0];
    string output = args[1];

    /// - Read an image data
    vtkImageData* outputData = vtkImageData::New();
    if (attrDim == 1) {
        ConvertImageT<ImageType>(input, outputData, scalarName.c_str(), 1, maskImage);
    } else if (attrDim == 3) {
        ConvertImageT<VectorImageType>(input, outputData, scalarName.c_str(), attrDim, maskImage);
    }

    vtkXMLImageDataWriter* w = vtkXMLImageDataWriter::New();
    w->SetFileName(output.c_str());
    w->SetDataModeToAppended();
    w->EncodeAppendedDataOff();
    w->SetCompressorTypeToZLib();
    w->SetDataModeToBinary();

    w->SetInput(outputData);
    w->Write();
}


/// @brief Convert an itk image file to vtkUnstructuredGrid
void runConvertITK2VTU(Options& opts, StringVector& args) {
    if (args.size() < 2) {
        cout << "requires input-image-file and output-vti-file" << endl;
        return;
    }

    int attrDim = opts.GetStringAsInt("-attrDim", 1);
    string scalarName = opts.GetString("-scalarName", "Intensity");

    string input = args[0];
    string output = args[1];

    string maskImageFile = opts.GetString("-maskImage");

    typedef itk::Image<ushort,3> MaskImageType;
    ImageIO<MaskImageType> maskIO;
    MaskImageType::Pointer maskImage = maskIO.ReadCastedImage(maskImageFile);

    /// - Read an image data
    vtkUnstructuredGrid* outputData = vtkUnstructuredGrid::New();
    if (attrDim == 1) {
        ConvertImageT<ImageType, MaskImageType>(input, outputData, maskImage, scalarName.c_str(), 1);
    } else if (attrDim == 3) {
        ConvertVectorImageT<VectorImageType, MaskImageType>(input, outputData, maskImage, scalarName.c_str(), attrDim);
    }

    vtkXMLUnstructuredGridWriter* w = vtkXMLUnstructuredGridWriter::New();
    w->SetFileName(output.c_str());
    w->SetDataModeToAppended();
    w->EncodeAppendedDataOff();
    w->SetCompressorTypeToZLib();
    w->SetDataModeToBinary();

    w->SetInput(outputData);
    w->Write();
}

bool endswith(std::string str, std::string substr) {
    size_t i = str.rfind(substr);
    return (i != string::npos) && (i == (str.length() - substr.length()));
}


/// @brief perform a line clipping to fit within the object
bool performLineClipping(vtkPolyData* streamLines, vtkModifiedBSPTree* tree, int lineId, vtkCell* lineToClip, vtkPolyData* object, vtkPoints* outputPoints, vtkCellArray* outputLines, double &length) {
    /// - Iterate over all points in a line
    vtkIdList* ids = lineToClip->GetPointIds();
    /// - Identify a line segment included in the line

    int nIntersections = 0;
    bool foundEndpoint = false;
    std::vector<vtkIdType> idList;
    for (int j = 2; j < ids->GetNumberOfIds(); j++) {
        double p1[3], p2[3];
        streamLines->GetPoint(ids->GetId(j-1), p1);
        streamLines->GetPoint(ids->GetId(j), p2);

        // handle initial condition
        if (j == 2) {
            double p0[3];
            streamLines->GetPoint(ids->GetId(0), p0);
            idList.push_back(outputPoints->GetNumberOfPoints());
            outputPoints->InsertNextPoint(p0);

            idList.push_back(outputPoints->GetNumberOfPoints());
            outputPoints->InsertNextPoint(p1);

            length = sqrt(vtkMath::Distance2BetweenPoints(p0, p1));
        }

        int subId;
        double x[3] = {-1,-1,-1};
        double t = 0;

        double pcoords[3] = { -1, };
        int testLine = tree->IntersectWithLine(p1, p2, 0.01, t, x, pcoords, subId);
        if (testLine) {
            nIntersections ++;
            if (nIntersections > 0) {
                idList.push_back(outputPoints->GetNumberOfPoints());
                outputPoints->InsertNextPoint(x);
                length += sqrt(vtkMath::Distance2BetweenPoints(p1, x));
                foundEndpoint = true;
                break;
            }
        }
//        cout << testLine << "; " << x[0] << "," << x[1] << "," << x[2] << endl;


        idList.push_back(outputPoints->GetNumberOfPoints());
        outputPoints->InsertNextPoint(p2);
        length += sqrt(vtkMath::Distance2BetweenPoints(p1, p2));
    }

    if (foundEndpoint) {
        outputLines->InsertNextCell(idList.size(), &idList[0]);
        return true;
    }
    return false;
}


/// @brief Perform a line clipping task
void runTraceClipping(Options& opts, StringVector& args) {
    string inputStreamsFile = args[0];
    string inputObjectFile = args[1];
    string outputStreamsFile = args[2];

    vtkIO vio;
    vtkPolyData* inputStream = vio.readFile(inputStreamsFile);
    vtkPolyData* inputObject = vio.readFile(inputObjectFile);
    vtkPolyData* outputObject = vtkPolyData::New();


    vtkCellArray* lines = inputStream->GetLines();
    vtkModifiedBSPTree* tree = vtkModifiedBSPTree::New();
    tree->SetDataSet(inputObject);
    tree->BuildLocator();

    vtkPoints* outputPoints = vtkPoints::New();
    vtkCellArray* outputLines = vtkCellArray::New();

    for (int i = 0; i < lines->GetNumberOfCells(); i++) {
        vtkCell* line = inputStream->GetCell(i);
        double length = 0;
        performLineClipping(inputStream, tree, i, line, inputObject, outputPoints, outputLines, length);
    }
    vio.writeFile("test.vtp", outputObject);
}

/// @brief Execute the stream tracer
void runStreamTracer(Options& opts, StringVector& args) {
    string inputVTUFile = args[0];
    string inputPointsFile = args[1];
    string outputStreamFile = args[2];
    string outputPointFile = args[3];
    bool zRotate = opts.GetBool("-zrotate", false);


    vtkDataSet* inputData;

    // FIXME - create a dataset reader
    if (endswith(inputVTUFile, string(".vtu"))) {
        vtkXMLUnstructuredGridReader* reader = vtkXMLUnstructuredGridReader::New();
        reader->SetFileName(inputVTUFile.c_str());
        reader->Update();

        vtkUnstructuredGrid* inputVTU = reader->GetOutput();
        inputData = inputVTU;
    } else if (endswith(inputVTUFile, ".vti")) {
        vtkXMLImageDataReader* reader = vtkXMLImageDataReader::New();
        reader->SetFileName(inputVTUFile.c_str());
        reader->Update();

        vtkImageData* inputVTI = reader->GetOutput();
        inputData = inputVTI;
    }

    vtkIO vio;
    vtkPolyData* inputPoints = vio.readFile(inputPointsFile);
    vtkPoints* points = inputPoints->GetPoints();

    /// - Converting the input points to the image coordinate
    const int nInputPoints = inputPoints->GetNumberOfPoints();
    for (int i = 0; i < nInputPoints; i++) {
        double p[3];
        points->GetPoint(i, p);
        // FixMe: Do not use a specific scaling factor
        if (zRotate) {
            p[0] = -p[0];
            p[1] = -p[1];
            p[2] = p[2];
        }
        points->SetPoint(i, p);
    }
    inputPoints->SetPoints(points);

    /// - Set up tracer (Use RK45, both direction, initial step 0.05, maximum propagation 500
    StreamTracer* tracer = StreamTracer::New();
    tracer->SetInput(inputData);
    tracer->SetSource(inputPoints);
//    double seedPoint[3];
 //   inputPoints->GetPoint(24745, seedPoint);
 //   tracer->SetStartPosition(seedPoint);
    tracer->SetIntegratorTypeToRungeKutta45();


    bool isBothDirection = false;
    if (opts.GetString("-traceDirection") == "both") {
        tracer->SetIntegrationDirectionToBoth();
        isBothDirection = true;
        cout << "Forward/Backward Integration" << endl;
    } else if (opts.GetString("-traceDirection") == "backward") {
        tracer->SetIntegrationDirectionToBackward();
        cout << "Backward Integration" << endl;
    } else {
        tracer->SetIntegrationDirectionToForward();
        cout << "Forward Integration" << endl;
    }

    tracer->SetInterpolatorTypeToDataSetPointLocator();
    tracer->SetMaximumPropagation(500);
    tracer->SetInitialIntegrationStep(0.05);
    tracer->Update();



    vtkPolyData* streamLines = tracer->GetOutput();

    // loop over the cell and compute the length
    int nCells = streamLines->GetNumberOfCells();
    cout << "# of cells: " << nCells << endl;


    /// - Prepare the output as a scalar array
//    vtkDataArray* streamLineLength = streamLines->GetCellData()->GetScalars("Length");

    /// - Prepare the output for the input points
    vtkDoubleArray* streamLineLengthPerPoint = vtkDoubleArray::New();
    streamLineLengthPerPoint->SetNumberOfTuples(nInputPoints);
    streamLineLengthPerPoint->SetName("Length");
    streamLineLengthPerPoint->SetNumberOfComponents(1);
    streamLineLengthPerPoint->FillComponent(0, 0);

    vtkIntArray* lineCorrect = vtkIntArray::New();
    lineCorrect->SetName("LineOK");
    lineCorrect->SetNumberOfValues(nInputPoints);
    lineCorrect->FillComponent(0, 0);

    inputPoints->GetPointData()->SetScalars(streamLineLengthPerPoint);
    inputPoints->GetPointData()->AddArray(lineCorrect);

    cout << "Assigning a length to each source vertex ..." << endl;
    vtkDataArray* seedIds = streamLines->GetCellData()->GetScalars("SeedId");
    if (seedIds) {
        // line clipping
        vtkPoints* outputPoints = vtkPoints::New();
        vtkCellArray* outputCells = vtkCellArray::New();

        /// construct a tree locator
        vtkModifiedBSPTree* tree = vtkModifiedBSPTree::New();
        tree->SetDataSet(inputPoints);
        tree->BuildLocator();


        vtkDoubleArray* lengthArray = vtkDoubleArray::New();
        lengthArray->SetName("Length");

        vtkIntArray* pointIds = vtkIntArray::New();
        pointIds->SetName("PointIds");

        for (int i = 0; i < nCells; i++) {
            int pid = seedIds->GetTuple1(i);
            double length = 0;
            if (pid > -1) {
                vtkCell* line = streamLines->GetCell(i);
                /// - Assume that a line starts from a point on the input mesh and must meet at the opposite surface of the starting point.
                bool lineAdded = performLineClipping(streamLines, tree, i, line, inputPoints, outputPoints, outputCells, length);

                if (lineAdded) {
                    pointIds->InsertNextValue(pid);
                    lengthArray->InsertNextValue(length);
                    streamLineLengthPerPoint->SetValue(pid, length);
                    lineCorrect->SetValue(pid, 1);
                } else {
                    lineCorrect->SetValue(pid, 2);
                }
            }
        }

        vtkPolyData* outputStreamLines = vtkPolyData::New();
        outputStreamLines->SetPoints(outputPoints);
        outputStreamLines->SetLines(outputCells);
        outputStreamLines->GetCellData()->AddArray(pointIds);
        outputStreamLines->GetCellData()->AddArray(lengthArray);


        vtkCleanPolyData* cleaner = vtkCleanPolyData::New();
        cleaner->SetInput(outputStreamLines);
        cleaner->Update();
        vio.writeFile(outputStreamFile, cleaner->GetOutput());
    } else {
        cout << "Can't find SeedId" << endl;
    }

    cout << lineCorrect->GetNumberOfTuples() << endl;
    cout << streamLineLengthPerPoint->GetNumberOfTuples() << endl;
    vio.writeFile(outputPointFile, inputPoints);
//    vio.writeXMLFile(outputVTKFile, streamLines);
}


/// @brief Copy a scalar list from a seed object to a stream line object
void runTraceScalarCombine(Options& opts, StringVector& args) {
    if (args.size() < 3) {
        cout << "requires input-seed input-stream output-stream-file" << endl;
        return;
    }

    string inputSeedFile = args[0];
    string inputStreamFile = args[1];
    string outputStreamFile = args[2];
    string scalarName = opts.GetString("-scalarName");

    if (scalarName == "") {
        cout << "requires -scalarName scalarName" << endl;
        return;
    }

    vtkIO vio;
    vtkPolyData* inputSeed = vio.readFile(inputSeedFile);
    vtkPolyData* inputStream = vio.readFile(inputStreamFile);

    vtkDataArray* pointIds = inputStream->GetCellData()->GetScalars("PointIds");
    if (pointIds == NULL) {
        cout << "Can't find PointIds" << endl;
        return;
    }
    vtkDataArray* scalars = inputSeed->GetPointData()->GetScalars(scalarName.c_str());
    if (scalars == NULL) {
        cout << "Can't find scalars: " << scalarName << endl;
        return;
    }

    vtkDoubleArray* outputScalars = vtkDoubleArray::New();
    outputScalars->SetName(scalarName.c_str());
    for (int i = 0; i < pointIds->GetNumberOfTuples(); i++) {
        int ptId = pointIds->GetTuple1(i);
        double value = scalars->GetTuple1(ptId);
        outputScalars->InsertNextTuple1(value);
    }
    inputStream->GetCellData()->AddArray(outputScalars);

    if (opts.GetBool("-zrotate")) {
        cout << "The output is rotated!" << endl;
        vio.zrotate(inputStream);
    }
    vio.writeFile(outputStreamFile, inputStream);
}


/// @brief Apply a filter to each stream line
void runFilterStream(Options& opts, StringVector& args) {
    string inputStream = args[0];
    string inputSeeds = args[1];
    string outputStreamFile = args[2];
    string scalarName = opts.GetString("-scalarName");

    double lowThreshold = opts.GetStringAsReal("-thresholdMin", itk::NumericTraits<float>::min());
    double highThreshold = opts.GetStringAsReal("-thresholdMax", itk::NumericTraits<float>::max());

    vtkIO vio;
    vtkPolyData* streamLines = vio.readFile(inputStream);
    vtkPolyData* streamSeeds = vio.readFile(inputSeeds);

    vtkPolyData* outputStream = vtkPolyData::New();
    outputStream->SetPoints(streamLines->GetPoints());

    /// - Lookup SeedId and a given scalar array
    vtkDataArray* seedList = streamLines->GetCellData()->GetScalars("SeedId");
    vtkDataArray* seedScalars = streamSeeds->GetPointData()->GetScalars(scalarName.c_str());
    vtkCellArray* lines = vtkCellArray::New();
    vtkDoubleArray* filteredScalars = vtkDoubleArray::New();
    filteredScalars->SetName(scalarName.c_str());
    filteredScalars->SetNumberOfComponents(1);

    /// - Lookup a corresponding point scalar, apply threashold filter, and add to the new object
    for (int i = 0; i < seedList->GetNumberOfTuples(); i++) {
        int seedId = seedList->GetTuple1(i);
        double value = seedScalars->GetTuple1(seedId);

        if (lowThreshold <= value && value <= highThreshold) {
            lines->InsertNextCell(streamLines->GetCell(i));
            filteredScalars->InsertNextValue(value);
        }
    }

    outputStream->SetLines(lines);
    outputStream->GetCellData()->AddArray(filteredScalars);
    outputStream->BuildCells();
    outputStream->BuildLinks();

    vio.writeFile(outputStreamFile, outputStream);
}


/// @brief Fit a model into a binary image
void runFittingModel(Options& opts, StringVector& args) {
    if (args.size() < 3) {
        cout << "requires input-model input-image output-model" << endl;
        return;
    }
    string inputModelFile = args[0];
    string inputImageFile = args[1];
    string outputModelFile = args[2];

    vtkIO vio;
    vtkPolyData* inputModel = vio.readFile(inputModelFile);
    const int nPoints = inputModel->GetNumberOfPoints();

    /// Apply z-rotation
    for (int i = 0; i < nPoints; i++) {
        double point[3];
        inputModel->GetPoint(i, point);
        if (opts.GetBool("-zrotate")) {
            point[0] = -point[0];
            point[1] = -point[1];
        }
        inputModel->GetPoints()->SetPoint(i, point);
    }

    ImageIO<MaskImageType> itkIO;
    MaskImageType::Pointer maskImage = itkIO.ReadCastedImage(inputImageFile);

    // for test, rescale 10 times
    VectorImageType::SpacingType spacing = maskImage->GetSpacing();
    for (int i = 0; i < 3; i++) {
//        spacing[i] *= 10;
    }
    maskImage->SetSpacing(spacing);


    // test
    VectorImageType::Pointer distImage = ComputeDistanceMap(maskImage);
    typedef itk::VectorNearestNeighborInterpolateImageFunction<VectorImageType> InterpolatorType;

    ImageIO<VectorImageType> distIO;
    distIO.WriteImage("dist.mha", distImage);

    InterpolatorType::Pointer distInterp = InterpolatorType::New();
    distInterp->SetInputImage(distImage);


    // iterate over the input model and project to the boundary
    const int nIters = 50;
    for (int i = 0; i < nIters; i++) {
        // Compute laplacian smoothing by taking iterative average

        vtkSmoothPolyDataFilter* filter = vtkSmoothPolyDataFilter::New();
        filter->SetInput(inputModel);
        filter->SetNumberOfIterations(1);
        filter->Update();
        inputModel = filter->GetOutput();

        // projection
        for (int j = 0; j < nPoints; j++) {
            VectorImageType::PointType point, nextPoint;
            inputModel->GetPoint(j, point.GetDataPointer());
            VectorType offset = distInterp->Evaluate(point);
            for (int k = 0; k < 3; k++) {
                nextPoint = point + offset[k] * spacing[k] * spacing[k];
            }
            cout << point << " => " << nextPoint << endl;
            inputModel->GetPoints()->SetPoint(j, nextPoint.GetDataPointer());
        }
    }

    vio.writeFile(outputModelFile, inputModel);
}


MaskImageType::Pointer Ellipse(int* outputSize, double *center, double *radius) {
    ImageType::SizeType size;    typedef itk::EllipseSpatialObject<3> EllipseType;
    typedef itk::SpatialObjectToImageFilter<EllipseType, MaskImageType> SpatialObjectToImageFilterType;

    SpatialObjectToImageFilterType::Pointer imageFilter = SpatialObjectToImageFilterType::New();

    for (int k = 0; k < 3; k++) {
        size[k] = outputSize[k];
    }
    imageFilter->SetSize(size);

    EllipseType::Pointer ellipse = EllipseType::New();
    ellipse->SetDefaultInsideValue(255);
    ellipse->SetDefaultOutsideValue(0);

    EllipseType::ArrayType axes;
    for (int k = 0; k < 3; k++) {
        axes[k] = radius[k];
    }
    ellipse->SetRadius(axes);

    EllipseType::TransformType::Pointer transform = EllipseType::TransformType::New();
    transform->SetIdentity();
    EllipseType::TransformType::OutputVectorType translation;
    for (int k = 0; k < 3; k++) {
        translation[k] = center[k];
    }
    transform->Translate(translation, false);

    ellipse->SetObjectToParentTransform(transform);
    imageFilter->SetInput(ellipse);
    imageFilter->SetUseObjectValue(true);
    imageFilter->SetOutsideValue(0);
    imageFilter->Update();
    return imageFilter->GetOutput();
}


/// @brief Create an ellipse binary image
void runEllipse(pi::Options &opts, StringVector &args) {
    if (!opts.GetBool("--ellipse")) {
        return;
    }

    if (args.size() < 3 * 3) {
        cout << "--ellipse output-image [image-size] [ellipse-center] [ellipse-radius] " << endl;
        exit(EXIT_FAILURE);
    }

    int size[3];
    double center[3], radius[3];
    for (int k = 0; k < 3; k++) {
        size[k] = atoi(args[k].c_str());
        center[k] = atof(args[3*1 + k].c_str());
        radius[k] = atof(args[3*2 + k].c_str());
    }

    MaskImageType::Pointer outputImage = Ellipse(size, center, radius);
    ImageIO<MaskImageType> io;
    string outputImageFile = opts.GetString("-o");
    io.WriteImage(outputImageFile, outputImage);

    exit(EXIT_SUCCESS);
}

/// @brief Sample pixel values from an image for a input model
void runSampleImage(Options& opts, StringVector& args) {
    vtkIO vio;

    string inputImageFile = args[0];
    string inputModelFile = args[1];
    string outputModelFile = args[2];

    ImageIO<ImageType> io;
    ImageType::Pointer inputImage = io.ReadCastedImage(inputImageFile);
    typedef itk::NearestNeighborInterpolateImageFunction<ImageType> ImageInterpolatorType;
    ImageInterpolatorType::Pointer interp = ImageInterpolatorType::New();
    interp->SetInputImage(inputImage);

    cout << "Building a image data ..." << flush;
    /// Create a vtu image
    /// - Create an instance for the output grid
    vtkUnstructuredGrid* imageData = vtkUnstructuredGrid::New();
    vtkPointData* pdata = imageData->GetPointData();

    /// - Create a point set to store valid points
    vtkPoints* pointSet = vtkPoints::New();

    /// - Create an array to store the pixel data
    vtkDoubleArray* attr = vtkDoubleArray::New();
    attr->SetNumberOfComponents(1);
    attr->SetName("Pixels");

    /// - Loop over the entire pixel of the mask
    itk::ImageRegionIteratorWithIndex<ImageType> iter(inputImage, inputImage->GetBufferedRegion());
    for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter) {
        double pixel = iter.Get();
        /// - Only sample non-negative pixels
        if (pixel > 0) {
            ImageType::PointType point;
            inputImage->TransformIndexToPhysicalPoint(iter.GetIndex(), point);

            /// - Add a point
            pointSet->InsertNextPoint(point[0], point[1], point[2]);

            /// - Add a pixel value (a scalar or a vector)
            attr->InsertNextValue(pixel);
        }
    }
    imageData->SetPoints(pointSet);
    pdata->AddArray(attr);

    cout << " done" << endl;

    /// FIXME: The number of points of the input should be greater than 0.
    vtkKdTreePointLocator* locator = vtkKdTreePointLocator::New();
    locator->SetDataSet(imageData);
    locator->BuildLocator();


    /// Now sample from the vtu
    vtkPolyData* inputModel = vio.readFile(inputModelFile);

    if (inputModel == NULL) {
        cout << "Can't read " << inputModelFile << endl;
        return;
    }

    vtkPoints* inputPoints = inputModel->GetPoints();
    const int nPoints = inputPoints->GetNumberOfPoints();

    vtkDoubleArray* pixels = vtkDoubleArray::New();
    pixels->SetNumberOfTuples(nPoints);
    pixels->SetName(opts.GetString("-outputScalarName", "PixelValue").c_str());
    pixels->SetNumberOfComponents(1);


    for (int i = 0; i < nPoints; i++) {
        ImageType::PointType p;
        inputPoints->GetPoint(i, p.GetDataPointer());

        if (opts.GetBool("-zrotate")) {
            p[0] = -p[0];
            p[1] = -p[1];
        }

        int pid = locator->FindClosestPoint(p.GetDataPointer());
        double pixel = attr->GetValue(pid);
        pixels->SetValue(i, pixel);

//        if (interp->IsInsideBuffer(p)) {
//            ImageType::PixelType pixel = interp->Evaluate(p);
//            pixels->SetValue(i, pixel);
//        }
    }

    inputModel->GetPointData()->AddArray(pixels);
    vio.writeFile(outputModelFile, inputModel);
}


/// @brief Compute the mean and Gaussian curvature for each point
void runComputeCurvature(Options& opts, StringVector& args) {
    vtkIO vio;

    string inputModelFile = args[0];
    string outputModelFile = args[1];

    /// Now sample from the vtu
    vtkPolyData* inputModel = vio.readFile(inputModelFile);

    if (inputModel == NULL) {
        cout << "Can't read " << inputModelFile << endl;
        return;
    }



    vtkCurvatures* curv1 = vtkCurvatures::New();
    curv1->SetInput(inputModel);
    curv1->SetCurvatureTypeToGaussian();
    curv1->Update();
    vtkPolyData* gx = curv1->GetOutput();
    inputModel->GetPointData()->AddArray(gx->GetPointData()->GetScalars("GAUSS_Curvature"));

    vtkCurvatures* curv2 = vtkCurvatures::New();
    curv2->SetInput(inputModel);
    curv2->SetCurvatureTypeToMean();
    curv2->Update();
    vtkPolyData* mx = curv2->GetOutput();
    inputModel->GetPointData()->AddArray(mx->GetPointData()->GetScalars("Mean_Curvature"));

    vtkCurvatures* curv3 = vtkCurvatures::New();
    curv3->SetInput(inputModel);
    curv3->SetCurvatureTypeToMaximum();
    curv3->Update();
    vtkPolyData* mx2 = curv3->GetOutput();
    inputModel->GetPointData()->AddArray(mx2->GetPointData()->GetScalars("Maximum_Curvature"));

    vtkCurvatures* curv4 = vtkCurvatures::New();
    curv4->SetInput(inputModel);
    curv4->SetCurvatureTypeToMinimum();
    curv4->Update();
    vtkPolyData* mx3 = curv4->GetOutput();
    inputModel->GetPointData()->AddArray(mx3->GetPointData()->GetScalars("Minimum_Curvature"));

    vio.writeFile(outputModelFile, mx3);
}


/// @brief Compute the average of scalars
void runAverageScalars(Options& opts, StringVector& args) {
    vtkIO vio;
    string outputFile = opts.GetString("-o");

    std::vector<vtkPolyData*> inputs;
    inputs.resize(args.size());

    inputs[0] = vio.readFile(args[0]);

    vtkDoubleArray* scalars = vtkDoubleArray::New();
    scalars->SetName(opts.GetString("-outputScalarName").c_str());
    scalars->SetNumberOfValues(inputs[0]->GetNumberOfPoints());

    for (int i = 0; i < args.size(); i++) {
        if (i == 0) {
            vtkDataArray* inputScalars = inputs[i]->GetPointData()->GetScalars(opts.GetString("-scalarName").c_str());
            for (int j = 0; j < scalars->GetNumberOfTuples(); j++) {
                scalars->SetValue(j, inputScalars->GetTuple1(j));
            }
        } else {
            inputs[i] = vio.readFile(args[i]);
            vtkDataArray* inputScalars = inputs[i]->GetPointData()->GetScalars(opts.GetString("-scalarName").c_str());
            for (int j = 0; j < scalars->GetNumberOfTuples(); j++) {
                scalars->SetValue(j, scalars->GetValue(j) + inputScalars->GetTuple1(j));
            }
        }
    }

    const double thresholdValue = opts.GetStringAsReal("-threshold", -1);
    const bool useThreshold =  thresholdValue != -1;

    cout << "Thresholding at " << thresholdValue << endl;
    for (int j = 0; j < scalars->GetNumberOfTuples(); j++) {
        double v = scalars->GetValue(j) / args.size();
        if (useThreshold) {
            scalars->SetValue(j, v < thresholdValue ? 1 : 2);
        } else {
            scalars->SetValue(j, v);
        }
    }

    for (int i = 0; i < args.size(); i++) {
        inputs[i]->GetPointData()->AddArray(scalars);
        cout << "Writing " << args[i] << endl;
        vio.writeFile(args[i], inputs[i]);
    }
}

/// @brief Compute the voronoi image from a surface model
void runVoronoiImage(Options& opts, StringVector& args) {
    string inputImageFile = args[0];
    string inputModelFile = args[1];
    string outputImageFile = args[2];

    ImageIO<MaskImageType> io;
    MaskImageType::Pointer maskImage = io.ReadCastedImage(inputImageFile);
    if (maskImage.IsNull()) {
        cout << "Can't read " << inputImageFile << endl;
        return;
    }
    itk::ImageRegionIteratorWithIndex<MaskImageType> iter(maskImage, maskImage->GetBufferedRegion());

    vtkIO vio;
    vtkPolyData* inputModel = vio.readFile(inputModelFile);
    vtkKdTreePointLocator* locator = vtkKdTreePointLocator::New();
    locator->SetDataSet(inputModel);
    locator->BuildLocator();

    vtkDataArray* scalars = inputModel->GetPointData()->GetScalars(opts.GetString("-scalarName").c_str());

    const bool zrotate = opts.GetBool("-zrotate");
    for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter) {
        MaskImageType::PointType voxelPoint;
        maskImage->TransformIndexToPhysicalPoint(iter.GetIndex(), voxelPoint);

        if (zrotate) {
            voxelPoint[0] = -voxelPoint[0];
            voxelPoint[1] = -voxelPoint[1];
        }

        int pid = locator->FindClosestPoint(voxelPoint.GetDataPointer());
        double value = scalars->GetTuple1(pid);
        iter.Set(value);
    }

    io.WriteImage(outputImageFile, maskImage);
}


/// @brief perform scan conversion
/// [input-vtk] [reference-image] [output-image]
///
int runScanConversion(pi::Options& opts, pi::StringVector& args) {
    vtkIO vio;
    string inputModelFile = args[0];
    string inputImageFile = args[1];
    string outputImageFile = args[2];

    vtkPolyData* pd = vio.readFile(inputModelFile);

    bool zrotate = opts.GetBool("-zrotate");

    // point flipping
    for (int i = 0; i < pd->GetNumberOfPoints(); i++) {
        double p[3];
        pd->GetPoint(i, p);
        if (zrotate) {
            p[0] = -p[0];
            p[1] = -p[1];
            pd->GetPoints()->SetPoint(i, p);
        }
    }

    vtkSmartPointer<vtkImageData> whiteImage = vtkSmartPointer<vtkImageData>::New();

    ImageIO<MaskImageType> imageIO;
    MaskImageType::Pointer refImage = imageIO.ReadImage(inputImageFile);


    // compute bounding box
    MaskImageType::RegionType region = refImage->GetBufferedRegion();
    MaskImageType::IndexType lowerIndex = region.GetIndex();
    MaskImageType::IndexType upperIndex = region.GetUpperIndex();

    MaskImageType::PointType lowerPoint, upperPoint;
    refImage->TransformIndexToPhysicalPoint(lowerIndex, lowerPoint);
    refImage->TransformIndexToPhysicalPoint(upperIndex, upperPoint);

    // mesh bounds
    double bounds[6];

    // image bounds
    bounds[0] = lowerPoint[0];
    bounds[1] = upperPoint[0];
    bounds[2] = lowerPoint[1];
    bounds[3] = upperPoint[1];
    bounds[4] = lowerPoint[2];
    bounds[5] = upperPoint[2];


    // make the same spacing as refImage
    double spacing[3]; // desired volume spacing
    double origin[3];
    for (int i = 0; i < 3; i++) {
        spacing[i] = refImage->GetSpacing()[i];
        origin[i] = refImage->GetOrigin()[i];
    }
    whiteImage->SetSpacing(spacing);
    whiteImage->SetOrigin(origin);

    // compute dimensions
    int dim[3];
    for (int i = 0; i < 3; i++)
    {
        dim[i] = static_cast<int>(ceil((bounds[i * 2 + 1] - bounds[i * 2]) / spacing[i])) + 1;
    }
    whiteImage->SetDimensions(dim);
    whiteImage->SetExtent(0, dim[0] - 1, 0, dim[1] - 1, 0, dim[2] - 1);

//    double origin[3];
//    origin[0] = bounds[0] + spacing[0] / 2;
//    origin[1] = bounds[2] + spacing[1] / 2;
//    origin[2] = bounds[4] + spacing[2] / 2;
//    whiteImage->SetOrigin(origin);

#if VTK_MAJOR_VERSION <= 5
    whiteImage->SetScalarTypeToUnsignedChar();
    whiteImage->AllocateScalars();
#else
    whiteImage->AllocateScalars(VTK_UNSIGNED_CHAR,1);
#endif
    // fill the image with foreground voxels:
    unsigned char inval = 255;
    unsigned char outval = 0;
    vtkIdType count = whiteImage->GetNumberOfPoints();
    for (vtkIdType i = 0; i < count; ++i)
    {
        whiteImage->GetPointData()->GetScalars()->SetTuple1(i, inval);
    }

    // polygonal data --> image stencil:
    vtkSmartPointer<vtkPolyDataToImageStencil> pol2stenc = vtkSmartPointer<vtkPolyDataToImageStencil>::New();
#if VTK_MAJOR_VERSION <= 5
    pol2stenc->SetInput(pd);
#else
    pol2stenc->SetInputData(pd);
#endif
    pol2stenc->SetOutputOrigin(origin);
    pol2stenc->SetOutputSpacing(spacing);
    pol2stenc->SetOutputWholeExtent(whiteImage->GetExtent());
    pol2stenc->Update();

    // cut the corresponding white image and set the background:
    vtkSmartPointer<vtkImageStencil> imgstenc = vtkSmartPointer<vtkImageStencil>::New();
#if VTK_MAJOR_VERSION <= 5
    imgstenc->SetInput(whiteImage);
    imgstenc->SetStencil(pol2stenc->GetOutput());
#else
    imgstenc->SetInputData(whiteImage);
    imgstenc->SetStencilConnection(pol2stenc->GetOutputPort());
#endif
    imgstenc->ReverseStencilOff();
    imgstenc->SetBackgroundValue(outval);
    imgstenc->Update();

    vtkSmartPointer<vtkMetaImageWriter> writer = vtkSmartPointer<vtkMetaImageWriter>::New();
    writer->SetFileName(outputImageFile.c_str());
#if VTK_MAJOR_VERSION <= 5
    writer->SetInput(imgstenc->GetOutput());
#else
    writer->SetInputData(imgstenc->GetOutput());
#endif
    writer->Write();

    return EXIT_SUCCESS;
}


int main(int argc, char * argv[])
{
    Options opts;
    // general options
    opts.addOption("-o", "Specify a filename for output; used with other options", "-o filename.nrrd", SO_REQ_SEP);
    opts.addOption("-scalarName", "scalar name [string]", SO_REQ_SEP);
    opts.addOption("-outputScalarName", "scalar name for output [string]", SO_REQ_SEP);
    opts.addOption("-sigma", "sigma value [double]", SO_REQ_SEP);
    opts.addOption("-threshold", "Threshold value [double]", SO_REQ_SEP);
    opts.addOption("-iter", "number of iterations [int]", SO_REQ_SEP);
    opts.addOption("-attrDim", "The number of components of attribute", "-attrDim 3 (vector)", SO_REQ_SEP);

    // scalar array handling
    opts.addOption("-exportScalars", "Export scalar values to a text file", "-exportScalars [in-mesh] [scalar.txt]", SO_NONE);
    opts.addOption("-importScalars", "Add scalar values to a mesh [in-mesh] [scalar.txt] [out-mesh]", SO_NONE);
    opts.addOption("-smoothScalars", "Gaussian smoothing of scalar values of a mesh. [in-mesh] [out-mesh]", SO_NONE);
    opts.addOption("-copyScalars", "Copy a scalar array of the input model to the output model", "-copyScalars input-model1 input-model2 output-model -scalarName name", SO_NONE);
    opts.addOption("-averageScalars", "Compute the average of scalars across given inputs", "-averageScalars -o output-vtk input1-vtk input2-vtk ... ", SO_NONE);

    // sampling from an image
    opts.addOption("-sampleImage", "Sample pixels for each point of a given model. Currently, only supported image type is a scalar", "-sampleImage image.nrrd model.vtp output.vtp -outputScalarName scalarName", SO_NONE);
    opts.addOption("-voronoiImage", "Compute the voronoi image from a given data set. A reference image should be given.", "-voronoiImage ref-image.nrrd input-dataset output-image.nrrd -scalarName voxelLabel", SO_NONE);
    opts.addOption("-scanConversion", "Compute a binary image from a surface model", "-scanConversion input-surface input-image.nrrd output-image.nrrd", SO_NONE);

    // mesh processing
    opts.addOption("-appendData", "Append input meshes into a single data [output-mesh]", SO_REQ_SEP);
    opts.addOption("-computeCurvature", "Compute curvature values for each point", "-computeCurvature input-vtk output-vtk", SO_NONE);

    opts.addOption("-vti", "Convert an ITK image to VTI format (VTKImageData)", "-vti imageFile outputFile [-attrDim 3] [-maskImage mask]", SO_NONE);
    opts.addOption("-vtu", "Convert an ITK image to VTU format (vtkUnstructuredGrid). This is useful when masking is needed.", "-vtu imageFile outputFile -maskImage maskImage", SO_NONE);
    opts.addOption("-maskImage", "A mask image for the use of -vtu", "-maskImage mask.nrrd", SO_REQ_SEP);
    opts.addOption("-traceStream", "Trace a stream line from a given point set", "-traceStream input-vtu-field input-vtk output-lines output-points", SO_NONE);
    opts.addOption("-traceDirection", "Choose the direction of stream tracing (both, forward, backward)", "-traceStream ... -traceDirection (both|forward|backward)", SO_REQ_SEP);
    opts.addOption("-zrotate", "Rotate all the points along the z-axis. Change the sign of x and y coordinate.", "-traceStream ... -zrotate", SO_NONE);
    opts.addOption("-traceClipping", "Clip stream lines to fit with an object", "-traceClipping stream_lines.vtp stream_object.vtp stream_lines_output.vtp", SO_NONE);
    opts.addOption("-traceScalarCombine", "Combine scalar values from a seed object to a stream line object. The stream line object must have PointIds for association. -zrotate option will produce the rotated output.", "-traceScalarCombine stream_seed.vtp stream_lines.vtp stream_lines_output.vtp -scalarName scalarToBeCopied", SO_NONE);
    opts.addOption("-filterStream", "Filter out stream lines which are lower than a given threshold", "-filterStream stream-line-input stream-seed-input stream-line-output -scalarName scalar -threshold xx", SO_NONE);
    opts.addOption("-thresholdMin", "Give a minimum threshold value for -filterStream", "-threshold 10 (select a cell whose attriubte is greater than 10)", SO_REQ_SEP);
    opts.addOption("-thresholdMax", "Give a maximum threshold value for -filterStream", "-threshold 10 (select a cell whose attriubte is lower than 10)", SO_REQ_SEP);
    opts.addOption("-fitting", "Fit a model into a binary image", "-fitting input-model binary-image output-model", SO_NONE);
    opts.addOption("-ellipse", "Create an ellipse with parameters []", "-ellipse 101 101 101 51 51 51 20 20 20 -o ellipse.nrrd", SO_NONE);

    opts.addOption("-h", "print help message", SO_NONE);
    StringVector args = opts.ParseOptions(argc, argv, NULL);


    if (opts.GetBool("-h")) {
        cout << "## *kmesh* Usage" << endl;
        opts.PrintUsage();
        return 0;
    } else if (opts.GetBool("-smoothScalars")) {
        runScalarSmoothing(opts, args);
    } else if (opts.GetBool("-importScalars")) {
        runImportScalars(opts, args);
    } else if (opts.GetBool("-exportScalars")) {
        runExportScalars(opts, args);
    } else if (opts.GetBool("-copyScalars")) {
        runCopyScalars(opts, args);
    } else if (opts.GetBool("-averageScalars")) {
        runAverageScalars(opts, args);
    } else if (opts.GetString("-appendData", "") != "") {
        runAppendData(opts, args);
    } else if (opts.GetBool("-sampleImage")) {
        runSampleImage(opts, args);
    } else if (opts.GetBool("-voronoiImage")) {
        runVoronoiImage(opts, args);
    } else if (opts.GetBool("-scanConversion")) {
        runScanConversion(opts, args);
    } else if (opts.GetBool("-vti")) {
        runConvertITK2VTI(opts, args);
    } else if (opts.GetBool("-vtu")) {
        runConvertITK2VTU(opts, args);
    } else if (opts.GetBool("-traceStream")) {
        runStreamTracer(opts, args);
    } else if (opts.GetBool("-filterStream")) {
        runFilterStream(opts, args);
    } else if (opts.GetBool("-fitting")) {
        runFittingModel(opts, args);
    } else if (opts.GetBool("-ellipse")) {
        runEllipse(opts, args);
    } else if (opts.GetBool("-traceClipping")) {
        runTraceClipping(opts, args);
    } else if (opts.GetBool("-computeCurvature")) {
        runComputeCurvature(opts, args);
    } else if (opts.GetBool("-traceScalarCombine")) {
        runTraceScalarCombine(opts, args);
    }
    return 0;
}
