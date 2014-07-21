
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
#include <algorithm>

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
#include <vtkTriangle.h>
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
#include <vtkPCAAnalysisFilter.h>
#include <vtkProcrustesAlignmentFilter.h>
#include <vtkLandmarkTransform.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkThreshold.h>
#include <vtkConnectivityFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkTransform.h>


#include <itkImage.h>
#include <itkVectorNearestNeighborInterpolateImageFunction.h>
#include <itkEllipseSpatialObject.h>
#include <itkSpatialObjectToImageFilter.h>
#include <itkNearestNeighborInterpolateImageFunction.h>

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_sparse_matrix.h>
#include <vnl/vnl_sparse_matrix_linear_system.h>
#include <vnl/algo/vnl_lsqr.h>
#include <vnl/algo/vnl_matrix_update.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>

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
static void cart2sph(const double *v, double *phi, double *theta) {
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


static void legendre(int n, float x, float* Y) {
    if (n < 0) return;

    float **P = new float*[n + 1];
    for (int i = 0; i <= n; i++) P[i] = new float[n + 1];
    float factor = -sqrt(1.0 - pow(x,2));

    // Init legendre
    P[0][0] = 1.0;        // P_0,0(x) = 1
    if (n == 0)
    {
        Y[0] = P[0][0];
        return;
    }

    // Easy values
    P[1][0] = x;      // P_1,0(x) = x
    P[1][1] = factor;     // P_1,1(x) = −sqrt(1 − x^2)
    if (n == 1)
    {
        Y[0] = P[1][0];
        Y[1] = P[1][1];
        return;
    }

    for (int l = 2; l <= n; l++)
    {
        for (int m = 0; m < l - 1 ; m++)
        {
            // P_l,m = (2l-1)*x*P_l-1,m - (l+m-1)*x*P_l-2,m / (l-k)
            P[l][m] = ((float)(2 * l - 1) * x * P[l - 1][m] - (float)(l + m - 1) * P[l - 2][m]) / (float)(l - m);
        }
        // P_l,l-1 = (2l-1)*x*P_l-1,l-1
        P[l][l-1] = (float)(2 * l - 1) * x * P[l-1][l-1];
        // P_l,l = (2l-1)*factor*P_l-1,l-1
        P[l][l] = (float)(2 * l - 1) * factor * P[l-1][l-1];
    }

    for (int i = 0; i <= n; i++) Y[i] = P[n][i];

    // release memory
    for (int i = 0; i <= n; i++) delete [] P[i];
    delete [] P;
}


/// @brief Compute the basis function for spherical harmonics
static void spharm_basis(int degree, double *p, double *Y) {
    // real spherical harmonics basis functions
    // polar coordinate
    double phi, theta;
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
        legendre(l, cos(theta), Pm);
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


std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}


void runExportPoints(Options& opts, StringVector& args) {
    string inputFile = args[0];
    vtkIO vio;
    vtkPolyData* input = vio.readFile(inputFile);
    const unsigned int nPoints = input->GetNumberOfPoints();
    
    string outputText = opts.GetString("-o");
    if (outputText == "") {
        cout << "output file must be specified" << endl;
        return;
    }
    
    ofstream fout(outputText);
    for (unsigned int j = 0; j < nPoints; j++) {
        double* p = input->GetPoint(j);
        for (int k = 0; k < 3; k++) {
            float nearest = floorf(p[k] * 1e8 + 0.5) / 1e8;
            fout << nearest << " ";
        }
        fout << endl;
    }
    fout.close();
}


void runImportPoints(Options& opts, StringVector& args) {
    string inputFile = args[0];
    string txtFile = args[1];
    string outputFile = opts.GetString("-o");

    if (outputFile == "") {
        cout << "output file is not set. check -o option." << endl;
        return;
    }
    
    cout << "reading mesh file: " << inputFile << endl;
    vtkIO vio;
    vtkPolyData* input = vio.readFile(inputFile);
    const unsigned int nPoints = input->GetNumberOfPoints();

    cout << "reading text file: " << txtFile << endl;
    ifstream infile(txtFile.c_str());
    std::string line;
    std::vector<double> txtPoints;
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        float a, b, c;
        if (!(iss >> a >> b >> c)) {
            cout << "Error while reading points!" << endl;
            break;
        } // error
        
        // process pair (a,b)
        txtPoints.push_back(a);
        txtPoints.push_back(b);
        txtPoints.push_back(c);
    }
    
    if (txtPoints.size() != nPoints * 3) {
        cout << "Insufficient number of points: ("
            << nPoints << " required, but " << txtPoints.size() << " found)" << endl;
        infile.close();
        return;
    }
    
    for (unsigned int j = 0; j < nPoints; j++) {
        double p[3];
        for (unsigned int k = 0; k < 3; k++) {
            p[k] = txtPoints[j*3+k];
        }
        input->GetPoints()->SetPoint(j, p);
    }
    

    vio.writeFile(outputFile, input);
}



void runTranslatePoints(Options& opts, StringVector& args) {
    string inputFile = args[0];
    string translate = opts.GetString("-translatePoints");
    cout << "translate by " << translate << endl;
    std::vector<string> strdata = split(translate, ',');
    std::vector<double> pntdata;
    for (unsigned int j = 0; j < strdata.size(); j++) {
        pntdata.push_back(atof(strdata[j].c_str()));
    }

    vtkIO vio;
    vtkPolyData* input = vio.readFile(inputFile);
    const unsigned int nPoints = input->GetNumberOfPoints();

    for (unsigned int j = 0; j < nPoints; j++) {
        double* p = input->GetPoint(j);
        for (unsigned int k = 0; k < pntdata.size(); k++) {
            p[k] += pntdata[k];
        }
        input->GetPoints()->SetPoint(j, p);
    }
    
    vio.writeFile(args[1], input);
}



void runSPHARMCoeff(Options& opts, StringVector& args) {
    string inputFile = args[0];

    vtkIO vio;
    vtkPolyData* input = vio.readFile(inputFile);
    const unsigned int nPoints = input->GetNumberOfPoints();

    vtkDataArray* scalars = input->GetPointData()->GetScalars(opts.GetString("-scalarName").c_str());
    if (scalars == NULL) {
        cout << "Can't read scalars: " << opts.GetString("-scalarName") << endl;
        return;
    }

    int degree = 20;


    /// Prepare values
    vnl_vector<double> values;
    values.set_size(scalars->GetNumberOfTuples());
    for (int i = 0; i < nPoints; i++) {
        values[i] = scalars->GetTuple1(i);
    }

    /// Prepare the vandermonde matrix from spharm basis
    vnl_matrix<double> bases;
    bases.set_size(nPoints, (degree-1)*(degree-1));

    for (uint i = 0; i < nPoints; i++) {
        double p[3];
        input->GetPoints()->GetPoint(i,p);
        spharm_basis(degree, p, bases[i]);
    }

    vnl_matrix_inverse<double> inv(bases);
    vnl_vector<double> coeffs = inv.pinverse() * values;

    vnl_vector<double> interpolated = bases * coeffs;
    vtkDoubleArray* newScalars = vtkDoubleArray::New();
    newScalars->SetName(opts.GetString("-outputScalarName", "NewScalars").c_str());
    newScalars->SetNumberOfComponents(1);
    newScalars->SetNumberOfTuples(nPoints);

    for (int i = 0; i < nPoints; i++) {
        newScalars->SetValue(i, interpolated[i]);
    }
    input->GetPointData()->AddArray(newScalars);

    vio.writeFile(args[1], input);

    return;
}

// append polydatas into one
void runAppendData(Options& opts, StringVector& args) {
    vtkAppendPolyData* appender = vtkAppendPolyData::New();
    vtkIO io;

    // read all files
    for (unsigned int i = 0; i < args.size(); i++) {
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


/// @brief Import vector values
void runImportVectors(Options& opts, StringVector& args) {
    string inputDataFile = args[0];
    string vectorDataFile = args[1];
    string outputDataFile = args[2];

    string attributeName = opts.GetString("-scalarName", "VectorValues");
    vtkIO vio;
    vtkPolyData* inputData = vio.readFile(inputDataFile);

    if (inputData == NULL) {
        cout << "Can't read " << inputDataFile << endl;
        return;
    }
    vtkFloatArray* vectorArray = vtkFloatArray::New();
    vectorArray->SetName(attributeName.c_str());

    char buff[1024];
    ifstream fin(vectorDataFile.c_str());

    int nCols = 0;
    int nRows = 0;
    while (fin.is_open() && !fin.eof()) {
        fin.getline(buff, sizeof(buff));
        stringstream ss(buff);
        if (!fin.eof()) {
            nCols = 0;
            std::vector<float> rows;
            while (ss.good()) {
                float value;
                ss >> value;
                rows.push_back(value);
                ++ nCols;
            }
            if (nRows == 0) {
                vectorArray->SetNumberOfComponents(nCols);
            }
            vectorArray->InsertNextTupleValue(&rows[0]);
            nRows ++;
        }
    }
    fin.close();

    cout << nRows << " x " << nCols << endl;
    if (nRows != inputData->GetNumberOfPoints()) {
        cout << "# of data is different from # of points." << endl;
        return;
    }


    /// Optionally, compute vector stats
    if (opts.GetBool("-computeVectorStats")) {
        vtkDoubleArray* meanValues = vtkDoubleArray::New();
        meanValues->SetName((attributeName + "_Avg").c_str());
        meanValues->SetNumberOfValues(nRows);

        vtkDoubleArray* stdValues = vtkDoubleArray::New();
        stdValues->SetName((attributeName + "_Std").c_str());
        stdValues->SetNumberOfValues(nRows);

        for (int i = 0; i < nRows; i++) {
            vector<double> values(nCols);
            vectorArray->GetTuple(i, &values[0]);
            MeanStd<double> meanStd;
            meanStd = for_each(values.begin(), values.end(), meanStd);

            meanValues->SetValue(i, meanStd.mean());
            stdValues->SetValue(i, meanStd.std());
        }

        inputData->GetPointData()->AddArray(meanValues);
        inputData->GetPointData()->AddArray(stdValues);
    }

    inputData->GetPointData()->AddArray(vectorArray);
    vio.writeFile(outputDataFile, inputData);
}


/// @brief Indirect import scalar values
void runIndirectImportScalars(Options& opts, StringVector& args) {

    /// input arguments
    string inputDataFile = args[0];
    string vectorDataFile = args[1];
    string outputDataFile = args[2];

    /// default attribute name to establish mapping
    string attributeName = opts.GetString("-scalarName", "MappingAttribute");

    /// read input vtk file
    vtkIO vio;
    vtkPolyData* inputData = vio.readFile(inputDataFile);
    if (inputData == NULL) {
        cout << "Can't read " << inputDataFile << endl;
        return;
    }

    const int nPoints = inputData->GetNumberOfPoints();
    vtkDataArray* mappingAttrs = inputData->GetPointData()->GetArray(attributeName.c_str());

    if (mappingAttrs == NULL) {
        cout << "Can't find the attribute: " << mappingAttrs << endl;
        return;
    }


    /// read data file into an array
    vtkFloatArray* vectorArray = vtkFloatArray::New();
    vectorArray->SetName(attributeName.c_str());

    char buff[1024];
    ifstream fin(vectorDataFile.c_str());

    int nCols = 0;
    int nRows = 0;
    while (fin.is_open() && !fin.eof()) {
        fin.getline(buff, sizeof(buff));
        stringstream ss(buff);
        if (!fin.eof()) {
            nCols = 0;
            std::vector<float> rows;
            while (ss.good()) {
                float value;
                ss >> value;
                rows.push_back(value);
                ++ nCols;
            }
            if (nRows == 0) {
                vectorArray->SetNumberOfComponents(nCols);
            }
            vectorArray->InsertNextTupleValue(&rows[0]);
            nRows ++;
        }
    }
    fin.close();

    cout << nRows << " x " << nCols << endl;

    /// assume the line no is an index to mapping
    vector<double> data;

    /// must have fixed elements per line
    data.resize(nCols);



    /// Copy the file data into an attribute array
    string importedAttributeName = opts.GetString("-outputScalarName", "ImportedAttribute");

    vtkFloatArray* importedArray = vtkFloatArray::New();
    importedArray->SetName(importedAttributeName.c_str());
    importedArray->SetNumberOfComponents(nCols);
    importedArray->SetNumberOfValues(nPoints);

    for (int i = 0; i < nPoints; i++) {
        int mappedId = (int) mappingAttrs->GetTuple1(i);
        vectorArray->GetTuple(mappedId, &data[0]);

        /// copy into new array
        importedArray->SetTuple(i, &data[0]);
    }

    inputData->GetPointData()->AddArray(importedArray);

    vio.writeFile(outputDataFile, inputData);

}

/// @brief Export vector values
void runExportVectors(Options& opts, StringVector& args) {
    string inputDataFile = args[0];
    string vectorDataFile = args[1];

    vtkIO vio;
    vtkPolyData* inputData = vio.readFile(inputDataFile);

    int nArr = inputData->GetPointData()->GetNumberOfArrays();
    vtkDataArray* vectorArray = NULL;

    for (int j = 0; j < nArr; j++) {
        if (opts.GetString("-scalarName", "VectorValues") == inputData->GetPointData()->GetArrayName(j)) {
            vectorArray = inputData->GetPointData()->GetArray(j);
            break;
        }
    }

    if (vectorArray == NULL) {
        cout << "Can't find the vector values : " << opts.GetString("-scalarName", "VectorValues") << endl;
        return;
    }

    const int nCols = vectorArray->GetNumberOfComponents();
    const int nRows = vectorArray->GetNumberOfTuples();

    ofstream out(vectorDataFile.c_str());
    std::vector<double> row;
    row.resize(nCols);
    for (int i = 0; i < nRows; i++) {
        vectorArray->GetTuple(i, &row[0]);
        for (int j = 0; j < nCols; j++) {
            if (j > 0) {
                out << "\t";
            }
            out << row[j];
        }
        out << endl;
    }
    out.close();

    return;
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

/// @brief Rescale the streamline with a given length
void runRescaleStream(Options& opts, StringVector& args) {
    string inputStreamFile = args[0];
    string inputDataFile = args[1];
    string scalarName = opts.GetString("-scalarName");

    vtkIO vio;
    vtkPolyData* inputStream = vio.readFile(inputStreamFile);
    vtkPolyData* inputData = vio.readFile(inputDataFile);
    vtkDataArray* inputScalars = inputData->GetPointData()->GetScalars(scalarName.c_str());

    for (int i = 0; i < inputStream->GetNumberOfLines(); i++) {
        double length = inputScalars->GetTuple1(i);
        vtkPolyLine* line = vtkPolyLine::SafeDownCast(inputStream->GetCell(i));
        if (line == NULL) {
            vtkLine* aline = vtkLine::SafeDownCast(inputStream->GetCell(i));
        } else {
            // reduce length
        }
    }
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

    cout << "Building an image data for sampling ..." << flush;
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

    if (pointSet->GetNumberOfPoints() < 1) {
        cout << "# of points should be greater than 0" << endl;
        return;
    }
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

    for (uint i = 0; i < args.size(); i++) {
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
    for (unsigned int j = 0; j < scalars->GetNumberOfTuples(); j++) {
        double v = scalars->GetValue(j) / args.size();
        if (useThreshold) {
            scalars->SetValue(j, v < thresholdValue ? 1 : 2);
        } else {
            scalars->SetValue(j, v);
        }
    }

    for (unsigned int i = 0; i < args.size(); i++) {
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



int runIndexToPoint(Options& opts, StringVector& args) {
    string inputIdx = opts.GetString("-indexToPoint");
    
    cout << "Input Index: " << inputIdx << endl;
    StringVector values = split(inputIdx, ',');
    ImageType::IndexType idx;
    for (int j = 0; j < values.size(); j++) {
        idx[j] = atoi(values[j].c_str());
    }


    ImageType::PointType point;
    ImageIO<ImageType> imgIO;
    ImageType::Pointer img = imgIO.ReadCastedImage(args[0]);
    img->TransformIndexToPhysicalPoint(idx, point);
    
    cout << point << endl;
    return EXIT_SUCCESS;
}


int runPointToIndex(Options& opts, StringVector& args) {
    string inputPoint = opts.GetString("-pointToIndex");
    cout << "Input Point: " << inputPoint << endl;
    StringVector values = split(inputPoint, ',');
    ImageType::PointType point;

    for (int j = 0; j < values.size(); j++) {
        point[j] = atof(values[j].c_str());
    }
    
    
    itk::ContinuousIndex<float,3> idx;
    ImageIO<ImageType> imgIO;
    ImageType::Pointer img = imgIO.ReadCastedImage(args[0]);
    img->TransformPhysicalPointToContinuousIndex(point, idx);
    cout << idx << endl;
    return EXIT_SUCCESS;

}


int runImageInfo(Options& opts, StringVector& args) {
    ImageInfo lastImageInfo;
    ImageIO<ImageType> imgIo;
    imgIo.ReadCastedImage(args[0], lastImageInfo);
    printf("Filename: %s\n", args[0].c_str());
    printf("Dims: %d %d %d\n", lastImageInfo.dimensions[0], lastImageInfo.dimensions[1], lastImageInfo.dimensions[2]);
    printf("Pixdims: %.2f %.2f %.2f\n", lastImageInfo.spacing[0], lastImageInfo.spacing[1], lastImageInfo.spacing[2]);
    printf("Origins: %.2f %.2f %.2f\n", lastImageInfo.origin[0], lastImageInfo.origin[1], lastImageInfo.origin[2]);
    return EXIT_SUCCESS;
}


/// @brief Run PCA analysis
void runPCA(Options& opts, StringVector& args) {
    vtkIO vio;
    std::vector<vtkPolyData*> inputs;
    inputs.resize(args.size());
    /// - Read a series of vtk files
    for (unsigned int i = 0; i < args.size(); i++) {
        inputs[i] = vio.readFile(args[i]);
    }

    vtkPCAAnalysisFilter* pcaFilter = vtkPCAAnalysisFilter::New();
    pcaFilter->SetNumberOfInputs(args.size());
    for (unsigned int i = 0; i < args.size(); i++) {
        pcaFilter->SetInput(i, inputs[i]);
    }
    pcaFilter->Update();
    vtkFloatArray* eigenValues = pcaFilter->GetEvals();
    for (unsigned int i = 0; i < eigenValues->GetNumberOfTuples(); i++) {
        cout << eigenValues->GetValue(i) << endl;
    }

    if (opts.GetString("-pcaMeanOut") != "") {
        vtkFloatArray* array = vtkFloatArray::New();
        pcaFilter->GetParameterisedShape(array, inputs[0]);
        vio.writeFile(opts.GetString("-pcaMeanOut"), inputs[0]);
    }
}

/// @brief Perform Procrustes alignment
void runProcrustes(Options& opts, StringVector& args) {
    int nInputs = args.size() / 2;

    vtkIO vio;
    vtkProcrustesAlignmentFilter* pros = vtkProcrustesAlignmentFilter::New();
    pros->SetNumberOfInputs(nInputs);
    for (int i = 0; i < nInputs; i++) {
        pros->SetInput(i, vio.readFile(args[i]));
    }
    pros->GetLandmarkTransform()->SetModeToSimilarity();
    pros->Update();

    for (unsigned int i = nInputs; i < args.size(); i++) {
        vio.writeFile(args[i], pros->GetOutput(i-nInputs));
    }

}


static bool sortSecond(const std::pair<int,double> &left, const std::pair<int,double> &right) {
    return left.second > right.second;
}



/// @brief Perform extraction of regions based on scalars. First, threshold with scalars.
/// Second, compute the connected components. Third, compute the area of each region.
/// Fourth, threshold each region with the computed area, exclude below 20%.
void runConnectScalars(Options& opts, StringVector& args) {
    string inputFile = args[0];
    string scalarName = opts.GetString("-scalarName");

    vtkIO vio;
    vtkPolyData* inputData = vio.readFile(inputFile.c_str());

    const int nPoints = inputData->GetNumberOfPoints();
    const double tmin = opts.GetStringAsReal("-thresholdMin", itk::NumericTraits<double>::min());
    const double tmax = opts.GetStringAsReal("-thresholdMax", itk::NumericTraits<double>::max());

    inputData->GetPointData()->SetScalars(inputData->GetPointData()->GetScalars(scalarName.c_str()));

    vtkIntArray* originalIds = vtkIntArray::New();
    originalIds->SetName("OriginalId");
    for (int i = 0; i < nPoints; i++) {
        originalIds->InsertNextValue(i);
    }
    inputData->GetPointData()->AddArray(originalIds);

    vtkThreshold* threshold = vtkThreshold::New();
    threshold->SetInput(inputData);
    threshold->ThresholdBetween(tmin, tmax);
    threshold->Update();

    vtkConnectivityFilter* conn = vtkConnectivityFilter::New();
    conn->SetInputConnection(threshold->GetOutputPort());
    conn->ColorRegionsOn();
    conn->SetExtractionModeToAllRegions();
    conn->Update();

    vtkUnstructuredGrid* connComp = conn->GetOutput();
    int nCells = connComp->GetNumberOfCells();

    vtkDataArray* regionIds = connComp->GetCellData()->GetScalars("RegionId");
    double regionIdRange[2];
    regionIds->GetRange(regionIdRange);
    int maxId = regionIdRange[1];

    /// Compute the area of each region
    typedef std::pair<int, double> IntPair;
    std::vector<IntPair> regionAreas;
    regionAreas.resize(maxId + 1);
    IntPair zero(0,0);
    std::fill(regionAreas.begin(), regionAreas.end(), zero);

    for (int i = 0; i < nCells; i++) {
        vtkCell* cell = connComp->GetCell(i);
        vtkTriangle* tri = vtkTriangle::SafeDownCast(cell);
        if (tri == NULL) {
            cout << "Non-triangle cell" << endl;
            continue;
        }
        double area = tri->ComputeArea();
        regionAreas[regionIds->GetTuple1(i)].first = regionIds->GetTuple1(i);
        regionAreas[regionIds->GetTuple1(i)].second += area;
    }

    /// Propagate areas to each cell
    vtkDoubleArray* cellArea = vtkDoubleArray::New();
    cellArea->SetName("RegionArea");
    cellArea->SetNumberOfTuples(nCells);

    for (int i = 0; i < nCells; i++) {
        int regionId = regionIds->GetTuple1(i);
        cellArea->SetTuple1(i, regionAreas[regionId].second);
    }
    connComp->GetCellData()->AddArray(cellArea);


    std::sort(regionAreas.begin(), regionAreas.end(), sortSecond);
    double totalArea = 0;
    for (int i = 0; i < maxId; i++) {
        totalArea += regionAreas[i].second;
        cout << totalArea << endl;
    }

    /// Re-label regionIds in the decreasing order of its area
    std::vector<int> regionRanking;
    regionRanking.resize(regionAreas.size());

    for (unsigned int i = 0; i < regionRanking.size(); i++) {
        regionRanking[regionAreas[i].first] = (i+1);
    }

    /// Update region Ids for every cell
    for (unsigned int i = 0; i < nCells; i++) {
        int regionId = regionIds->GetTuple1(i);
        regionIds->SetTuple1(i, regionRanking[regionId]);
    }

    /// Update the RegionId in point data
    vtkDataArray* regionIds2 = connComp->GetPointData()->GetScalars("RegionId");
    for (unsigned int i = 0; i < regionIds2->GetNumberOfTuples(); i++) {
        int regionId = regionIds2->GetTuple1(i);
        regionIds2->SetTuple1(i, regionRanking[regionId]);
    }


    /// Compute the 90% percentile
    double totalPer = 0, areaMin = 0;
    for (int i = 0; i < maxId; i++) {
        totalPer += (regionAreas[i].second / totalArea);
        if (totalPer > 0.90) {
            areaMin = regionAreas[i].second;
            break;
        }
    }

    cout << "Minimum Region Area for 90%: " << areaMin << endl;

    /// Threshold each cell with area
    vtkThreshold* areaThresholder = vtkThreshold::New();
    areaThresholder->SetInput(connComp);
    areaThresholder->ThresholdByUpper(areaMin);
    areaThresholder->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "RegionArea");
    areaThresholder->Update();
    vtkUnstructuredGrid* areaOutput = areaThresholder->GetOutput();



    /// New RegonId array
    vtkIntArray* pointRegionIds = vtkIntArray::New();
    pointRegionIds->SetName("RegionIds");
    pointRegionIds->SetNumberOfTuples(nPoints);
    for (int i = 0; i < pointRegionIds->GetNumberOfTuples(); i++) {
        pointRegionIds->SetTuple1(i, 0);
    }

    vtkDataArray* areaRegionIds = areaOutput->GetPointData()->GetScalars("RegionId");
    vtkDataArray* pointIds = areaOutput->GetPointData()->GetScalars("OriginalId");
    for (int i = 0; i < pointIds->GetNumberOfTuples(); i++) {
        int ptid = pointIds->GetTuple1(i);
        int regionId = areaRegionIds->GetTuple1(i);

        pointRegionIds->SetTuple1(ptid, regionId);
    }

    inputData->GetPointData()->AddArray(pointRegionIds);
    vio.writeFile(args[1], inputData);


//    vtkDataArray* scalars = inputData->GetPointData()->GetScalars(scalarName.c_str());
//    for (int i = 0; i < nPoints; i++) {
//        double value = scalars->GetTuple1(i);
//        if (value < tmin || value > tmax) {
//            value = std::numeric_limits<double>::quiet_NaN();
//        }
//    }


//    vtkPolyDataConnectivityFilter* conf = vtkPolyDataConnectivityFilter::New();
//    conf->SetInput(inputData);
//    conf->ScalarConnectivityOn();
//    conf->SetExtractionModeToAllRegions();
//    conf->FullScalarConnectivityOn();
//    conf->ColorRegionsOn();
//    conf->SetScalarRange(0, 0.05);
//    conf->Update();

//    vtkPolyData* output = conf->GetOutput();
//    vio.writeFile(args[1], output);
}



/// @brief Find neighbors of a point
void extractNeighbors(std::vector<int>& ids, vtkPolyData* data, std::set<int>& neighbors, int nRing) {
    /// for every id in ids
    for (unsigned int i = 0; i < ids.size(); i++) {
        int id = ids[i];

        /// if id is found in neighbors
        if (neighbors.find(id) == neighbors.end()) {
            /// add to neighbors
            neighbors.insert(id);
        }

        if (nRing > 0) {
            /// find every neighboring cells
            vtkIdList* cellIds = vtkIdList::New();
            data->GetPointCells(id, cellIds);


            std::vector<int> idSet;
            for (int j = 0; j < cellIds->GetNumberOfIds(); j++) {
                const int cellId = cellIds->GetId(j);
                vtkIdList* pointIds = vtkIdList::New();
                data->GetCellPoints(cellId, pointIds);
                for (int k = 0; k < pointIds->GetNumberOfIds(); k++) {
                    int newId = pointIds->GetId(k);
                    if (id == 110) {
                        //                        cout << newId << endl;
                    }
                    /// if id is found in neighbors, proceed to next point
                    if (neighbors.find(newId) == neighbors.end()) {
                        idSet.push_back(newId);
                    }
                }
                pointIds->Delete();
            }
            cellIds->Delete();

            extractNeighbors(idSet, data, neighbors, nRing - 1);
        }
    }
}

template <class T1, class T2>
struct SecondComp {
    bool operator()(const pair<T1,T2>& p1, const pair<T1,T2>& p2) {
        return p1.second < p2.second;
    }
};

/// @brief Perform correlation-based clustering
void runCorrelationClustering(Options& opts, StringVector& args) {
    string inputDataFile = args[0];
    string inputVectorName = opts.GetString("-scalarName", "VectorValues");
    string outputScalarName = opts.GetString("-outputScalarName", "CorrClusters");
    string outputDataFile = args[1];

    vtkIO vio;
    vtkPolyData* inputData = vio.readFile(inputDataFile);
    const int nPoints = inputData->GetNumberOfPoints();

    int nArr = inputData->GetPointData()->GetNumberOfArrays();
    vtkDataArray* vectorArray = NULL;

    for (int j = 0; j < nArr; j++) {
        if (opts.GetString("-scalarName", "VectorValues") == inputData->GetPointData()->GetArrayName(j)) {
            vectorArray = inputData->GetPointData()->GetArray(j);
            break;
        }
    }

    /// # of tuples
    const int nCols = vectorArray->GetNumberOfComponents();

    /// storage for vector data
    std::vector<double> idata;
    idata.resize(vectorArray->GetNumberOfComponents());

    std::vector<double> jdata;
    jdata.resize(vectorArray->GetNumberOfComponents());


    /// Compute the mean and stdev for each vector
    std::vector<double> dataMean, dataStd;
    for (int i = 0; i < nPoints; i++) {
        vectorArray->GetTuple(i, &idata[0]);
        double sum = 0;
        double sum2 = 0;
        for (int j = 0; j < nCols; j++) {
            double val = idata[j];
//            cout << val << "\t";
            sum += val;
            sum2 += val * val;
        }
//        cout << endl;
        double avg = sum / nCols;
        double var = sum2/nCols - avg*avg;
        double std = sqrt(var);

        dataMean.push_back(avg);
        dataStd.push_back(std);

//        cout << std << endl;
    }



    /// Iterate over all points and compute the maximum correlation between neighbors
    vtkIntArray* maxCorrNeighbor = vtkIntArray::New();
    maxCorrNeighbor->SetName("CorrNeighbors");
    maxCorrNeighbor->SetNumberOfValues(nPoints);

    for (int i = 0; i < nPoints; i++) {
        // set up for neighbor finding
        std::set<int> neighbors;
        std::vector<int> ids;
        ids.push_back(i);

        vectorArray->GetTuple(i, &idata[0]);
        double imean = dataMean[i];
        double istd = dataStd[i];

        extractNeighbors(ids, inputData, neighbors, 1);
        std::set<int>::iterator iter = neighbors.begin();

        int maxNeighbor = -1;
        double maxCorr = -10;

        for (; iter != neighbors.end(); iter++) {
            int j = *iter;
            if (j == i) {
                continue;
            }
            vectorArray->GetTuple(j, &jdata[0]);

            /// Compute the cross correlation
            double jmean = dataMean[j];
            double jstd = dataStd[j];

            double corr = 0;
            for (int k = 0; k < nCols; k++) {
                corr += ((idata[k] - imean) * (jdata[k] - jmean));
            }
            corr /= (istd * jstd);
//            cout << corr << ",";
            if (corr > maxCorr) {
                maxCorr = corr;
                maxNeighbor = j;
            }
        }
//        cout << endl;
        maxCorrNeighbor->SetValue(i, maxNeighbor);
    }

    /// connected components with maxCorrNeighbors
    /// Loop over all points and mark its region
    /// If it encounters previously marked region, then stop and mark its region to the previous region
    std::vector<int> markups;
    markups.resize(nPoints);
    std::fill(markups.begin(), markups.end(), -1);

    std::vector<int> curRegion;

    cout << nPoints << endl;
    int regionCounter = 0;
    for (int i = 0; i < nPoints; i++) {
        if (markups[i] >= 0) {
            continue;
        }

        int j = i;
        curRegion.clear();
        curRegion.push_back(j);
        markups[j] = regionCounter;

        while (true) {
            j = maxCorrNeighbor->GetValue(j);
//            cout << j << " ";
            int jMark = markups[j];

            bool keepConnecting = (jMark < 0);
            if (keepConnecting) {
                markups[j] = regionCounter;
                curRegion.push_back(j);
            } else if (jMark == regionCounter) {
                regionCounter++;
                break;
            } else if (jMark < regionCounter) {
                // replace all regionCounter into jMark
                for (unsigned int k = 0; k < curRegion.size(); k++) {
                    markups[curRegion[k]] = jMark;
                }
                break;
            }
        }
//        cout << ":" << regionCounter << endl;
    }

    vtkIntArray* regionArray = vtkIntArray::New();
    regionArray->SetName(outputScalarName.c_str());
    regionArray->SetNumberOfValues(markups.size());
    vtkIntArray* regionModArray = vtkIntArray::New();
    regionModArray->SetName((outputScalarName + "Mod").c_str());
    regionModArray->SetNumberOfValues(markups.size());
    for (unsigned int i = 0; i < markups.size(); i++) {
        regionArray->SetValue(i, markups[i]);
        regionModArray->SetValue(i, markups[i] / 10);
//        cout << markups[i] << " ";
    }
//    cout << endl;

    inputData->GetPointData()->AddArray(regionArray);
    inputData->GetPointData()->AddArray(maxCorrNeighbor);
    inputData->GetPointData()->AddArray(regionModArray);
    vio.writeFile(outputDataFile, inputData);


//    if (true) return;


    /// output separation
    /// will create two files: control & experimental
    /// Assume that the offset is given
    typedef pair<int, int> IdValuePair;
    typedef std::vector<IdValuePair> IdValueVector;

    IdValueVector idValues;
    for (int i = 0; i < nPoints; i++) {
        idValues.push_back(make_pair(i, regionArray->GetValue(i)));
    }

    /// sort idValues by the second element
    sort(idValues.begin(), idValues.end(), SecondComp<int, int>());



    /// print data files
    string data1 = args[2];
    string data2 = args[3];

    /// file handle for two data files
    ofstream of1(data1.c_str()), of2(data2.c_str());

    /// need change to use the option
    const int nControls = 8;
    int pointPerRegion = 0;
    int prevRegionId = -1;

    for (int i = 0; i < nPoints; i++) {
        int pointId = idValues[i].first;
        int regionId = idValues[i].second;
        pointPerRegion++;

        if (i > 0 && regionId != prevRegionId) {
            // handle previous region
            of1 << endl;
            of2 << endl;

            cout << prevRegionId << "\t" << (pointPerRegion - 1) << endl;
            
            /// Prepare for the next region
            pointPerRegion = 1;

            /// process next region
        }


        double* data = vectorArray->GetTuple(pointId);
        for (int j = 0; j < nCols; j++) {
            if (j < nControls) {
                of1 << data[j] << " ";
            } else {
                of2 << data[j] << " ";
            }
        }

        prevRegionId = regionId;
    }
    cout << prevRegionId << "\t" << (pointPerRegion) << endl;

    of1 << endl;
    of2 << endl;

    of1.close();
    of2.close();
}

/// @brief perform ridge detection
void runDetectRidge(Options& opts, StringVector& args) {
    vtkIO vio;
    string inputFile = args[0];
    const int nRing = opts.GetStringAsInt("-n", 2);

    vtkPolyData* inputData = vio.readFile(inputFile);
    vtkPoints* inputPoints = inputData->GetPoints();
    const int nPoints = inputData->GetNumberOfPoints();

    vtkPolyDataNormals* normalFilter = vtkPolyDataNormals::New();
    normalFilter->SetInput(inputData);
    normalFilter->ComputePointNormalsOn();
    normalFilter->ConsistencyOn();
    normalFilter->GetAutoOrientNormals();
    normalFilter->Update();
    vtkFloatArray* normals = vtkFloatArray::SafeDownCast(normalFilter->GetOutput()->GetPointData()->GetScalars("Normals"));
    if (normals == NULL) {
        cout << "Can't compute point normals" << endl;
        return;
    }

    cout << "# points: " << nPoints << endl;
    cout << "# rings (neighbor degrees): " << nRing << endl;

    vtkDoubleArray* gaussianCurv = vtkDoubleArray::New();
    gaussianCurv->SetNumberOfTuples(nPoints);
    gaussianCurv->SetName("GaussianCurvature");

    vtkDoubleArray* e1 = vtkDoubleArray::New();
    e1->SetNumberOfTuples(nPoints);
    e1->SetName("CurvatureExtreme1");

    vtkDoubleArray* e2 = vtkDoubleArray::New();
    e2->SetNumberOfTuples(nPoints);
    e2->SetName("CurvatureExtreme2");

    vtkDoubleArray* shapeOperator = vtkDoubleArray::New();
    shapeOperator->SetNumberOfComponents(3);
    shapeOperator->SetNumberOfTuples(nPoints);
    shapeOperator->SetName("ShapeOperator");

    vtkDoubleArray* majorDirection = vtkDoubleArray::New();
    majorDirection->SetNumberOfComponents(3);
    majorDirection->SetNumberOfTuples(nPoints);
    majorDirection->SetName("MajorDirection");

    vtkDoubleArray* minorDirection = vtkDoubleArray::New();
    minorDirection->SetNumberOfComponents(3);
    minorDirection->SetNumberOfTuples(nPoints);
    minorDirection->SetName("MinorDirection");

    vtkDoubleArray* kappa1 = vtkDoubleArray::New();
    kappa1->SetNumberOfComponents(1);
    kappa1->SetNumberOfTuples(nPoints);
    kappa1->SetName("Kappa1");

    vtkDoubleArray* kappa2 = vtkDoubleArray::New();
    kappa2->SetNumberOfComponents(1);
    kappa2->SetNumberOfTuples(nPoints);
    kappa2->SetName("Kappa2");


    /// Loop over all points

    for (int i = 0; i < nPoints; i++) {
        /// Find the center
        double x[3];
        inputPoints->GetPoint(i, x);

        /// Extract neighbor points
        std::set<int> neighbors;
        std::vector<int> ids;
        ids.push_back(i);
        extractNeighbors(ids, inputData, neighbors, nRing);


        /// Compute the point normal
        double normal[3];
        normals->GetTuple(i, normal);

        /// Compute the rotation matrix
        const double northPole[3] = { 0, 0, 1 };

        vnl_matrix<double> rotation(3,3);
        vio.rotateVector(normal, northPole, rotation);

        int nNeighbors = neighbors.size();

        vnl_matrix<double> U(3*nNeighbors, 7);
        vnl_vector<double> d(3*nNeighbors);

        std::set<int>::iterator iter = neighbors.begin();
        for (int j = 0; iter != neighbors.end(); j++, iter++) {
            int nbrId = *iter;
            double p[3];

            inputPoints->GetPoint(nbrId, p);

            vnl_vector<double> px(3), pn(3);
            vtkMath::Subtract(p, x, px.data_block());
            normals->GetTuple(nbrId, pn.data_block());

            /// Rotate the neighbor to the north pole
            vnl_vector<double> q = rotation * px;
            vnl_vector<double> qn = rotation * pn;

//            double k2 = sqrt(2/q[0]*q[0] + q[1]*q[1]);
//
//            if (isnan(k2)) {
//                k2 = 0;
//            }
//
//            q[0] *= k2;
//            q[1] *= k2;
//            q[2] *= k2;

            U[3*j][0] = 0.5 * q[0]*q[0];
            U[3*j][1] = q[0]*q[1];
            U[3*j][2] = 0.5 * q[1]*q[1];
            U[3*j][3] = q[0]*q[0]*q[0];
            U[3*j][4] = q[0]*q[0]*q[1];
            U[3*j][5] = q[0]*q[1]*q[1];
            U[3*j][6] = q[1]*q[1]*q[1];
            d[3*j] = q[2];

            U[3*j+1][0] = q[0];
            U[3*j+1][1] = q[1];
            U[3*j+1][2] = 0;
            U[3*j+1][3] = 3*q[0]*q[0];
            U[3*j+1][4] = 2*q[0]*q[1];
            U[3*j+1][5] = q[1]*q[1];
            U[3*j+1][6] = 0;

            if (qn[2] == 0) {
                d[3*j+1] = 0;
            } else {
                d[3*j+1] = -qn[0] / qn[2];
            }

            U[3*j+2][0] = 0;
            U[3*j+2][1] = q[0];
            U[3*j+2][2] = q[1];
            U[3*j+2][3] = 0;
            U[3*j+2][4] = q[0]*q[0];
            U[3*j+2][5] = 2*q[0]*q[1];
            U[3*j+2][6] = 3*q[1]*q[1];

            if (qn[2] == 0) {
                d[3*j+2] = 0;
            } else {
                d[3*j+2] = -qn[1] / qn[2];
            }
        }

        vnl_matrix_inverse<double> Uinv(U);
        vnl_vector<double> coeff = Uinv.pinverse() * d;
        vnl_matrix<double> W(2,2);

        W[0][0] = coeff[0];
        W[0][1] = W[1][0] = coeff[1];
        W[1][1] = coeff[2];

        vnl_symmetric_eigensystem<double> Weigen(W);


        double k1 = Weigen.get_eigenvalue(1);
        double k2 = Weigen.get_eigenvalue(0);
        double gk = k1 * k2;

        double e[2] = { 0, 0 };
        for (int k = 0; k < 2; k++) {
            double t1 = Weigen.get_eigenvector(k)[1];
            double t2 = Weigen.get_eigenvector(k)[0];
            e[k] = t1*(t1*t1*coeff[3] + t2*t2*coeff[5]/3.0) + t2*(t1*t1*coeff[4]/3.0 + t2*t2*coeff[6]);
        }
        e[0] = e[0]*e[0] + e[1]*e[1];
//        cout << e[0] << "," << e[1] << endl;

        kappa1->SetValue(i, k1);
        kappa2->SetValue(i, k2);

        vnl_matrix<double> rotation_inv = rotation.transpose();
        vnl_vector<double> majorDir(3);
        majorDir[0] = Weigen.get_eigenvector(1)[0];
        majorDir[1] = Weigen.get_eigenvector(1)[1];
        majorDir[2] = 0;
        majorDir = rotation_inv * majorDir;

        vnl_vector<double> minorDir(3);
        minorDir[0] = Weigen.get_eigenvector(0)[0];
        minorDir[1] = Weigen.get_eigenvector(0)[1];
        minorDir[2] = 0;
        minorDir = rotation_inv * minorDir;


//        cout << k1 << "," << k2 << endl;
//        cout << Weigen.get_eigenvector(1) << " => " << majorDir << ";" << endl;
//        cout << Weigen.get_eigenvector(0) << " => " << minorDir << ";" << endl;

        majorDirection->SetTuple(i, majorDir.data_block());
        minorDirection->SetTuple(i, minorDir.data_block());

        gaussianCurv->SetValue(i, gk);
        shapeOperator->SetTuple3(i, W[0][0], W[0][1], W[1][1]);
        e1->SetValue(i, e[0]);
        e2->SetValue(i, e[1]);
    }

    inputData->GetPointData()->AddArray(gaussianCurv);
    inputData->GetPointData()->AddArray(shapeOperator);
    inputData->GetPointData()->AddArray(e1);
    inputData->GetPointData()->AddArray(e2);
    inputData->GetPointData()->AddArray(kappa1);
    inputData->GetPointData()->AddArray(kappa2);
    inputData->GetPointData()->AddArray(majorDirection);
    inputData->GetPointData()->AddArray(minorDirection);

    vio.writeFile(args[1], inputData);
}



/// @brief run TPS warping
void runTPSWarp(Options& opts, StringVector& args) {
    
}


/// @brief run a test code
void runTest(Options& opts, StringVector& args) {
    vtkIO io;

    double v1[3] = {0, 0, 5 };
    double v2[3] = { 0, 1, 0 };
    vnl_matrix<double> rotation;

    io.rotateVector(v1, v2, rotation);

    vnl_vector<double> vv(3);
    vv.set_size(3);
    vv.copy_in(v1);
    cout << vv << endl;
    cout << rotation << endl;
    cout << rotation * vv << endl;
}

int main(int argc, char * argv[])
{
    Options opts;
    // general options
    opts.addOption("-o", "Specify a filename for output; used with other options", "-o filename.nrrd", SO_REQ_SEP);
    opts.addOption("-n", "Specify n (integer number) for an operation. Refer related options", "-n integer",SO_REQ_SEP);
    opts.addOption("-scalarName", "scalar name [string]", SO_REQ_SEP);
    opts.addOption("-outputScalarName", "scalar name for output [string]", SO_REQ_SEP);
    opts.addOption("-sigma", "sigma value [double]", SO_REQ_SEP);
    opts.addOption("-threshold", "Threshold value [double]", SO_REQ_SEP);
    opts.addOption("-iter", "number of iterations [int]", SO_REQ_SEP);
    opts.addOption("-attrDim", "The number of components of attribute", "-attrDim 3 (vector)", SO_REQ_SEP);
    opts.addOption("-thresholdMin", "Give a minimum threshold value for -filterStream, -connectScalars", "-thresholdMin 10 (select a cell whose attriubte is greater than 10)", SO_REQ_SEP);
    opts.addOption("-thresholdMax", "Give a maximum threshold value for -filterStream -connectScalars", "-thresholdMax 10 (select a cell whose attriubte is lower than 10)", SO_REQ_SEP);
    opts.addOption("-test", "Run test code", SO_NONE);

    // points handling
    opts.addOption("-exportPoints", "Save points into a text file", "-exportPoints [in-mesh] -o [out-txt]", SO_NONE);
    opts.addOption("-importPoints", "Read points in a text file into a vtk file", "-importPoints [in-mesh] [txt] -o [out-vtk]", SO_NONE);
    opts.addOption("-translatePoints", "Translate points by adding given tuples", "-translatePoints=x,y,z [in-vtk] [out-vtk]", SO_REQ_SEP);
    
    // scalar array handling
    opts.addOption("-exportScalars", "Export scalar values to a text file", "-exportScalars [in-mesh] [scalar.txt]", SO_NONE);
    opts.addOption("-importScalars", "Add scalar values to a mesh [in-mesh] [scalar.txt] [out-mesh]", SO_NONE);
    opts.addOption("-indirectImportScalars", "Import scalar values indirectly via another mapping attribute", "-indirectImportScalars in-vtk in-txt out-vtk -scalarName mapping-attribute -outputScalarName imported-attribute", SO_NONE);
    opts.addOption("-smoothScalars", "Gaussian smoothing of scalar values of a mesh. [in-mesh] [out-mesh]", SO_NONE);
    opts.addOption("-importVectors", "Add vector values to a mesh [in-mesh] [scalar.txt] [out-mesh] [-computeVectorStats]", SO_NONE);
    opts.addOption("-exportVectors", "Export vector values to a mesh [in-mesh] [scalar.txt]", SO_NONE);
    opts.addOption("-computeVectorStats", "Compute mean and std for a vector attribute", "-importVectors ... [-computeVectorStats]", SO_NONE);

    opts.addOption("-copyScalars", "Copy a scalar array of the input model to the output model", "-copyScalars input-model1 input-model2 output-model -scalarName name", SO_NONE);
    opts.addOption("-averageScalars", "Compute the average of scalars across given inputs", "-averageScalars -o output-vtk input1-vtk input2-vtk ... ", SO_NONE);
    opts.addOption("-connectScalars", "Compute the connected components based on scalars and assign region ids", "-connectScalars input.vtk output.vtk -scalarName scalar -thresholdMin min -thresholdMax max", SO_NONE);
    opts.addOption("-corrClustering", "Compute correlational clusters -corrClustering input-vtk -scalarName values-to-compute-correlation -outputScalarName clusterId", SO_NONE);

    // sampling from an image
    opts.addOption("-sampleImage", "Sample pixels for each point of a given model. Currently, only supported image type is a scalar", "-sampleImage image.nrrd model.vtp output.vtp -outputScalarName scalarName", SO_NONE);
    opts.addOption("-voronoiImage", "Compute the voronoi image from a given data set. A reference image should be given.", "-voronoiImage ref-image.nrrd input-dataset output-image.nrrd -scalarName voxelLabel", SO_NONE);
    opts.addOption("-scanConversion", "Compute a binary image from a surface model", "-scanConversion input-surface input-image.nrrd output-image.nrrd", SO_NONE);
    opts.addOption("-indexToPoint", "Transform an index to a physical point with a given image", "-indexToPoint=i,j,k [image-file]", SO_REQ_SEP);
    opts.addOption("-pointToIndex", "Transform a physical point to an index", "-pointToIndex=x,y,z [image-file]", SO_REQ_SEP);
    opts.addOption("-imageInfo", "Print out basic image inforation like ImageStat", "-imageInfo [image1] [image2] ...", SO_NONE);


    // mesh processing
    opts.addOption("-appendData", "Append input meshes into a single data [output-mesh]", SO_REQ_SEP);
    opts.addOption("-computeCurvature", "Compute curvature values for each point", "-computeCurvature input-vtk output-vtk", SO_NONE);
    opts.addOption("-pca", "Perform PCA analysis", "-pca input1-vtk input2-vtk ... -o output.txt", SO_NONE);
    opts.addOption("-pcaMeanOut", "A filename for PCA mean output", "-pca ... -pcaMeanOut file.vtk", SO_REQ_SEP);
    opts.addOption("-procrustes", "Perform Procrustes alignment", "-procrustes input1.vtk input2.vtk ... output1.vtk output2.vtk ...", SO_NONE);

    opts.addOption("-vti", "Convert an ITK image to VTI format (VTKImageData)", "-vti imageFile outputFile [-attrDim 3] [-maskImage mask]", SO_NONE);
    opts.addOption("-vtu", "Convert an ITK image to VTU format (vtkUnstructuredGrid). This is useful when masking is needed.", "-vtu imageFile outputFile -maskImage maskImage", SO_NONE);
    opts.addOption("-maskImage", "A mask image for the use of -vtu", "-maskImage mask.nrrd", SO_REQ_SEP);
    opts.addOption("-traceStream", "Trace a stream line from a given point set", "-traceStream input-vtu-field input-vtk output-lines output-points", SO_NONE);
    opts.addOption("-traceDirection", "Choose the direction of stream tracing (both, forward, backward)", "-traceStream ... -traceDirection (both|forward|backward)", SO_REQ_SEP);
    opts.addOption("-zrotate", "Rotate all the points along the z-axis. Change the sign of x and y coordinate.", "-traceStream ... -zrotate", SO_NONE);
    opts.addOption("-traceClipping", "Clip stream lines to fit with an object", "-traceClipping stream_lines.vtp stream_object.vtp stream_lines_output.vtp", SO_NONE);
    opts.addOption("-traceScalarCombine", "Combine scalar values from a seed object to a stream line object. The stream line object must have PointIds for association. -zrotate option will produce the rotated output.", "-traceScalarCombine stream_seed.vtp stream_lines.vtp stream_lines_output.vtp -scalarName scalarToBeCopied", SO_NONE);
    opts.addOption("-rescaleStream", "Rescale streamlines to fit with given lengths", "-rescaleStream input-stream-lines length.txt or input.vtp -scalarName scalarname", SO_NONE);
    
    
    /// TPS warping
    opts.addOption("-tpsWarp", "Run TPS warping", "-tpsWarp [source-landmark] [target-landmark] [input0] [output0] [input1] [output1] ...", SO_NONE);


    /// Compute SPHARM smoothing
    opts.addOption("-spharmCoeff", "Compute SPHARM coefficients", "-spharmCoeff input-vtk output-txt -scalarName scalarValueToEvaluate", SO_NONE);

    /// Mesh operation
    opts.addOption("-detectRidge", "Run ridge detection", "-n nRings -detectRidge input.vtk output.vtk", SO_NONE);

    opts.addOption("-filterStream", "Filter out stream lines which are lower than a given threshold", "-filterStream stream-line-input stream-seed-input stream-line-output -scalarName scalar -threshold xx", SO_NONE);
    opts.addOption("-fitting", "Fit a model into a binary image", "-fitting input-model binary-image output-model", SO_NONE);
    opts.addOption("-ellipse", "Create an ellipse with parameters []", "-ellipse 101 101 101 51 51 51 20 20 20 -o ellipse.nrrd", SO_NONE);

    opts.addOption("-h", "print help message", SO_NONE);
    StringVector args = opts.ParseOptions(argc, argv, NULL);


    if (opts.GetBool("-h")) {
        cout << "## *kmesh* Usage" << endl;
        opts.PrintUsage();
        return 0;
    } else if (opts.GetBool("-importPoints")) {
        runImportPoints(opts, args);
    } else if (opts.GetBool("-exportPoints")) {
        runExportPoints(opts, args);
    } else if (opts.GetString("-translatePoints") != "") {
        runTranslatePoints(opts, args);
    } else if (opts.GetBool("-smoothScalars")) {
        runScalarSmoothing(opts, args);
    } else if (opts.GetBool("-importScalars")) {
        runImportScalars(opts, args);
    } else if (opts.GetBool("-indirectImportScalars")) {
        runIndirectImportScalars(opts, args);
    } else if (opts.GetBool("-exportScalars")) {
        runExportScalars(opts, args);
    } else if (opts.GetBool("-importVectors")) {
        runImportVectors(opts, args);
    } else if (opts.GetBool("-exportVectors")) {
        runExportVectors(opts, args);
    } else if (opts.GetBool("-copyScalars")) {
        runCopyScalars(opts, args);
    } else if (opts.GetBool("-averageScalars")) {
        runAverageScalars(opts, args);
    } else if (opts.GetBool("-connectScalars")) {
        runConnectScalars(opts, args);
    } else if (opts.GetBool("-corrClustering")) {
        runCorrelationClustering(opts, args);
    } else if (opts.GetString("-appendData", "") != "") {
        runAppendData(opts, args);
    } else if (opts.GetBool("-sampleImage")) {
        runSampleImage(opts, args);
    } else if (opts.GetBool("-voronoiImage")) {
        runVoronoiImage(opts, args);
    } else if (opts.GetBool("-scanConversion")) {
        runScanConversion(opts, args);
    } else if (opts.GetString("-indexToPoint") != "") {
        return runIndexToPoint(opts, args);
    } else if (opts.GetString("-pointToIndex") != "") {
        return runPointToIndex(opts, args);
    } else if (opts.GetBool("-imageInfo")) {
        return runImageInfo(opts, args);
    } else if (opts.GetBool("-vti")) {
        runConvertITK2VTI(opts, args);
    } else if (opts.GetBool("-vtu")) {
        runConvertITK2VTU(opts, args);
    } else if (opts.GetBool("-traceStream")) {
        runStreamTracer(opts, args);
    } else if (opts.GetBool("-filterStream")) {
        runFilterStream(opts, args);
    } else if (opts.GetBool("-rescaleStream")) {
        runRescaleStream(opts, args);
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
    } else if (opts.GetBool("-pca")) {
        runPCA(opts, args);
    } else if (opts.GetBool("-procrustes")) {
        runProcrustes(opts, args);
    } else if (opts.GetBool("-spharmCoeff")) {
        runSPHARMCoeff(opts, args);
    } else if (opts.GetBool("-detectRidge")) {
        runDetectRidge(opts, args);
    } else if (opts.GetBool("-tpsWarp")) {
        runTPSWarp(opts, args);
    } else if (opts.GetBool("-test")) {
        runTest(opts, args);
    }
    return 0;
}
