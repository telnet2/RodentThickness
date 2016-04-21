#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkProcrustesAlignmentFilter.h>
#include <vtkLandmarkTransform.h>
#include <vtkThinPlateSplineTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkTimerLog.h>

#include <iostream>
#include <fstream>

#include "piOptions.h"
#include "vtkSmartPointer.h"
using namespace std;

vtkPolyData* readMesh(string file) {
    vtkPolyDataReader* reader = vtkPolyDataReader::New();
    cout << "reading [" << file << "] ..." << flush;
    reader->SetFileName(file.c_str());
    reader->Update();
    cout << " done" << endl;
    vtkPolyData* vtk = reader->GetOutput();
    return vtk;
}

void writeMesh(vtkPolyData* mesh, string file) {
  vtkPolyDataWriter* writer = vtkPolyDataWriter::New();
  writer->SetInputData(mesh);
  writer->SetFileName(file.c_str());
  writer->Update();
}

void computeAverageMesh(string outputMesh, pi::StringVector& args) {
  int argc = args.size();
  vtkProcrustesAlignmentFilter* filter = vtkProcrustesAlignmentFilter::New();
  //filter->SetNumberOfInputs(argc);
  vtkSmartPointer<vtkPolyData> input0;
  for (int i = 0; i < argc; i++) {
    vtkPolyDataReader* reader = vtkPolyDataReader::New();
    cout << "Reading " << args[i] << endl;
    reader->SetFileName(args[i].c_str());
    reader->Update();
    vtkPolyData* vtk = reader->GetOutput();
    cout << "# points: " << vtk->GetNumberOfPoints() << endl;
    //filter->SetInput(i, vtk);
    filter->AddInputDataObject(vtk);

    if(i == 0){
      input0 = vtk;
    }
  }
  cout << "Computing mean ..." << flush;
  filter->GetLandmarkTransform()->SetModeToRigidBody();
  filter->Update();
  cout << " done" << endl;
  vtkPoints* averagePoints = filter->GetMeanPoints();
  vtkPolyData* firstInput = input0;
  firstInput->SetPoints(averagePoints);
  vtkPolyDataWriter* writer = vtkPolyDataWriter::New();
  writer->SetInputData(firstInput);
  writer->SetFileName(outputMesh.c_str());
  writer->Update();
  filter->Delete();
}

void computeTPSWarp(pi::StringVector& args, bool inputFlip) {
  string sourceFile = args[0];
  string targetFile = args[1];
  string inputFile = args[2];
  string outputFile = args[3];

  vtkPolyData* sourceMesh = readMesh(sourceFile);
  vtkPolyData* targetMesh = readMesh(targetFile);
  vtkPolyData* inputMesh  = readMesh(inputFile);

//    for (int i = 0; i < sourceMesh->GetNumberOfPoints(); i++) {
//      double p[3];
//      sourceMesh->GetPoint(i, p);
//      p[0] = -p[0]; p[1] = -p[1];
//      sourceMesh->GetPoints()->SetPoint(i, p);
//    sourceMesh->GetPoint(i, p);
//    cout << i << ", " << p[0] << ", " << p[1] << endl;

  if (inputFlip) {
    for (int i = 0; i < inputMesh->GetNumberOfPoints(); i++) {
      double p[3];
      inputMesh->GetPoint(i, p);
      p[0] = -p[0]; p[1] = -p[1];
      inputMesh->GetPoints()->SetPoint(i, p);
    }
  }

  vtkTimerLog* timer = vtkTimerLog::New();
  timer->StartTimer();

  vtkThinPlateSplineTransform* tpsWarp = vtkThinPlateSplineTransform::New();
  tpsWarp->SetSourceLandmarks(sourceMesh->GetPoints());
  tpsWarp->SetTargetLandmarks(targetMesh->GetPoints());
  tpsWarp->SetBasisToR();
  cout << "computing tps ... " << flush;
  tpsWarp->Update();
  timer->StopTimer();
  double elapsedSecs = timer->GetElapsedTime();
  cout << "done (" << elapsedSecs << ") secs" << endl;

  vtkTransformPolyDataFilter* transformer = vtkTransformPolyDataFilter::New();
  transformer->SetTransform(tpsWarp);
  transformer->SetInputData(inputMesh);
  transformer->Update();

  vtkPolyData* warpedMesh = transformer->GetOutput();
  writeMesh(warpedMesh, outputFile);
}

void computePointFlip(pi::StringVector& args) {
  vtkPolyData* mesh = readMesh(args[0]);
  for (int i = 0; i < mesh->GetNumberOfPoints(); i++) {
    double p[3];
    mesh->GetPoint(i, p);
    p[0] = -p[0]; p[1] = -p[1];
    mesh->GetPoints()->SetPoint(i, p);
  }
  writeMesh(mesh, args[1]);
}


void computePointMatch(pi::StringVector& args) {
  vtkPolyData* source = readMesh(args[0]);
  vtkPolyData* target = readMesh(args[1]);

  if (source->GetNumberOfPoints() != target->GetNumberOfPoints()) {
    cout << "# of points are different" << endl;
    return;
  }

  vtkDoubleArray* array = vtkDoubleArray::New();
  array->SetName("Point Displacement");
  array->SetNumberOfComponents(3);
  array->SetNumberOfTuples(source->GetNumberOfPoints());

  for (int i = 0; i < source->GetNumberOfPoints(); i++) {
    double p[3], q[3];
    source->GetPoint(i, p);
    target->GetPoint(i, q);
    double d[3];
    for (int j = 0; j < 3; j++) {
      d[j] = q[j] - p[j];
    }
    array->SetTuple3(i, d[0], d[1], d[2]);
  }

  source->GetPointData()->AddArray(array);
  writeMesh(source, args[2]);
}

void computePointList(pi::StringVector& args, bool inputflip) {
  vtkPolyData* source = readMesh(args[0]);
  for (int i = 0; i < source->GetNumberOfPoints(); i++) {
    double p[3];
    source->GetPoint(i, p);
    if (inputflip) {
      p[0] = -p[0];
      p[1] = -p[1];
    }
    cout << p[0] << " " << p[1] << " " << p[2] << endl;
  }
}

void computePointAttributeById(pi::StringVector& args) {
  vtkPolyData* source = readMesh(args[0]);

  vtkIntArray* array = vtkIntArray::New();
  array->SetName("PointId");
  array->SetNumberOfComponents(1);
  array->SetNumberOfValues(source->GetNumberOfPoints());

  for (int i = 0; i < source->GetNumberOfPoints(); i++) {
    array->SetValue(i, i);
  }

  source->GetPointData()->AddArray(array);
  writeMesh(source, args[1]);
}

void addPointAttr(pi::StringVector& args) {
  if (args.size() == 0) {
    cout << "--addPointAttr input attribute attribute_name output" << endl;
    return;
  }
  vtkPolyData* source = readMesh(args[0]);

  vtkDoubleArray* array = vtkDoubleArray::New();
  array->SetName(args[2].c_str());
  array->SetNumberOfComponents(1);

  ifstream in(args[1].c_str());
  string line;
  while (!in.eof()) {
    getline(in, line);
    if (in.good()) {
      array->InsertNextValue(atof(line.c_str()));     
    }
  }
  cout << "# of values read: " << array->GetNumberOfTuples() << endl;

  source->GetPointData()->AddArray(array);
  writeMesh(source, args[3]);
}

int main(int argc, char** argv) {
  CSimpleOpt::SOption options[] ={
    { 1, "-o", SO_REQ_SEP },
    { 2, "--average", SO_NONE },
    { 3, "--tps", SO_NONE },
    { 4, "--inputflip", SO_NONE },
    { 5, "--pointflip", SO_NONE },
    { 6, "--pointmatch", SO_NONE },
    { 7, "--pointattrbyid", SO_NONE },
    { 8, "--pointlist", SO_NONE },
    { 9, "--addPointAttr", SO_NONE },
    SO_END_OF_OPTIONS
  };

  pi::Options parser;
  pi::StringVector& args = parser.ParseOptions(argc, argv, options);

  // compute average mesh
  bool _computeAverage = parser.GetBool("--average");
  bool _computeTPS = parser.GetBool("--tps");

  if (_computeAverage) {
    computeAverageMesh(parser.GetString("-o"), args);
  } else if (_computeTPS) {
    computeTPSWarp(args, parser.GetBool("--inputflip"));
  } else if (parser.GetBool("--pointflip")) {
    computePointFlip(args);
  } else if (parser.GetBool("--pointmatch")) {
    computePointMatch(args);
  } else if (parser.GetBool("--pointattrbyid")) {
    computePointAttributeById(args);
  } else if (parser.GetBool("--pointlist")) {
    computePointList(args, parser.GetBool("--inputflip"));
  } else if (parser.GetBool("--addPointAttr")) {
    addPointAttr(args);
  } else {
    cout << "available options: --average, --tps, --pointflip, --pointmatch, --pointattrbyid, --pointlist, --addPointAttr" << endl;
  }
}
