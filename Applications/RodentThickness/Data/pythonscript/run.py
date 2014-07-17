#!/tools/Python/Current/bin/python

'''
This python script accepts .csv data set, .bms configurations, and output directory then produces statistical analysis results from the processed data
'''


from optparse import OptionParser
import csv
import os,os.path,re,sys
import system

'''
concatToColumns.py
'''
def countLines(files):
    lines = []
    for file in files:
        fin = open(file)
        lines.append(len(fin.readlines()))
        fin.close()
    return lines

def readAllFiles(files):
    lines = []
    for file in files:
        fin = open(file)
        lines.append(fin.readlines())
        fin.close()
    return lines

def file2string(f):
    fin = open(f)
    lines = fin.readlines()
    lines = map(lambda x: x.strip(), lines)
    return lines


def find_executables(exes):
    for exe in exes:
        if system.find_program(exe) is None:
            raise RuntimeError("can't find %s in PATH variable" % exe)


def concatColumnsToFile(input, output, inputIsList=False):
    if (inputIsList):
        args = input
    else:
        if (input != ""):
            args = file2string(input)
        else:
            return
    lines = countLines(args)
    if (len(lines) == 0):
        print "Empty files"
        return
    if (min(lines) != max(lines)):
        print "Some of input files have different line numbers: ", lines
        return
    lineStrings = readAllFiles(args)
    if (output != ""):
        fout = open(output, "w")
    else:
        fout = sys.stdout
    for i in range(0,lines[0]):
        tempOutputLine = [ x[i].strip() for x in lineStrings]
        outputLineString = "\t".join(tempOutputLine)
        fout.write(outputLineString)
        fout.write("\n")
    fout.close()

'''
execute MeshPointsIntensitySampling
'''
def performAnalysis(cvsdata, config, outputDir, opt="initial_dense"):
    ids = []
    groups = []
    for (id, labelmap, group) in cvsdata:
        ids.append(id)
        groups.append(group)

    meshDir = "%s/Processing/1.MeasurementandSPHARM" % (outputDir)
    shapeworksDir = "%s/Processing/2.shapeworks" % (outputDir)
    sampleDir = "%s/Processing/3.meshintensity" % (outputDir)
    statDir = "%s/Statistical" % (outputDir)
    pathKmesh = config["kmeshPath"]

    if (opt == "initial_dense"):
        samplingCmd = config["MeshPointIntensitySamplingPath"] + " --workDir %s --inputAsPhysicalCoord --inputPointFlip --distanceVector %s -i nn -a %s -m %s --smoothedAttributeOutput %s --originalMeshOutput %s %s %s"
    else:
        samplingCmd = config["MeshPointIntensitySamplingPath"] + " --workDir %s --inputAsPhysicalCoord --distanceVector %s -i nn -a %s -m %s --smoothedAttributeOutput %s --originalMeshOutput %s %s %s"



    # depending on options, choose different meshes
    # loop over ids
    for id in ids:
        distanceVector = "%s/%s.distanceVector.nrrd" % (sampleDir, id)
        inputMeasurement = "%s/%s_measurementoutput.nrrd" % (meshDir, id)
        #    inputMeasurement = "%s/%s.thickness.mha" % (meshDir, id)
        if (opt == "initial_dense"):
            samplingTxt = "%s/%s.initialDenseSampling.txt" % (sampleDir, id)
            smoothSamplingTxt = "%s/%s.smoothInitialDenseSampling.txt" % (sampleDir, id)
            samplingMeshOutput = "%s/%s.initialDenseSampling.vtk" % (sampleDir, id)
            originalMeshOutput = "%s/%s.initialDenseOriginal.vtk" % (sampleDir, id)
            inputMesh = "%s/%s.subj.SPHARM.vtk" % (meshDir, id)
        elif (opt == "correspondence"):
            samplingTxt = "%s/%s.sampling.txt" % (sampleDir, id)
            smoothSamplingTxt = "%s/%s.smoothSampling.txt" % (sampleDir, id)
            samplingMeshOutput = "%s/%s.sampling.vtk" % (sampleDir, id)
            originalMeshOutput = "%s/%s.original.vtk" % (sampleDir, id)
            inputMesh = "%s/%s.correspondence.vtk" % (shapeworksDir, id)
        elif (opt == "dense_correspondence"):
            samplingTxt = "%s/%s.denseSampling.txt" % (sampleDir, id)
            smoothSamplingTxt = "%s/%s.smoothDenseSampling.txt" % (sampleDir, id)
            samplingMeshOutput = "%s/%s.denseSampling.vtk" % (sampleDir, id)
            originalMeshOutput = "%s/%s.denseOriginal.vtk" % (sampleDir, id)
            inputMesh = "%s/%s.dense_correspondence.vtk" % (shapeworksDir, id)
        elif (opt == "spharm_sampling"):
            samplingTxt = "%s/%s.spharmSampling.txt" % (sampleDir, id)
            smoothSamplingTxt = "%s/%s.smoothSpharmSampling.txt" % (sampleDir, id)
            samplingMeshOutput = "%s/%s.spharmSampling.vtk" % (sampleDir, id)
            originalMeshOutput = "%s/%s.spharmOriginal.vtk" % (sampleDir, id)
            inputMesh = "%s/%s.subj.SPHARM.vtk" % (meshDir, id)

        # execute sampling
        # exeCmd = samplingCmd % (outputDir, distanceVector, samplingTxt, samplingMeshOutput, smoothSamplingTxt, originalMeshOutput, inputMeasurement, inputMesh)
        # new sampling command
        useStreamLine = True and opt == "spharm_sampling"
        if (not useStreamLine):
            exeCmd = "%s -sampleImage %s %s %s -outputScalarName Thickness -zrotate" % (pathKmesh, inputMeasurement, inputMesh, samplingMeshOutput)
            system.run_process(exeCmd,verbose=True)

            exeCmd = "%s -smoothScalars %s %s -scalarName Thickness -outputScalarName smoothThickness -iter 3" % (pathKmesh, samplingMeshOutput, samplingMeshOutput)
            system.run_process(exeCmd,verbose=True)

            exeCmd = "%s -exportScalars %s %s -scalarName Thickness " % (pathKmesh, samplingMeshOutput, samplingTxt)
            system.run_process(exeCmd,verbose=True)

            exeCmd = "%s -exportScalars %s %s -scalarName smoothThickness " % (pathKmesh, samplingMeshOutput, smoothSamplingTxt)
            system.run_process(exeCmd,verbose=True)
        else:
            gradientFile = "%s/Processing/1.MeasurementandSPHARM/%s/laplacianGradient.nrrd" % (outputDir, id)
            gradientVTIFile = "%s/Processing/1.MeasurementandSPHARM/%s/laplacianGradient.vti" % (outputDir, id)
            streamLineFile = "%s/Processing/1.MeasurementandSPHARM/%s/stream_lines.vtp" % (outputDir, id)

            # vti image creation
            exeCmd = "%s -vti %s %s -attrDim 3" % (pathKmesh, gradientFile, gradientVTIFile)
            system.run_process(exeCmd,verbose=True)

            # RK45 integration
            exeCmd = "%s -traceStream %s %s %s %s -zrotate -traceDirection backward" % (pathKmesh, gradientVTIFile, inputMesh, streamLineFile, samplingMeshOutput)
            system.run_process(exeCmd,verbose=True)

            # scalar export
            exeCmd = "%s -exportScalars %s %s -scalarName Length " % (pathKmesh, samplingMeshOutput, smoothSamplingTxt)
            system.run_process(exeCmd,verbose=True)



    # aggregate sampling files into a single file
    groupSet = set(groups)
    datafiles = []
    tag = ""
    for group in groupSet:
        # setup output filename template
        if (opt == "initial_dense"):
            tag = "_initialDenseSampling"
            datafilename = "%s/data_initialDenseSampling_%s.txt" % (statDir, group)
            dataregionfilename = "%s/data_region_initialDenseSampling_%s.txt" % (statDir, group)
            outfilename = "%s/list_initialDenseSampling_%s.txt" % (statDir, group)
        elif (opt == "correspondence"):
            tag = ""
            datafilename = "%s/data_%s.txt" % (statDir, group)
            dataregionfilename = "%s/data_region_%s.txt" % (statDir, group)
            outfilename = "%s/list_%s.txt" % (statDir, group)
        elif (opt == "dense_correspondence"):
            tag = "_dense"
            datafilename = "%s/data_dense_%s.txt" % (statDir, group)
            dataregionfilename = "%s/data_region_dense_%s.txt" % (statDir, group)
            outfilename = "%s/list_dense_%s.txt" % (statDir, group)
        elif (opt == "spharm_sampling"):
            tag = "_spharm"
            datafilename = "%s/data_spharm_%s.txt" % (statDir, group)
            dataregionfilename = "%s/data_region_spharm_%s.txt" % (statDir, group)
            outfilename = "%s/list_spharm_%s.txt" % (statDir, group)


        # process and write the aggregated data
        outfile = open(outfilename, 'w')
        for (id, g) in zip(ids, groups):
            if (g != group):
                continue
            if (opt == "initial_dense"):
                samplingTxt = "%s/%s.smoothInitialDenseSampling.txt" % (sampleDir, id)
            elif (opt == "correspondence"):
                samplingTxt = "%s/%s.smoothSampling.txt" % (sampleDir, id)
            elif (opt == "dense_correspondence"):
                samplingTxt = "%s/%s.smoothDenseSampling.txt" % (sampleDir, id)
            elif (opt == "spharm_sampling"):
                samplingTxt = "%s/%s.smoothSpharmSampling.txt" % (sampleDir, id)
            if (os.path.isfile(samplingTxt)):
                outfile.write(samplingTxt)
                outfile.write("\n")
            else:
                print "can't find", samplingTxt
        outfile.close()

        # concatenate those attribute files into a file
        concatColumnsToFile(outfilename, datafilename)
        datafiles.append((group, datafilename))

    # perform statistical analysis
    rscriptExePath = config["RScriptPath"]
    rscriptPath = config["RScriptRodentThicknessPath"]
    if (opt == "initial_dense"):
        outputfile = "%s/stats.initialDenseSampling.txt" % statDir
        outputVTK = "%s/stats.initialDenseSampling.vtk" % (statDir)
        inputVTK = "%s/%s.subj.SPHARM.vtk" % (meshDir, ids[0])
    elif (opt == "correspondence"):
        outputfile = "%s/stats.correspondence.txt" % statDir
        outputVTK = "%s/stats.correspondence.vtk" % (statDir)
        inputVTK = "%s/%s.correspondence.vtk" % (shapeworksDir, ids[0])
    elif (opt == "dense_correspondence"):
        outputfile = "%s/stats.denseCorrespondence.txt" % statDir
        outputVTK = "%s/stats.denseCorrespondence.vtk" % (statDir)
        inputVTK = "%s/%s.dense_correspondence.vtk" % (shapeworksDir, ids[0])
    elif (opt == "spharm_sampling"):
        outputfile = "%s/stats.spharmCorrespondence.txt" % statDir
        outputVTK = "%s/stats.spharmCorrespondence.vtk" % (statDir)
        inputVTK = "%s/%s.subj.SPHARM.vtk" % (meshDir, ids[0])

    pythonPath = config["PythonPath"]
    pythonScriptPath = config["vtkPointAttributesPythonScriptPath"]

    datafilelist = ""
    datagroups = ""
    for (group, file) in datafiles:
        datafilelist = datafilelist + " " + file
        datagroups = datagroups + " " + group

        # compute per-group statistics
        cmd = "%s -computeVectorStats -scalarName %s -importVectors %s %s %s/%s_thickness.vtk" % (pathKmesh, group + "_Thickness", inputVTK, file, statDir, group)
        print cmd
        system.run_process(cmd,verbose=True)


    if (len(datafiles) > 1):
        runCmd = "%s %s %s %s %s" % (rscriptExePath, rscriptPath, datafilelist, outputfile, datagroups)
        system.run_process(runCmd,verbose=True)
        visCmd = "%s %s %s %s -i %s -t" % (pythonPath, pythonScriptPath, inputVTK, outputVTK, outputfile)
        system.run_process(visCmd,verbose=True)

        # connected components for p-value
        print "computing connected components with t.pvalue..."
        runCmd = "%s -connectScalars %s -scalarName t.pvalue %s -thresholdMin 0 -thresholdMax 0.05" % (pathKmesh, outputVTK, outputVTK)
        system.run_process(runCmd,verbose=True)

        print "export scalars..."
        runCmd = "%s -exportScalars -scalarName RegionIds %s %s/regionIds.txt" % (pathKmesh, outputVTK, statDir)
        system.run_process(runCmd,verbose=True)

        # create a combined data file with regions
        regionIdFile = "%s/regionIds.txt" % (statDir)
        if os.path.exists(regionIdFile):
            for (group, file) in datafiles:
                dataregionfilename = "%s/data_region%s_%s.txt" % (statDir, tag, group)
                concatColumnsToFile([regionIdFile, "%s/data%s_%s.txt" % (statDir, tag, group)], dataregionfilename, inputIsList = True)
                # perform splitting
                performSplit([dataregionfilename], "%s/data_%s_region_" % (statDir, group) + "%02d.txt")

            # boxplot
            groupList = list(groupSet)
            runCmd = "%s %s '%s' '%s' %s/%s_%s_regions.pdf --label1 %s --label2 %s" % (pythonPath, pythonScriptPath.replace("vtkPointAttributes.py", "boxplot.py"), "%s/data_%s_region_*.txt" % (statDir, groupList[0]), "%s/data_%s_region_*.txt" % (statDir, groupList[1]), statDir, groupList[0], groupList[1], groupList[0], groupList[1])
            print runCmd
            system.run_process(runCmd,verbose=True)






'''
Label map processing for thickness measurement
def labelprocessing(config, csvdata, outputdir, labelNCleft, labelNCright, labelOB):
  workdir = "%s/Processing/1.MeasurementandSPHARM"
  for (id, labelmap, group):
    brain_nc_left = "%s/%s_brain_left.nrrd" % (workdir, id)
    brain_nc_left = "%s/%s_brain_left.nrrd" % (workdir, id)
    brain_nc_left = "%s/%s_brain_left.nrrd" % (workdir, id)
    brain_nc_left = "%s/%s_brain_left.nrrd" % (workdir, id)
    brain_nc_left = "%s/%s_brain_left.nrrd" % (workdir, id)
'''

# run thickness measurement
# if the measurement file exists, skip it
def compute_thickness(config, data, outputdir, ids, idl, idh):
    pathTool = config["measureThicknessFilterPath"]
    for (id, labelmap, group) in data:
        workdir = "%s/Processing/1.MeasurementandSPHARM/%s" % (outputdir, id)
        measurementoutput = "%s/Processing/1.MeasurementandSPHARM/%s_measurementoutput.nrrd" % (outputdir, id)

        # temporary
        #    workdir = "%s/Processing/1.MeasurementandSPHARM/%s_2" % (outputdir, id)
        labelmap = "%s/Processing/1.MeasurementandSPHARM/%s.boundaryMap.mha" % (outputdir, id)
        #    measurementoutput = "%s/Processing/1.MeasurementandSPHARM/%s.thickness.mha" % (outputdir, id)

        if (not os.path.exists(workdir)):
            os.makedirs(workdir)
        if (not os.path.isfile(measurementoutput)):
            cmd = "%s --mr --sbt --workdir %s --ids %s --idl %s --idh %s --ttrns 500 %s %s" % (pathTool, workdir, ids, idl, idh, labelmap, measurementoutput)
            print cmd
            system.run_process(cmd,verbose=True)
        else:
            print "Skipping", measurementoutput




def genparamesh(opts, config, data, outputDir):
    pathGenParaMesh = config["GenParaMeshCLPPath"]
    pathImageMath = config["ImageMathPath"]
    pathSegPostProcess = "SegPostProcess" #config["SegPostProcess"]
    workdir = "%s/Processing/1.MeasurementandSPHARM" % (outputDir)

    # check if executables exist
    find_executables([pathGenParaMesh, pathImageMath, pathSegPostProcess])

    def file_exists(f):
        return not opts.overwrite and os.path.isfile(f) and os.stat(f).st_size > 0

    for (id, labelmap, group) in data:
        eulerName = "%s/%s.euler.txt" % (workdir, id)
        paraFile = "%s/%s.para.vtk" % (workdir, id)
        surfFile = "%s/%s.surf.vtk" % (workdir, id)
        surfaceLabel = "%s/%s.seg.nii.gz" % (workdir, id)
        segpostLabel = "%s/%s.segpost.nii.gz" % (workdir, id)

        if (not file_exists(surfaceLabel)):
            cmd = "%s %s -extractLabel %s -outfile %s" \
                  % (pathImageMath, labelmap, opts.ids, surfaceLabel)
            system.run_process(cmd,verbose=opts.verbose)

        if (not file_exists(segpostLabel)):
            cmd = "%s %s -o %s" \
                  % (pathSegPostProcess, surfaceLabel, segpostLabel)
            system.run_process(cmd,verbose=opts.verbose)

        if (not file_exists(paraFile)):
            num_iters = opts.genparamesh_iter
            cmd = "%s --EulerFile --outEulerName %s %s --iter %d --label %s %s %s" \
                  % (pathGenParaMesh, eulerName, segpostLabel, num_iters, 1, paraFile, surfFile)
            system.run_process(cmd,verbose=opts.verbose)


def paratospharm(opts, config, data, outputDir):
    workdir = "%s/Processing/1.MeasurementandSPHARM" % (outputDir)
    pathParaToSPHARM = config["ParaToSPHARMMeshCLPPath"]

    idlist = [id for id,l,g in data]
    parafiles = ["%s/%s.para.vtk"%(workdir,id) for id,l,g in data]
    surffiles = ["%s/%s.surf.vtk"%(workdir,id) for id,l,g in data]
    lowresfiles = ["%s/%s.ip.SPHARM.vtk"%(workdir,id) for id,l,g in data]
    highresfiles = ["%s/%s.subj.SPHARM.vtk"%(workdir,id) for id,l,g in data]
    lowresparams = ["%s/%s.ip.Param.vtk"%(workdir,id) for id,l,g in data]
    highresparams = ["%s/%s.subj.Param.vtk"%(workdir,id) for id,l,g in data]
    dataset = zip(parafiles, surffiles, lowresfiles, lowresparams, highresfiles, highresparams)

    def file_exists(f):
        return not opts.overwrite and os.path.isfile(f) and os.stat(f).st_size > 0

    for idx,(pf,sf,lm,lp,hm,hp) in enumerate(dataset):
        if not file_exists(lm):
            if idx == 0:
                cmd = "%s %s %s --subdivLevel 10 --spharmDegree 20  %s/%s.ip. --paraOut" % (pathParaToSPHARM, pf, sf, workdir, idlist[idx])
            else:
                cmd = "%s %s %s --subdivLevel 10 --spharmDegree 20  %s/%s.ip. --flipTemplateOn --flipTemplate %s/%s.ip.SPHARM.coef --paraOut" % (pathParaToSPHARM, pf, sf, workdir, idlist[idx], workdir, idlist[idx-1])
            system.run_process(cmd,verbose=opts.verbose)

        if not file_exists(hm):
            if idx == 0:
                cmd = "%s %s %s --subdivLevel 50 --spharmDegree 20  %s/%s.subj. --paraOut" % (pathParaToSPHARM, pf, sf, workdir, idlist[idx])
            else:
                cmd = "%s %s %s --subdivLevel 50 --spharmDegree 20  %s/%s.subj. --flipTemplateOn --flipTemplate %s/%s.subj.SPHARM.coef --paraOut" % (pathParaToSPHARM, pf, sf, workdir, idlist[idx], workdir, idlist[idx-1])
            system.run_process(cmd,verbose=opts.verbose)


def resample_segmentations(config,data,outputdir):
    workdir = "%s/Processing/1.MeasurementandSPHARM" % (outputDir)
    reference_image = "%s/%s.segpost.nii.gz" % (workdir,data[0][0])

    pathResampleVolume = config["ResampleVolume2Path"]
    for id,l,g in data:
        segpost_image = "%s/%s.segpost.nii.gz" % (workdir,id)
        labelMapNCResampled = "%s/%s.labelMapNCsegmentation.nrrd" % (workdir,id)
        cmd = "%s %s %s -R %s -i nn" % (pathResampleVolume, \
                                        segpost_image, labelMapNCResampled, reference_image)
        system.run_process(cmd,verbose=True)


def setup_particle_tools(config, data, outputDir, ids):
    inputSurfaceFiles = open("%s/Processing/2.shapeworks/inputsurfacemodels.txt" % (outputDir), 'w')
    inputImageFiles = open("%s/Processing/2.shapeworks/inputimage.txt" % (outputDir), 'w')
    outputCorrespondenceFiles = open("%s/Processing/2.shapeworks/outputCorrespondence.txt" % (outputDir), 'w')
    outputWarpedFiles = open("%s/Processing/2.shapeworks/outputWarped.txt" % (outputDir), 'w')

    workdir = "%s/Processing/1.MeasurementandSPHARM" % (outputDir)
    lowresfiles = ["%s/%s.ip.SPHARM.vtk"%(workdir,id) for id,l,g in data]
    highresfiles = ["%s/%s.subj.SPHARM.vtk"%(workdir,id) for id,l,g in data]
    shapedir = "%s/Processing/2.shapeworks" % (outputDir)
    correspondencefiles = ["%s/%s.correspondence.vtk"%(shapedir,id) for id,l,g in data]
    warpedfiles = ["%s/%s.warped.vtk"%(shapedir,id) for id,l,g in data]

    # file generation
    inputSurfaceFiles.write("\n".join(lowresfiles))
    inputImageFiles.write("\n".join(segmentationfiles))
    outputCorrespondenceFiles.write("\n".join(correspondencefiles))
    outputWarpedFiles.write("\n".join(warpedfiles))

    inputSurfaceFiles.close()
    inputImageFiles.close()
    outputCorrespondenceFiles.close()
    outputWarpedFiles.close()



# @brief compute consistent boundary conditions from surface maps
def regenerate_segmentations(data, config, outputdir):
    pathKmesh = config["kmeshPath"]
    pathKcalc = config["kcalcPath"]
    workdir = "%s/Processing/1.MeasurementandSPHARM" % (outputdir)


    # before executing the boundary correction
    # check if none of input is not generated.
    labelMapList = [labelMap for (tag, labelMap, group) in data]
    testLabelMapsFail = False in [os.path.exists(x) for x in labelMapList]
    if testLabelMapsFail:
        raise RuntimeError("There are missing label maps")

    surfaceInputList = ["%s/%s.subj.SPHARM.vtk" % (workdir, tag) for (tag,labelMap,group) in data]
    testSurfaceInputsFail = False in [os.path.exists(x) for x in surfaceInputList]
    if testSurfaceInputsFail:
        raise RuntimeError("There are missing surface inputs. Check SPHARM results!")

    for (idx, item) in enumerate(data):
        tag = item[0]
        labelMap = item[1]
        group = item[2]

        labelOutput = workdir + "/" + tag + "/labelmap_zeroSolution.nrrd"
        surfaceInput = workdir + "/" + tag + ".subj.SPHARM.vtk"
        surfaceLabels = workdir + "/" + tag + ".labels.vtp"
        voronoiImage = workdir + "/" + tag + ".voronoi.mha"
        surfaceImage = workdir + "/" + tag + ".solution.mha"
        boundaryMap = workdir + "/" + tag + ".boundaryMap.mha"

        cmd = "%s -e 'A==3?0:A' -o %s %s" % (pathKcalc, labelOutput, labelMap)
        system.run_process(cmd,verbose=True)

        cmd = "%s -sampleImage -zrotate %s %s %s -outputScalarName labels" % (pathKmesh, labelOutput, surfaceInput, surfaceLabels)
        system.run_process(cmd,verbose=True)

    cmd = "%s -averageScalars -threshold 1.8 -scalarName labels -outputScalarName meanLabels" % (pathKmesh)
    for (tag, labelmap, group) in data:
        cmd += " " + workdir + "/" + tag + ".labels.vtp"
    print cmd
    system.run_process(cmd,verbose=True)

    for (idx, item) in enumerate(data):
        tag = item[0]
        labelMap = item[1]
        group = item[2]

        labelOutput = workdir + "/" + tag + "/labelmap_zeroSolution.nrrd"
        surfaceInput = workdir + "/" + tag + ".subj.SPHARM.vtk"
        surfaceLabels = workdir + "/" + tag + ".labels.vtp"
        voronoiImage = workdir + "/" + tag + ".voronoi.mha"
        surfaceImage = workdir + "/" + tag + ".solution.mha"
        boundaryMap = workdir + "/" + tag + ".boundaryMap.mha"

        cmd = "%s -voronoiImage -zrotate %s %s %s -scalarName meanLabels" % (pathKmesh, labelOutput, surfaceLabels, voronoiImage)
        print cmd
        system.run_process(cmd,verbose=True)

        cmd = "%s -scanConversion -zrotate %s %s %s" % (pathKmesh, surfaceLabels, labelOutput, surfaceImage)
        print cmd
        system.run_process(cmd,verbose=True)

        cmd = "%s -e 'B>0?3:A' -o %s %s %s" % (pathKcalc, boundaryMap, voronoiImage, surfaceImage)
        print cmd
        system.run_process(cmd,verbose=True)

def run_particle_tools(config, data, outputdir, ids):
    pathPython = config["PythonPath"]
    pathShapeWorksScript = config["ShapeWorksPythonScriptPath"]
    pathShapeWorksRun = config["ShapeWorksRunPath"]
    pathImageMath = config["ImageMathPath"]
    pathShapeWorksGroom = config["ShapeWorksGroomPath"]
    pathBinaryToDistanceMap = config["BinaryToDistanceMapPath"]

    workdir = "%s/Processing/2.shapeworks" % (outputdir)

    # generate binary distance map for each label map
    # cmd = "%s %s %s/inputimage.txt %s/inputsurfacemodels.txt -c %s/outputCorrespondence.txt -w %s/outputWarped.txt --workingDir %s --pathShapeWorksRun %s --pathShapeWorksGroom %s --pathImageMath %s --pathBinaryToDistanceMap %s" % (pathPython, pathShapeWorksScript, workdir, workdir, workdir, workdir, workdir, pathShapeWorksRun, pathShapeWorksGroom, pathImageMath, pathBinaryToDistanceMap)
    cmd = "%s %s %s/inputimage.txt %s/inputsurfacemodels.txt " + \
            "-c %s/outputCorrespondence.txt -w %s/outputWarped.txt " +\
            "--workingDir %s" % (pathPython, pathShapeWorksScript, \
            workdir, workdir, workdir, workdir, workdir)
    system.run_process(cmd,verbose=True)


def performSplit(args, outputPattern):
    # read each argument file
    for inputFile in args:
        fin = open(inputFile)
        lines = fin.readlines()
        outputLines = {}
        # iterate over the entire lines
        for line in lines:
            words = line.split("\t")
            regionId = int(words[0])
            # skip if the region is 0
            if (regionId == 0):
                continue
            # otherwise, add to outputLines separately
            if (regionId in outputLines):
                outputLines[regionId].append(line)
            else:
                outputLines[regionId] = [ line ]
        fin.close()
        # write all region lines
        for region in outputLines:
            ofile = outputPattern % (region)
            fout = open(ofile, "w")
            for line in outputLines[region]:
                fout.write(line)
            fout.close()


# @brief Main function
if (__name__ == "__main__"):
    parser = OptionParser(usage="usage: %prog dataset.csv config.bms output-directory")
    #  parser.add_option("--labelprocessing", help="preprocessing for label map generation", action="store_true", dest="labelprocessing")
    parser.add_option("--path", help="add directories to PATH environment (e.g. --path=/usr/bin:/home/bin)", destion="dirs", metavar="PATH", default="")
    parser.add_option("--genparamesh", help="generate SPHARM parameter and surface mesh", action="store_true", dest="genparamesh")
    parser.add_option("--genparamesh-iter", metavar="NUM", help="set the number of iteration for GenParaMesh (default: 100)", type="int", dest="genparamesh_iter", default=100)
    parser.add_option("--paratospharm", help="reconstruct surface mesh from SPHARM parameters", action="store_true", dest="paratospharm")
    parser.add_option("--resample-segmenations", help="resample segmentation images to have consistent dimensions", action="store_true", dest="resample_segmentations")
    parser.add_option("--run-particle-tools", help="run particle correspondence", action="store_true", dest="run_particle_tools")
    parser.add_option("--regenerate-segmentations", dest="regenerate_segmentations", help="Regenerate segmentation images based on reconstructed mesh", action="store_true")
    parser.add_option("--compute-thickness", help="compute thickness for a given label", action="store_true", dest="compute_thickness")
    parser.add_option("--compute-statistics", help="compute statistics for hypothesis testing", action="store_true", dest="compute_statistics")
    parser.add_option("--run-preprocessing", help="generate surface meshes and regenerate consistent labelmaps", action="store_true", dest="run_processing")
    parser.add_option("--run-all", help="run all pipeline steps", action="store_true", dest="run_all")

    parser.add_option("--ids", metavar="3", help="solution label", dest="ids", default="3")
    parser.add_option("--idl", metavar="2", help="low boundary label", dest="idl", default="2")
    parser.add_option("--idh", metavar="1", help="high boundary label", dest="idh", default="1")

    parser.add_option("--labelNCleft", help="label id for NC left", dest="labelNCleft", default="30")
    parser.add_option("--labelNCright", help="label id for NC right", dest="labelNCright", default="14")
    parser.add_option("--labelOB", help="label id for Olfactory Bulb", dest="labelOB", default="16")

    parser.add_option("--regionSplit", help="split the data file with its region (1st column)", dest="regionSplit", action="store_true")
    parser.add_option("--outputPattern", help="Specify the output pattern ex) --outputPattern region_control_%02d.txt", dest="outputPattern")
    parser.add_option("--overwrite", help="specify if the script overwrites previous results", action="store_true", dest="overwrite")
    parser.add_option("--log-file", help="Specify the filename for logging stdout and stderr outputs", dest="logfileName", default=None)
    parser.add_option("--log-commands", help="Specify if the script logs only command lines not executing those", action="store_true", dest="log_commands")
    parser.add_option("--verbose", help="Print out logs to cosnole", action="store_true", dest="verbose", default=False)

    (opts, args) = parser.parse_args()



    # set up the logger for external programs
    if opts.log_commands:
        logformat = "%(message)s"
    else:
        logformat = "%(levelname)s %(asctime)s %(message)s"

    if opts.logfileName is not None:
        system.rename_logfile(opts.logfileName)
        system.setup_logger(logfilename=opts.logfileName, logformat=logformat)
    else:
        system.setup_logger(logformat=logformat)


    # argument processing
    if (len(args) < 2):
        parser.print_help()
    else:
        csvfile = args[0]
        outputdir = args[1]
        #bmsfile = args[1]

        if (not os.path.isdir("%s/Processing/1.MeasurementandSPHARM" % (outputdir))):
            os.makedirs("%s/Processing/1.MeasurementandSPHARM" % (outputdir))

        if (not os.path.isdir("%s/Processing/2.shapeworks" % (outputdir))):
            os.makedirs("%s/Processing/2.shapeworks" % (outputdir))

        if (not os.path.isdir("%s/Processing/3.meshintensity" % (outputdir))):
            os.makedirs("%s/Processing/3.meshintensity" % (outputdir))

        if (not os.path.isdir("%s/Statistical" % (outputdir))):
            os.makedirs("%s/Statistical" % (outputdir))

        programs_list = ["ResampleVolume2", "kmesh", "ImageMath", \
                  "MeshPointsIntensitySampling", "SegPostProcessCLP", \
                  "measureThicknessFilter", "ParaToSPHARMMeshCLP", "GenParaMeshCLP", \
                  "Rscript", "BinaryToDistanceMap", "kcalc", \
                  "ShapeWorksRun", "ShapeWorksGroom"]
 
        newpath = sys.path[0]+"/../Rscript-build/bin"+((":"+opts.dirs) if opts.dirs else "")
        config = configure_paths(newpath, programs_list, opts.verbose)

        csvreader = csv.reader(open(csvfile, "r"))
        csvheader = csvreader.next()

        csvdata = []
        dataIds = []
        dataGroups = []

        for csvline in csvreader:
            dataIds.append(csvline[0])
            dataGroups.append(csvline[2])
            csvdata.append(csvline)

        if (opts.verbose):
            print "Verbose=On"

        if (opts.log_commands):
            system.set_log_commands(True)

        if (opts.regionSplit):
            performSplit(args, opts.outputPattern)
            # just finish the script
            sys.exit(0)


        if (opts.genparamesh):
            genparamesh(opts, config, csvdata, outputdir)

        if (opts.paratospharm):
            paratospharm(opts, config, csvdata, outputdir)

        if (opts.resample_segmentations):
            resample_segmentations(opts,config,csvdata,outputdir)

        if (opts.regenerate_segmentations):
            regenerate_segmentations(csvdata, config, outputdir)


          
        if (opts.compute_thickness):
            compute_thickness(config, csvdata, outputdir, opts.ids, opts.idl, opts.idh)

        if (opts.run_particle_tools):
            setup_particle_tools(config, csvdata, outputdir, opts.ids)
            run_particle_tools(config, csvdata, outputdir, opts.ids)

        if (opts.compute_statistics):
            performAnalysis(csvdata, config, outputdir, "initial_dense")
            # performAnalysis(csvdata, config, outputdir, "spharm_sampling")
            # performAnalysis(csvdata, config, outputdir, "correspondence")

        if (opts.preprocessing):
            genparamesh(opts, config, csvdata, outputdir)
            paratospharm(opts, config, csvdata, outputdir)
            resample_segmentations(opts,config,csvdata,outputdir)
            regenerate_segmentations(csvdata, config, outputdir)

        if (opts.run_all):
            genparamesh(opts, config, csvdata, outputdir)
            paratospharm(opts, config, csvdata, outputdir)
            resample_segmentations(opts,config,csvdata,outputdir)
            regenerate_segmentations(csvdata, config, outputdir)
            compute_thickness(config, csvdata, outputdir, opts.ids, opts.idl, opts.idh)
            setup_particle_tools(config, csvdata, outputdir, opts.ids)
            run_particle_tools(config, csvdata, outputdir, opts.ids)
            performAnalysis(csvdata, config, outputdir, "initial_dense")
